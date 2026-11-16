"""
MODULE:     t.graminr (graminr_lib.grasprog)

AUTHOR(S):  Yann Chemin

PURPOSE:    GrasProg (FOPROQ) grassland model, ported from the GraMinR R
            package (R/grasprog.R) and vectorised over pixels.

SPDX-License-Identifier: Unlicense
"""

import numpy as np

from .core import as_array, limit01, safe_div

OUTPUTS = (
    "w",
    "age",
    "rs",
    "sw",
    "lai",
    "ri",
    "ti",
    "wi",
    "gi",
    "eta",
    "etp",
    "wt",
    "yield",
    "tempsum",
    "di",
    "irrig",
)

REQUIRED = ("tavg", "rad_par", "rain", "et_pot")
OPTIONAL = {"gcut": 0.0, "growth_flag": 1.0}

ALIASES = {
    "Tmin": "tmin",
    "Topt": "topt",
    "Tmax": "tmax",
    "LAI50": "lai50",
    "a": "a_lai",
    "b": "b_lai",
    "c": "c_lai",
    "Rlow": "rlow",
    "Rhigh": "rhigh",
    "LAIlow": "lailow",
    "LAIhigh": "laihigh",
    "PAW": "paw",
}

DEFAULTS = {
    "begsum": 250.0,
    "begbas": 0.0,
    "swthr": 0.8,
    "droughtindex": "yes",
    "waterstress": "yes",
}

#: Parameters collected into per-regrowth-compartment vectors.
COMPARTMENT_PREFIXES = ("rs", "ak", "w0")


def _estimate_rs(week):
    """Initial relative growth rate from the calendar week (Eq. 17)."""
    return np.minimum(1.0, 1.34 * np.exp(-0.0558 * week))


def _estimate_w0(week):
    """Start-of-growth-period biomass from the calendar week (Eq. 18)."""
    w0 = -257.0 + 40.3 * week - 2.21 * week**2 + 0.0508 * week**3 - 0.000413 * week**4
    return 10.0 * np.maximum(10.0, w0)


def _compartment(table, cutnum, fallback):
    """Value of a per-compartment parameter for the current cut count.

    A cut count beyond the last supplied compartment reuses the last one,
    mirroring the original's replaceParams(): a DMn block that does not
    redefine a value simply leaves the previous compartment's in place.
    """
    if table is None:
        return fallback
    tbl = np.asarray(table, dtype=np.float64).ravel()
    idx = np.clip(np.asarray(cutnum, dtype=np.int64) - 1, 0, tbl.size - 1)
    return tbl[idx]


def _calculate_ri(radiation, lai, lailow, laihigh, dr, cr, rlow, rhigh):
    """Radiation index (Eqs 8-10)."""
    c_lai = np.where(
        lai <= lailow,
        0.0,
        np.where(
            lai >= laihigh,
            1.0,
            1.0
            - np.exp(-dr * safe_div(lai - lailow, laihigh, 0.0)) / (1.0 - np.exp(-dr)),
        ),
    )
    r_opt = rlow + c_lai * (rhigh - rlow)
    ratio = safe_div(
        1.0 - np.exp(-cr * safe_div(radiation, r_opt, 0.0)), 1.0 - np.exp(-cr), 0.0
    )
    ri = np.where(radiation >= rhigh, 1.0, np.where(r_opt != 0, ratio, 0.0))
    return np.where(np.isfinite(ri), ri, 0.0)


def _calculate_ti(temp, tmin, topt, tmax, ct):
    """Temperature index (Eqs 11-12)."""
    z = np.where(
        temp < topt,
        np.abs(temp - topt) * safe_div(1.0, topt - tmin, 0.0),
        np.abs(temp - topt) * safe_div(1.0, topt, 0.0),
    )
    ti = np.where(
        z < 0.5,
        1.0 - 0.5 * (2.0 * z) ** ct,
        np.where(z < 1.0, 0.5 * (2.0 * (1.0 - z)) ** ct, 0.0),
    )
    inside = np.isfinite(temp) & (temp >= tmin) & (temp <= tmax)
    return np.where(inside, ti, 0.0)


def run(params, weather, npixels, outputs=None, progress=None):
    """Run GrasProg over `npixels` pixels."""
    p = dict(params)
    n = npixels
    year = np.asarray(weather["year"])
    doy = np.asarray(weather["doy"])
    ndays = len(doy)
    wanted = tuple(OUTPUTS) if outputs is None else tuple(outputs)

    tavg = weather["tavg"]
    rad_par = weather["rad_par"]
    rain = weather["rain"]
    et_pot = weather["et_pot"]
    gcut = weather["gcut"]
    growth_flag = weather.get("growth_flag")

    def g(key):
        """Parameter `key` broadcast over the pixel axis."""
        return as_array(p[key], n)

    tmin, topt, tmax = g("tmin"), g("topt"), g("tmax")
    lai50 = g("lai50")
    a_lai, b_lai, c_lai = g("a_lai"), g("b_lai"), g("c_lai")
    rlow, rhigh = g("rlow"), g("rhigh")
    lailow, laihigh = g("lailow"), g("laihigh")
    dr, cr, ct = g("dr"), g("cr"), g("ct")
    paw = g("paw")
    begsum, begbas, swthr = g("begsum"), g("begbas"), g("swthr")

    rs_table = p.get("rs_compartments")
    ak_table = p.get("ak_compartments")
    w0_table = p.get("w0_compartments")

    drought_on = str(p.get("droughtindex", "yes")).lower() == "yes"
    waterstress_on = str(p.get("waterstress", "yes")).lower() == "yes"

    cutnum = np.ones(n, dtype=np.int64)
    w_carry = np.zeros(n)
    sw_carry = paw.copy()
    nfk_counter = np.zeros(n)
    di_carry = np.ones(n)
    tempsum = np.zeros(n)

    out = {name: np.zeros((ndays, n), dtype=np.float32) for name in wanted}

    def store(name, value):
        if name in out:
            out[name][i] = value

    for i in range(ndays):
        week = float(round(doy[i] / 7.0))

        # A fresh start: the first simulated day, or a new calendar year.
        fresh = (i == 0) or (doy[i] == 1)
        if fresh:
            cutnum = np.ones(n, dtype=np.int64)
            w_carry = _compartment(w0_table, cutnum, _estimate_w0(week)) * np.ones(n)
            sw_carry = paw.copy()
            nfk_counter = np.zeros(n)
            di_carry = np.ones(n)
            tempsum = np.maximum(0.0, tavg[i] - begbas)
        else:
            tempsum = tempsum + np.maximum(0.0, tavg[i] - begbas)

        started = tempsum >= begsum

        ak_current = _compartment(ak_table, cutnum, a_lai)
        rs_cut = _compartment(rs_table, cutnum, _estimate_rs(week))

        lai = b_lai * (1.0 - np.exp(-c_lai * (w_carry / 10.0)))
        lai = np.where(np.isfinite(lai), lai, 0.0)

        etp = et_pot[i] * np.minimum(1.0, (9.0 + 7.0 * lai) / 30.0)
        eta = etp * np.minimum(1.0, safe_div(sw_carry, swthr * paw, 0.0))

        ri = limit01(
            _calculate_ri(rad_par[i], lai, lailow, laihigh, dr, cr, rlow, rhigh)
        )
        ti = limit01(_calculate_ti(tavg[i], tmin, topt, tmax, ct))
        wi = limit01(np.where(etp == 0, 1.0, safe_div(eta, etp, 0.0)))

        # The original lets the water balance go negative under a sustained
        # deficit: there is deliberately no floor at 0 here.
        sw_new = np.minimum(paw, sw_carry + rain[i] - eta)
        irrig = np.zeros(n)
        if not waterstress_on:
            irrig = paw - sw_new
            sw_new = paw.copy()

        # Drought index: a lockout of at least 14 days once soil water has
        # been below 30% of usable capacity on two consecutive days.
        if drought_on:
            low_twice = (safe_div(sw_new, paw, 1.0) < 0.3) & (
                safe_div(sw_carry, paw, 1.0) < 0.3
            )
            nfk_counter = np.where(
                low_twice, np.maximum(14.0, nfk_counter + 1.0), nfk_counter
            )
            locked = (nfk_counter >= 14) | ((nfk_counter > 0) & (di_carry < 1))
            di = np.where(locked, 0.0, 1.0)
            nfk_counter = np.where(locked, nfk_counter - 1.0, nfk_counter)
            di_carry = di
        else:
            di = np.ones(n)
            nfk_counter = np.zeros(n)

        gi = limit01(ri) * limit01(ti) * limit01(wi) * di
        age = np.maximum(
            0.0, safe_div(1.0, 1.0 + (safe_div(lai, lai50, 0.0)) ** ak_current, 1.0)
        )
        rs = rs_cut * age
        rs = np.where(np.isfinite(rs), rs, 0.0)
        wt = w_carry * rs * gi
        wt = np.where(np.isfinite(wt), wt, 0.0)

        if growth_flag is None:
            flag_on = np.ones(n, dtype=bool)
        else:
            gf = growth_flag[i]
            flag_on = ~np.isfinite(gf) | (gf == 1)

        w0_now = _compartment(w0_table, cutnum, _estimate_w0(week)) * np.ones(n)
        w_grown = np.where(started & flag_on, w_carry + wt, w0_now)

        # On a cut day, today's row holds the grown, pre-cut biomass; the new
        # compartment's floor only applies from tomorrow.
        cut_today = started & (gcut[i] == 1)
        cutnum_next = np.where(cut_today, cutnum + 1, cutnum)
        w0_next = _compartment(w0_table, cutnum_next, _estimate_w0(week)) * np.ones(n)
        yield_today = np.where(cut_today, np.maximum(0.0, w_grown - w0_next), 0.0)
        w_carry = np.where(cut_today, w0_next, w_grown)
        cutnum = cutnum_next
        sw_carry = sw_new

        store("w", w_grown)
        store("age", age)
        store("rs", rs)
        store("sw", sw_new)
        store("lai", lai)
        store("ri", ri)
        store("ti", ti)
        store("wi", wi)
        store("gi", gi)
        store("eta", eta)
        store("etp", etp)
        store("wt", wt)
        store("yield", yield_today)
        store("tempsum", tempsum)
        store("di", di)
        store("irrig", irrig)

        if progress is not None:
            progress(i, ndays)

    return out
