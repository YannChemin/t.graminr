"""
MODULE:     t.graminr (graminr_lib.modvege)

AUTHOR(S):  Yann Chemin

PURPOSE:    ModVege grassland model (Jouven et al. 2006), ported from the
            GraMinR R package (R/modvege.R, R/modvege_lib.R) and vectorised
            over pixels.

SPDX-License-Identifier: Unlicense
"""

import numpy as np

from .core import as_array, safe_div, trailing_mean, year_cumsum

#: Output variables produced by this model, in reporting order.
OUTPUTS = (
    "gvb",
    "dvb",
    "grb",
    "drb",
    "hb",
    "ib",
    "g",
    "abc",
    "stp",
    "gva",
    "gra",
    "dva",
    "dra",
    "sea",
    "ftm",
    "env",
    "pgr",
    "atr",
    "lai",
    "omdgv",
    "omdgr",
    "omddv",
    "omddr",
)

#: Weather variables this model consumes, with the default used when the
#: driver cannot supply one.
REQUIRED = ("tavg", "rad_par", "rain", "et_pot")
OPTIONAL = {
    "et_act_rs": 0.0,
    "lai_rs": 0.0,
    "gcut": 0.0,
    "graze_count": 0.0,
    "graze_weight": 0.0,
}

#: Parameter aliases of read_config(model = "modvege") in graminr_conf.R.
ALIASES = {
    "ST1": "st1",
    "ST2": "st2",
    "NI": "ni",
    "WHC": "waterHoldingCapacity",
    "WR": "waterReserve",
    "minSEA": "minsea",
    "maxSEA": "maxsea",
    "W_GV": "wgv",
    "alpha_PAR": "alphapar",
    "T0": "t0",
    "T1": "t1",
    "T2": "t2",
    "beta_T": "betaT",
    "SLA": "sla",
    "LLS": "lls",
    "rho_GV": "rhogv",
    "percentLAM": "pctlam",
    "W_GR": "wgr",
    "a_NI": "allocNI",
    "max_fNI": "maxFNI",
    "rho_GR": "rhogr",
    "W_DV": "wdv",
    "K_DV": "kdv",
    "Kl_DV": "kldv",
    "rho_DV": "rhodv",
    "W_DR": "wdr",
    "K_DR": "kdr",
    "Kl_DR": "kldr",
    "rho_DR": "rhodr",
    "init_AGE_GV": "gv_init_age",
    "init_AGE_GR": "gr_init_age",
    "init_AGE_DV": "dv_init_age",
    "init_AGE_DR": "dr_init_age",
    "RUEmax": "ruemax",
    "gammaGV": "gv_gamma",
    "gammaGR": "gr_gamma",
}

#: Parameters with a default, mirroring read_config()'s fallbacks.
DEFAULTS = {
    "DEFAULT_CUT_HEIGHT": 0.05,
    "REP_base": 0.2,
    "REP_coef": 0.8 / 0.65,
    "ET_METHOD": "split",
    "cellSurface": 0.01,
    "maxAmountToIngest": 0.0,
}


def f_water_stress(water_reserve, whc, pet):
    """Water-limitation factor fW, Jouven et al. (2006) Eq. 5."""
    ws = np.minimum(safe_div(water_reserve, whc, 0.0), 1.0)
    low = np.where(
        ws <= 0.2,
        4.0 * ws,
        np.where(
            ws <= 0.4,
            0.75 * ws + 0.65,
            np.where(ws <= 0.6, 0.25 * ws + 0.85, 1.0),
        ),
    )
    mid = np.where(
        ws <= 0.2,
        2.0 * ws,
        np.where(
            ws <= 0.4,
            1.5 * ws + 0.1,
            np.where(
                ws <= 0.6,
                ws + 0.3,
                np.where(ws <= 0.8, 0.5 * ws + 0.6, 1.0),
            ),
        ),
    )
    return np.where(pet <= 3.8, low, np.where(pet <= 6.5, mid, ws))


def f_temperature(mean_ten_days_t, today_t, t0, t1, t2):
    """Temperature-limitation factor fT, gated by the trailing 10-day mean."""
    shape = np.where(
        today_t < t1,
        safe_div(today_t - t0, t1 - t0, 0.0),
        np.where(today_t < t2, 1.0, np.maximum(0.0, 1.0 - 0.05 * (today_t - t2))),
    )
    gate = (
        np.isfinite(mean_ten_days_t)
        & np.isfinite(today_t)
        & (mean_ten_days_t > t0)
        & (today_t > t0)
    )
    return np.where(gate, shape, 0.0)


def f_sea(maxsea, minsea, sum_t, st2, st1):
    """Seasonal-effect factor fSEA."""
    return np.where(
        (sum_t < 200) | (sum_t >= st2),
        minsea,
        np.where(
            sum_t < st1 - 200,
            minsea + (maxsea - minsea) * safe_div(sum_t - 200, st1 - 200, 0.0),
            np.where(
                sum_t < st1 - 100,
                maxsea,
                maxsea
                + (minsea - maxsea) * safe_div(sum_t - st1 + 100, st2 - st1 + 100, 0.0),
            ),
        ),
    )


def f_par_i(par_i, alphapar):
    """Radiation-saturation factor fPARi."""
    return np.where(par_i < 5, 1.0, np.maximum(1.0 - alphapar * (par_i - 5), 0.0))


def f_clai(pctlam, sla, gv_biomass):
    """LAI derived from green vegetative biomass."""
    return sla * (gv_biomass / 10.0) * pctlam


def potential_growth(par_i, ruemax, pctlam, sla, gv_biomass, lai):
    """Potential growth PGRO, light interception by a Beer-Lambert canopy."""
    lai = np.where(lai <= 0, sla * pctlam * (gv_biomass / 10.0), lai)
    return par_i * ruemax * (1.0 - np.exp(-0.6 * lai)) * 10.0


def actual_et(pet, pctlam, sla, gv_biomass, water_reserve, whc, lai, method):
    """Actual evapotranspiration, 'split' (default) or legacy 'linear'."""
    lai = np.where(lai <= 0, sla * pctlam * (gv_biomass / 10.0), lai)
    if method == "linear":
        threshold = safe_div(3.0 * 10.0, sla * pctlam, np.inf)
        return np.where(
            gv_biomass > threshold, pet, pet * safe_div(gv_biomass, threshold, 0.0)
        )
    light = 1.0 - np.exp(-0.6 * lai)
    pt = pet * light
    pe = pet - pt
    ta = pt * f_water_stress(water_reserve, whc, pet)
    ea = pe * np.minimum(safe_div(water_reserve, whc, 0.0), 1.0)
    return ta + ea


def _age_class_three(ratio):
    """Discrete 1/2/3 age class used by the dead compartments."""
    return np.where(ratio < 1.0 / 3.0, 1.0, np.where(ratio < 2.0 / 3.0, 2.0, 3.0))


def _age_class_ramp(ratio):
    """Continuous age factor used by the green compartments."""
    return np.where(ratio < 1.0 / 3.0, 1.0, np.where(ratio < 1.0, 3.0 * ratio, 3.0))


def _exe_cut(cut_height, bulk_density, biomass):
    """Cut down to cut_height; returns (harvested, remaining)."""
    after = cut_height * bulk_density * 10.0
    remaining = np.minimum(biomass, after)
    return biomass - remaining, remaining


def _available_above(biomass, cut_height, bulk_density):
    return np.maximum(0.0, biomass - cut_height * bulk_density * 10.0)


def run(params, weather, npixels, outputs=None, progress=None):
    """Run ModVege over `npixels` pixels.

    :param params: dict of parameters (scalars or per-pixel arrays)
    :param weather: dict of (ndays, npixels) float arrays, plus 'year'/'doy'
    :param npixels: number of pixels
    :param outputs: iterable of wanted output names (default: all)
    :param progress: callable(day_index, ndays) for progress reporting
    :return: dict name -> (ndays, npixels) float32 array
    """
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
    et_act_rs = weather.get("et_act_rs")
    lai_rs = weather.get("lai_rs")
    gcut = weather["gcut"]
    graze_count = weather.get("graze_count")
    graze_weight = weather.get("graze_weight")

    # Parameters, broadcast to the pixel axis.
    def g(key):
        """Parameter `key` broadcast over the pixel axis."""
        return as_array(p[key], n)

    st1, st2 = g("st1"), g("st2")
    t0, t1, t2 = g("t0"), g("t1"), g("t2")
    sla, pctlam = g("sla"), g("pctlam")
    lls = g("lls")
    kdv, kldv, kdr, kldr = g("kdv"), g("kldv"), g("kdr"), g("kldr")
    rhogv, rhodv, rhogr, rhodr = g("rhogv"), g("rhodv"), g("rhogr"), g("rhodr")
    ruemax, alphapar = g("ruemax"), g("alphapar")
    minsea, maxsea = g("minsea"), g("maxsea")
    whc = g("waterHoldingCapacity")
    gv_gamma, gr_gamma = g("gv_gamma"), g("gr_gamma")
    rep_base, rep_coef = g("REP_base"), g("REP_coef")
    default_cut = g("DEFAULT_CUT_HEIGHT")
    cell_surface = g("cellSurface")
    max_ingest = g("maxAmountToIngest")
    min_omd_gv, max_omd_gv = g("minOMDgv"), g("maxOMDgv")
    min_omd_gr, max_omd_gr = g("minOMDgr"), g("maxOMDgr")
    mean_omd_dv, mean_omd_dr = g("meanOMDdv"), g("meanOMDdr")
    et_method = str(p.get("ET_METHOD", "split"))

    # Nitrogen nutrition index, floored at 0.35 (Belanger et al. 1994).
    ni = np.maximum(g("ni"), 0.35)

    # State variables.
    gv_biomass = g("wgv")
    gr_biomass = g("wgr")
    dv_biomass = g("wdv")
    dr_biomass = g("wdr")
    gv_avg_age = g("gv_init_age")
    gr_avg_age = g("gr_init_age")
    dv_avg_age = g("dv_init_age")
    dr_avg_age = g("dr_init_age")
    water_reserve = g("waterReserve")
    harvested = np.zeros(n)
    ingested = np.zeros(n)
    is_cut = np.zeros(n, dtype=bool)
    a2r_flag = np.zeros(n, dtype=bool)

    # Degree-day sum: cumulated max(tavg, 0), reset every calendar year.
    sum_t_all = year_cumsum(np.maximum(tavg, 0.0), year)
    mean10_all = trailing_mean(tavg, 10)

    out = {name: np.zeros((ndays, n), dtype=np.float32) for name in wanted}

    def store(name, value):
        if name in out:
            out[name][i] = value

    for i in range(ndays):
        temperature = tavg[i]
        mean10 = mean10_all[i]
        sum_t = sum_t_all[i]
        cut_height = gcut[i]

        sea = f_sea(maxsea, minsea, sum_t, st2, st1)
        ftm = f_temperature(mean10, temperature, t0, t1, t2)

        is_harvested = cut_height != 0
        if graze_count is None or graze_weight is None:
            is_grazed = np.zeros(n, dtype=bool)
        else:
            is_grazed = (graze_count[i] != 0) & (graze_weight[i] != 0)
        # isCut latches until a day with neither a harvest nor grazing.
        is_cut = np.where(is_grazed | is_harvested, is_cut, False)

        # Evapotranspiration and LAI: remote sensing when available.
        lai = (
            np.zeros(n)
            if lai_rs is None
            else np.asarray(lai_rs[i], dtype=np.float64).copy()
        )
        eta = (
            np.zeros(n)
            if et_act_rs is None
            else np.asarray(et_act_rs[i], dtype=np.float64).copy()
        )
        lai = np.where(~np.isfinite(lai), 0.0, lai)
        eta = np.where(~np.isfinite(eta), 0.0, eta)
        need_lai = (eta <= 0) & (lai <= 0)
        lai = np.where(need_lai, f_clai(pctlam, sla, gv_biomass), lai)
        eta = np.where(
            eta <= 0,
            actual_et(
                et_pot[i], pctlam, sla, gv_biomass, water_reserve, whc, lai, et_method
            ),
            eta,
        )

        water_reserve = np.minimum(np.maximum(0.0, water_reserve + rain[i] - eta), whc)

        # Senescence and abscission are taken from the pre-cut biomass.
        gv_precut, dv_precut = gv_biomass, dv_biomass
        gr_precut, dr_precut = gr_biomass, dr_biomass

        harvested_today = np.zeros(n)
        ingested_today = np.zeros(n)

        if np.any(is_harvested):
            gv_h, gv_rem = _exe_cut(cut_height, rhogv, gv_biomass)
            dv_h, dv_rem = _exe_cut(cut_height, rhodv, dv_biomass)
            gr_h, gr_rem = _exe_cut(cut_height, rhogr, gr_biomass)
            dr_h, dr_rem = _exe_cut(cut_height, rhodr, dr_biomass)
            total_h = (gv_h + dv_h + gr_h + dr_h) * cell_surface
            harvested_today = np.where(is_harvested, total_h, 0.0)
            gv_biomass = np.where(is_harvested, gv_rem, gv_biomass)
            dv_biomass = np.where(is_harvested, dv_rem, dv_biomass)
            gr_biomass = np.where(is_harvested, gr_rem, gr_biomass)
            dr_biomass = np.where(is_harvested, dr_rem, dr_biomass)
            is_cut = is_cut | is_harvested

        if np.any(is_grazed):
            # NB: GraMinR's R defoliation() is unreachable (it references an
            # undefined 'maxAmountToIngest') and never subtracts the ingested
            # biomass from the compartments. The intended behaviour is
            # implemented here: ingest up to maxAmountToIngest, shared out in
            # proportion to what each compartment offers above cut_height.
            av_gv = _available_above(gv_biomass, cut_height, rhogv)
            av_dv = _available_above(dv_biomass, cut_height, rhodv)
            av_gr = _available_above(gr_biomass, cut_height, rhogr)
            av_dr = _available_above(dr_biomass, cut_height, rhodr)
            available = av_gv + av_dv + av_gr + av_dr
            take = np.where(is_grazed, np.minimum(available, max_ingest), 0.0)
            share = safe_div(take, available, 0.0)
            gv_biomass = gv_biomass - av_gv * share
            dv_biomass = dv_biomass - av_dv * share
            gr_biomass = gr_biomass - av_gr * share
            dr_biomass = dr_biomass - av_dr * share
            ingested_today = take
            is_cut = is_cut | is_grazed

        # Allocation to the reproductive compartment.
        a2r = np.where(
            (sum_t > st1) & (sum_t < st2), rep_base + rep_coef * (ni - 0.35), 0.0
        )
        a2r_flag = a2r_flag | (is_cut & (sum_t > st1))
        a2r = np.where(a2r_flag, 0.0, a2r)

        env = (
            ftm
            * ni
            * f_par_i(rad_par[i], alphapar)
            * f_water_stress(water_reserve, whc, et_pot[i])
        )
        pgr = potential_growth(rad_par[i], ruemax, pctlam, sla, gv_biomass, lai)
        gro = env * pgr * sea

        # --- Green vegetative ---
        gv_age_ratio = safe_div(gv_avg_age, lls, 0.0)
        sen_gv = np.where(
            temperature > t0,
            kdv * gv_precut * temperature * _age_class_ramp(gv_age_ratio),
            np.where(temperature < 0, kdv * gv_precut * np.abs(temperature), 0.0),
        )
        gv_biomass = gv_biomass - sen_gv
        grow_gv = np.where(temperature > t0, gro * (1.0 - a2r), 0.0)
        gv_avg_age = np.where(
            is_cut,
            gv_avg_age,
            np.where(
                gv_biomass + grow_gv > 0,
                (np.maximum(0.0, temperature) + gv_avg_age)
                * safe_div(gv_biomass, gv_biomass + grow_gv, 0.0),
                0.0,
            ),
        )
        gv_biomass = gv_biomass + grow_gv

        # --- Dead vegetative ---
        abs_dv = np.where(
            temperature > 0,
            kldv
            * dv_precut
            * temperature
            * _age_class_three(safe_div(dv_avg_age, lls, 0.0)),
            0.0,
        )
        dv_biomass = dv_biomass - abs_dv
        grow_dv = (1.0 - gv_gamma) * sen_gv
        dv_avg_age = np.where(
            is_cut,
            dv_avg_age,
            np.where(
                dv_biomass + grow_dv > 0,
                (np.maximum(0.0, temperature) + dv_avg_age)
                * safe_div(dv_biomass, dv_biomass + grow_dv, 0.0),
                0.0,
            ),
        )
        dv_biomass = dv_biomass + grow_dv

        # --- Green reproductive ---
        gr_span = safe_div(gr_avg_age, st2 - st1, 0.0)
        sen_gr = np.where(
            temperature > t0,
            kdr * gr_precut * temperature * _age_class_ramp(gr_span),
            np.where(temperature < 0, kdr * gr_precut * np.abs(temperature), 0.0),
        )
        gr_biomass = gr_biomass - sen_gr
        grow_gr = np.where(temperature > t0, gro * a2r, 0.0)
        gr_avg_age = np.where(
            is_cut,
            gr_avg_age,
            np.where(
                gr_biomass + grow_gr > 0,
                (np.maximum(0.0, temperature) + gr_avg_age)
                * safe_div(gr_biomass, gr_biomass + grow_gr, 0.0),
                0.0,
            ),
        )
        gr_biomass = gr_biomass + grow_gr

        # --- Dead reproductive ---
        abs_dr = np.where(
            temperature > 0,
            kldr
            * dr_precut
            * temperature
            * _age_class_three(safe_div(dr_avg_age, st2 - st1, 0.0)),
            0.0,
        )
        dr_biomass = dr_biomass - abs_dr
        grow_dr = (1.0 - gr_gamma) * sen_gr
        dr_avg_age = np.where(
            is_cut,
            dr_avg_age,
            np.where(
                dr_biomass + grow_dr > 0,
                (np.maximum(0.0, temperature) + dr_avg_age)
                * safe_div(dr_biomass, dr_biomass + grow_dr, 0.0),
                0.0,
            ),
        )
        dr_biomass = dr_biomass + grow_dr

        harvested = harvested + np.where(is_cut, harvested_today, 0.0)
        ingested = ingested + ingested_today

        height = np.where(is_cut, cut_height, default_cut)
        abc = (
            _available_above(gv_biomass, height, rhogv)
            + _available_above(dv_biomass, height, rhodv)
            + _available_above(gr_biomass, height, rhogr)
            + _available_above(dr_biomass, height, rhodr)
        )

        store("gvb", gv_biomass)
        store("dvb", dv_biomass)
        store("grb", gr_biomass)
        store("drb", dr_biomass)
        store("hb", harvested)
        store("ib", ingested)
        store("g", gro)
        store("abc", abc)
        store("stp", sum_t)
        store("gva", gv_avg_age)
        store("gra", gr_avg_age)
        store("dva", dv_avg_age)
        store("dra", dr_avg_age)
        store("sea", sea)
        store("ftm", ftm)
        store("env", env)
        store("pgr", pgr)
        store("atr", a2r)
        store("lai", lai)
        store(
            "omdgv",
            np.maximum(
                min_omd_gv,
                max_omd_gv - gv_avg_age * safe_div(max_omd_gv - min_omd_gv, lls, 0.0),
            ),
        )
        store(
            "omdgr",
            np.maximum(
                min_omd_gr,
                max_omd_gr
                - gr_avg_age * safe_div(max_omd_gr - min_omd_gr, st2 - st1, 0.0),
            ),
        )
        store("omddv", mean_omd_dv)
        store("omddr", mean_omd_dr)

        if progress is not None:
            progress(i, ndays)

    return out
