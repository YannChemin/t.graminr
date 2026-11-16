"""
MODULE:     t.graminr (graminr_lib.grassim)

AUTHOR(S):  Yann Chemin

PURPOSE:    GrasSim grassland model, ported from the GraMinR R package
            (R/grassim.R, R/grassim_lib.R) and vectorised over pixels.

SPDX-License-Identifier: Unlicense
"""

import numpy as np

from .core import as_array, safe_div

OUTPUTS = (
    "gvb",
    "grb",
    "dvb",
    "drb",
    "g",
    "pgr",
    "env",
    "sea",
    "ftm",
    "lai",
    "hb",
    "stp",
    "gva",
    "gra",
    "dva",
    "dra",
    "fw",
    "fn",
    "rnc",
    "nmin",
    "norg",
    "sward_height",
    "tsb",
    "abc",
    "omdgv",
    "omdgr",
)

REQUIRED = ("tavg", "rad_par", "rain", "et_pot")
OPTIONAL = {"gcut": 0.0, "kc": np.nan, "st_grass": np.nan}

ALIASES = {}

DEFAULTS = {
    "DEFAULT_CUT_HEIGHT": 0.05,
    "initial_reproductive_biomass": 0.0,
    "apex_grazed": 0.0,
    "notRunoff": 0.0,
    "ST_Tmin": 0.0,
    "ST_Tmax": 18.0,
    "ST_start": 32.0,
    "RNC_fixed": np.nan,
}

REQUIRED_PARAMS = (
    "SLA",
    "percentageLAM",
    "ST1",
    "ST2",
    "maxSEA",
    "minSEA",
    "LLS",
    "maxOMDGV",
    "minOMDGV",
    "maxOMDGR",
    "minOMDGR",
    "BDGV",
    "BDDV",
    "BDGR",
    "BDDR",
    "RUEmax",
    "sigmaGV",
    "sigmaGR",
    "T0",
    "T1",
    "T2",
    "Tlimit",
    "STmin",
    "KGV",
    "KGR",
    "KlDV",
    "KlDR",
    "Tref",
    "K",
    "percentageofNmin",
    "NH3volatfactor",
    "a_Ncrit",
    "b_Ncrit",
    "a_Nmax",
    "b_Nmax",
    "FNH_coef1",
    "FNH_coef2",
    "FNAmax",
    "NSc",
    "RNCmax",
    "RNCmin",
    "FNHmax",
    "Initial_height",
    "AGEGV",
    "AGEGR",
    "AGEDV",
    "AGEDR",
    "Nmin",
    "Norg",
)

#: Upper day-of-year of each calendar month, as used by the original's
#: monthly crop-coefficient lookup.
MONTH_LIMITS = (30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330, 366)


def soil_hydro(sand, clay, org):
    """Pedotransfer for the soil water limits [mm], as in the original."""
    wc = (0.2576 - 0.002 * sand + 0.0036 * clay + 0.0299 * org) * 1000.0
    return {
        "waterCapacity": wc,
        "waterSaturation": 100.0 / 88.0 * wc,
        "Wiltingpoint": (0.026 + 0.005 * clay + 0.0158 * org) * 1000.0,
    }


def daily_kc(doy, kc_monthly):
    """Daily crop coefficient from the 12 monthly values."""
    kc = np.asarray(kc_monthly, dtype=np.float64).ravel()
    if kc.size != 12:
        raise ValueError("daily_kc: kc_monthly must hold 12 values (January..December)")
    idx = np.searchsorted(np.asarray(MONTH_LIMITS), np.asarray(doy), side="left")
    return kc[np.clip(idx, 0, 11)]


def thermal_sum(tavg, doy, t_min, t_max, st_start):
    """Grass thermal sum [degree-days], accumulated as in the original.

    Note that day j uses day j-1's temperature, and that accumulation only
    begins after day-of-year `st_start`.
    """
    tavg = np.asarray(tavg, dtype=np.float64)
    ndays = tavg.shape[0]
    st = np.zeros_like(tavg)
    if ndays < 2:
        return st
    doy = np.asarray(doy)
    for j in range(1, ndays):
        prev = tavg[j - 1]
        after_start = doy[j] > st_start
        increment = np.where(
            after_start & (prev >= t_min) & (prev <= t_max),
            prev - t_min,
            np.where(after_start & (prev > t_max), t_max - t_min, 0.0),
        )
        st[j] = st[j - 1] + increment
    return st


def f_w(w, pet):
    """Water-limitation factor fW."""
    low = np.where(
        w < 0.2,
        4.0 * w,
        np.where(w < 0.4, 0.75 * w + 0.65, np.where(w < 0.6, 0.25 * w + 0.85, 1.0)),
    )
    mid = np.where(
        w < 0.2,
        2.0 * w,
        np.where(
            w < 0.4,
            1.5 * w + 0.1,
            np.where(w < 0.6, w + 0.3, np.where(w < 0.8, 0.5 * w + 0.6, 1.0)),
        ),
    )
    return np.where(pet < 3.81, low, np.where(pet < 6.35, mid, w))


def f_t(temp, t0, t1, t2, tlimit):
    """Temperature-limitation factor fT."""
    return np.where(
        temp <= t0,
        0.0,
        np.where(
            temp <= t1,
            safe_div(temp - t0, t1 - t0, 0.0),
            np.where(
                temp <= t2,
                1.0,
                np.where(
                    temp <= tlimit, safe_div(tlimit - temp, tlimit - t2, 0.0), 0.0
                ),
            ),
        ),
    )


def f_par(par_i):
    """Radiation-saturation factor fPARi."""
    return np.where(par_i <= 5, 1.0, -par_i / 22.0 + 27.0 / 22.0)


def f_sea(st, st_min, st1, st2, min_sea, max_sea):
    """Seasonal-effect factor SEA."""
    return np.where(
        st < st_min,
        min_sea,
        np.where(
            st <= st1 - 200,
            safe_div(max_sea - min_sea, st1 - 100 - st_min, 0.0) * (st - st_min)
            + min_sea,
            np.where(
                st <= st1 - 100,
                max_sea,
                np.where(
                    st <= st2,
                    safe_div(min_sea - max_sea, st2 - st1, 0.0) * (st - st1) + max_sea,
                    min_sea,
                ),
            ),
        ),
    )


def f_n(rnc):
    """Nitrogen-limitation factor fN."""
    return np.clip(0.99 * (1.0 - 3.78 * np.exp(-5.36 * rnc)), 0.0, 1.0)


def run(params, weather, npixels, outputs=None, progress=None):
    """Run GrasSim over `npixels` pixels."""
    p = dict(params)
    n = npixels
    doy = np.asarray(weather["doy"])
    ndays = len(doy)
    wanted = tuple(OUTPUTS) if outputs is None else tuple(outputs)

    tavg = weather["tavg"]
    rad_par = weather["rad_par"]
    rain = weather["rain"]
    et_pot = weather["et_pot"]
    gcut = weather["gcut"]

    def g(key):
        """Parameter `key` broadcast over the pixel axis."""
        return as_array(p[key], n)

    sla, pct_lam = g("SLA"), g("percentageLAM")
    st1, st2 = g("ST1"), g("ST2")
    max_sea, min_sea = g("maxSEA"), g("minSEA")
    lls = g("LLS")
    max_omd_gv, min_omd_gv = g("maxOMDGV"), g("minOMDGV")
    max_omd_gr, min_omd_gr = g("maxOMDGR"), g("minOMDGR")
    bd_gv, bd_dv, bd_gr, bd_dr = g("BDGV"), g("BDDV"), g("BDGR"), g("BDDR")
    rue_max = g("RUEmax")
    sigma_gv, sigma_gr = g("sigmaGV"), g("sigmaGR")
    t0, t1, t2, tlimit = g("T0"), g("T1"), g("T2"), g("Tlimit")
    st_min = g("STmin")
    k_gv, k_gr, kl_dv, kl_dr = g("KGV"), g("KGR"), g("KlDV"), g("KlDR")
    t_ref, k_nitro = g("Tref"), g("K")
    pct_of_nmin, nh3_volat = g("percentageofNmin"), g("NH3volatfactor")
    a_ncrit, b_ncrit = g("a_Ncrit"), g("b_Ncrit")
    a_nmax, b_nmax = g("a_Nmax"), g("b_Nmax")
    fnh_c1, fnh_c2, fna_max = g("FNH_coef1"), g("FNH_coef2"), g("FNAmax")
    nsc, rnc_max, rnc_min, fnh_max = g("NSc"), g("RNCmax"), g("RNCmin"), g("FNHmax")
    water_capacity = g("waterCapacity")
    water_saturation = g("waterSaturation")
    wilting_point = g("Wiltingpoint")
    default_cut = g("DEFAULT_CUT_HEIGHT")
    rnc_fixed = g("RNC_fixed")

    # Daily crop coefficient: the weather column where supplied, else the
    # monthly Kc parameters.
    kc = weather.get("kc")
    if kc is None:
        kc = np.full((ndays, n), np.nan)
    kc = np.asarray(kc, dtype=np.float64).copy()
    if np.any(~np.isfinite(kc)):
        if "kc_monthly" not in p:
            raise ValueError(
                "GrasSim needs a daily 'kc' forcing or the monthly Kc_1..Kc_12 "
                "parameters."
            )
        monthly = daily_kc(doy, p["kc_monthly"])
        kc = np.where(np.isfinite(kc), kc, monthly[:, None])

    # Grass thermal sum: the weather column when complete, else computed.
    st_supplied = weather.get("st_grass")
    if st_supplied is not None and np.all(np.isfinite(st_supplied)):
        st_all = np.asarray(st_supplied, dtype=np.float64)
    else:
        if st_supplied is not None and np.any(np.isfinite(st_supplied)):
            raise ValueError(
                "'st_grass' is cumulative and must be complete, or absent "
                "altogether so that it can be computed."
            )
        st_all = thermal_sum(
            tavg, doy, g("ST_Tmin"), g("ST_Tmax"), as_array(p["ST_start"], n)
        )

    # State variables.
    initial_height = g("Initial_height")
    gv_biomass = initial_height * 10.0 * bd_gv
    gr_biomass = g("initial_reproductive_biomass")
    dv_biomass = initial_height * 10.0 * bd_dv
    dr_biomass = g("initial_reproductive_biomass")
    gv_avg_age, gr_avg_age = g("AGEGV"), g("AGEGR")
    dv_avg_age, dr_avg_age = g("AGEDV"), g("AGEDR")
    water = water_capacity.copy()
    n_min, n_org = g("Nmin"), g("Norg")
    n_conc = a_ncrit * 0.01
    qn_gv = n_conc * gv_biomass
    qn_gr = n_conc * gr_biomass
    qn_dv = 0.008 * dv_biomass
    qn_dr = 0.008 * dr_biomass
    not_runoff = g("notRunoff")
    apex_grazed_ever = np.zeros(n)
    sward_height = initial_height.copy()
    harvested = np.zeros(n)

    fert_min = weather.get("fert_min")
    fert_org = weather.get("fert_org")

    out = {name: np.zeros((ndays, n), dtype=np.float32) for name in wanted}

    def store(name, value):
        if name in out:
            out[name][i] = value

    for i in range(ndays):
        pp = rain[i]
        temp = tavg[i]
        pet = et_pot[i]
        par_i = rad_par[i]
        st = st_all[i]
        kc_day = kc[i]
        cut_today = gcut[i] > 0

        lai = sla * gv_biomass / 10.0 * pct_lam

        ngv = safe_div(qn_gv, gv_biomass, 0.0)
        ngr = safe_div(qn_gr, gr_biomass, 0.0)
        ndv = safe_div(qn_dv, dv_biomass, 0.0)
        ndr = safe_div(qn_dr, dr_biomass, 0.0)

        # --- Water balance ---
        aet = pet * kc_day
        water_t = water + (pp - aet) + not_runoff
        leached = np.where(
            (water_t >= water_capacity) & (water_t < water_saturation),
            0.2 * (water_t - water_capacity),
            np.where(water_t > water_saturation, 0.2 * (water_t - water_capacity), 0.0),
        )
        not_runoff = np.where(
            water_t > water_saturation, 0.2 * (water_t - water_saturation), 0.0
        )
        water_t = water_t - leached
        water_t = np.minimum(water_t, water_saturation)
        water_t = np.maximum(water_t, wilting_point)
        w_rel = np.clip(
            safe_div(water_t - wilting_point, water_capacity - wilting_point, 0.0),
            0.0,
            1.0,
        )

        fw_day = f_w(w_rel, pet)
        ft_day = f_t(temp, t0, t1, t2, tlimit)
        fpar = f_par(par_i)
        sea = f_sea(st, st_min, st1, st2, min_sea, max_sea)

        # --- Age factors, abscission, senescence ---
        gv_ratio = safe_div(gv_avg_age, lls, 0.0)
        gr_ratio = safe_div(gr_avg_age, st2 - st1, 0.0)
        dv_ratio = safe_div(dv_avg_age, lls, 0.0)
        dr_ratio = safe_div(dr_avg_age, st2 - st1, 0.0)
        f_age_gv = np.where(
            gv_ratio < 1 / 3, 1.0, np.where(gv_ratio < 1, 3 * gv_ratio, 3.0)
        )
        f_age_gr = np.where(
            gr_ratio < 1 / 3, 1.0, np.where(gr_ratio < 1, 3 * gr_ratio, 3.0)
        )
        f_age_dv = np.where(dv_ratio < 1 / 3, 1.0, np.where(dv_ratio < 2 / 3, 2.0, 3.0))
        f_age_dr = np.where(dr_ratio < 1 / 3, 1.0, np.where(dr_ratio < 2 / 3, 2.0, 3.0))

        abs_dv = np.where(temp > 0, kl_dv * dv_biomass * temp * f_age_dv, 0.0)
        abs_dr = np.where(temp > 0, kl_dr * dr_biomass * temp * f_age_dr, 0.0)
        sen_gv = np.where(
            temp > t0,
            k_gv * gv_biomass * temp * f_age_gv,
            np.where(temp < 0, k_gv * gv_biomass * np.abs(temp), 0.0),
        )
        sen_gr = np.where(
            temp > t0,
            k_gr * gr_biomass * temp * f_age_gr,
            np.where(temp < 0, k_gr * gr_biomass * np.abs(temp), 0.0),
        )

        # --- Soil nitrogen ---
        g0 = 0.8 * w_rel + 0.2
        f_t_nitro = np.exp(k_nitro * (temp - t_ref))
        vp = (
            (0.0929 + (0.1833 - 0.0929) * np.exp(-0.2173 * n_org / 1000.0))
            * n_org
            / 1000.0
        )
        mineralisation = vp * f_t_nitro * g0
        ip = np.maximum(4.0 / 1000.0 * n_min, 0.0)
        immobilisation = ip * f_t_nitro * g0
        fna = np.minimum(fna_max, fna_max * safe_div(n_min, nsc, 0.0))
        n_supply = np.where(n_min < 0, 0.0, fna * n_min)

        # --- Plant nitrogen status ---
        bm = gv_biomass + gr_biomass + dv_biomass + dr_biomass
        # The 0.05 m reference height is hardcoded in the original,
        # independently of the actual residual cut height.
        cut_off_gv = 0.05 * 10.0 * bd_gv
        cut_off_gr = 0.05 * 10.0 * bd_gr
        cut_off_dv = 0.05 * 10.0 * bd_dv
        cut_off_dr = 0.05 * 10.0 * bd_dr
        gv_over5 = np.maximum(gv_biomass - cut_off_gv, 0.0)
        gr_over5 = np.maximum(gr_biomass - cut_off_gr, 0.0)
        dv_over5 = np.maximum(dv_biomass - cut_off_dv, 0.0)
        dr_over5 = np.maximum(dr_biomass - cut_off_dr, 0.0)
        bm_over5 = np.maximum(
            bm - (cut_off_gv + cut_off_gr + cut_off_dv + cut_off_dr), 0.0
        )

        ncrit = np.where(
            bm_over5 <= 1000,
            a_ncrit * 0.01,
            a_ncrit * 0.01 * np.power(np.maximum(bm_over5, 1e-9) / 1000.0, -b_ncrit),
        )
        nmax = np.where(
            bm_over5 <= 1000,
            a_nmax * 0.01,
            a_nmax * 0.01 * np.power(np.maximum(bm_over5, 1e-9) / 1000.0, -b_nmax),
        )
        # BMover5 == 0 right after a cut at 0.05 m forces RNC to RNCmin,
        # whatever the nitrogen supply: original behaviour, kept.
        n_act = np.where(
            bm_over5 <= 0,
            rnc_min * ncrit,
            safe_div(
                ngv * gv_over5 + ngr * gr_over5 + ndv * dv_over5 + ndr * dr_over5,
                bm_over5,
                0.0,
            ),
        )
        n_act = np.minimum(ncrit * rnc_max, n_act)
        rnc_day = np.clip(safe_div(n_act, ncrit, 0.0), 0.25, 1.0)
        rnc_day = np.where(np.isfinite(rnc_fixed), rnc_fixed, rnc_day)
        rnc_day = np.clip(rnc_day, 0.25, 1.0)
        fn_day = f_n(rnc_day)

        # --- Growth ---
        env = fpar * ft_day * np.minimum(fw_day, fn_day)
        pgro = par_i * rue_max * (1.0 - np.exp(-0.6 * lai)) * 10.0
        gro = pgro * env * sea

        apex_grazed = ((st > st1) & (st < st2) & cut_today).astype(np.float64)
        apex_grazed_ever = np.maximum(apex_grazed_ever, apex_grazed)

        rep = np.where(
            st < st1,
            0.0,
            np.where(
                (st < st2)
                & (apex_grazed == 0)
                & (~cut_today)
                & (rnc_day > 0.35)
                & (apex_grazed_ever == 0),
                0.25 + ((1.0 - 0.25) * (rnc_day - 0.35)) / (1.0 - 0.35),
                0.0,
            ),
        )
        gro_gv = gro * (1.0 - rep)
        gro_gr = gro * rep

        gv_new = np.maximum(gv_biomass + gro_gv - sen_gv, 0.0)
        gr_new = np.maximum(gr_biomass + gro_gr - sen_gr, 0.0)
        dv_new = np.maximum(dv_biomass + (1.0 - sigma_gv) * sen_gv - abs_dv, 0.0)
        dr_new = np.maximum(dr_biomass + (1.0 - sigma_gr) * sen_gr - abs_dr, 0.0)
        gv_biomass, gr_biomass = gv_new, gr_new
        dv_biomass, dr_biomass = dv_new, dr_new

        sward_height = np.maximum.reduce(
            [
                safe_div(gv_biomass / 10.0, bd_gv, 0.0),
                safe_div(gr_biomass / 10.0, bd_gr, 0.0),
                safe_div(dv_biomass / 10.0, bd_dv, 0.0),
                safe_div(dr_biomass / 10.0, bd_dr, 0.0),
            ]
        )

        # --- Ages (weights use the already-updated biomass, as in the original) ---
        warm = temp > 0
        den_gv = gv_biomass - sen_gv + gro_gv
        den_gr = gr_biomass - sen_gr + gro_gr
        den_dv = dv_biomass - abs_dv + (1.0 - sigma_gv) * sen_gv
        den_dr = dr_biomass - abs_dr + (1.0 - sigma_gr) * sen_gr
        gv_avg_age = np.where(
            warm & (den_gv > 0),
            safe_div(gv_biomass - sen_gv, den_gv, 0.0) * (gv_avg_age + temp),
            gv_avg_age,
        )
        gr_avg_age = np.where(
            warm & (den_gr > 0),
            safe_div(gr_biomass - sen_gr, den_gr, 0.0) * (gr_avg_age + temp),
            gr_avg_age,
        )
        dv_avg_age = np.where(
            warm & (den_dv > 0),
            safe_div(dv_biomass - abs_dv, den_dv, 0.0) * (dv_avg_age + temp),
            dv_avg_age,
        )
        dr_avg_age = np.where(
            warm & (den_dr > 0),
            safe_div(dr_biomass - abs_dr, den_dr, 0.0) * (dr_avg_age + temp),
            dr_avg_age,
        )

        # --- Digestibility ---
        omd_gv = max_omd_gv - safe_div(gv_avg_age * (max_omd_gv - min_omd_gv), lls, 0.0)
        omd_gr = np.where(
            st < st1,
            max_omd_gr,
            np.where(
                st > st2,
                min_omd_gr,
                max_omd_gr
                - safe_div(gr_avg_age * (max_omd_gr - min_omd_gr), st2 - st1, 0.0),
            ),
        )

        # --- Plant nitrogen demand and uptake ---
        fnh = fnh_max * np.maximum(fnh_c1, fnh_c2 * (1.0 + rnc_day))
        n_demand = safe_div(bm_over5 * (nmax - n_act), fnh, 0.0)
        n_uptake = np.minimum(n_demand, n_supply)
        n_plant_litter = (abs_dv + abs_dr) * 0.008
        qn_dv = qn_dv + ((1.0 - sigma_gv) * sen_gv - abs_dv) * 0.008
        qn_dr = qn_dr + ((1.0 - sigma_gr) * sen_gr - abs_dr) * 0.008
        growing = gro > 0
        qn_gv = np.where(
            growing,
            qn_gv + (safe_div(n_uptake * fnh * gro_gv, gro, 0.0) - sen_gv * 0.008),
            qn_gv - sen_gv * 0.008,
        )
        qn_gr = np.where(
            growing,
            qn_gr + (safe_div(n_uptake * fnh * gro_gr, gro, 0.0) - sen_gr * 0.008),
            qn_gr - sen_gr * 0.008,
        )
        qn_gv = np.maximum(qn_gv, 0.0)
        qn_gr = np.maximum(qn_gr, 0.0)
        qn_dv = np.maximum(qn_dv, 0.0)
        qn_dr = np.maximum(qn_dr, 0.0)

        # --- Soil nitrogen pools ---
        n_leached = np.where(water_t > 0, safe_div(n_min, water_t, 0.0) * leached, 0.0)
        n_from_rain = 0.009 * pp
        min_fert = np.zeros(n) if fert_min is None else fert_min[i]
        org_fert = np.zeros(n) if fert_org is None else fert_org[i]
        n_org = (
            n_org
            + immobilisation
            - mineralisation
            + (1.0 - pct_of_nmin) * org_fert
            + n_plant_litter
        )
        n_min = (
            n_min
            + n_from_rain
            + mineralisation
            + min_fert
            + pct_of_nmin * (1.0 - nh3_volat) * org_fert
            - immobilisation
            - n_uptake
            - n_leached
        )

        # --- Cut (plant N pools are not reduced by the harvest, as in the original) ---
        if np.any(cut_today):
            height = gcut[i]
            res_gv = np.minimum(gv_biomass, height * 10.0 * bd_gv)
            res_gr = np.minimum(gr_biomass, height * 10.0 * bd_gr)
            res_dv = np.minimum(dv_biomass, height * 10.0 * bd_dv)
            res_dr = np.minimum(dr_biomass, height * 10.0 * bd_dr)
            taken = (
                gv_biomass
                + gr_biomass
                + dv_biomass
                + dr_biomass
                - res_gv
                - res_gr
                - res_dv
                - res_dr
            )
            harvested = harvested + np.where(cut_today, taken, 0.0)
            gv_biomass = np.where(cut_today, res_gv, gv_biomass)
            gr_biomass = np.where(cut_today, res_gr, gr_biomass)
            dv_biomass = np.where(cut_today, res_dv, dv_biomass)
            dr_biomass = np.where(cut_today, res_dr, dr_biomass)
            sward_height = np.where(cut_today, height, sward_height)

        tsb = gv_biomass + gr_biomass + dv_biomass + dr_biomass
        water = water_t

        height_abc = np.where(gcut[i] > 0, gcut[i], default_cut)
        abc = (
            np.maximum(gv_biomass - height_abc * 10.0 * bd_gv, 0.0)
            + np.maximum(gr_biomass - height_abc * 10.0 * bd_gr, 0.0)
            + np.maximum(dv_biomass - height_abc * 10.0 * bd_dv, 0.0)
            + np.maximum(dr_biomass - height_abc * 10.0 * bd_dr, 0.0)
        )

        store("gvb", gv_biomass)
        store("grb", gr_biomass)
        store("dvb", dv_biomass)
        store("drb", dr_biomass)
        store("g", gro)
        store("pgr", pgro)
        store("env", env)
        store("sea", sea)
        store("ftm", ft_day)
        store("lai", lai)
        store("hb", harvested)
        store("stp", st)
        store("gva", gv_avg_age)
        store("gra", gr_avg_age)
        store("dva", dv_avg_age)
        store("dra", dr_avg_age)
        store("fw", fw_day)
        store("fn", fn_day)
        store("rnc", rnc_day)
        store("nmin", n_min)
        store("norg", n_org)
        store("sward_height", sward_height)
        store("tsb", tsb)
        store("abc", abc)
        store("omdgv", omd_gv)
        store("omdgr", omd_gr)

        if progress is not None:
            progress(i, ndays)

    return out
