"""
MODULE:     t.graminr (graminr_lib.lingra)

AUTHOR(S):  Yann Chemin

PURPOSE:    LINGRA grassland model with the PCSE water balance, ported from
            the GraMinR R package (R/lingra.R, R/lingra_lib.R) and vectorised
            over pixels. Latitude and elevation are taken per pixel, so the
            solar geometry and the reference evapotranspiration vary across
            the region instead of being fixed at a single site.

SPDX-License-Identifier: Unlicense
"""

import numpy as np

from .core import afgen, as_array, clamp, safe_div

OUTPUTS = (
    "lai",
    "dvs",
    "tsum",
    "tiller_number",
    "weight_lv_green",
    "weight_lv_dead",
    "weight_harv",
    "weight_harv_annual",
    "weight_abg",
    "weight_re",
    "weight_rt",
    "leaf_length",
    "slaint",
    "rd",
    "sm",
    "w",
    "wlow",
    "ss",
    "rftra",
    "tra",
    "tramx",
    "et0",
    "e0",
    "es0",
    "evs",
    "evw",
    "perc",
    "loss",
    "runoff",
    "mowing",
)

REQUIRED = ("rad_glob", "tmin", "tmax", "vp", "wind", "rain")
OPTIONAL = {
    "et_act_rs": np.nan,
    "t_act_rs": np.nan,
    "lai_rs": np.nan,
    "gcut": 0.0,
    "rirr": 0.0,
}

ALIASES = {}

DEFAULTS = {
    "ETMODEL": "PM",
    # PCSE's LINGRA starts the soil temperature at a fixed 5 degree C; a
    # parameter set may say otherwise.
    "TemperatureSoilinit": 5.0,
    # Only the biomass-threshold mowing rule (imopt = 1) uses cwght.
    "cwght": float("inf"),
    # What suspends photosynthesis on a cutting day. "harvest" (the default)
    # suspends it only when biomass is actually removed, reproducing PCSE's
    # LINGRA to machine precision (see tests/test_lingra_pcse.py). "mowflag"
    # suspends it on every scheduled cutting day, which is what GraMinR's R
    # does: on a day whose cut removes nothing, because green leaf weight is
    # already at or below cwlvg, a whole day of assimilation is still lost.
    # Keep "mowflag" only to reproduce GraMinR's R output exactly.
    "HarvestGate": "harvest",
    "ANGSTA": 0.18,
    "ANGSTB": 0.55,
    "ELEV": 0.0,
    "IFUNRN": 0.0,
    "NOTINF": 0.0,
    "NINFTB": (0.0, 0.0, 0.5, 0.0, 1.5, 1.0),
}

#: Tabulated parameters, given as flat x1,y1,x2,y2,... sequences.
TABLE_PARAMS = (
    "KDIFTB",
    "LUEreductionSoilTempTB",
    "LUEreductionRadiationTB",
    "PartitioningRootsTB",
    "NINFTB",
)

#: Parameters run() reads and cannot default; checked after the derivations
#: so that a set given in either naming scheme is validated the same way.
REQUIRED_PARAMS = (
    "LAT",
    "KDIFTB",
    "LUEreductionSoilTempTB",
    "LUEreductionRadiationTB",
    "PartitioningRootsTB",
    "SLA",
    "LUEmax",
    "CO2A",
    "LAIcrit",
    "SiteFillingMax",
    "TillerFormRateA0",
    "TillerFormRateB0",
    "TillerFormRateA8",
    "TillerFormRateB8",
    "TSUMmax",
    "TempBase",
    "RDRbase",
    "RDRshading",
    "RDRdrought",
    "CFET",
    "DEPNR",
    "IAIRDU",
    "IOX",
    "CRAIRC",
    "SM0",
    "SMW",
    "SMFCF",
    "SMLIM",
    "SOPE",
    "KSUB",
    "RDMSOL",
    "SSMAX",
    "SSI",
    "WAV",
    "imopt",
    "cwlvg",
    "RDI",
    "RRI",
    "RDMCR",
    "LAIinit",
    "TillerNumberinit",
    "WeightREinit",
    "WeightRTinit",
)

#: Derived parameters computed from the raw LINGRA names when absent.
DERIVED = {
    "wai": lambda p: 1000.0 * p["rootdi"] * p["wci"],
    "wlvgi": lambda p: p["LAIi"] / p["sla"],
    "cwlvg": lambda p: p["cLAI"] / p["sla"],
    "fsmax": lambda p: p["nitr"] / p["nitmax"] * 0.693,
    # rootdi/rootdm/rrdmax are in metres; the PCSE water balance and root
    # dynamics work in centimetres throughout, hence the x100.
    "RDI": lambda p: 100.0 * p["rootdi"],
    "RRI": lambda p: 100.0 * p["rrdmax"],
    "RDMCR": lambda p: 100.0 * p["rootdm"],
    "TempBase": lambda p: p["tmbas1"],
    "TillerNumberinit": lambda p: p["tilli"],
    "WeightREinit": lambda p: p["wrei"],
    "WeightRTinit": lambda p: p["wrti"],
    "LAIinit": lambda p: p["LAIi"],
    "SLA": lambda p: p["sla"],
    "CO2A": lambda p: p["co2a"],
    "LAIcrit": lambda p: p["LAIcr"],
    "SiteFillingMax": lambda p: p["fsmax"],
}


def check_angstrom(angsta, angstb):
    """Reject implausible Angstrom coefficients, as the R port does."""
    a, b = abs(float(angsta)), abs(float(angstb))
    if not 0.1 <= a <= 0.4:
        raise ValueError(
            "Implausible Angstrom coefficient A ({:.3f}); expected within "
            "[0.1, 0.4]. Calibrate locally (see FAO-56 or WOFOST/CABO).".format(a)
        )
    if not 0.3 <= b <= 0.7:
        raise ValueError(
            "Implausible Angstrom coefficient B ({:.3f}); expected within "
            "[0.3, 0.7].".format(b)
        )
    if not 0.6 <= a + b <= 0.9:
        raise ValueError(
            "Implausible Angstrom sum A+B ({:.3f}); expected within [0.6, 0.9].".format(
                a + b
            )
        )


def sat_vapour_pressure(temp):
    """Saturated vapour pressure [kPa] from temperature [degree C]."""
    return 0.6108 * np.exp((17.27 * temp) / (237.3 + temp))


def sweaf(et0, depnr):
    """Soil water easily available fraction (Doorenbos & Kassam 1979)."""
    value = 1.0 / (0.76 + 1.5 * et0) - (5.0 - depnr) * 0.10
    if float(np.min(depnr)) < 3:
        value = value + (et0 - 0.6) / (depnr * (depnr + 3.0))
    return np.clip(value, 0.10, 0.95)


def astro(doy, latitude, radiation):
    """Solar geometry after van Kraalingen (1991), vectorised over pixels.

    :param doy: scalar day of year
    :param latitude: per-pixel latitude [degree]
    :param radiation: per-pixel global radiation [J m-2 d-1]
    """
    rad = np.pi / 180.0
    angle = -4.0
    dec = -np.arcsin(np.sin(23.45 * rad) * np.cos(2.0 * np.pi * (doy + 10.0) / 365.0))
    sc = 1370.0 * (1.0 + 0.033 * np.cos(2.0 * np.pi * doy / 365.0))

    sinld = np.sin(rad * latitude) * np.sin(dec)
    cosld = np.cos(rad * latitude) * np.cos(dec)
    aob = safe_div(sinld, cosld, 0.0)

    inside = np.abs(aob) <= 1
    aob_c = np.clip(aob, -1.0, 1.0)
    root = np.sqrt(np.maximum(0.0, 1.0 - aob_c**2))

    dayl_in = 12.0 * (1.0 + 2.0 * np.arcsin(aob_c) / np.pi)
    dayl_out = np.where(aob > 1, 24.0, 0.0)
    dayl = np.where(inside, dayl_in, dayl_out)

    dsinb = np.where(
        inside,
        3600.0 * (dayl * sinld + 24.0 * cosld * root / np.pi),
        3600.0 * (dayl * sinld),
    )

    aob_corr = safe_div(-np.sin(angle * rad) + sinld, cosld, 0.0)
    aob_corr_c = np.clip(aob_corr, -1.0, 1.0)
    daylp = np.where(
        np.abs(aob_corr) <= 1,
        12.0 * (1.0 + 2.0 * np.arcsin(aob_corr_c) / np.pi),
        np.where(aob_corr > 1, 24.0, 0.0),
    )

    angot = sc * dsinb
    atmtr = np.where(dayl > 0, safe_div(radiation, angot, 0.0), 0.0)

    frdif = np.where(
        atmtr > 0.75,
        0.23,
        np.where(
            atmtr > 0.35,
            1.33 - 1.46 * atmtr,
            np.where(atmtr > 0.07, 1.0 - 2.3 * (atmtr - 0.07) ** 2, 1.0),
        ),
    )
    return {
        "dayl": dayl,
        "daylp": daylp,
        "atmtr": atmtr,
        "angot": angot,
        "difpp": frdif * atmtr * 0.5 * sc,
    }


def penman(doy, lat, elev, tmin, tmax, avrad, vap, wind2, angsta, angstb):
    """Modified Penman open-water, bare-soil and crop evaporation [mm/d]."""
    psycon, lhvap = 0.67, 2.45e6
    stbc = 5.670373e-8 * 24 * 60 * 60
    refcfw, refcfs, refcfc = 0.05, 0.15, 0.25

    tmpa = (tmin + tmax) / 2.0
    bu = 0.54 + 0.35 * np.clip((tmax - tmin - 12.0) / 4.0, 0.0, 1.0)
    pbar = 1013.0 * np.exp(-0.034 * elev / (tmpa + 273.0))
    gamma = psycon * pbar / 1013.0
    svap = 6.10588 * np.exp(17.32491 * tmpa / (tmpa + 238.102))
    delta = 238.102 * 17.32491 * svap / (tmpa + 238.102) ** 2
    vap = np.minimum(vap, svap)

    a = astro(doy, lat, avrad)
    relssd = np.clip(safe_div(a["atmtr"] - abs(angsta), abs(angstb), 0.0), 0.0, 1.0)
    rb = (
        stbc
        * (tmpa + 273.0) ** 4
        * (0.56 - 0.079 * np.sqrt(np.maximum(vap, 0.0)))
        * (0.1 + 0.9 * relssd)
    )
    rnw = (avrad * (1.0 - refcfw) - rb) / lhvap
    rns = (avrad * (1.0 - refcfs) - rb) / lhvap
    rnc = (avrad * (1.0 - refcfc) - rb) / lhvap
    deficit = np.maximum(0.0, svap - vap)
    ea = 0.26 * deficit * (0.5 + bu * wind2)
    eac = 0.26 * deficit * (1.0 + bu * wind2)
    return {
        "e0": np.maximum(0.0, (delta * rnw + gamma * ea) / (delta + gamma)),
        "es0": np.maximum(0.0, (delta * rns + gamma * ea) / (delta + gamma)),
        "et0": np.maximum(0.0, (delta * rnc + gamma * eac) / (delta + gamma)),
        "angot": a["angot"],
    }


def penman_monteith(doy, lat, elev, tmin, tmax, avrad, vap_hpa, wind2):
    """FAO-56 reference crop evapotranspiration ET0 [mm/d]."""
    psycon, lhvap = 0.665, 2.45e6
    refcfc, cres, stbc = 0.23, 70.0, 4.903e-3

    tmpa = (tmin + tmax) / 2.0
    vap = vap_hpa / 10.0
    t_std = 293.0
    patm = 101.3 * ((t_std - 0.0065 * elev) / t_std) ** 5.26
    gamma = psycon * patm * 1.0e-3
    svap_tmpa = sat_vapour_pressure(tmpa)
    delta = (4098.0 * svap_tmpa) / (tmpa + 237.3) ** 2
    svap = (sat_vapour_pressure(tmax) + sat_vapour_pressure(tmin)) / 2.0
    vap = np.minimum(vap, svap)

    stb_tmax = stbc * (tmax + 273.16) ** 4
    stb_tmin = stbc * (tmin + 273.16) ** 4
    rnl_tmp = ((stb_tmax + stb_tmin) / 2.0) * (
        0.34 - 0.14 * np.sqrt(np.maximum(vap, 0.0))
    )

    a = astro(doy, lat, avrad)
    csky = (0.75 + 2e-05 * elev) * a["angot"]
    rnl = rnl_tmp * (1.35 * safe_div(avrad, csky, 0.0) - 0.35)
    rn = ((1.0 - refcfc) * avrad - rnl) / lhvap
    ea = (900.0 / (tmpa + 273.0)) * wind2 * (svap - vap)
    mgamma = gamma * (1.0 + (cres / 208.0 * wind2))
    et0 = (delta * rn) / (delta + mgamma) + (gamma * ea) / (delta + mgamma)
    return np.where(csky > 0, np.maximum(0.0, et0), 0.0)


def reference_et(
    doy, lat, elev, tmin, tmax, avrad, vap_hpa, wind2, angsta, angstb, etmodel="PM"
):
    """Reference evapotranspiration, in cm/day as LINGRA expects."""
    if etmodel not in ("PM", "P"):
        raise ValueError("etmodel must be 'PM' or 'P'")
    out = penman(doy, lat, elev, tmin, tmax, avrad, vap_hpa, wind2, angsta, angstb)
    et0 = out["et0"]
    if etmodel == "PM":
        et0 = penman_monteith(doy, lat, elev, tmin, tmax, avrad, vap_hpa, wind2)
    return out["e0"] / 10.0, out["es0"] / 10.0, et0 / 10.0


def run(params, weather, npixels, outputs=None, progress=None):
    """Run LINGRA over `npixels` pixels."""
    p = dict(params)
    n = npixels
    year = np.asarray(weather["year"])
    doy = np.asarray(weather["doy"])
    ndays = len(doy)
    wanted = tuple(OUTPUTS) if outputs is None else tuple(outputs)

    check_angstrom(p["ANGSTA"], p["ANGSTB"])

    def g(key):
        """Parameter `key` broadcast over the pixel axis."""
        return as_array(p[key], n)

    lat = g("LAT")
    elev = g("ELEV")
    angsta, angstb = float(p["ANGSTA"]), float(p["ANGSTB"])
    etmodel = str(p.get("ETMODEL", "PM"))

    kdiftb = p["KDIFTB"]
    lue_soiltemp_tb = p["LUEreductionSoilTempTB"]
    lue_rad_tb = p["LUEreductionRadiationTB"]
    part_roots_tb = p["PartitioningRootsTB"]
    ninftb = p.get("NINFTB", DEFAULTS["NINFTB"])

    sla = g("SLA")
    luemax, co2a = g("LUEmax"), g("CO2A")
    lai_crit = g("LAIcrit")
    site_filling_max = g("SiteFillingMax")
    tfa0, tfb0 = g("TillerFormRateA0"), g("TillerFormRateB0")
    tfa8, tfb8 = g("TillerFormRateA8"), g("TillerFormRateB8")
    tsum_max = g("TSUMmax")
    temp_base = g("TempBase")
    rdr_base, rdr_shading, rdr_drought = g("RDRbase"), g("RDRshading"), g("RDRdrought")
    cfet, depnr = g("CFET"), g("DEPNR")
    iairdu, iox = int(p["IAIRDU"]), int(p["IOX"])
    crairc = g("CRAIRC")
    sm0, smw, smfcf = g("SM0"), g("SMW"), g("SMFCF")
    sope, ksub = g("SOPE"), g("KSUB")
    rdmsol, ssmax = g("RDMSOL"), g("SSMAX")
    harvest_gate = str(p.get("HarvestGate", "harvest")).lower()
    if harvest_gate not in ("harvest", "mowflag"):
        raise ValueError(
            "LINGRA parameter 'HarvestGate' must be 'harvest' or 'mowflag', "
            "not '{}'.".format(harvest_gate)
        )
    imopt = int(p["imopt"])
    if imopt == 1 and not np.all(np.isfinite(as_array(p["cwght"], n))):
        raise ValueError(
            "LINGRA with imopt = 1 mows on a green-leaf-weight threshold, so "
            "the parameter 'cwght' must be given."
        )
    cwght, cwlvg = g("cwght"), g("cwlvg")
    rri = g("RRI")
    notinf = g("NOTINF")
    ifunrn = int(float(p.get("IFUNRN", 0))) == 1

    rad_glob = weather["rad_glob"]
    tmin, tmax = weather["tmin"], weather["tmax"]
    rain = weather["rain"]
    wind = weather["wind"]
    vp = weather["vp"]
    gcut = weather["gcut"]
    rirr = weather.get("rirr")
    et_act_rs = weather.get("et_act_rs")
    t_act_rs = weather.get("t_act_rs")

    # Unit conversions of model_init() for the "lingra" model.
    davtmp_all = (tmin + tmax) / 2.0
    irrad_all = rad_glob * 1.0e6  # MJ m-2 d-1 -> J m-2 d-1
    vap_hpa_all = vp * 10.0  # kPa -> hPa
    rain_cm_all = rain / 10.0  # mm -> cm

    # --- Initial state (PCSE LINGRA.initialize) ---
    rdi = g("RDI")
    tsum = np.zeros(n)
    dvs = np.zeros(n)
    lai = g("LAIinit").copy()
    days_after_harvest = np.zeros(n)
    cutting_number = np.zeros(n)
    tiller_number = g("TillerNumberinit").copy()
    weight_lv_green = g("LAIinit") / sla
    weight_lv_dead = np.zeros(n)
    weight_harv = np.zeros(n)
    weight_re = g("WeightREinit").copy()
    weight_rt = g("WeightRTinit").copy()
    leaf_length = np.zeros(n)
    slaint = sla.copy()
    rd = rdi.copy()
    rdm = np.maximum(rdi, np.minimum(g("RDMCR"), rdmsol))
    temperature_soil = g("TemperatureSoilinit").copy()

    # Water balance: PCSE always uses a fixed 10 cm rooting depth for its
    # initial state, whatever the crop's own RDI; the redistribution towards
    # the crop's actual RD happens on day 1. It also tracks its own maximum
    # depth, generally different from the crop's.
    rd_wb_init = 10.0
    rdm_wb = np.maximum(rd_wb_init, rdmsol)
    smlim = clamp(smw, sm0, g("SMLIM"))
    sm = clamp(smw, smlim, smw + g("WAV") / rd_wb_init)
    w = sm * rd_wb_init
    wlow = clamp(np.zeros(n), sm0 * (rdm_wb - rd_wb_init), g("WAV") + rdm_wb * smw - w)
    ss = g("SSI").copy()
    dsos = np.zeros(n)
    dslr = np.where(sm >= (smw + 0.5 * (smfcf - smw)), 1.0, 5.0)
    rinold = np.zeros(n)
    wb_rd_old = np.full(n, rd_wb_init)

    out = {name: np.zeros((ndays, n), dtype=np.float32) for name in wanted}
    harv_year_base = np.zeros(n)
    previous_year = None

    def store(name, value):
        if name in out:
            out[name][i] = value

    for i in range(ndays):
        if previous_year is not None and year[i] != previous_year:
            harv_year_base = weight_harv.copy()
        previous_year = year[i]

        davtmp = davtmp_all[i]
        irrad = irrad_all[i]

        # --- Reference evapotranspiration ---
        e0, es0, et0 = reference_et(
            float(doy[i]),
            lat,
            elev,
            tmin[i],
            tmax[i],
            irrad,
            vap_hpa_all[i],
            wind[i],
            angsta,
            angstb,
            etmodel,
        )

        # --- Crop evapotranspiration ---
        kglob = 0.75 * afgen(kdiftb, dvs)
        et0_crop = np.maximum(0.0, cfet * et0)
        ekl = np.exp(-kglob * lai)
        evwmx = e0 * ekl
        evsmx = np.maximum(0.0, es0 * ekl)
        tramx = et0_crop * (1.0 - ekl)
        swdep = sweaf(et0_crop, depnr)
        smcr = (1.0 - swdep) * (smfcf - smw) + smw
        rfws = np.clip(safe_div(sm - smw, smcr - smw, 0.0), 0.0, 1.0)
        if iairdu == 0 and iox == 1:
            rfosmx = np.clip(safe_div(sm0 - sm, crairc, 0.0), 0.0, 1.0)
            rfos = rfosmx + (1.0 - np.minimum(dsos, 4.0) / 4.0) * (1.0 - rfosmx)
        else:
            rfos = np.ones(n)
        rftra = rfos * rfws
        tra = tramx * rftra
        if t_act_rs is not None:
            forced = t_act_rs[i]
            tra = np.where(np.isfinite(forced), forced, tra)

        d_temperature_soil = (davtmp - temperature_soil) / 10.0
        rr = np.minimum(rdm - rd, rri)

        # --- Mowing ---
        if imopt == 1:
            mowing = weight_lv_green > cwght
        else:
            mowing = gcut[i] == 1
        d_weight_harv = np.where(mowing, np.maximum(0.0, weight_lv_green - cwlvg), 0.0)
        d_days_after_harvest = np.where(
            mowing, -days_after_harvest, np.where(cutting_number > 0, 1.0, 0.0)
        )
        d_cutting_number = np.where(mowing, 1.0, 0.0)

        # --- Death rates ---
        rdrdrought = clamp(np.zeros(n), rdr_drought, rdr_drought * (1.0 - rftra))
        rdrshading = clamp(
            np.zeros(n),
            rdr_shading,
            rdr_shading * safe_div(lai - lai_crit, lai_crit, 0.0),
        )
        rdrtotal = rdr_base + np.maximum(rdrshading, rdrdrought)
        lai_death = np.where(
            mowing, d_weight_harv * slaint, lai * (1.0 - np.exp(-rdrtotal))
        )

        # --- Partitioning ---
        rtfraction = afgen(part_roots_tb, rftra)
        lvfraction = 1.0 - rtfraction

        # --- Source-limited growth ---
        dtr = irrad / 1.0e6
        par = dtr * 0.5
        rf_temperature = afgen(lue_soiltemp_tb, temperature_soil)
        rf_radiation = afgen(lue_rad_tb, dtr)
        fint = 1.0 - np.exp(-afgen(kdiftb, dvs) * lai)
        parint = fint * par
        lue_act = luemax * rf_temperature * rf_radiation * rftra
        suspended = mowing if harvest_gate == "mowflag" else (d_weight_harv > 0)
        growth_source = np.where(
            suspended,
            0.0,
            lue_act * parint * (1.0 + 0.8 * np.log(co2a / 360.0)) * 10.0,
        )

        # --- Sink-limited growth ---
        leaf_app_rate = np.where(rf_temperature > 0, temperature_soil * 0.01, 0.0)
        tiller_form_rate = np.where(
            days_after_harvest < 8,
            np.maximum(0.0, tfa0 - tfb0 * lai) * rf_temperature,
            clamp(np.zeros(n), site_filling_max, tfa8 - tfb8 * lai) * rf_temperature,
        )
        tiller_death_rate = np.maximum(
            0.01 * (1.0 + safe_div(tsum, tsum_max, 0.0)),
            0.05 * safe_div(lai - lai_crit, lai_crit, 0.0),
        )
        d_tiller_number = np.where(
            tiller_number <= 14000,
            (tiller_form_rate - tiller_death_rate) * leaf_app_rate * tiller_number,
            -tiller_death_rate * leaf_app_rate * tiller_number,
        )
        d_leaf_length_pot = np.where(
            (davtmp - temp_base) > 0,
            0.83 * np.log(np.maximum(davtmp, 2.0)) - 0.8924,
            0.0,
        )
        lai_growth_sink = (tiller_number * 1e4 * (d_leaf_length_pot * 0.3)) * 1e-8
        growth_sink = np.where(
            d_weight_harv <= 0,
            lai_growth_sink * safe_div(1.0, sla, 0.0) * safe_div(1.0, lvfraction, 0.0),
            0.0,
        )

        # --- Growth switch, with capped reserve mobilisation ---
        gap = growth_sink - growth_source
        used = np.minimum(weight_re, np.maximum(gap, 0.0))
        source_limited = growth_source < growth_sink
        growth_act = np.where(source_limited, growth_source + used, growth_sink)
        d_weight_re = np.where(source_limited, -used, growth_source - growth_sink)

        lai_growth_act = growth_act * lvfraction * sla
        lv_growth = np.where(d_weight_harv <= 0, growth_act * lvfraction, 0.0)
        lv_death = safe_div(lai_death, slaint, 0.0)
        d_lai = lai_growth_act - lai_death
        d_weight_lv = lv_growth - lv_death
        d_weight_rt = growth_act * rtfraction
        d_leaf_length_act = np.where(
            d_weight_harv > 0,
            -leaf_length,
            np.where(
                growth_sink > 0,
                d_leaf_length_pot * safe_div(growth_act, growth_sink, 0.0),
                0.0,
            ),
        )
        d_tsum = np.maximum(davtmp - temp_base, 0.0)

        # --- Water balance (PCSE WaterbalanceFD) ---
        rain_cm = rain_cm_all[i]
        irr = np.zeros(n) if rirr is None else rirr[i]
        wet = ss > 1
        evw = np.where(wet, evwmx, 0.0)
        dry_fresh = rinold >= 1
        evsmxt = evsmx * (np.sqrt(dslr + 1.0) - np.sqrt(dslr))
        evs_dry = np.where(dry_fresh, evsmx, np.minimum(evsmx, evsmxt + rinold))
        dslr_dry = np.where(dry_fresh, 1.0, dslr + 1.0)
        evs = np.where(wet, 0.0, evs_dry)
        dslr = np.where(wet, dslr, dslr_dry)

        if ifunrn:
            rinpre = (1.0 - notinf * afgen(ninftb, rain_cm)) * rain_cm
        else:
            rinpre = (1.0 - notinf) * rain_cm
        rinpre = rinpre + irr + ss
        rinpre = np.where(ss > 0.1, np.minimum(sope, rinpre - evw), rinpre)

        we = smfcf * rd
        perc1 = clamp(np.zeros(n), sope, (w - we) - tra - evs)
        welow = smfcf * (rdm_wb - rd)
        loss = clamp(np.zeros(n), ksub, wlow - welow + perc1)
        perc2 = ((rdm_wb - rd) * sm0 - wlow) + loss
        perc = np.minimum(perc1, perc2)
        rin = np.minimum(rinpre, (sm0 - sm) * rd + tra + evs + perc)
        rinold = rin

        if et_act_rs is not None:
            forced = et_act_rs[i]
            evs = np.where(np.isfinite(forced), forced, evs)
        dw = rin - tra - evs - perc
        dwlow = perc - loss
        wtmp = w + dw
        negative = wtmp < 0
        evs = np.where(negative, evs + wtmp, evs)
        dw = np.where(negative, -w, dw)
        if np.any(evs < -1e-9):
            raise ValueError(
                "LINGRA water balance: negative soil evaporation; check the "
                "soil parameters (SM0/SMW/SMFCF/RDMSOL) and the forcings."
            )
        evs = np.maximum(evs, 0.0)

        sstmp = rain_cm + irr - evw - rin
        dss = np.minimum(sstmp, ssmax - ss)
        dtsr = sstmp - dss

        # --- Record the day, then integrate ---
        store("lai", lai)
        store("dvs", dvs)
        store("tsum", tsum)
        store("tiller_number", tiller_number)
        store("weight_lv_green", weight_lv_green)
        store("weight_lv_dead", weight_lv_dead)
        store("weight_harv", weight_harv)
        store("weight_harv_annual", weight_harv - harv_year_base)
        store("weight_abg", weight_harv + weight_lv_green)
        store("weight_re", weight_re)
        store("weight_rt", weight_rt)
        store("leaf_length", leaf_length)
        store("slaint", slaint)
        store("rd", rd)
        store("sm", sm)
        store("w", w)
        store("wlow", wlow)
        store("ss", ss)
        store("rftra", rftra)
        store("tra", tra)
        store("tramx", tramx)
        store("et0", et0)
        store("e0", e0)
        store("es0", es0)
        store("evs", evs)
        store("evw", evw)
        store("perc", perc)
        store("loss", loss)
        store("runoff", dtsr)
        store("mowing", mowing.astype(np.float64))

        tsum = tsum + d_tsum
        lai = lai + d_lai
        days_after_harvest = days_after_harvest + d_days_after_harvest
        cutting_number = cutting_number + d_cutting_number
        tiller_number = tiller_number + d_tiller_number
        weight_lv_green = weight_lv_green + d_weight_lv
        weight_lv_dead = weight_lv_dead + lv_death
        weight_harv = weight_harv + d_weight_harv
        weight_re = weight_re + d_weight_re
        weight_rt = weight_rt + d_weight_rt
        leaf_length = leaf_length + d_leaf_length_act
        slaint = safe_div(lai, weight_lv_green, 0.0)
        dvs = safe_div(tsum, tsum_max, 0.0)
        temperature_soil = temperature_soil + d_temperature_soil
        rd_new = rd + rr

        w_int = w + dw
        wlow_int = wlow + dwlow
        ss = ss + dss

        # Redistribution between the root zone and the layer below it.
        rdchange = rd_new - wb_rd_old
        wdr = np.where(
            rdchange > 0.001,
            np.minimum(
                wlow_int, wlow_int * safe_div(rdchange, rdmsol - wb_rd_old, 0.0)
            ),
            w_int * safe_div(rdchange, wb_rd_old, 0.0),
        )
        wlow = wlow_int - wdr
        w = w_int + wdr
        sm = safe_div(w, rd_new, 0.0)
        dsos = np.where(sm >= (sm0 - crairc), dsos + 1.0, 0.0)
        rd = rd_new
        wb_rd_old = rd_new

        # Remote-sensing LAI overrides the integrated value when present.
        if weather.get("lai_rs") is not None:
            forced = weather["lai_rs"][i]
            lai = np.where(np.isfinite(forced), forced, lai)

        if progress is not None:
            progress(i, ndays)

    return out
