"""
MODULE:     t.graminr (graminr_lib.basgra)

AUTHOR(S):  Yann Chemin

PURPOSE:    BASGRA grassland model (carbon, nitrogen, snow, soil frost and
            winter hardening), ported from the GraMinR R package
            (R/basgra.R, R/basgra_lib.R) and vectorised over pixels.

SPDX-License-Identifier: Unlicense
"""

import numpy as np

from .core import as_array, safe_div

OUTPUTS = (
    "nell_vg",
    "lai",
    "phen",
    "evap",
    "tran",
    "tranrf",
    "tiller_v",
    "tiller_g1",
    "tiller_g2",
    "gv_biomass",
    "dv_biomass",
    "harvested",
    "dry_matter",
    "nsh",
    "nrt",
    "nmin",
    "clitt",
    "csomf",
    "csoms",
    "nlitt",
    "nsomf",
    "nsoms",
    "fngrowth",
    "nfert_tot",
    "f_wall_dm",
    "f_wall_dmsh",
    "f_wall_lv",
    "f_wall_st",
    "f_digest_dm",
    "f_digest_dmsh",
    "f_digest_lv",
    "f_digest_st",
    "f_digest_wall",
    # Diagnostics: the daily rates and indices behind the state variables,
    # in the same names BASGRA.f90 reports them under.
    "dayl",
    "parav",
    "parint",
    "sla_new",
    "resnor",
    "rleaf",
    "lerg",
    "tv2til",
    "rdrt",
    "rdrs",
    "fspot",
    "rgrtv",
    "rwa",
    "wfps",
    "luemxq",
    "nsource",
    "nsink",
    "nupt",
    "nleaching",
    "nmineralisation",
    "nshmob",
    "dmlv",
    "dmst",
    "yield_last",
    # BASGRA's own state variables, in the names BASGRA.f90 reports them
    # under: the carbon pools, the water pools, the frost and anoxia state.
    "clv",
    "clvd",
    "cres",
    "crt",
    "cst",
    "cstub",
    "tilg",
    "tiltot",
    "rootd",
    "wal",
    "wapl",
    "waps",
    "was",
    "drystor",
    "wetstor",
    "sdepth",
    "fdepth",
    "lt50",
    "o2",
    "tanaer",
    "vern",
)

REQUIRED = ("rad_glob", "tmin", "tmax", "vp", "rain", "wind", "et_pot")
# model_basgra() uses two different "no data" sentinels: et_act_rs/t_act_rs
# are tested with != 0, so 0 means "not available", while lai_rs is tested
# with is.na(), so NaN means "not available".
OPTIONAL = {"et_act_rs": 0.0, "t_act_rs": 0.0, "lai_rs": np.nan, "gcut": 0.0}

ALIASES = {}

DEFAULTS = {
    "nfertmult": 1.0,
    "rdrroot": 0.02,
    "rdrstub": 0.2,
    "fncgshmin": 0.0,
    "tcnshmob": 8.0,
    "tcnupt": 8.0,
    # Day-of-year offset applied when the day length is computed. GraMinR's R
    # uses the current day (0), which is what this port does by default. The
    # BASGRA Fortran reports, and uses, the previous day's day length: its
    # fortran_reference_output.csv DAYL column equals this formula evaluated
    # at doy - 1 on every one of its 254 days, and feeding that lagged series
    # in reproduces its whole state trajectory to single precision (see
    # tests/test_basgra_fortran.py). Set -1 to reproduce the Fortran.
    "dayl_doy_offset": 0.0,
}

#: Site constants hardcoded in parameters_site.f90 / model_init().
SITE_CONSTANTS = {
    "delt": 1.0,
    "co2a": 350.0,
    "drate": 50.0,
    "lambdaice": 1.9354e5,
    "latentheat": 335000.0,
    "poolinfillimit": 0.2,
    "rhowater": 1000.0,
    "irrigf": 0.0,
    "kmin": 4.0,
    "ampl": 0.625,
    "knfix": 0.0,
    "rrunbulk": 0.05,
}


def day_length(doy, lat):
    """Day length as a fraction of the day, BASGRA's DDAYL()."""
    rad = np.pi / 180.0
    dec = -np.arcsin(np.sin(23.45 * rad) * np.cos(2.0 * np.pi * (doy + 10.0) / 365.0))
    # tan(rad*lat) is singular at the equator; nudge it off zero so that a
    # pixel exactly on the equator stays finite (its day length is 0.5).
    tan_lat = np.tan(rad * lat)
    tan_lat = np.where(np.abs(tan_lat) < 1e-12, 1e-12, tan_lat)
    bound = np.arctan(1.0 / tan_lat)
    decc = np.maximum(np.arctan(-1.0 / tan_lat), np.minimum(bound, dec))
    inner = np.clip(tan_lat * np.tan(decc), -1.0, 1.0)
    return 0.5 * (1.0 + 2.0 * np.arcsin(inner) / np.pi)


def penman(dtr, davtmp, wind, vp, lai, rnintc):
    """BASGRA's own Penman potential soil evaporation and transpiration."""
    dtrjm2 = dtr * 1e6
    boltzm = 5.668e-8
    lhvap = 2.4e6
    psych = 0.067
    bbrad = boltzm * (davtmp + 273.0) ** 4 * 86400.0
    svp = 0.611 * np.exp(17.4 * davtmp / (davtmp + 239.0))
    slope = 4158.6 * svp / (davtmp + 239.0) ** 2
    rlwn = bbrad * np.maximum(0.0, 0.55 * (1.0 - safe_div(vp, svp, 0.0)))
    nrads = dtrjm2 * (1.0 - 0.15) - rlwn
    nradc = dtrjm2 * (1.0 - 0.25) - rlwn
    penmrs = nrads * slope / (slope + psych)
    penmrc = nradc * slope / (slope + psych)
    wdf = 2.63 * (1.0 + 0.54 * wind)
    penmd = lhvap * wdf * (svp - vp) * psych / (slope + psych)
    pevap = np.exp(-0.5 * lai) * (penmrs + penmd) / lhvap
    ptran = (1.0 - np.exp(-0.5 * lai)) * (penmrc + penmd) / lhvap
    return pevap, np.maximum(0.0, ptran - 0.5 * rnintc)


def lue_co2_tm(davtmp, rubisc, co2a, kluetilg, fractv, k, parav):
    """Light-use efficiency, temperature- and CO2-dependent (LUECO2TM)."""
    t = davtmp
    rubiscn = rubisc * (1e6 / 550000.0)
    r = 8.314
    co2i = 0.7 * co2a
    vcmax = rubiscn * 20.0 * np.exp((1.0 / 298.0 - 1.0 / (t + 273.0)) * 68000.0 / r)
    kmc = 460.0 * np.exp((1.0 / 298.0 - 1.0 / (t + 273.0)) * 65800.0 / r)
    kmo = 33.0 * np.exp((1.0 / 298.0 - 1.0 / (t + 273.0)) * 1400.0 / r)
    o2 = 21.0
    gammax = 0.5 * 0.21 * kmc * o2 / kmo
    pmax = vcmax * (co2i - gammax) / (co2i + kmc * (1.0 + o2 / kmo))
    tmpfac = np.clip((t + 4.0) / 5.0, 0.0, 1.0)
    eff = tmpfac * (1.0 / 2.1) * (co2i - gammax) / (4.5 * co2i + 10.5 * gammax)
    return safe_div(
        eff * pmax * (1.0 + kluetilg * (1.0 - fractv)), eff * k * parav + pmax, 0.0
    )


def run(params, weather, npixels, outputs=None, progress=None):
    """Run BASGRA over `npixels` pixels."""
    p = dict(params)
    n = npixels
    year = np.asarray(weather["year"])
    doy = np.asarray(weather["doy"])
    ndays = len(doy)
    wanted = tuple(OUTPUTS) if outputs is None else tuple(outputs)

    def g(key):
        """Parameter `key` broadcast over the pixel axis."""
        return as_array(p[key], n)

    c = SITE_CONSTANTS
    delt = c["delt"]

    lat = g("lat")
    k = g("k")
    ksnow = g("ksnow")
    gamma_p = g("gamma")
    lambdasoil = g("lambdasoil")
    trainsnow = g("trainsnow")
    tmeltfreeze = g("tmeltfreeze")
    swret, swrf = g("swret"), g("swrf")
    rhopack, rhonewsnow = g("rhopack"), g("rhonewsnow")
    wpoolmax = g("wpoolmax")
    wcst = g("wcst")
    wcad = g("fwcad") * wcst
    wcwp = g("fwcwp") * wcst
    wcfc = g("fwcfc") * wcst
    wcwet = g("fwcwet") * wcst
    tranco = g("tranco")
    rootdm, rrdmax = g("rootdm"), g("rrdmax")
    fgas, fo2mx = g("fgas"), g("fo2mx")
    claiv, hagere = g("claiv"), g("hagere")
    cocresmx = g("cocresmx")
    daylb, daylp, dlmxge, phencr = g("daylb"), g("daylp"), g("dlmxge"), g("phencr")
    tbase, slamax, fslamin = g("tbase"), g("slamax"), g("fslamin")
    rubisc, kluetilg = g("rubisc"), g("kluetilg")
    rehardredday = g("rehardredday")
    thardmx, lt50mn, lt50mx = g("thardmx"), g("lt50mn"), g("lt50mx")
    hparam, kresphard = g("hparam"), g("kresphard")
    tcres, cstavm, simax1t = g("tcres"), g("cstavm"), g("simax1t")
    nellvm = g("nellvm")
    lfwidv, lfwidg, shape, yg = g("lfwidv"), g("lfwidg"), g("shape"), g("yg")
    ncshmax, ncr = g("ncshmax"), g("ncr")
    tcnshmob, tcnupt = g("tcnshmob"), g("tcnupt")
    fncgshmin = g("fncgshmin")
    krtotaer = g("krtotaer")
    ldt50a, ldt50b, krdranaer = g("ldt50a"), g("ldt50b"), g("krdranaer")
    krsr3h, dparam, tsurfdiff, ratedmx = (
        g("krsr3h"),
        g("dparam"),
        g("tsurfdiff"),
        g("ratedmx"),
    )
    laicr, rdrsco, rdrsmx = g("laicr"), g("rdrsco"), g("rdrsmx")
    dayl_doy_offset = float(np.asarray(g("dayl_doy_offset")).reshape(-1)[0])
    # An externally supplied day length, in fraction of a day, overrides the
    # computed one. The BASGRA Fortran's own series is needed to reproduce
    # its reference run, because the value it uses on the first day is left
    # over from initialisation and is not any day's day length.
    dayl_forcing = weather.get("dayl")
    if dayl_forcing is not None:
        dayl_forcing = np.asarray(dayl_forcing, dtype=np.float64).reshape(len(doy), -1)
    rdrtmin, rdrtem = g("rdrtmin"), g("rdrtem")
    rdrstub, rdrroot = g("rdrstub"), g("rdrroot")
    phy, laitil, laieft, fsmax = g("phy"), g("laitil"), g("laieft"), g("fsmax")
    toptge, rgenmx = g("toptge"), g("rgenmx")
    daylg1g2, rgrtg1g2 = g("daylg1g2"), g("rgrtg1g2")
    tvern = g("tvern")
    tmaxf, tsigmaf = g("tmaxf"), g("tsigmaf")
    tclitt, tcsomf, tcsoms = g("tclitt"), g("tcsomf"), g("tcsoms")
    flittsomf, fsomfsoms = g("flittsomf"), g("fsomfsoms")
    rnleach, knemit = g("rnleach"), g("knemit")
    rfn2o, wfps50n2o = g("rfn2o"), g("wfps50n2o")
    f_wall_lv_fmin, f_wall_lv_max = g("f_wall_lv_fmin"), g("f_wall_lv_max")
    f_wall_st_fmin, f_wall_st_max = g("f_wall_st_fmin"), g("f_wall_st_max")
    f_digest_wall_fmin, f_digest_wall_max = (
        g("f_digest_wall_fmin"),
        g("f_digest_wall_max"),
    )
    bias = c["kmin"] + c["ampl"]

    # --- Initial state (BASGRA.f90 initial constants) ---
    clv = 10.0 ** g("log10clvi")
    cres = 10.0 ** g("log10cresi")
    crt = 10.0 ** g("log10crti")
    laii = 10.0 ** g("log10laii")
    cst = g("csti").copy()
    lai = laii.copy()
    lt50 = g("lt50i").copy()
    o2 = fgas * rootdm * fo2mx * 1000.0 / 22.4
    phen = g("pheni").copy()
    rootd = rootdm.copy()
    tiltoti, frtilgi, frtilgg1i = g("tiltoti"), g("frtilgi"), g("frtilgg1i")
    tilg1 = tiltoti * frtilgi * frtilgg1i
    tilg2 = tiltoti * frtilgi * (1.0 - frtilgg1i)
    tilv = tiltoti * (1.0 - frtilgi)
    wal = 1000.0 * rootdm * g("wci")

    nrt = ncr * crt
    ncshi = safe_div(ncshmax * (1.0 - np.exp(-k * laii)), k * laii, 0.0)
    nsh = ncshi * (clv + cst)

    clitt = g("clitt0").copy()
    csom0, fcsomf0 = g("csom0"), g("fcsomf0")
    csomf = csom0 * fcsomf0
    csoms = csom0 * (1.0 - fcsomf0)
    nlitt = safe_div(g("clitt0"), g("cnlitt0"), 0.0)
    nsomf = safe_div(csom0 * fcsomf0, g("cnsomf0"), 0.0)
    nsoms = safe_div(csom0 * (1.0 - fcsomf0), g("cnsoms0"), 0.0)
    nmin = g("nmin0").copy()
    nfert_tot = np.zeros(n)

    clvd = np.zeros(n)
    cstub = np.zeros(n)
    drystor = np.zeros(n)
    fdepth = np.zeros(n)
    sdepth = np.zeros(n)
    tanaer = np.zeros(n)
    vern = np.ones(n)
    wapl = np.zeros(n)
    waps = np.zeros(n)
    was = np.zeros(n)
    wetstor = np.zeros(n)
    tsurf = np.zeros(n)
    hv_biomass = np.zeros(n)
    yield_last = np.zeros(n)
    evap = np.zeros(n)
    tran = np.zeros(n)

    rad_glob = weather["rad_glob"]
    tmin, tmax = weather["tmin"], weather["tmax"]
    vp_all, wind_all = weather["vp"], weather["wind"]
    rain_all = weather["rain"]
    et_pot_all = weather["et_pot"]
    gcut = weather["gcut"]
    rs_evap = weather.get("et_act_rs")
    rs_tran = weather.get("t_act_rs")
    rs_lai = weather.get("lai_rs")
    nfert_day = weather.get("nfert")
    ndep_day = weather.get("ndep")

    out = {name: np.zeros((ndays, n), dtype=np.float32) for name in wanted}

    def store(name, value):
        if name in out:
            out[name][i] = value

    for i in range(ndays):
        rain = rain_all[i]
        gr = rad_glob[i]
        vp = vp_all[i]
        wind = wind_all[i]
        rsetpot = et_pot_all[i]
        rsevap = np.zeros(n) if rs_evap is None else rs_evap[i]
        rstran = np.zeros(n) if rs_tran is None else rs_tran[i]
        rslai = np.full(n, np.nan) if rs_lai is None else rs_lai[i]
        rsgcut = gcut[i]

        # Remote-sensing overrides, with the sentinels model_basgra() uses.
        evap = np.where(rsevap != 0, rsevap, evap)
        tran = np.where(rstran != 0, rstran, tran)
        lai = np.where(np.isfinite(rslai), rslai, lai)

        davtmp = (tmin[i] + tmax[i]) / 2.0
        dtr = gr * np.exp(-ksnow * drystor)
        par = 0.5 * 4.56 * dtr
        if dayl_forcing is None:
            dayl = day_length(float(doy[i]) + dayl_doy_offset, lat)
        else:
            dayl = np.broadcast_to(dayl_forcing[i], lat.shape)

        # --- Physics: surface temperature and soil frost ---
        wcl = np.where(fdepth < rootd, safe_div(wal * 0.001, rootd - fdepth, 0.0), 0.0)
        frozen = fdepth > 0.0
        tsurf = np.where(
            frozen,
            safe_div(davtmp, 1.0 + 10.0 * safe_div(sdepth, fdepth, 0.0), 0.0),
            davtmp * np.exp(-gamma_p * sdepth),
        )
        fperm = np.where(frozen, 0.0, 1.0)

        wceff = np.where(
            fdepth > rootd,
            wcfc,
            np.where(frozen, safe_div(0.001 * was, fdepth, 0.0), wcl),
        )
        no_frost = ((fdepth == 0.0) & (tsurf > 0.0)) | (wceff == 0.0)
        alpha = safe_div(lambdasoil, c["rhowater"] * wceff * c["latentheat"], 0.0)
        pfrate = np.where(
            no_frost,
            0.0,
            np.sqrt(np.maximum(0.0, fdepth**2 - 2.0 * alpha * tsurf)) - fdepth,
        )
        frate = np.where(
            (pfrate >= 0.0) & (fdepth > 0.0) & (fdepth < rootd),
            safe_div(pfrate * safe_div(0.001 * was, fdepth, 0.0), wcfc, 0.0),
            np.where(pfrate + fdepth / delt < 0.0, -fdepth / delt, pfrate),
        )

        # --- Microclimate: snow and the surface pool ---
        warm = davtmp > trainsnow
        pwater = np.where(warm, rain, 0.0)
        psnow = np.where(warm, 0.0, rain)

        melt = bias + c["ampl"] * dayl
        snowmelt = np.where(
            davtmp > tmeltfreeze,
            np.maximum(0.0, np.minimum(drystor / delt, melt * (davtmp - tmeltfreeze))),
            0.0,
        )
        wmaxstore = drystor * swret

        refreezemax = swrf * (tmeltfreeze - davtmp)
        refreeze = np.where(
            (wetstor > 0) & (davtmp < tmeltfreeze),
            np.minimum(wetstor / delt, refreezemax),
            0.0,
        )
        staywet = wetstor / delt - refreeze

        wavail = staywet + snowmelt + pwater
        wremain = np.minimum(wavail, wmaxstore)
        wsupply = wavail - wremain

        swe = drystor + wetstor
        density = np.where(
            sdepth > 0, np.minimum(480.0, safe_div(swe, sdepth, 0.0)), 0.0
        )
        packmelt = np.where(
            (sdepth > 0) & (density > 0),
            np.maximum(
                0.0,
                np.minimum(
                    sdepth / delt, sdepth * rhopack - safe_div(snowmelt, density, 0.0)
                ),
            ),
            0.0,
        )

        rnintc = np.minimum(wsupply, 0.25 * lai)
        pinfil = wsupply - rnintc
        infil = np.where(fdepth <= c["poolinfillimit"], pinfil, 0.0)
        runon = pinfil - infil

        poolvolremain = np.maximum(0.0, wpoolmax - wapl - waps)
        poolinfil = np.minimum(runon, poolvolremain)
        poolwavail = poolinfil + wapl / delt
        pooldrain = np.where(
            poolwavail == 0,
            0.0,
            np.where(
                fdepth <= c["poolinfillimit"],
                poolwavail,
                np.maximum(0.0, np.minimum(-frate * 1000.0, poolwavail)),
            ),
        )
        eta_ice = c["lambdaice"] / (c["rhowater"] * c["latentheat"])
        pirate = np.where(
            (tsurf > 0) & (wapl == 0) & (waps == 0),
            0.0,
            np.sqrt(np.maximum(0.0, (0.001 * waps) ** 2 - 2.0 * eta_ice * tsurf * delt))
            / delt
            - (0.001 * waps) / delt,
        )
        freezepl = np.where(
            pirate < 0,
            0.0,
            np.maximum(
                0.0,
                np.minimum(poolinfil + wapl / delt - pooldrain * delt, pirate * 1000.0),
            ),
        )
        thawps = np.where(pirate < 0, np.minimum(waps / delt, -pirate * 1000.0), 0.0)

        permgas = np.where(waps == 0, 1.0, 0.0)

        # --- Evapotranspiration ---
        rs_all = (rsevap != 0) & (rstran != 0) & (rsetpot != 0)
        pevap, ptran = penman(dtr, davtmp, wind, vp, lai, rnintc)
        rwa = np.clip(safe_div(wcl - wcad, wcfc - wcad, 0.0), 0.0, 1.0)
        wfps = np.clip(safe_div(wcl - wcad, wcst - wcad, 0.0), 0.0, 1.0)
        waad = 1000.0 * wcad * (rootd - fdepth)
        evap_c = pevap * rwa
        wccr = wcwp + np.maximum(
            0.01, safe_div(ptran, ptran + tranco, 0.0) * (wcfc - wcwp)
        )
        fr = np.where(
            wcl > wccr,
            np.where(
                (wcst - wcwet) != 0,
                np.clip(safe_div(wcst - wcl, wcst - wcwet, 0.0), 0.0, 1.0),
                1.0,
            ),
            np.clip(safe_div(wcl - wcwp, wccr - wcwp, 0.0), 0.0, 1.0),
        )
        tran_c = ptran * fr
        availf = np.where(
            (evap_c + tran_c) > 0,
            np.minimum(1.0, safe_div((wal - waad) / delt, evap_c + tran_c, 0.0)),
            0.0,
        )
        evap_c = evap_c * availf
        tran_c = tran_c * availf
        tranrf_c = np.where(ptran > 0, safe_div(tran_c, ptran, 0.0), 1.0)

        tranrf_rs = np.where(
            rsetpot > (rstran + rsevap), safe_div(rstran + rsevap, rsetpot, 1.0), 1.0
        )
        evap = np.where(rs_all, rsevap, evap_c)
        tran = np.where(rs_all, rstran, tran_c)
        tranrf = np.where(rs_all, tranrf_rs, tranrf_c)
        # rwa/wfps have no remote-sensing counterpart (this override is a
        # GraMinR extension with no equivalent in BASGRA.f90).
        rwa = np.where(rs_all, 0.0, rwa)
        wfps = np.where(rs_all, 0.0, wfps)

        # --- Light interception ---
        parav = np.where(dayl > 0, safe_div(par * (1e6 / (24 * 3600)), dayl, 0.0), 0.0)
        parint = par * (1.0 - np.exp(-k * lai))

        # --- Root growth ---
        rrootd = np.where(
            (rootd < rootdm) & (wcl > wcwp),
            np.minimum(rrdmax, (rootdm - rootd) / delt),
            0.0,
        )
        explor = 1000.0 * rrootd * wcfc

        # --- Soil water fluxes ---
        wafc = 1000.0 * wcfc * np.maximum(0.0, rootd - fdepth)
        wast = 1000.0 * wcst * np.maximum(0.0, rootd - fdepth)
        infiltot = infil + pooldrain
        freezel = np.where(
            fdepth < rootd,
            np.maximum(
                0.0,
                np.minimum(
                    wal / delt + (infiltot - evap - tran),
                    safe_div(frate, rootd - fdepth, 0.0) * wal,
                ),
            ),
            0.0,
        )
        thaws = np.where(
            (fdepth > 0) & (fdepth <= rootd),
            np.maximum(
                0.0, np.minimum(was / delt, safe_div(-frate * was, fdepth, 0.0))
            ),
            0.0,
        )
        net = infiltot - evap - tran - freezel + thaws
        drain = np.maximum(0.0, np.minimum(c["drate"], (wal - wafc) / delt + net))
        runoff = np.maximum(0.0, (wal - wast) / delt + (net - drain))
        irrig = c["irrigf"] * ((wafc - wal) / delt - (net - drain - runoff))

        fo2 = safe_div(o2, rootd * fgas * 1000.0 / 22.4, 0.0)

        # --- Harvest ---
        harv = np.where(rsgcut != 0, 1.0, 0.0)
        noharv = 1.0 - harv
        fractv = safe_div(tilv + tilg1, tilv + tilg1 + tilg2, 0.0)
        clai = fractv * claiv
        harvfr = np.where(lai <= clai, 0.0, 1.0 - safe_div(clai, lai, 0.0))
        harvla = harv * lai * harvfr / delt
        harvlv = harv * clv * harvfr / delt
        harvph = harv * phen / delt
        harvre = harv * ((harvfr * fractv) + (1.0 - fractv) * hagere) * cres / delt
        harvst = harv * cst / delt
        gstub = harvst * (1.0 - hagere)
        harvtilg2 = harv * tilg2 / delt

        # --- Reserves ---
        cresmx = cocresmx * (clv + cres + cst)
        resnor = np.clip(safe_div(cres, cresmx, 0.0), 0.0, 1.0)

        # --- Phenology ---
        gphen = np.maximum(
            0.0, (davtmp - 0.01) * 0.000144 * 24.0 * (np.minimum(daylp, dayl) - 0.24)
        )
        dphen = np.where(dayl < daylb, phen / delt, 0.0)
        phenrf = np.clip(safe_div(1.0 - phen, 1.0 - phencr, 0.0), 0.0, 1.0)
        daylge = np.clip(safe_div(dayl - daylb, dlmxge - daylb, 0.0), 0.0, 1.0)

        # --- Leaf growth potential ---
        efftmp = np.maximum(tbase, davtmp)
        lerv = np.maximum(0.0, (-0.76 + 0.52 * efftmp) / 1000.0)
        lerg = daylge * np.maximum(0.0, (-5.46 + 2.80 * efftmp) / 1000.0)
        slamin = slamax * fslamin
        slanew = slamax - resnor * (slamax - slamin)

        luemxq = lue_co2_tm(davtmp, rubisc, c["co2a"], kluetilg, fractv, k, parav)

        # --- Hardening sink ---
        rehardredend = np.where(lat > 0, 91.0, 91.0 + 183.0)
        rehardredstart = np.mod(rehardredend - rehardredday, 365.0)
        doysincestart = np.mod(float(doy[i]) - rehardredstart, 365.0)
        rehardperiod = np.where(
            doysincestart < (rehardredday + 0.5 * (365.0 - rehardredday)),
            np.maximum(0.0, 1.0 - safe_div(doysincestart, rehardredday, 0.0)),
            1.0,
        )
        rateh = np.where(
            (tsurf > thardmx) | (lt50 < lt50mn),
            0.0,
            rehardperiod * hparam * (thardmx - tsurf) * (lt50 - lt50mn),
        )
        resphardsi = rateh * clv * kresphard * np.clip(resnor * 5.0, 0.0, 1.0)

        # --- Growth ---
        phot = parint * tranrf * 12.0 * luemxq * noharv
        resmob = (cres * noharv / tcres) * np.clip(davtmp / 5.0, 0.0, 1.0)
        source = resmob + phot
        resphard = np.minimum(source, resphardsi)
        allotot = source - resphard
        gressi = 0.5 * (resmob + np.maximum(0.0, (cresmx - cres) / delt))
        cstav = np.where(tilg2 != 0.0, safe_div(cst, tilg2, 0.0), 0.0)
        sink1t = np.maximum(0.0, 1.0 - safe_div(cstav, cstavm, 0.0)) * simax1t
        nellvg = np.minimum(phenrf * nellvm, nellvm)
        glaisi = (
            (
                (lerv * (tilv + tilg1) * nellvm * lfwidv)
                + (lerg * tilg2 * nellvg * lfwidg)
            )
            * shape
            * tranrf
        )
        glvsi = np.maximum(0.0, safe_div(glaisi * noharv, slanew, 0.0) / yg)
        gstsi = np.maximum(0.0, (sink1t * tilg2 * tranrf * noharv) / yg)

        # Nitrogen limitation of growth
        nshk = (clv + cst) * ncshmax * safe_div(1.0 - np.exp(-k * lai), k * lai, 0.0)
        nshexcess = np.maximum(0.0, nsh - nshk)
        nshmob = noharv * nshexcess / tcnshmob
        nsource = nmin / tcnupt + nshmob
        nsink = np.maximum(0.0, (glvsi + gstsi) * ncshmax)
        fngrowth = np.minimum(1.0, safe_div(nsource, nsink, 1.0))
        glaisi = glaisi * fngrowth
        glvsi = glvsi * fngrowth
        gstsi = gstsi * fngrowth

        # Allocation
        gshsi = glvsi + gstsi
        spring = daylge >= 0.1
        allosh_spring = np.minimum(allotot, gshsi)
        gres_spring = np.minimum(allotot - allosh_spring, gressi)
        gres_autumn = np.minimum(allotot, gressi)
        allosh_autumn = np.minimum(allotot - gres_autumn, gshsi)
        allosh = np.where(spring, allosh_spring, allosh_autumn)
        gres = np.where(spring, gres_spring, gres_autumn)
        allort = allotot - allosh - gres
        gshsi_safe = np.where(gshsi == 0.0, 1.0, gshsi)
        allolv = glvsi * (allosh / gshsi_safe)
        allost = gstsi * (allosh / gshsi_safe)
        glv = allolv * yg
        gst = allost * yg
        grt = allort * yg
        respgsh = (allolv + allost) * (1.0 - yg)
        respgrt = allort * (1.0 - yg)

        faer = np.clip(safe_div(fo2, fo2mx, 0.0), 0.0, 1.0)
        rplantaer = faer * (respgrt + respgsh + resphard)

        # --- Senescence ---
        dtanaer = np.where(permgas == 0, 1.0, -tanaer / delt)
        ld50 = ldt50a + ldt50b * lt50
        rdrtox = np.where(
            tanaer > 0,
            krdranaer / (1.0 + np.exp(-krdranaer * (tanaer - ld50))),
            0.0,
        )
        rsr3h = 1.0 / (1.0 + np.exp(-krsr3h * (tsurf - lt50)))
        rdrfrost = np.minimum(0.5, 1.0 - rsr3h)
        rated = np.minimum(
            dparam * (lt50mx - lt50) * (tsurf + tsurfdiff), (lt50mx - lt50) / delt
        )
        dehardrate = np.clip(rated, 0.0, ratedmx)
        hardrate = safe_div(resphard, clv * kresphard, 0.0)

        tv1 = np.where(lai < laicr, 0.0, rdrsco * safe_div(lai - laicr, laicr, 0.0))
        rdrs = np.minimum(tv1, rdrsmx)
        rdrt = np.maximum(rdrtmin, rdrtem * tsurf)
        tv2 = noharv * np.maximum.reduce([rdrs, rdrt, rdrfrost, rdrtox])
        tv2til = noharv * np.maximum.reduce([rdrs, rdrfrost, rdrtox])
        dlai = lai * tv2
        dlv = clv * tv2
        dstub = cstub * rdrstub
        dtilv = tilv * tv2til
        drt = crt * rdrroot

        # --- Foliage 2: tiller dynamics ---
        glai = slanew * glv
        tv1f = np.where(tsurf < tbase, 0.0, safe_div(tsurf, phy, 0.0))
        rleaf = (
            tv1f
            * noharv
            * tranrf
            * daylge
            * (fractv + phenrf * (1.0 - fractv))
            * fngrowth
        )
        fspot = np.clip(laitil - laieft * lai, 0.0, fsmax)
        rgrtv = np.maximum(0.0, fspot * resnor * rleaf)
        gtilv = tilv * rgrtv
        tge = np.maximum(
            0.0, 1.0 - safe_div(np.abs(davtmp - toptge), toptge - tbase, 0.0)
        )
        rgrtvg1 = np.minimum(1.0 - tv2til, noharv * daylge * tge * rgenmx) * vern
        tilvg1 = tilv * rgrtvg1
        tilg1g2 = np.where(dayl > daylg1g2, tilg1 * rgrtg1g2, 0.0)

        # --- Oxygen fluxes ---
        o2mx = fo2mx * rootd * fgas * 1000.0 / 22.4
        o2out = rplantaer * krtotaer / 12.0
        o2in = permgas * ((o2mx - o2) + o2out * delt)

        nfert = np.zeros(n) if nfert_day is None else nfert_day[i]
        ndep = np.zeros(n) if ndep_day is None else ndep_day[i]

        # --- Soil carbon and nitrogen ---
        ftsoil = np.exp(
            (tsurf - 10.0) * (2.0 * tmaxf - tsurf - 10.0) / (2.0 * tsigmaf**2)
        )
        rclitt = safe_div(clitt * runoff, rootd, 0.0) * c["rrunbulk"] * 0.001
        dclitt = clitt * ftsoil / tclitt
        dclittsomf = flittsomf * dclitt
        dclittrsoil = dclitt - dclittsomf
        rcsomf = safe_div(csomf * runoff, rootd, 0.0) * c["rrunbulk"] * 0.001
        dcsomf = csomf * ftsoil / tcsomf
        dcsomfsoms = fsomfsoms * dcsomf
        dcsomfrsoil = dcsomf - dcsomfsoms
        dcsoms = csoms * ftsoil / tcsoms
        rnlitt = safe_div(nlitt * runoff, rootd, 0.0) * c["rrunbulk"] * 0.001
        dnlitt = safe_div(nlitt * dclitt, clitt, 0.0)
        nlittsomf = dnlitt * flittsomf
        nlittnmin = dnlitt - nlittsomf
        rnsomf = safe_div(nsomf * runoff, rootd, 0.0) * c["rrunbulk"] * 0.001
        dnsomf = safe_div(nsomf * dcsomf, csomf, 0.0)
        nsomfsoms = dnsomf * fsomfsoms
        nsomfnmin = dnsomf - nsomfsoms
        dnsoms = safe_div(nsoms * dcsoms, csoms, 0.0)
        nmineralisation = nlittnmin + nsomfnmin + dnsoms
        nfixation = grt * c["knfix"]
        nleaching = np.where(
            (wal > 0) & (nmin > 0), safe_div(nmin * rnleach * drain, wal, 0.0), 0.0
        )
        nemission = nmin * knemit * rwa
        fn2o = 1.0 / (1.0 + np.exp(-rfn2o * (wfps - wfps50n2o)))
        nemission_tot = nemission

        # --- Plant nitrogen ---
        nshnor = safe_div(nsh, (clv + cst) * ncshmax, 0.0)
        knmax = np.where(nshnor * lai > 0, safe_div(1.0, nshnor * lai, 1e10), 1e10)
        kn_hi = safe_div(
            1.5 - 3.0 * np.sqrt(np.maximum(0.0, 0.25 - (1.0 - nshnor) * 2.0 / 3.0)),
            lai,
            0.0,
        )
        kn = np.clip(np.where(nshnor < 5.0 / 8.0, knmax, kn_hi), 0.0, knmax)
        fncgrowth = fncgshmin + (1.0 - fncgshmin) * fngrowth
        gnsh = np.where(
            glai > 0,
            np.where(
                kn > 0,
                fncgrowth
                * ncshmax
                * (glv + gst)
                * safe_div(1.0 - np.exp(-kn * glai), kn * glai, 0.0),
                fncgrowth * ncshmax * (glv + gst),
            ),
            0.0,
        )
        dnsh_raw = np.where(
            kn > 0,
            nsh
            * safe_div(clv, clv + cst, 0.0)
            * (
                1.0
                - safe_div(
                    1.0 - np.exp(-kn * (lai - dlai)), 1.0 - np.exp(-kn * lai), 0.0
                )
            ),
            nsh * safe_div(clv, clv + cst, 0.0) * safe_div(dlai, lai, 0.0),
        )
        dnsh = np.where(
            dlai > 0,
            np.maximum(dnsh_raw, nsh - ncshmax * ((clv + cst) - dlv)),
            0.0,
        )
        harvnsh_raw = np.where(
            kn > 0,
            nsh * safe_div(1.0 - np.exp(-kn * harvla), 1.0 - np.exp(-kn * lai), 0.0),
            nsh * safe_div(harvla, lai, 0.0),
        )
        harvnsh = np.where(
            harvla > 0,
            np.maximum(harvnsh_raw, nsh - ncshmax * ((clv + cst) - (harvlv + harvst))),
            0.0,
        )
        gnrt = ncr * grt * fngrowth
        dnrt = drt * safe_div(nrt, crt, 0.0)
        gnmob = np.minimum(nshmob, gnsh + gnrt)
        nshmobsoil = nshmob - gnmob
        nupt = gnsh + gnrt - gnmob

        # --- Dry matter and digestibility ---
        dmlv = clv / 0.45
        dmst = cst / 0.45
        dmstub = cstub / 0.45
        dmres = cres / 0.40
        dmsh = dmlv + dmst + dmres
        dm_today = dmsh + dmstub
        f_wall_lv_min = f_wall_lv_fmin * f_wall_lv_max
        f_wall_st_min = f_wall_st_fmin * f_wall_st_max
        f_wall_lv = f_wall_lv_min + phen * (f_wall_lv_max - f_wall_lv_min)
        f_wall_st = f_wall_st_min + phen * (f_wall_st_max - f_wall_st_min)
        f_wall_dm = safe_div(
            f_wall_lv * dmlv + f_wall_st * dmst + dmstub, dm_today, 0.0
        )
        f_wall_dmsh = safe_div(f_wall_lv * dmlv + f_wall_st * dmst, dmsh, 0.0)
        f_digest_wall_min = f_digest_wall_fmin * f_digest_wall_max
        f_digest_wall = f_digest_wall_max - phen * (
            f_digest_wall_max - f_digest_wall_min
        )
        f_digest_lv = 1.0 - f_wall_lv + f_digest_wall * f_wall_lv
        f_digest_st = 1.0 - f_wall_st + f_digest_wall * f_wall_st
        f_digest_dmsh = safe_div(
            f_digest_lv * dmlv + f_digest_st * dmst + dmres, dmsh, 0.0
        )
        f_digest_dm = safe_div(
            f_digest_lv * dmlv + f_digest_st * dmst + dmres + f_digest_wall * dmstub,
            dm_today,
            0.0,
        )

        # --- Outputs: recorded before the state increments, as in BASGRA.f90 ---
        store("nell_vg", nellvg)
        store("lai", lai)
        store("phen", phen)
        store("evap", evap)
        store("tran", tran)
        store("tranrf", tranrf)
        store("tiller_v", tilv)
        store("tiller_g1", tilg1)
        store("tiller_g2", tilg2)
        store("gv_biomass", clv + cst)
        store("dv_biomass", clvd)
        store("harvested", hv_biomass)
        store("dry_matter", dm_today)
        store("nsh", nsh)
        store("nrt", nrt)
        store("nmin", nmin)
        store("clitt", clitt)
        store("csomf", csomf)
        store("csoms", csoms)
        store("nlitt", nlitt)
        store("nsomf", nsomf)
        store("nsoms", nsoms)
        store("fngrowth", fngrowth)
        store("nfert_tot", nfert_tot)
        store("f_wall_dm", f_wall_dm)
        store("f_wall_dmsh", f_wall_dmsh)
        store("f_wall_lv", f_wall_lv)
        store("f_wall_st", f_wall_st)
        store("f_digest_dm", f_digest_dm)
        store("f_digest_dmsh", f_digest_dmsh)
        store("f_digest_lv", f_digest_lv)
        store("f_digest_st", f_digest_st)
        store("f_digest_wall", f_digest_wall)
        store("dayl", dayl)
        store("parav", parav)
        store("parint", parint)
        store("sla_new", slanew)
        store("resnor", resnor)
        store("rleaf", rleaf)
        store("lerg", lerg)
        store("tv2til", tv2til)
        store("rdrt", rdrt)
        store("rdrs", rdrs)
        store("fspot", fspot)
        store("rgrtv", rgrtv)
        store("rwa", rwa)
        store("wfps", wfps)
        store("luemxq", luemxq)
        store("nsource", nsource)
        store("nsink", nsink)
        store("nupt", nupt)
        store("nleaching", nleaching)
        store("nmineralisation", nmineralisation)
        store("nshmob", nshmob)
        store("dmlv", dmlv)
        store("dmst", dmst)
        store("yield_last", yield_last)
        store("clv", clv)
        store("clvd", clvd)
        store("cres", cres)
        store("crt", crt)
        store("cst", cst)
        store("cstub", cstub)
        store("tilg", tilg1 + tilg2)
        store("tiltot", tilv + tilg1 + tilg2)
        store("rootd", rootd)
        store("wal", wal)
        store("wapl", wapl)
        store("waps", waps)
        store("was", was)
        store("drystor", drystor)
        store("wetstor", wetstor)
        store("sdepth", sdepth)
        store("fdepth", fdepth)
        store("lt50", lt50)
        store("o2", o2)
        store("tanaer", tanaer)
        store("vern", vern)

        # --- State equations ---
        cres = cres + gres - resmob - harvre
        clv = clv + glv - dlv - harvlv
        clvd = clvd + dlv
        dm_hv = (harvlv + harvst * hagere) / 0.45 + harvre / 0.40
        hv_biomass = hv_biomass + np.where(dm_hv > 0, dm_hv, 0.0)
        yield_last = np.where(dm_hv > 0, dm_hv, yield_last)
        crt = crt + grt - drt
        cst = cst + gst - harvst
        cstub = cstub + gstub - dstub
        drystor = drystor + refreeze + psnow - snowmelt
        fdepth = fdepth + frate
        lai_int = lai + glai - dlai - harvla
        lai = np.where(np.isfinite(rslai), rslai, lai_int)
        lt50 = lt50 + dehardrate - hardrate
        o2 = o2 + o2in - o2out
        phen = np.minimum(1.0, phen + gphen - dphen - harvph)
        rootd = rootd + rrootd
        sdepth = sdepth + safe_div(psnow, rhonewsnow, 0.0) - packmelt
        tanaer = tanaer + dtanaer
        tilg1 = tilg1 + tilvg1 - tilg1g2
        tilg2 = tilg2 + tilg1g2 - harvtilg2
        tilv = tilv + gtilv - tilvg1 - dtilv

        # Vernalisation: cleared on a fixed hemisphere-dependent day, then
        # re-armed as soon as it gets cold enough.
        reset_north = (lat > 0) & (doy[i] == 305)
        reset_south = (lat < 0) & (doy[i] == 122)
        vern = np.where(reset_north | reset_south, 0.0, vern)
        vern = np.where(davtmp < tvern, 1.0, vern)

        # NB: `poolinfil` is deliberately absent from the wapl balance below.
        # GraMinR's surfacepool() computes it but does not return it, so the
        # R model always adds 0 here; preserved so that results match.
        wal = (
            wal
            + thaws
            - freezel
            + pooldrain
            + infil
            + explor
            + irrig
            - drain
            - runoff
            - evap
            - tran
        )
        wapl = wapl + thawps - freezepl - pooldrain
        waps = waps - thawps + freezepl
        was = was - thaws + freezel
        # BASGRA.f90 has WETSTOR = WETSTOR + Wremain - WETSTOR, i.e. a reset.
        wetstor = wremain

        nrt = nrt + gnrt - dnrt
        nsh = nsh + gnsh - dnsh - harvnsh - nshmob
        clitt = clitt + dlv + dstub - rclitt - dclittrsoil - dclittsomf
        csomf = csomf + drt + dclittsomf - rcsomf - dcsomfrsoil - dcsomfsoms
        csoms = csoms + dcsomfsoms - dcsoms
        nlitt = nlitt + dnsh - rnlitt - dnlitt
        nsomf = nsomf + dnrt + nlittsomf - rnsomf - dnsomf
        nsoms = nsoms + nsomfsoms - dnsoms
        nmin = np.maximum(
            0.0,
            nmin
            + ndep
            + nfert
            + nmineralisation
            + nfixation
            + nshmobsoil
            - nupt
            - nleaching
            - nemission_tot,
        )
        nfert_tot = nfert_tot + nfert

        if progress is not None:
            progress(i, ndays)

    return out
