#!/usr/bin/env python3

"""
MODULE:     t.graminr

AUTHOR(S):  Yann Chemin

PURPOSE:    Spatially distributed grassland growth modelling: runs the
            GraMinR grassland models (LINGRA, ModVege, BASGRA, GrasProg,
            GrasSim) on every grassland pixel of the current region and
            writes one space-time raster dataset per output variable.

            The module can acquire everything it needs by itself: the
            climate forcings through t.in.era5 (ERA5-Land reanalysis) or
            t.in.ecmwf (IFS/AIFS open-data forecasts), the grassland mask
            through r.in.landcover, and the optional remotely sensed LAI
            and evapotranspiration forcings through r.in.modis.

            Progress is reported through gs.percent(), which MetaMonitor
            picks up; see https://github.com/YannChemin/MetaMonitor

SPDX-License-Identifier: Unlicense
"""

# %module
# % description: Runs the GraMinR grassland growth models over every grassland pixel of a region and writes one space-time raster dataset per output variable.
# % keyword: temporal
# % keyword: raster
# % keyword: grassland
# % keyword: biomass
# % keyword: yield
# % keyword: agriculture
# % keyword: model
# % keyword: strds
# %end

# %option
# % key: model
# % type: string
# % required: yes
# % multiple: no
# % options: lingra,modvege,basgra,grasprog,grassim
# % answer: modvege
# % descriptions: lingra;LINGRA, light-interception driven, with the PCSE water balance;modvege;ModVege (Jouven et al. 2006), four biomass compartments;basgra;BASGRA, carbon and nitrogen with snow, soil frost and winter hardening;grasprog;GrasProg (FOPROQ), growth indices per regrowth compartment;grassim;GrasSim, four compartments with a full nitrogen balance
# % description: Grassland model to run
# % guisection: Model
# %end

# %option
# % key: start
# % type: string
# % required: yes
# % description: First simulated day, YYYY-MM-DD
# % guisection: Model
# %end

# %option
# % key: end
# % type: string
# % required: yes
# % description: Last simulated day, YYYY-MM-DD (inclusive)
# % guisection: Model
# %end

# %option G_OPT_R_BASENAME_OUTPUT
# % key: output
# % required: yes
# % description: Prefix of the outputs; each variable becomes the STRDS <output>_<variable>
# % guisection: Model
# %end

# %option
# % key: variables
# % type: string
# % required: no
# % multiple: yes
# % description: Output variables to write (default: a small, model-specific selection; use -l to list every variable a model offers, or variables=all for all of them)
# % guisection: Model
# %end

# %option G_OPT_F_INPUT
# % key: parameters
# % required: no
# % description: Two-column name,value parameter CSV, in GraMinR's read_config() format (default: the set shipped for the chosen model)
# % guisection: Model
# %end

# %option
# % key: parameter_maps
# % type: string
# % required: no
# % multiple: yes
# % key_desc: name=raster
# % description: Parameters that vary per pixel, as name=raster pairs (for example WHC=soil_whc,NI=nutrition_index)
# % guisection: Model
# %end

# %option
# % key: parameter_values
# % type: string
# % required: no
# % multiple: yes
# % key_desc: name=value
# % description: Parameters overridden on the command line, as name=value pairs (for example ST1=550,NI=0.8)
# % guisection: Model
# %end

# %option G_OPT_R_INPUT
# % key: grassland
# % required: no
# % description: Raster restricting the simulation to grassland. Either a land-cover map, with its classes given by grassland_classes=, or a ready-made mask whose non-null, non-zero cells are grassland. Downloaded through r.in.landcover when left empty
# % guisection: Grassland
# %end

# %option
# % key: grassland_classes
# % type: string
# % required: no
# % multiple: yes
# % description: Class codes of grassland in grassland=, when it is a land-cover map (default: the grassland classes of the chosen landcover_collection)
# % guisection: Grassland
# %end

# %option
# % key: landcover_collection
# % type: string
# % required: no
# % multiple: no
# % options: esa-worldcover,io-lulc-annual-v02,io-lulc-9-class,io-lulc,nrcan-landcover,drcog-lulc,cop-lc,modis-mcd12q1-061,dynamicworld
# % answer: esa-worldcover
# % description: Land-cover collection r.in.landcover downloads when grassland= is empty
# % guisection: Grassland
# %end

# %option
# % key: landcover_year
# % type: string
# % required: no
# % description: Year of the downloaded land cover (default: the year the simulation starts, or the nearest year the collection covers)
# % guisection: Grassland
# %end

# %option
# % key: climate
# % type: string
# % required: no
# % multiple: no
# % options: era5,ecmwf,strds
# % answer: era5
# % descriptions: era5;ERA5-Land reanalysis through t.in.era5, for a past or current period;ecmwf;IFS/AIFS open-data forecast through t.in.ecmwf, for the coming days;strds;no download, every forcing is given as an existing space-time raster dataset
# % description: Where the climate forcings come from
# % guisection: Climate
# %end

# %option
# % key: climate_prefix
# % type: string
# % required: no
# % answer: graminr_climate
# % description: STRDS name prefix of the downloaded climate data
# % guisection: Climate
# %end

# %option G_OPT_M_DIR
# % key: cache_dir
# % required: no
# % description: Directory to cache downloaded climate files in, reused by later runs
# % guisection: Climate
# %end

# %option
# % key: par_fraction
# % type: double
# % required: no
# % answer: 0.48
# % description: Photosynthetically active radiation as a fraction of global radiation, used to derive rad_par
# % guisection: Climate
# %end

# %option
# % key: et_pot_method
# % type: string
# % required: no
# % multiple: no
# % options: dataset,hargreaves
# % answer: dataset
# % descriptions: dataset;take potential evapotranspiration from the climate dataset (ERA5 potential_evaporation);hargreaves;estimate it from temperature and latitude (Hargreaves and Samani 1985), which is the only option for the ECMWF open-data forecasts
# % description: How potential evapotranspiration is obtained
# % guisection: Climate
# %end

# %option
# % key: ecmwf_dataset
# % type: string
# % required: no
# % answer: ifs_oper
# % options: ifs_oper,aifs_single
# % description: ECMWF open-data forecast dataset used with climate=ecmwf
# % guisection: Climate
# %end

# %option G_OPT_STRDS_INPUT
# % key: tavg
# % required: no
# % description: STRDS of mean air temperature [degree C]
# % guisection: Forcings
# %end

# %option G_OPT_STRDS_INPUT
# % key: tmin
# % required: no
# % description: STRDS of minimum air temperature [degree C]
# % guisection: Forcings
# %end

# %option G_OPT_STRDS_INPUT
# % key: tmax
# % required: no
# % description: STRDS of maximum air temperature [degree C]
# % guisection: Forcings
# %end

# %option G_OPT_STRDS_INPUT
# % key: rain
# % required: no
# % description: STRDS of precipitation [mm/day]
# % guisection: Forcings
# %end

# %option G_OPT_STRDS_INPUT
# % key: rad_glob
# % required: no
# % description: STRDS of global radiation [MJ/m2/day]
# % guisection: Forcings
# %end

# %option G_OPT_STRDS_INPUT
# % key: rad_par
# % required: no
# % description: STRDS of photosynthetically active radiation [MJ/m2/day]; derived from rad_glob through par_fraction= when absent
# % guisection: Forcings
# %end

# %option G_OPT_STRDS_INPUT
# % key: et_pot
# % required: no
# % description: STRDS of potential evapotranspiration [mm/day]
# % guisection: Forcings
# %end

# %option G_OPT_STRDS_INPUT
# % key: vp
# % required: no
# % description: STRDS of actual vapour pressure [kPa]
# % guisection: Forcings
# %end

# %option G_OPT_STRDS_INPUT
# % key: wind
# % required: no
# % description: STRDS of wind speed at 2 m [m/s]
# % guisection: Forcings
# %end

# %option G_OPT_STRDS_INPUT
# % key: lai_rs
# % required: no
# % description: STRDS of remotely sensed leaf area index [m2/m2], overriding the modelled LAI where it has a value
# % guisection: Remote sensing
# %end

# %option G_OPT_STRDS_INPUT
# % key: et_act_rs
# % required: no
# % description: STRDS of remotely sensed actual evapotranspiration [mm/day], overriding the modelled soil evaporation
# % guisection: Remote sensing
# %end

# %option G_OPT_STRDS_INPUT
# % key: t_act_rs
# % required: no
# % description: STRDS of remotely sensed actual transpiration [mm/day], overriding the modelled transpiration
# % guisection: Remote sensing
# %end

# %option
# % key: modis_prefix
# % type: string
# % required: no
# % answer: graminr_modis
# % description: Raster name prefix of the MODIS products downloaded with -r
# % guisection: Remote sensing
# %end

# %option
# % key: gap_fill
# % type: integer
# % required: no
# % answer: 0
# % description: Days a map may be carried forward when no map covers a simulated day; 8 or 16 lets a MODIS composite force every day of its period
# % guisection: Remote sensing
# %end

# %option G_OPT_STRDS_INPUT
# % key: cut_strds
# % required: no
# % description: STRDS of cut/harvest days, one value per pixel per day: a flag for LINGRA, BASGRA and GrasProg, a residual cut height [m] for ModVege and GrasSim
# % guisection: Management
# %end

# %option
# % key: cut_dates
# % type: string
# % required: no
# % multiple: yes
# % description: Cut/harvest calendar shared by every pixel: a day of year (1-366, repeated every year) and/or a date (YYYY-MM-DD). Mutually exclusive with cut_strds=
# % guisection: Management
# %end

# %option
# % key: cut_height
# % type: double
# % required: no
# % description: Residual sward height [m] assigned to the days of cut_dates= for ModVege and GrasSim (default: the parameter set's DEFAULT_CUT_HEIGHT)
# % guisection: Management
# %end

# %option
# % key: fert_dates
# % type: string
# % required: no
# % multiple: yes
# % key_desc: date/amount
# % description: Mineral nitrogen fertilisation of BASGRA and GrasSim, as date/amount pairs in kg N/ha (for example 2024-04-05/60,2024-06-20/40)
# % guisection: Management
# %end

# %option
# % key: fert_org_dates
# % type: string
# % required: no
# % multiple: yes
# % key_desc: date/amount
# % description: Organic nitrogen fertilisation of GrasSim, as date/amount pairs in kg N/ha
# % guisection: Management
# %end

# %option
# % key: ndep_dates
# % type: string
# % required: no
# % multiple: yes
# % key_desc: date/rate
# % description: Background nitrogen deposition of BASGRA, as date/rate pairs in kg N/ha/day, interpolated linearly between the points
# % guisection: Management
# %end

# %option G_OPT_R_INPUT
# % key: elevation
# % required: no
# % description: Elevation raster [m] used by the LINGRA reference evapotranspiration (default: the ELEV parameter, the same everywhere)
# % guisection: Model
# %end

# %option
# % key: chunk
# % type: integer
# % required: no
# % answer: 200000
# % description: Number of grassland pixels held in memory at once; lower it on a memory-tight machine, raise it for speed
# % guisection: Model
# %end

# %flag
# % key: l
# % description: List the output variables of the chosen model and exit
# % suppress_required: yes
# %end

# %flag
# % key: r
# % description: Download the remotely sensed leaf area index (MODIS MOD15A2H, 8-day, 500 m) through r.in.modis
# % guisection: Remote sensing
# %end

# %flag
# % key: g
# % description: Use Google Earth Engine instead of Microsoft Planetary Computer for the land-cover download
# % guisection: Grassland
# %end

# %flag
# % key: k
# % description: Keep the downloaded and intermediate maps instead of leaving only the model outputs
# %end

# %rules
# % exclusive: cut_strds,cut_dates
# % requires: grassland_classes,grassland
# %end

import atexit
import datetime as dt
import os
import subprocess
import sys

import numpy as np

import grass.script as gs

#: Variables written when the user names none, one meaningful set per model.
DEFAULT_VARIABLES = {
    "modvege": ("gvb", "dvb", "grb", "drb", "hb", "lai", "abc", "g"),
    "lingra": (
        "lai",
        "weight_lv_green",
        "weight_harv",
        "weight_harv_annual",
        "tsum",
        "tiller_number",
        "sm",
    ),
    "basgra": (
        "lai",
        "gv_biomass",
        "harvested",
        "dry_matter",
        "tiller_v",
        "nmin",
        "f_digest_dm",
    ),
    "grasprog": ("w", "yield", "lai", "gi", "sw", "tempsum"),
    "grassim": ("gvb", "grb", "dvb", "drb", "hb", "lai", "tsb", "abc"),
}

#: Unit of each output variable, for t.support.
#: Units of the output variables, per model. They cannot be shared: GraMinR's
#: models do not work in the same units, and the same name means different
#: things in different ones. BASGRA follows the Fortran in g C/m2, g N/m2 and
#: g DM/m2, while ModVege, GrasProg and GrasSim work in kg DM/ha and
#: kg N/ha -- so `nmin`, for instance, is g N/m2 in BASGRA and kg N/ha in
#: GrasSim. An output absent from its model's table is dimensionless (an
#: index, a fraction, a status) and gets no unit.
COMMON_UNITS = {
    "lai": "m2/m2",
    "rootd": "m",
    "sward_height": "m",
    "leaf_length": "m",
    "dayl": "fraction of day",
    "parav": "MJ/m2/day",
    "parint": "MJ/m2/day",
}

MODEL_UNITS = {
    "lingra": {
        "tsum": "degree-day",
        "tiller_number": "tillers/m2",
        "weight_harv": "kg DM/ha",
        "weight_harv_annual": "kg DM/ha",
        "weight_lv_green": "kg DM/ha",
        "weight_lv_dead": "kg DM/ha",
        "weight_abg": "kg DM/ha",
        "weight_re": "kg DM/ha",
        "weight_rt": "kg DM/ha",
        "evs": "cm/day",
        "evw": "cm/day",
        "tra": "cm/day",
        "tramx": "cm/day",
        "et0": "cm/day",
        "e0": "cm/day",
        "es0": "cm/day",
        "perc": "cm/day",
        "loss": "cm/day",
        "runoff": "cm/day",
        "sm": "m3/m3",
        "w_balance": "cm",
        "wlow": "cm",
        "ss": "cm",
        "rd": "cm",
        # w is the water in the rooted zone, so that sm = w / rd.
        "w": "cm",
        "slaint": "ha/kg",
    },
    "modvege": {
        "gvb": "kg DM/ha",
        "dvb": "kg DM/ha",
        "grb": "kg DM/ha",
        "drb": "kg DM/ha",
        "hb": "kg DM/ha",
        "ib": "kg DM/ha",
        "abc": "kg DM/ha",
        "tsb": "kg DM/ha",
        "g": "kg DM/ha/day",
        "pgr": "kg DM/ha/day",
        "gva": "degree-day",
        "gra": "degree-day",
        "dva": "degree-day",
        "dra": "degree-day",
        "stp": "degree-day",
    },
    "grasprog": {
        "w": "kg DM/ha",
        "wt": "kg DM/ha/day",
        "yield": "kg DM/ha",
        "tempsum": "degree-day",
        "sw": "mm",
        "eta": "mm/day",
        "etp": "mm/day",
        "irrig": "mm/day",
        "age": "day",
    },
    "grassim": {
        "gvb": "kg DM/ha",
        "grb": "kg DM/ha",
        "dvb": "kg DM/ha",
        "drb": "kg DM/ha",
        "hb": "kg DM/ha",
        "abc": "kg DM/ha",
        "tsb": "kg DM/ha",
        "g": "kg DM/ha/day",
        "pgr": "kg DM/ha/day",
        "gva": "degree-day",
        "gra": "degree-day",
        "dva": "degree-day",
        "dra": "degree-day",
        "stp": "degree-day",
        "nmin": "kg N/ha",
        "norg": "kg N/ha",
    },
    # BASGRA keeps the Fortran's units throughout.
    "basgra": {
        "clv": "g C/m2",
        "clvd": "g C/m2",
        "cres": "g C/m2",
        "crt": "g C/m2",
        "cst": "g C/m2",
        "cstub": "g C/m2",
        "clitt": "g C/m2",
        "csomf": "g C/m2",
        "csoms": "g C/m2",
        "nmin": "g N/m2",
        "nsh": "g N/m2",
        "nrt": "g N/m2",
        "nlitt": "g N/m2",
        "nsomf": "g N/m2",
        "nsoms": "g N/m2",
        "nsource": "g N/m2/day",
        "nsink": "g N/m2/day",
        "nupt": "g N/m2/day",
        "nleaching": "g N/m2/day",
        "nmineralisation": "g N/m2/day",
        "nshmob": "g N/m2/day",
        "nfert_tot": "g N/m2",
        "dry_matter": "g DM/m2",
        "dmlv": "g DM/m2",
        "dmst": "g DM/m2",
        "yield_last": "g DM/m2",
        "wal": "mm",
        "wapl": "mm",
        "waps": "mm",
        "was": "mm",
        "drystor": "mm",
        "wetstor": "mm",
        "evap": "mm/day",
        "tran": "mm/day",
        "sdepth": "m",
        "fdepth": "m",
        "lt50": "degree C",
        "o2": "mol/m2",
        "tanaer": "day",
        "tiller_v": "tillers/m2",
        "tiller_g1": "tillers/m2",
        "tiller_g2": "tillers/m2",
        "tilg": "tillers/m2",
        "tiltot": "tillers/m2",
        "sla_new": "m2/g",
        "lerg": "m/day",
        "luemxq": "mol CO2/mol PAR",
        "rleaf": "1/day",
        "rgrtv": "1/day",
        "rdrt": "1/day",
        "rdrs": "1/day",
        "tv2til": "1/day",
        "gv_biomass": "g C/m2",
        "dv_biomass": "g C/m2",
        "harvested": "g DM/m2",
    },
}


def unit_of(model, variable):
    """Unit of `variable` in `model`, or an empty string when it has none."""
    return MODEL_UNITS.get(model, {}).get(variable, COMMON_UNITS.get(variable, ""))


#: Variables that are a running total and so must not be aggregated as a mean.
CUMULATIVE = {
    "hb",
    "ib",
    "harvested",
    "weight_harv",
    "weight_harv_annual",
    "nfert_tot",
    "yield",
}

#: kg/ha in one g/m2. BASGRA works in g N/m2 and g C/m2, the other models
#: and every user-facing nitrogen option in kg N/ha.
KG_HA_PER_G_M2 = 10.0

CLEANUP = []


def cleanup():
    """Remove the intermediate maps, unless -k was given.

    This runs from atexit, including after a fatal error, so it must not
    raise: a failure to remove an intermediate map would otherwise mask
    whatever actually went wrong.
    """
    for kind, name in CLEANUP:
        try:
            gs.run_command(
                "g.remove",
                type=kind,
                name=name,
                flags="f",
                quiet=True,
                stderr=subprocess.DEVNULL,
            )
        except Exception:
            pass


def parse_date(text, label):
    try:
        return dt.datetime.strptime(text.strip(), "%Y-%m-%d").date()
    except (ValueError, AttributeError):
        gs.fatal(
            _("{label}= must be a date, YYYY-MM-DD, not '{text}'.").format(
                label=label, text=text
            )
        )


#: The forcing options a user may point straight at an STRDS. A forcing
#: given this way is used as it stands, in the units the manual documents.
FORCING_OPTIONS = (
    "tavg",
    "tmin",
    "tmax",
    "rain",
    "rad_glob",
    "rad_par",
    "et_pot",
    "vp",
    "wind",
    "lai_rs",
    "et_act_rs",
    "t_act_rs",
)

#: Forcings this module derives from another one, so that a set supplying
#: only the source is complete: rad_par comes out of rad_glob through
#: par_fraction= (see build_forcings).
DERIVABLE_FROM = {"rad_par": "rad_glob"}


def preflight(options, flags, model, module):
    """Reject an unusable request before anything is downloaded or read.

    Everything here is local and costs nothing: the output variable names,
    the completeness of a supplied forcing set, and the existence of the
    maps named on the command line. It runs before build_mask(), because
    that one may download a land-cover product -- and a typo in
    `variables=` or a forgotten forcing should not cost the user a download
    first. The return value is the resolved variable list, so it is not
    computed twice.
    """
    variables = resolve_variables(
        model, [value for value in options["variables"].split(",") if value], module
    )

    for option in ("grassland", "elevation"):
        name = options.get(option)
        if name and not gs.find_file(name, element="cell")["file"]:
            gs.fatal(
                _("Raster map <{name}> given as {option}= was not found.").format(
                    name=name, option=option
                )
            )

    # parameter_maps= is a list of name=raster pairs, so both halves are
    # checked: a malformed pair and a missing raster are equally fatal.
    for pair in options["parameter_maps"].split(","):
        if not pair:
            continue
        if "=" not in pair:
            gs.fatal(
                _(
                    "parameter_maps= takes name=raster pairs; '{pair}' is not one."
                ).format(pair=pair)
            )
        parameter, _sep, raster = pair.partition("=")
        if not parameter.strip() or not raster.strip():
            gs.fatal(
                _(
                    "parameter_maps= takes name=raster pairs; '{pair}' leaves "
                    "one side empty."
                ).format(pair=pair)
            )
        if not gs.find_file(raster.strip(), element="cell")["file"]:
            gs.fatal(
                _(
                    "Raster map <{raster}>, given in parameter_maps= for the "
                    "'{parameter}' parameter, was not found."
                ).format(raster=raster.strip(), parameter=parameter.strip())
            )

    if options["climate"] == "strds":
        # climate=strds promises that every forcing is already in the
        # mapset, so an incomplete set is an error now rather than a
        # failure after the mask has been built.
        supplied = {name for name in FORCING_OPTIONS if options.get(name)}
        missing = [
            name
            for name in module.REQUIRED
            if name not in supplied and DERIVABLE_FROM.get(name) not in supplied
        ]
        if options["et_pot_method"] == "hargreaves":
            missing = [name for name in missing if name != "et_pot"]
            missing += [name for name in ("tmin", "tmax") if name not in supplied]
        if missing:
            gs.fatal(
                _(
                    "climate=strds needs every forcing of the '{model}' model "
                    "to be given, and {names} {verb} missing. Either supply "
                    "{them} as a space-time raster dataset, or choose "
                    "climate=era5 or climate=ecmwf to have them downloaded."
                ).format(
                    model=model,
                    names=", ".join(sorted(set(missing))),
                    verb="is" if len(set(missing)) == 1 else "are",
                    them="it" if len(set(missing)) == 1 else "them",
                )
            )

    return variables


def resolve_variables(model, requested, module):
    """Output variables to write, validated against the model's own list."""
    available = tuple(module.OUTPUTS)
    if not requested:
        return tuple(
            name
            for name in DEFAULT_VARIABLES.get(model, available[:6])
            if name in available
        )
    if len(requested) == 1 and requested[0] == "all":
        return available
    unknown = [name for name in requested if name not in available]
    if unknown:
        gs.fatal(
            _(
                "The '{model}' model has no output variable {names}. Its "
                "variables are: {available}"
            ).format(
                model=model, names=", ".join(unknown), available=", ".join(available)
            )
        )
    return tuple(requested)


def build_mask(options, flags, start):
    """Flat indices of the grassland pixels, and the mask raster used.

    <output>_grassland and <output>_landcover are this module's own
    intermediates: it derives their names from output=, writes them itself
    and removes them again unless -k was given. They are always rewritten,
    rather than honouring --overwrite like a requested output, because a
    run that died before its cleanup would otherwise leave a map behind
    that blocks every later run with the same output= prefix.
    """
    from graminr_lib import fetch
    from graminr_lib.gio import read_raster, region_shape

    rows, cols = region_shape()
    given = options["grassland"]
    classes = [value for value in options["grassland_classes"].split(",") if value]

    if given:
        if classes:
            mask_map = fetch.grassland_mask(
                given,
                options["landcover_collection"],
                classes,
                output="{}_grassland".format(options["output"]),
                overwrite=True,
            )
            if not flags["k"]:
                CLEANUP.append(("raster", mask_map))
        else:
            mask_map = given
            gs.message(
                _(
                    "Grassland is taken from <{}>: every non-null, non-zero cell "
                    "counts as grassland."
                ).format(mask_map)
            )
    else:
        year = options["landcover_year"] or str(start.year)
        prefix = "{}_landcover".format(options["output"])
        landcover = fetch.download_landcover(
            options["landcover_collection"],
            year,
            year,
            prefix,
            use_gee=flags["g"],
            overwrite=True,
        )
        if not flags["k"]:
            CLEANUP.append(("raster", landcover))
        mask_map = fetch.grassland_mask(
            landcover,
            options["landcover_collection"],
            classes or None,
            output="{}_grassland".format(options["output"]),
            overwrite=True,
        )
        if not flags["k"]:
            CLEANUP.append(("raster", mask_map))

    data = read_raster(mask_map).ravel()
    index = np.nonzero(np.isfinite(data) & (data != 0))[0]
    if index.size == 0:
        gs.fatal(
            _(
                "No grassland pixel in the current region. Check <{mask}>, the "
                "grassland classes, and the region itself."
            ).format(mask=mask_map)
        )
    gs.message(
        _(
            "{count} grassland pixel(s) out of {total} in the region ({percent:.1f}%)."
        ).format(
            count=index.size,
            total=rows * cols,
            percent=100.0 * index.size / (rows * cols),
        )
    )
    return index, mask_map, rows, cols


def collect_sources(options, flags, model_module, start, end, mask_index):
    """STRDS available for each forcing, downloading what is missing."""
    from graminr_lib import fetch

    needed = list(model_module.REQUIRED)
    optional = dict(model_module.OPTIONAL)
    sources = {}

    # Forcings the user supplied directly always win.
    for name in FORCING_OPTIONS:
        if options.get(name):
            sources[name] = options[name]

    kind = options["climate"]

    missing = [
        name
        for name in needed
        if name not in sources and sources.get(DERIVABLE_FROM.get(name)) is None
    ]
    if options["et_pot_method"] == "hargreaves" and "et_pot" in missing:
        # Hargreaves needs tmin and tmax, even for a model that does not.
        for name in ("tmin", "tmax"):
            if name not in needed and name not in sources:
                missing.append(name)
        missing = [name for name in missing if name != "et_pot"]

    if missing and kind == "era5":
        downloaded = fetch.download_era5(
            missing,
            start,
            end,
            options["climate_prefix"],
            cache_dir=options["cache_dir"] or None,
            overwrite=gs.overwrite(),
        )
        sources.update(downloaded)
        if not flags["k"]:
            for strds in downloaded.values():
                CLEANUP.append(("strds", strds))
    elif missing and kind == "ecmwf":
        if options["et_pot_method"] != "hargreaves" and "et_pot" in needed:
            gs.fatal(
                _(
                    "The ECMWF open-data forecasts publish no potential "
                    "evapotranspiration, which the '{model}' model needs. Run "
                    "with et_pot_method=hargreaves, or supply et_pot= "
                    "yourself."
                ).format(model=options["model"])
            )
        downloaded = fetch.download_ecmwf(
            missing,
            options["climate_prefix"],
            dataset=options["ecmwf_dataset"],
            cache_dir=options["cache_dir"] or None,
            overwrite=gs.overwrite(),
        )
        sources.update(downloaded)
        if not flags["k"]:
            for strds in downloaded.values():
                CLEANUP.append(("strds", strds))
    elif missing:
        gs.fatal(
            _(
                "With climate=strds every forcing must be given; missing: {names}"
            ).format(names=", ".join(sorted(set(missing))))
        )

    # Optional remotely sensed forcings.
    if flags["r"]:
        from graminr_lib.gio import register_strds

        # Only what r.in.modis can actually deliver: see
        # fetch.MODIS_PRODUCTS and fetch.MODIS_UNAVAILABLE.
        wanted = [
            name
            for name in fetch.MODIS_PRODUCTS
            if name in optional and name not in sources
        ]
        skipped = [
            name
            for name in fetch.MODIS_UNAVAILABLE
            if name in optional and name not in sources
        ]
        if skipped:
            gs.warning(
                _(
                    "-r cannot download {names}: no MOD16 evapotranspiration "
                    "collection is published per 8-day composite, only per "
                    "year. The {model} model will compute them itself; supply "
                    "them as an STRDS to override that."
                ).format(names=", ".join(sorted(skipped)), model=model)
            )
        for name in wanted:
            maps, product = fetch.download_modis(
                name,
                start,
                end,
                options["modis_prefix"],
                overwrite=gs.overwrite(),
            )
            scaled = []
            for raster in maps:
                plain = raster.split("@")[0]
                target = "{}_scaled".format(plain)
                gs.mapcalc(
                    "{out} = if({src} > {vmax}, null(), {src} * {scale})".format(
                        out=target,
                        src=raster,
                        vmax=product["valid_max"],
                        scale=product["scale"],
                    ),
                    overwrite=True,
                    quiet=True,
                )
                scaled.append(target)
                if not flags["k"]:
                    CLEANUP.append(("raster", raster))
                    CLEANUP.append(("raster", target))
            dates = [_modis_date(raster) for raster in maps]
            strds = "{}_{}".format(options["modis_prefix"], name)
            register_strds(
                strds,
                list(zip(scaled, dates)),
                title="MODIS {}".format(product["collection"]),
                description="t.graminr forcing '{}' from {}".format(
                    name, product["collection"]
                ),
                unit=product["unit"],
                overwrite=True,
                temporary_dir=gs.tempdir(),
            )
            sources[name] = strds
            if not flags["k"]:
                CLEANUP.append(("strds", strds))
            if int(options["gap_fill"]) < product["gap"]:
                gs.message(
                    _(
                        "The {collection} composites step every {gap} days; "
                        "gap_fill={gap} would let each one force all of its "
                        "days."
                    ).format(collection=product["collection"], gap=product["gap"])
                )

    return sources, needed, optional


def _modis_date(raster):
    """Date encoded in an r.in.modis map name, <prefix>_<YYYYMMDD>_<band>."""
    for part in raster.split("@")[0].split("_"):
        if len(part) == 8 and part.isdigit():
            return dt.datetime.strptime(part, "%Y%m%d").date()
    gs.fatal(_("Cannot read a date out of the map name <{}>.").format(raster))


def main():
    options, flags = gs.parser()

    # graminr_lib/ and params/ sit beside this script when it is run from
    # the source tree, and in <base>/etc/t.graminr/ once installed. Which
    # base that is depends on how the module was installed: an addon built
    # into the GRASS source tree lands under GISBASE, one installed with
    # g.extension under GRASS_ADDON_BASE, and the script itself always ends
    # up in <base>/scripts/, so ../etc/t.graminr covers both without
    # consulting the environment at all.
    here = os.path.dirname(os.path.abspath(__file__))
    candidates = [here, os.path.join(os.path.dirname(here), "etc", "t.graminr")]
    for base in (os.environ.get("GRASS_ADDON_BASE"), os.environ.get("GISBASE")):
        if base:
            candidates.append(os.path.join(base, "etc", "t.graminr"))
    for candidate in candidates:
        if candidate and candidate not in sys.path and os.path.isdir(candidate):
            sys.path.insert(0, candidate)

    from graminr_lib import cutdates, fetch, params
    from graminr_lib.gio import (
        latitude_longitude_arrays,
        register_strds,
        set_raster_metadata,
        simulated_days,
        write_raster,
    )

    model = options["model"]
    module = params.MODELS[model]

    if flags["l"]:
        gs.message(_("Output variables of the '{model}' model:").format(model=model))
        for name in module.OUTPUTS:
            unit = unit_of(model, name)
            sys.stdout.write(
                "{}{}\n".format(name, " [{}]".format(unit) if unit else "")
            )
        return 0

    start = parse_date(options["start"], "start")
    end = parse_date(options["end"], "end")
    days = simulated_days(start, end)
    ndays = len(days)
    year = np.array([day.year for day in days])
    doy = np.array([day.timetuple().tm_yday for day in days])
    gs.message(
        _(
            "Simulating {days} day(s), {start} to {end}, with the '{model}' model."
        ).format(days=ndays, start=start, end=end, model=model)
    )

    # Everything local and free is checked first, so that a typo costs
    # nothing: build_mask() below may download a land-cover product.
    variables = preflight(options, flags, model, module)

    atexit.register(cleanup)

    mask_index, mask_map, rows, cols = build_mask(options, flags, start)
    npixels = int(mask_index.size)

    gs.message(
        _("Writing {count} output variable(s): {names}").format(
            count=len(variables), names=", ".join(variables)
        )
    )

    parameters = params.load(
        model,
        path=options["parameters"] or None,
        parameter_maps=[v for v in options["parameter_maps"].split(",") if v],
        mask_index=mask_index,
        overrides=[v for v in options["parameter_values"].split(",") if v],
    )

    # Latitude and elevation per pixel: LINGRA's and BASGRA's solar geometry
    # and reference evapotranspiration both depend on them, so a region
    # spanning several degrees is handled correctly rather than being pinned
    # to the single site of the parameter file.
    needs_latitude = model in ("lingra", "basgra")
    latitude = None
    if needs_latitude:
        lat_grid, _lon_grid = latitude_longitude_arrays()
        latitude = lat_grid.ravel()[mask_index]
        parameters["LAT" if model == "lingra" else "lat"] = latitude
        gs.verbose(
            _(
                "Latitude over the grassland pixels: {low:.3f} to {high:.3f} degree."
            ).format(low=float(np.min(latitude)), high=float(np.max(latitude)))
        )
        if model == "lingra" and options["elevation"]:
            from graminr_lib.gio import read_raster

            elevation = read_raster(options["elevation"]).ravel()[mask_index]
            parameters["ELEV"] = np.where(
                np.isfinite(elevation), elevation, float(parameters["ELEV"])
            )

    sources, needed, optional = collect_sources(
        options, flags, module, start, end, mask_index
    )
    if options["et_pot_method"] == "hargreaves":
        sources.pop("et_pot", None)

    forcing_names = list(needed)
    for name in optional:
        if name in sources:
            forcing_names.append(name)
    defaults = {
        name: value for name, value in optional.items() if not isinstance(value, str)
    }

    forcings = fetch.build_forcings(
        forcing_names,
        sources,
        mask_index,
        par_fraction=float(options["par_fraction"]),
        gap_fill=int(options["gap_fill"]),
        source_kind=options["climate"],
        defaults=defaults,
    )

    if options["et_pot_method"] == "hargreaves" and "et_pot" in needed:
        if latitude is None:
            lat_grid, _lon_grid = latitude_longitude_arrays()
            latitude = lat_grid.ravel()[mask_index]
        helper = fetch.build_forcings(
            ["tmin", "tmax"],
            sources,
            mask_index,
            par_fraction=float(options["par_fraction"]),
            gap_fill=int(options["gap_fill"]),
            source_kind=options["climate"],
            defaults=defaults,
        )
        forcings["et_pot"] = fetch.HargreavesET0(
            helper["tmin"], helper["tmax"], latitude
        )
        gs.message(
            _(
                "Potential evapotranspiration is estimated with Hargreaves and "
                "Samani (1985), from temperature and latitude."
            )
        )

    # Cut and fertilisation calendars.
    cut_value = 1.0
    if model in ("modvege", "grassim"):
        cut_value = float(
            options["cut_height"] or parameters.get("DEFAULT_CUT_HEIGHT", 0.05)
        )
    cuts = cutdates.CutCalendar(
        days,
        mask_index,
        strds=options["cut_strds"] or None,
        cut_dates=options["cut_dates"] or None,
        cut_value=cut_value,
        gap_fill=int(options["gap_fill"]),
    )
    fert = cutdates.EventCalendar(days, options["fert_dates"], "fert_dates")
    fert_org = cutdates.EventCalendar(days, options["fert_org_dates"], "fert_org_dates")
    ndep = cutdates.EventCalendar(
        days, options["ndep_dates"], "ndep_dates", interpolate=True
    )

    # --- Read the forcings once, into a per-pixel, per-day store ---
    # The models are independent from pixel to pixel, so the run is chunked
    # over pixels; reading each map once into a memory-mapped store keeps a
    # multi-year run from re-reading every raster for every chunk.
    workdir = gs.tempdir()
    store = {}
    gs.message(_("Reading the forcings over {} day(s)...").format(ndays))
    for name, forcing in forcings.items():
        path = os.path.join(workdir, "forcing_{}.dat".format(name))
        store[name] = np.memmap(
            path, dtype=np.float32, mode="w+", shape=(ndays, npixels)
        )
    for position, day in enumerate(days):
        for name, forcing in forcings.items():
            store[name][position] = forcing.values(day)
        gs.percent(position, ndays, 1)
    gs.percent(1, 1, 1)

    # Management arrays, per day and per pixel.
    gcut = np.memmap(
        os.path.join(workdir, "gcut.dat"),
        dtype=np.float32,
        mode="w+",
        shape=(ndays, npixels),
    )
    for position, day in enumerate(days):
        gcut[position] = cuts.values(position, day)
    store["gcut"] = gcut

    # Management amounts are the same on every pixel, so they are kept as one
    # value per day and broadcast inside the model loop; materialising them
    # over every pixel would cost as much memory as a whole forcing.
    uniform = {}
    if model == "basgra":
        # fert_dates= and ndep_dates= are in kg N/ha and kg N/ha/day, which
        # is how the amounts are quoted agronomically and what GrasSim uses.
        # BASGRA, following the Fortran, works in g N/m2 throughout, so the
        # amounts are converted here: 1 g/m2 = 10 kg/ha. Without this a
        # 60 kg N/ha dressing would be applied as 600.
        uniform["nfert"] = fert.amounts / KG_HA_PER_G_M2
        uniform["ndep"] = ndep.amounts / KG_HA_PER_G_M2
    elif model == "grassim":
        # GrasSim already works in kg N/ha.
        uniform["fert_min"] = fert.amounts
        uniform["fert_org"] = fert_org.amounts

    # --- Run the model, chunk of pixels by chunk of pixels ---
    chunk_size = max(1, int(options["chunk"]))
    results = {}
    for name in variables:
        path = os.path.join(workdir, "out_{}.dat".format(name))
        results[name] = np.memmap(
            path, dtype=np.float32, mode="w+", shape=(ndays, npixels)
        )

    nchunks = (npixels + chunk_size - 1) // chunk_size
    gs.message(
        _("Running '{model}' over {pixels} pixel(s) in {chunks} chunk(s)...").format(
            model=model, pixels=npixels, chunks=nchunks
        )
    )
    done = [0]
    total_steps = nchunks * ndays

    for chunk in range(nchunks):
        low = chunk * chunk_size
        high = min(npixels, low + chunk_size)
        weather = {"year": year, "doy": doy}
        for name, data in store.items():
            weather[name] = np.asarray(data[:, low:high], dtype=np.float64)
        for name, amounts in uniform.items():
            weather[name] = amounts.astype(np.float64)[:, None]

        def progress(_day_index, _ndays, _base=done[0]):
            gs.percent(_base + _day_index, total_steps, 1)

        chunk_out = module.run(
            parameters
            if high - low == npixels
            else _slice_parameters(parameters, low, high),
            weather,
            high - low,
            outputs=variables,
            progress=progress,
        )
        for name in variables:
            results[name][:, low:high] = chunk_out[name]
        done[0] += ndays
    gs.percent(1, 1, 1)

    # --- Write one STRDS per output variable ---
    # r.support records this in each output raster's history, as
    # gs.raster_history() would, but without its four extra subprocesses.
    command_line = os.environ.get("CMDLINE", "t.graminr")
    gs.message(_("Writing {} space-time raster dataset(s)...").format(len(variables)))
    canvas = np.full(rows * cols, np.nan, dtype=np.float32)
    for name in variables:
        registration = []
        for position, day in enumerate(days):
            raster = "{}_{}_{}".format(options["output"], name, day.strftime("%Y%m%d"))
            canvas[:] = np.nan
            canvas[mask_index] = results[name][position]
            write_raster(
                canvas,
                raster,
                rows,
                cols,
                title="{} {} {}".format(model, name, day),
                overwrite=True,
            )
            set_raster_metadata(raster, unit=unit_of(model, name), history=command_line)
            registration.append((raster, day))
            gs.percent(position, ndays, 2)
        strds = "{}_{}".format(options["output"], name)
        register_strds(
            strds,
            registration,
            title="{} {}".format(model, name),
            # No apostrophes here: GRASS's temporal framework interpolates
            # the title and description straight into SQL, so a single
            # quote makes t.create fail with an sqlite3 syntax error.
            description=(
                "t.graminr: {name} of the {model} grassland model, "
                "{start} to {end}".format(name=name, model=model, start=start, end=end)
            ),
            unit=unit_of(model, name),
            semantic="sum" if name in CUMULATIVE else "mean",
            overwrite=gs.overwrite(),
            temporary_dir=workdir,
        )
        gs.message(_("<{strds}> written.").format(strds=strds))

    gs.message(
        _("Done. {count} dataset(s) with the prefix <{prefix}_>.").format(
            count=len(variables), prefix=options["output"]
        )
    )
    return 0


def _slice_parameters(parameters, low, high):
    """Restrict the per-pixel parameters to one chunk of pixels."""
    sliced = {}
    for name, value in parameters.items():
        if isinstance(value, np.ndarray) and value.ndim == 1 and value.size > 1:
            sliced[name] = value[low:high]
        else:
            sliced[name] = value
    return sliced


if __name__ == "__main__":
    sys.exit(main())
