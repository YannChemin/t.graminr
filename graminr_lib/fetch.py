"""
MODULE:     t.graminr (graminr_lib.fetch)

AUTHOR(S):  Yann Chemin

PURPOSE:    Acquisition of every dataset t.graminr needs, so that a run can
            start from nothing but a region and a date range:

              * climate   -- t.in.era5 (ERA5-Land reanalysis, the default)
                             or t.in.ecmwf (IFS/AIFS open-data forecasts)
              * land cover -- r.in.landcover, reduced to a grassland mask
              * remote sensing -- r.in.modis (MOD15A2H LAI, MOD16 ET)

            Nothing is derived into new maps: the unit conversions and the
            composite variables (mean temperature, PAR, vapour pressure,
            wind speed) are applied to the pixel values as each day is read,
            which keeps a multi-year run from littering the mapset with tens
            of thousands of intermediate rasters.

SPDX-License-Identifier: Unlicense
"""

import numpy as np

import grass.script as gs

from .gio import Forcing

#: ERA5 variables of t.in.era5 needed by each logical forcing.
ERA5_SOURCES = {
    "tavg": ("temperature",),
    "tmin": ("temperature_min",),
    "tmax": ("temperature_max",),
    "rain": ("precipitation",),
    "rad_glob": ("solar_radiation",),
    "rad_par": ("solar_radiation",),
    "et_pot": ("potential_evaporation",),
    "vp": ("dewpoint_temperature",),
    "wind": ("wind_u", "wind_v"),
}

#: ECMWF open-data parameters of t.in.ecmwf needed by each logical forcing.
ECMWF_SOURCES = {
    "tavg": ("2t",),
    "tmin": ("2t",),
    "tmax": ("2t",),
    "rain": ("tp",),
    "rad_glob": ("ssrd",),
    "rad_par": ("ssrd",),
    "vp": ("2d",),
    "wind": ("10u", "10v"),
}

#: Grassland classes of each r.in.landcover collection.
GRASSLAND_CLASSES = {
    "esa-worldcover": (30,),
    "io-lulc-annual-v02": (11,),
    "io-lulc-9-class": (11,),
    "io-lulc": (9,),
    "nrcan-landcover": (9, 10, 12),
    "drcog-lulc": (2,),
    "cop-lc": (30,),
    "modis-mcd12q1-061": (10,),
    "dynamicworld": (2,),
}

#: Band holding the class code, per r.in.landcover collection.
LANDCOVER_BANDS = {
    "esa-worldcover": "map",
    "cop-lc": "discrete_classification",
    "modis-mcd12q1-061": "LC_Type1",
    "dynamicworld": "label",
}

#: MODIS collections and bands for the optional remote-sensing forcings.
#: The MODIS products r.in.modis can supply for a t.graminr forcing.
#:
#: Only leaf area index is here. The obvious companion, MOD16
#: evapotranspiration, is not available at a useful step: Planetary
#: Computer carries modis-16A3GF-061, which is *yearly* gap-filled ET, and
#: not the 8-day modis-16A2GF-061. One value per year cannot force a daily
#: model, and spreading it evenly over 365 days would quietly replace a
#: measurement with a constant, so et_act_rs and t_act_rs are not
#: downloaded at all; they have to be supplied as an STRDS.
MODIS_PRODUCTS = {
    "lai_rs": {
        "collection": "modis-15A2H-061",
        "band": "Lai_500m",
        "scale": 0.1,
        "valid_max": 100.0,
        "gap": 8,
        "unit": "m2/m2",
    },
}

#: Forcings a user may reasonably expect -r to fetch, and why it cannot.
MODIS_UNAVAILABLE = {
    "et_act_rs": "actual evapotranspiration",
    "et_pot_rs": "potential evapotranspiration",
    "t_act_rs": "actual transpiration",
}


def kelvin_to_celsius(values):
    return values - 273.15


def dewpoint_to_vapour_pressure(dewpoint_c):
    """Actual vapour pressure [kPa] from the dewpoint [degree C] (FAO-56)."""
    return 0.6108 * np.exp((17.27 * dewpoint_c) / (237.3 + dewpoint_c))


def wind_10m_to_2m(speed_10m):
    """FAO-56 logarithmic reduction of a 10 m wind speed to 2 m."""
    return speed_10m * 4.87 / np.log(67.8 * 10.0 - 5.42)


def wind_speed(u_component, v_component):
    return np.hypot(u_component, v_component)


def era5_variables_for(forcings):
    """ERA5 variables of t.in.era5 needed to build `forcings`."""
    needed = []
    for name in forcings:
        for variable in ERA5_SOURCES.get(name, ()):
            if variable not in needed:
                needed.append(variable)
    return needed


def ecmwf_parameters_for(forcings):
    """ECMWF open-data parameters needed to build `forcings`."""
    needed = []
    for name in forcings:
        for parameter in ECMWF_SOURCES.get(name, ()):
            if parameter not in needed:
                needed.append(parameter)
    return needed


def download_era5(forcings, start, end, prefix, cache_dir=None, overwrite=False):
    """Import the ERA5 variables needed by `forcings` as one STRDS each."""
    variables = era5_variables_for(forcings)
    if not variables:
        return {}
    gs.message(
        _(
            "Downloading {count} ERA5-Land variable(s) through t.in.era5: {names}"
        ).format(count=len(variables), names=", ".join(variables))
    )
    kwargs = dict(
        variables=",".join(variables),
        start=start.strftime("%Y-%m-%d"),
        end=end.strftime("%Y-%m-%d"),
        output_prefix=prefix,
        overwrite=overwrite,
    )
    if cache_dir:
        kwargs["cache_dir"] = cache_dir
    gs.run_command("t.in.era5", **kwargs)
    return {variable: "{}_{}".format(prefix, variable) for variable in variables}


def download_ecmwf(
    forcings,
    prefix,
    dataset="ifs_oper",
    steps="0-360",
    date=None,
    time=None,
    cache_dir=None,
    overwrite=False,
):
    """Import the ECMWF open-data forecast parameters, aggregated to daily.

    t.in.ecmwf delivers sub-daily forecast steps in native ECMWF units, so
    each parameter is aggregated to one map per day here: a sum for the
    accumulated fields and a mean, minimum or maximum for the rest.
    """
    parameters = ecmwf_parameters_for(forcings)
    if not parameters:
        return {}
    gs.message(
        _(
            "Downloading {count} ECMWF open-data parameter(s) through t.in.ecmwf: "
            "{names}"
        ).format(count=len(parameters), names=", ".join(parameters))
    )
    kwargs = dict(
        datasets=dataset,
        parameters=",".join(parameters),
        steps=steps,
        output_prefix=prefix,
        overwrite=overwrite,
    )
    if date:
        kwargs["date"] = date
    if time:
        kwargs["time"] = time
    if cache_dir:
        kwargs["cache_dir"] = cache_dir
    # -d de-accumulates the accumulated fields into per-interval amounts,
    # which is what a daily sum must be built from.
    gs.run_command("t.in.ecmwf", flags="d", **kwargs)

    sources = {}
    accumulated = {"tp", "ssrd"}
    for parameter in parameters:
        raw = "{}_{}_{}".format(prefix, dataset, parameter)
        methods = {"sum" if parameter in accumulated else "average"}
        if "tmin" in forcings and parameter == "2t":
            methods.add("minimum")
        if "tmax" in forcings and parameter == "2t":
            methods.add("maximum")
        for method in methods:
            daily = "{}_daily_{}".format(raw, method)
            gs.run_command(
                "t.rast.aggregate",
                input=raw,
                output=daily,
                basename="{}_{}_{}_d".format(prefix, parameter, method),
                granularity="1 day",
                method=method,
                overwrite=overwrite,
                quiet=True,
            )
            sources["{}:{}".format(parameter, method)] = daily
    return sources


def download_landcover(
    collection, start, end, prefix, resolution=0, use_gee=False, overwrite=False
):
    """Import a land-cover map and return the raster holding the classes."""
    gs.message(
        _("Downloading land cover '{collection}' through r.in.landcover").format(
            collection=collection
        )
    )
    kwargs = dict(
        collection=collection,
        start=start,
        end=end,
        output=prefix,
        overwrite=overwrite,
    )
    if resolution:
        kwargs["resolution"] = int(resolution)
    band = LANDCOVER_BANDS.get(collection)
    if band:
        kwargs["bands"] = band
    gs.run_command("r.in.landcover", flags="g" if use_gee else "", **kwargs)

    candidates = gs.list_strings("raster", pattern="{}*".format(prefix))
    if not candidates:
        gs.fatal(
            _("r.in.landcover produced no raster with the prefix <{}>.").format(prefix)
        )
    # Prefer the class band when several were imported.
    if band:
        for name in candidates:
            if name.split("@")[0].endswith(band):
                return name
    return sorted(candidates)[0]


def grassland_mask(
    landcover_map, collection, classes=None, output="graminr_grassland", overwrite=False
):
    """Reduce a land-cover map to a 1/NULL grassland mask."""
    if classes:
        codes = [int(value) for value in classes]
    else:
        codes = list(GRASSLAND_CLASSES.get(collection, ()))
        if not codes:
            gs.fatal(
                _(
                    "No grassland class is known for the collection "
                    "'{collection}'; name the class code(s) with "
                    "grassland_classes=."
                ).format(collection=collection)
            )
    gs.message(
        _("Grassland classes of '{collection}': {codes}").format(
            collection=collection, codes=", ".join(str(code) for code in codes)
        )
    )
    test = " || ".join("{} == {}".format(landcover_map, code) for code in codes)
    gs.mapcalc(
        "{out} = if({test}, 1, null())".format(out=output, test=test),
        overwrite=overwrite,
        quiet=True,
    )
    return output


def download_modis(variable, start, end, prefix, resolution=500, overwrite=False):
    """Import a MODIS product and register it as a daily-resolved STRDS.

    The composites keep their native 8-day step; the day index of
    graminr_lib.gio carries each one forward over its period, so no
    resampling into daily maps is needed.
    """
    if variable in MODIS_UNAVAILABLE:
        gs.fatal(
            _(
                "No MODIS product can supply '{variable}' ({what}). The only "
                "MOD16 evapotranspiration collection on Planetary Computer is "
                "modis-16A3GF-061, which is gap-filled per year, not per 8-day "
                "composite, so it cannot force a daily model. Supply "
                "{variable}= as an STRDS instead, or let the model compute it."
            ).format(variable=variable, what=MODIS_UNAVAILABLE[variable])
        )
    if variable not in MODIS_PRODUCTS:
        gs.fatal(_("No MODIS product is configured for '{}'.").format(variable))
    product = MODIS_PRODUCTS[variable]
    gs.message(
        _("Downloading MODIS {collection} ({band}) through r.in.modis").format(
            collection=product["collection"], band=product["band"]
        )
    )
    gs.run_command(
        "r.in.modis",
        collection=product["collection"],
        bands=product["band"],
        start=start.strftime("%Y-%m-%d"),
        end=end.strftime("%Y-%m-%d"),
        resolution=int(resolution),
        output=prefix,
        overwrite=overwrite,
    )
    maps = gs.list_strings("raster", pattern="{}_*_{}".format(prefix, product["band"]))
    if not maps:
        gs.fatal(
            _(
                "r.in.modis imported no map for {collection}/{band} over this "
                "region and period."
            ).format(collection=product["collection"], band=product["band"])
        )
    return sorted(maps), product


def build_forcings(
    names,
    sources,
    mask_index,
    par_fraction=0.48,
    gap_fill=0,
    source_kind="era5",
    defaults=None,
):
    """Build the Forcing objects the models read, one per logical variable.

    :param names: logical forcing names wanted by the model
    :param sources: mapping built by download_era5()/download_ecmwf(), plus
        any STRDS the user supplied directly under its logical name
    :param mask_index: flat indices of the grassland pixels
    :param par_fraction: PAR as a fraction of global radiation
    :param gap_fill: days a map may be carried forward
    :param source_kind: 'era5', 'ecmwf' or 'strds'
    :param defaults: fallback constant per forcing when no dataset has it
    """
    defaults = defaults or {}
    forcings = {}

    if source_kind == "strds":
        # The Forcings options are documented in the same units t.in.era5
        # delivers, so an externally supplied set goes through the same
        # derivations: rad_par out of rad_glob, for one.
        sources = dict(sources)
        for logical, era5_name in (
            ("rad_glob", "solar_radiation"),
            ("tavg", "temperature"),
            ("tmin", "temperature_min"),
            ("tmax", "temperature_max"),
            ("rain", "precipitation"),
            ("et_pot", "potential_evaporation"),
        ):
            if logical in sources:
                sources.setdefault(era5_name, sources[logical])
        source_kind = "era5"

    def pick(*keys):
        for key in keys:
            if key in sources:
                return sources[key]
        return None

    for name in names:
        # A directly supplied STRDS always wins over a downloaded one.
        direct = sources.get(name)
        if direct and name not in ERA5_SOURCES and name not in ECMWF_SOURCES:
            forcings[name] = Forcing(
                name,
                [direct],
                mask_index,
                max_gap=gap_fill,
                fill=defaults.get(name),
            )
            continue

        if source_kind == "era5":
            forcing = _era5_forcing(
                name, pick, mask_index, par_fraction, gap_fill, defaults, direct
            )
        elif source_kind == "ecmwf":
            forcing = _ecmwf_forcing(
                name, pick, mask_index, par_fraction, gap_fill, defaults, direct
            )
        else:
            forcing = None
            if direct:
                forcing = Forcing(
                    name,
                    [direct],
                    mask_index,
                    max_gap=gap_fill,
                    fill=defaults.get(name),
                )

        if forcing is None:
            if name in defaults:
                forcings[name] = _constant(name, mask_index, defaults[name])
                continue
            gs.fatal(
                _(
                    "No dataset supplies the '{var}' forcing. Pass it as an "
                    "STRDS with {var}=, or let t.graminr download it."
                ).format(var=name)
            )
        forcings[name] = forcing
    return forcings


def _constant(name, mask_index, value):
    """A forcing that is the same everywhere, every day."""

    class Constant(object):
        unit = None

        def __init__(self):
            self.name = name

        def values(self, day):
            return np.full(mask_index.size, float(value))

    return Constant()


def _era5_forcing(name, pick, mask_index, par_fraction, gap_fill, defaults, direct):
    """Build one forcing out of the ERA5 STRDS, already in model units."""
    fill = defaults.get(name)

    def make(strds_list, **kw):
        """This forcing, read from `strds_list`."""
        return Forcing(name, strds_list, mask_index, max_gap=gap_fill, fill=fill, **kw)

    # A directly supplied STRDS is already in the units the t.graminr options
    # document (degree C, mm/day, MJ/m2/day, kPa, m/s at 2 m), so it is taken
    # as it stands. Only a downloaded source needs the derivations below --
    # the PAR fraction, dewpoint to vapour pressure, the wind components and
    # the 10 m to 2 m reduction -- and applying them to a supplied dataset
    # would corrupt it.
    if direct:
        return make([direct])

    if name in ("tavg", "tmin", "tmax", "rain", "rad_glob", "et_pot"):
        # t.in.era5 already converts to degree C, mm/d and MJ m-2 d-1.
        strds = pick(*ERA5_SOURCES[name])
        return make([strds]) if strds else None
    if name == "rad_par":
        strds = pick("solar_radiation")
        return make([strds], transform=lambda v: v * par_fraction) if strds else None
    if name == "vp":
        strds = pick("dewpoint_temperature")
        return make([strds], transform=dewpoint_to_vapour_pressure) if strds else None
    if name == "wind":
        u_strds, v_strds = pick("wind_u"), pick("wind_v")
        if not (u_strds and v_strds):
            return None
        return make(
            [u_strds, v_strds],
            combine=wind_speed,
            transform=wind_10m_to_2m,
        )
    return None


def _ecmwf_forcing(name, pick, mask_index, par_fraction, gap_fill, defaults, direct):
    """Build one forcing out of the daily-aggregated ECMWF STRDS.

    ECMWF open-data is delivered in native units, so the conversions that
    t.in.era5 already applies are done here instead.
    """
    fill = defaults.get(name)

    def make(strds_list, **kw):
        """This forcing, read from `strds_list`."""
        return Forcing(name, strds_list, mask_index, max_gap=gap_fill, fill=fill, **kw)

    # A directly supplied STRDS is already in the units the t.graminr options
    # document (degree C, mm/day, MJ/m2/day, kPa, m/s at 2 m), so it is taken
    # as it stands. Only a downloaded source needs the derivations below --
    # the PAR fraction, dewpoint to vapour pressure, the wind components and
    # the 10 m to 2 m reduction -- and applying them to a supplied dataset
    # would corrupt it.
    if direct:
        return make([direct])

    if name == "tavg":
        strds = pick("2t:average")
        return make([strds], transform=kelvin_to_celsius) if strds else None
    if name == "tmin":
        strds = pick("2t:minimum", "2t:average")
        return make([strds], transform=kelvin_to_celsius) if strds else None
    if name == "tmax":
        strds = pick("2t:maximum", "2t:average")
        return make([strds], transform=kelvin_to_celsius) if strds else None
    if name == "rain":
        strds = pick("tp:sum")
        # IFS total precipitation is in metres of water.
        return make([strds], transform=lambda v: v * 1000.0) if strds else None
    if name in ("rad_glob", "rad_par"):
        strds = pick("ssrd:sum")
        if not strds:
            return None
        factor = 1.0e-6 * (par_fraction if name == "rad_par" else 1.0)
        return make([strds], transform=lambda v: v * factor)
    if name == "vp":
        strds = pick("2d:average")
        if not strds:
            return None
        return make(
            [strds],
            transform=lambda v: dewpoint_to_vapour_pressure(kelvin_to_celsius(v)),
        )
    if name == "wind":
        u_strds, v_strds = pick("10u:average"), pick("10v:average")
        if not (u_strds and v_strds):
            return None
        return make([u_strds, v_strds], combine=wind_speed, transform=wind_10m_to_2m)
    if name == "et_pot":
        # ECMWF open-data publishes no potential evaporation; it has to be
        # estimated, so the caller is told rather than silently given zeros.
        return None
    return None


def extraterrestrial_radiation(doy, latitude_deg):
    """Extraterrestrial radiation Ra [MJ m-2 d-1], FAO-56 Eq. 21."""
    phi = np.radians(latitude_deg)
    dr = 1.0 + 0.033 * np.cos(2.0 * np.pi * doy / 365.0)
    declination = 0.409 * np.sin(2.0 * np.pi * doy / 365.0 - 1.39)
    sunset = np.arccos(np.clip(-np.tan(phi) * np.tan(declination), -1.0, 1.0))
    return (
        (24.0 * 60.0 / np.pi)
        * 0.0820
        * dr
        * (
            sunset * np.sin(phi) * np.sin(declination)
            + np.cos(phi) * np.cos(declination) * np.sin(sunset)
        )
    )


class HargreavesET0(object):
    """Reference evapotranspiration estimated from temperature alone.

    Hargreaves and Samani (1985), as recommended by FAO-56 when humidity,
    wind and radiation are missing. This is what supplies the 'et_pot'
    forcing for the ECMWF open-data path, which publishes no potential
    evaporation of its own, and it is available as a fallback whenever an
    ERA5 potential-evaporation series cannot be used.
    """

    unit = "mm/day"

    def __init__(self, tmin_forcing, tmax_forcing, latitude):
        self.name = "et_pot"
        self.tmin = tmin_forcing
        self.tmax = tmax_forcing
        self.latitude = np.asarray(latitude, dtype=np.float64)

    def values(self, day):
        tmin = self.tmin.values(day)
        tmax = self.tmax.values(day)
        tmean = (tmin + tmax) / 2.0
        ra = extraterrestrial_radiation(day.timetuple().tm_yday, self.latitude)
        spread = np.sqrt(np.maximum(0.0, tmax - tmin))
        return np.maximum(0.0, 0.0023 * (tmean + 17.8) * spread * ra * 0.408)
