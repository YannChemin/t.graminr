"""
MODULE:     t.graminr (graminr_lib.gio)

AUTHOR(S):  Yann Chemin

PURPOSE:    Space-time raster dataset access for t.graminr: resolve the map
            covering each simulated day, read only the grassland pixels of
            it, and write model outputs back as a registered STRDS.

            Progress is reported to MetaMonitor through gs.percent(), which
            the libmmgrass.so shim credits to the calling script.
            See https://github.com/YannChemin/MetaMonitor

SPDX-License-Identifier: Unlicense
"""

import datetime as dt
import os
import tempfile

import numpy as np

import grass.script as gs
import grass.script.array as garray


def region_shape():
    """(rows, cols) of the current computational region."""
    region = gs.region()
    return int(region["rows"]), int(region["cols"])


def latitude_longitude_arrays():
    """Per-cell latitude and longitude [degree] of the current region.

    In a latitude/longitude project the cell centres are already degrees;
    in any other projection they are reprojected with m.proj, so that the
    solar geometry of LINGRA and BASGRA is correct across a wide region.
    """
    region = gs.region()
    rows, cols = int(region["rows"]), int(region["cols"])
    east = np.linspace(
        float(region["w"]) + float(region["ewres"]) / 2.0,
        float(region["e"]) - float(region["ewres"]) / 2.0,
        cols,
    )
    north = np.linspace(
        float(region["n"]) - float(region["nsres"]) / 2.0,
        float(region["s"]) + float(region["nsres"]) / 2.0,
        rows,
    )
    xx, yy = np.meshgrid(east, north)

    projection = gs.parse_command("g.proj", flags="g")
    if projection.get("unit", "").lower() in ("degree", "degrees") or projection.get(
        "proj"
    ) in ("ll", "longlat"):
        return yy, xx

    # Reproject every cell centre to WGS84 in one m.proj call: -o takes the
    # current project as input and lon/lat WGS84 as output, -d asks for
    # decimal degrees. m.proj writes "lon lat z" per line.
    #
    # The coordinates go through a file rather than a pipe. grass.script's
    # read_command() cannot take a string on stdin at all, and feeding a
    # whole region through a pipe while not draining stdout would deadlock
    # as soon as the region is larger than the pipe buffer -- which is most
    # of them.
    handle, path = tempfile.mkstemp(prefix="graminr_mproj_", suffix=".txt")
    try:
        with os.fdopen(handle, "w") as stream:
            stream.writelines(
                "{} {}\n".format(x, y) for x, y in zip(xx.ravel(), yy.ravel())
            )
        out = gs.read_command(
            "m.proj", flags="od", input=path, separator="space", quiet=True
        )
    finally:
        try:
            os.remove(path)
        except OSError:
            pass
    rendered = [
        line.replace("|", " ").split() for line in out.strip().splitlines() if line
    ]
    if len(rendered) != rows * cols:
        gs.fatal(
            _(
                "m.proj returned {got} reprojected coordinate(s) for the "
                "{want} cell(s) of the region."
            ).format(got=len(rendered), want=rows * cols)
        )
    values = np.array([[float(field) for field in fields[:2]] for fields in rendered])
    lon = values[:, 0].reshape(rows, cols)
    lat = values[:, 1].reshape(rows, cols)
    return lat, lon


def read_raster(name, dtype=np.float64):
    """Read a raster over the current region, NULL becoming NaN."""
    data = garray.array(mapname=name, null=np.nan, dtype=np.float32)
    return np.asarray(data, dtype=dtype)


def write_raster(array, name, rows, cols, title=None, overwrite=False):
    """Write a 2-D array as a raster, NaN becoming NULL."""
    buffer = garray.array(dtype=np.float32)
    buffer[...] = np.asarray(array, dtype=np.float32).reshape(rows, cols)
    buffer.write(
        mapname=name, title=title, null=np.nan, overwrite=overwrite, quiet=True
    )


def strds_maps(strds):
    """Maps of an STRDS as a list of (name, start_date, end_date).

    Dates are dt.date; end_date is None for an instant-in-time map.
    """
    text = gs.read_command(
        "t.rast.list",
        input=strds,
        columns="name,start_time,end_time",
        separator="|",
        flags="u",
        quiet=True,
    )
    maps = []
    for line in text.strip().splitlines():
        parts = line.split("|")
        if len(parts) < 2 or not parts[0]:
            continue
        name = parts[0]
        start = _parse_time(parts[1])
        end = _parse_time(parts[2]) if len(parts) > 2 else None
        if start is None:
            continue
        maps.append((name, start, end))
    maps.sort(key=lambda item: item[1])
    if not maps:
        gs.fatal(_("Space-time raster dataset <{}> holds no map.").format(strds))
    return maps


def _parse_time(text):
    """Parse a GRASS absolute time stamp into a date."""
    if text is None:
        return None
    text = text.strip()
    if not text or text.lower() == "none":
        return None
    for pattern in ("%Y-%m-%d %H:%M:%S", "%Y-%m-%d %H:%M", "%Y-%m-%d"):
        try:
            return dt.datetime.strptime(text, pattern).date()
        except ValueError:
            continue
    return None


class DayIndex:
    """Maps a simulated day to the map of an STRDS that covers it.

    A map whose interval contains the day is used directly. Otherwise the
    most recent earlier map is carried forward, up to `max_gap` days: this
    is what turns an 8-day or 16-day remote-sensing composite into a daily
    forcing without resampling it into thousands of new maps.
    """

    def __init__(self, strds, max_gap=0):
        self.strds = strds
        self.max_gap = int(max_gap)
        self.maps = strds_maps(strds)
        self._starts = [item[1] for item in self.maps]

    def map_for(self, day):
        """Name of the map covering `day`, or None."""
        best = None
        for name, start, end in self.maps:
            if start > day:
                break
            if end is not None and end <= day and (day - end).days > self.max_gap:
                continue
            if end is None and (day - start).days > self.max_gap:
                continue
            best = name
        return best


class Forcing:
    """One model forcing variable, read day by day over the chosen pixels.

    `transform` is applied to the raw values of each map; `combine` builds a
    variable out of several maps of the same day (the wind speed out of its
    two components, for instance).
    """

    def __init__(
        self,
        name,
        sources,
        mask_index,
        transform=None,
        combine=None,
        max_gap=0,
        fill=None,
        unit=None,
    ):
        self.name = name
        self.sources = list(sources)
        self.mask_index = mask_index
        self.transform = transform
        self.combine = combine
        self.fill = fill
        self.unit = unit
        self.indexes = [DayIndex(s, max_gap=max_gap) for s in self.sources]
        self._cache = {}

    def values(self, day):
        """Values of this forcing for `day`, over the masked pixels."""
        parts = []
        for index in self.indexes:
            name = index.map_for(day)
            if name is None:
                if self.fill is None:
                    gs.fatal(
                        _(
                            "No map of <{strds}> covers {day}, needed for the "
                            "'{var}' forcing. Extend the dataset, shorten the "
                            "simulated period, or raise gap_fill=."
                        ).format(strds=index.strds, day=day, var=self.name)
                    )
                parts.append(None)
                continue
            if name not in self._cache:
                # Keep only as many maps as this forcing reads per day, so a
                # multi-component variable does not evict its own first part.
                if len(self._cache) >= len(self.indexes):
                    self._cache.clear()
                self._cache[name] = read_raster(name).ravel()[self.mask_index]
            parts.append(self._cache[name])

        if any(part is None for part in parts):
            return np.full(self.mask_index.size, float(self.fill))

        if self.combine is not None:
            value = self.combine(*parts)
        else:
            value = parts[0]
        if self.transform is not None:
            value = self.transform(value)
        value = np.asarray(value, dtype=np.float64)
        if self.fill is not None:
            value = np.where(np.isfinite(value), value, float(self.fill))
        return value


def simulated_days(start, end):
    """List of dt.date from `start` to `end`, inclusive."""
    if end < start:
        gs.fatal(_("The end date is before the start date."))
    count = (end - start).days + 1
    return [start + dt.timedelta(days=offset) for offset in range(count)]


def set_raster_metadata(raster, unit=None, history=None):
    """Record the data unit and the command history of one output raster.

    Both go through a single r.support call. That matters: a run writes one
    raster per variable per day, and the obvious pair -- set the unit, then
    call grass.script.raster_history() -- costs about five subprocesses per
    raster, because raster_history() itself shells out to g.gisenv,
    g.findfile, g.tempfile and r.support. Over a year and ten variables
    that is tens of thousands of process spawns, which dominates the whole
    run; one call per raster is five times cheaper.

    An STRDS has no data unit of its own (t.support's `unit` is the unit of
    *relative* time, not of the data), so the unit belongs here, on each
    registered raster, where r.info and the GUI show it.
    """
    options = {}
    if unit:
        options["units"] = unit
    if history:
        options["history"] = history
    if not options:
        return
    gs.run_command("r.support", map=raster, quiet=True, **options)


def register_strds(
    name,
    maps,
    title,
    description,
    unit=None,
    semantic="mean",
    overwrite=False,
    temporary_dir=None,
):
    """Create an STRDS and register `maps`, a list of (map name, date)."""
    # GRASS's temporal framework interpolates the title and the description
    # straight into SQL, so an apostrophe in either makes t.create fail with
    # "sqlite3.OperationalError: near ...: syntax error". Callers are not
    # expected to know that, and a parameter-set or output name could carry
    # one, so they are replaced here.
    title = str(title).replace("'", "")
    description = str(description).replace("'", "")
    # An STRDS carries no data unit, so it is named in the description, where
    # t.info shows it; set_raster_metadata() puts it on the rasters themselves.
    if unit:
        description = "{} [{}]".format(description, unit)
    gs.run_command(
        "t.create",
        output=name,
        type="strds",
        temporaltype="absolute",
        semantictype=semantic,
        title=title,
        description=description,
        overwrite=overwrite,
        quiet=True,
    )
    listing = os.path.join(
        temporary_dir, "register_{}.txt".format(name.replace("@", "_"))
    )
    # Explicit start and end give each map a one-day interval, which is what
    # t.rast.aggregate and t.rast.univar expect of a daily series.
    with open(listing, "w") as handle:
        for map_name, day in maps:
            handle.write(
                "{}|{}|{}\n".format(
                    map_name,
                    day.strftime("%Y-%m-%d"),
                    (day + dt.timedelta(days=1)).strftime("%Y-%m-%d"),
                )
            )
    gs.run_command(
        "t.register",
        input=name,
        type="raster",
        file=listing,
        separator="|",
        quiet=True,
    )
    # The unit is not an STRDS property: t.support has no unit= option (the
    # STRDS "unit" in GRASS's temporal framework is the unit of *relative*
    # time, not of the data). It belongs on the rasters, and set_raster_metadata()
    # above puts it there as each one is written.
