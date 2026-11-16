# t.graminr

A GRASS GIS temporal addon that runs the [GraMinR](https://github.com/YannChemin/graminr)
grassland growth models over every grassland pixel of a region, writing one
space-time raster dataset per output variable.

Five of GraMinR's models are ported: **LINGRA**, **ModVege**, **BASGRA**,
**GrasProg/FOPROQ** and **GrasSim**. PaSim is not included.

The module can provision its own inputs: climate from *t.in.era5* or
*t.in.ecmwf*, the grassland mask from *r.in.landcover*, and optional
remotely sensed LAI and evapotranspiration from *r.in.modis*.

See `t.graminr.md` for the manual.

## Layout

```
t.graminr.py        the GRASS module: parser, data acquisition, raster I/O
graminr_lib/        the models, with no GRASS dependency
  core.py           shared numerics (afgen, clamp, year_cumsum, ...)
  lingra.py         LINGRA, with FAO-56 Penman / Penman-Monteith
  modvege.py        ModVege
  basgra.py         BASGRA
  grasprog.py       GrasProg / FOPROQ
  grassim.py        GrasSim
  params.py         parameter reading, renames, derived values, rasters
  cutdates.py       cutting and fertilisation calendars
  fetch.py          t.in.era5 / t.in.ecmwf / r.in.landcover / r.in.modis
  gio.py            raster and STRDS I/O, day indexing, forcings
  msg.py            messaging shim, so the models run outside GRASS
params/             the default parameter set of each model
tests/              validation against GraMinR's reference outputs
```

Each model is a daily loop whose single array axis is the pixel: every
scalar branch of the R code is a `numpy.where()`. Pixels are independent,
so a run is chunked over them, with the forcings and outputs in memory-mapped
files, which keeps a multi-year run over a million pixels within ordinary
memory.

The models, the parameter handling and the management calendars import no
GRASS module, so they can be run and tested without a mapset; `msg.py`
stands in for `grass.script`'s messaging and raises instead of exiting.
Only `gio.py` and `fetch.py` need a session, since they do the raster and
STRDS I/O and call the other modules.

## Build

Standalone, against the GRASS source tree, through the symlink in
`grass-addons`:

```sh
make -C ~/dev/grass-addons/src/temporal/t.graminr MODULE_TOPDIR=~/dev/grass
```

To build *and* install into every GRASS on this machine:

```sh
./install_local.sh
```

That second step is not optional housekeeping. GRASS resolves a module from
`GRASS_ADDON_BASE` (`~/.grass8/addons`) **before** the source tree's own
`dist` directory, so a stale copy there silently shadows a freshly built
one -- a rebuilt module appears to change nothing. The script copies the
script, `graminr_lib/`, `params/`, the HTML page and the man page into
`~/.grass8/addons` and `/usr/local/grass86`, then compares each installed
file with the build and fails if any differs.


## Tests

```sh
python3 -m pytest tests
```

53 tests. The reference comparisons need a GraMinR checkout for its data;
they look in `~/dev/graminr` or at `$GRAMINR_DIR`, and skip if it is
absent. The calendar tests need nothing.

| Model | Reference | Agreement |
|---|---|---|
| LINGRA | PCSE, Zegveld 1986, 366 days | all 14 published variables to 6e-8 |
| BASGRA | BASGRA Fortran, Saerheim 2000, 254 days | 42 of 48 quantities to 1e-5, the rest at the single-precision floor |
| GrasSim | the original GrasSim workbook, 2019 | all 23 daily quantities within 0.42%, most within 0.01% |
| ModVege | none shipped by GraMinR | behavioural tests only |
| GrasProg | none shipped by GraMinR | behavioural tests only |

ModVege and GrasProg get behavioural tests because GraMinR ships no reference
output for them, only inputs.

The deviations from GraMinR's R, and two discrepancies found in GraMinR's own
reference material, are documented under NOTES in the manual and in the test
docstrings.

The models were compared against the reference outputs the original
implementations produced, not against GraMinR's R: R is not installed on this
machine, so GraMinR itself could not be run.

## Licence

Unlicense; see `LICENSE`.
