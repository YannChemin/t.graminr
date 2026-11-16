## DESCRIPTION

*t.graminr* runs the grassland growth models of the
[GraMinR](https://github.com/YannChemin/graminr) package over every
grassland pixel of the current computational region and writes one
space-time raster dataset (STRDS) per output variable.

Five models are available through **model**:

| **model** | Description |
|---|---|
| `lingra` | LINGRA: light interception and use, with the PCSE water balance |
| `modvege` | ModVege: four biomass compartments (Jouven et al., 2006) |
| `basgra` | BASGRA: carbon and nitrogen, snow, soil frost, winter hardening |
| `grasprog` | GrasProg/FOPROQ: growth indices per regrowth compartment |
| `grassim` | GrasSim: four compartments with a full nitrogen balance |

Each model is a daily time step applied to every pixel independently.
The simulation runs from **start** to **end** inclusive, and the results
are written as `<output>_<variable>`, one STRDS per variable, registered
at one-day intervals.

Use the **-l** flag to list the variables a model offers. By default a
small, model-specific selection is written; **variables** selects others,
and `variables=all` writes every one.

### Self-provisioning

The module can acquire everything it needs by itself, so that a run needs
no prepared data:

- the climate forcings through *t.in.era5* (ERA5-Land reanalysis) or
  *t.in.ecmwf* (IFS/AIFS open-data forecasts), selected with **climate**;
- the grassland mask through *r.in.landcover*, when **grassland** is left
  empty;
- the optional remotely sensed leaf area index through *r.in.modis*, with
  the **-r** flag.

ERA5-Land is a reanalysis and so suits a past or current period; the
ECMWF open data are forecasts and suit the coming days. Both are reduced
to the same daily forcings: temperature from kelvin to degrees Celsius,
dewpoint to actual vapour pressure, wind from 10 m to 2 m, and radiation
to MJ/m²/day. The conversions are applied to the pixel values as they are
read, so no derived maps are left in the mapset.

The ECMWF open data publish no potential evaporation, so a model that
needs **et_pot** must be run with `et_pot_method=hargreaves`, which
estimates reference evapotranspiration after Hargreaves and Samani
(1985) from temperature and latitude. Asking for `climate=ecmwf` with
`et_pot_method=dataset` fails rather than substituting anything silently.

With `climate=strds` the forcings are taken from existing datasets named
through **tavg**, **tmin**, **tmax**, **rain**, **rad_glob**,
**rad_par**, **et_pot**, **vp** and **wind**. Any of them may also be
given alongside `climate=era5` or `climate=ecmwf`, in which case the
named dataset is used and that variable is not downloaded.

### Restricting the run to grassland

**grassland** takes either a land-cover map, whose grassland classes are
named by **grassland_classes**, or a ready-made mask whose non-null and
non-zero cells are grassland. Left empty, the map is downloaded with
*r.in.landcover* from the collection given by **landcover_collection**,
and the grassland classes of that collection are used.

Only the masked pixels are simulated, and the outputs are NULL elsewhere.

**-r** downloads MODIS leaf area index (`modis-15A2H-061`, `Lai_500m`,
8-day, 500 m). It does not download evapotranspiration: the only MOD16
collection published on Planetary Computer is `modis-16A3GF-061`, which is
gap-filled per *year*, and one value a year cannot force a daily model.
**et_act_rs** and **t_act_rs** therefore have to be supplied as datasets if
they are wanted at all; left empty, the model computes them itself.

### Composite forcings and gap filling

An 8-day or 16-day MODIS composite does not cover every simulated day.
**gap_fill** is the number of days a map may be carried forward when no
map covers a day, so `gap_fill=8` lets an 8-day composite force each day
of its period without resampling it into thousands of maps.

### Cuts, grazing and fertilisation

A cutting calendar is given either per pixel through **cut_strds** or,
shared by every pixel, through **cut_dates**; the two are mutually
exclusive. **cut_dates** takes days of year (1-366, repeated every year)
and dates (`YYYY-MM-DD`). With neither, the model falls back on its own
cutting rule where it has one, as LINGRA does with `imopt=1`.

GraMinR's `gcut` does not mean the same thing in every model, and
*t.graminr* feeds each model the right quantity:

| Model | What `gcut` carries |
|---|---|
| LINGRA, BASGRA, GrasProg | a flag: the day is a cutting day |
| ModVege, GrasSim | the residual cutting height, in metres |

For ModVege and GrasSim the height comes from **cut_height**, defaulting
to the parameter set's `DEFAULT_CUT_HEIGHT` (0.05 m). Passing a flag
where a height is expected is silently wrong rather than an error -- a
height of 1 m is above the sward, so almost nothing is removed -- which
is why the height is a separate option.

Nitrogen is supplied as date/amount pairs: **fert_dates** (mineral
nitrogen, BASGRA and GrasSim), **fert_org_dates** (organic nitrogen,
GrasSim) and **ndep_dates** (background deposition, BASGRA, interpolated
linearly between the given points).

### Parameters

**parameters** takes a two-column `name,value` CSV in the format of
GraMinR's `read_config()`. Without it, the set shipped for the chosen
model is used. Individual values are overridden on the command line with
**parameter_values**, and any parameter can be made to vary across the
region by naming a raster for it in **parameter_maps**, for instance
`parameter_maps=WHC=soil_whc,NI=nutrition_index`.

### Memory

The forcings are read once into memory-mapped per-day, per-pixel stores
and the run is then chunked over pixels, **chunk** at a time, so a
multi-year run over a million pixels stays within the memory of an
ordinary machine. Lower **chunk** on a memory-tight machine and raise it
for speed.

Progress is reported through `G_percent()`, which
[MetaMonitor](https://github.com/YannChemin/MetaMonitor) picks up.

### Units

The models do not share units, and the same name can mean different things
in different ones. BASGRA keeps the Fortran's g C/m², g N/m² and g DM/m²,
while ModVege, GrasProg and GrasSim work in kg DM/ha and kg N/ha -- so
`nmin` is g N/m² under `model=basgra` and kg N/ha under `model=grassim`.

The **-l** flag prints each variable with its unit, every output raster
carries it (`r.info`), and it is named in the STRDS description (`t.info`).
A variable with no unit printed is dimensionless: an index, a fraction or a
status.

The *inputs* are not asymmetric in the same way. **fert_dates**,
**fert_org_dates** and **ndep_dates** are in kg N/ha and kg N/ha/day for
every model, because that is how a dressing is quoted; for BASGRA they are
converted to g N/m² on the way in. So a 60 kg N/ha dressing is written
`2024-04-10/60` whichever model is run, while BASGRA then reports the
resulting `nmin` as 6 g N/m² and GrasSim as 60 kg N/ha.

## NOTES

### PaSim is not included

GraMinR's sixth model, PaSim, is deliberately absent.

### Deviations from GraMinR's R code

Two places where this port does not reproduce GraMinR's R exactly. Both
are deliberate and both are covered by the tests.

**LINGRA, photosynthesis on a cutting day.** GraMinR suspends
assimilation on every scheduled cutting day; PCSE's LINGRA, which GraMinR
was validated against, suspends it only when biomass is actually removed.
The two agree until a scheduled cut lands on a day when green leaf weight
is already at or below `cwlvg`: nothing is harvested, yet GraMinR still
loses that day's assimilation. On the Zegveld 1986 reference run this
happens on day of year 240 and costs about 11% of green leaf weight for
the rest of the season. The LINGRA parameter `HarvestGate` selects
between the two rules and defaults to `harvest`, PCSE's. Set
`HarvestGate=mowflag` to reproduce GraMinR's R output.

**ModVege, grazing.** The grazing branch of GraMinR's `modvege` cannot be
reached, so grazing is reimplemented here from the published model, as
`graze_count` animals of `graze_weight` each, ingesting at most
`maxAmountToIngest` and taking from each compartment in proportion to
what it offers above the cutting height.

One place where a quirk is preserved rather than corrected: BASGRA's
`surfacepool()` computes `poolinfil` but does not return it, so GraMinR's
water balance always adds zero there. That branch *is* reachable, so the
omission is kept, with a comment, so that results match.

### Installing over an older copy

GRASS resolves an addon from `GRASS_ADDON_BASE` (`~/.grass8/addons` by
default) before the source tree's own `dist` directory, so an older copy
installed there shadows a newer build and the module appears not to have
changed. `install_local.sh` in the source directory builds and then copies
the script, the `graminr_lib` package, the parameter sets and the
documentation into every GRASS installation on the machine, comparing each
installed file with the build afterwards.

### Validation

The ported models are compared against the reference outputs that
GraMinR ships, in `tests/`:

| Model | Reference | Agreement |
|---|---|---|
| LINGRA | PCSE, Zegveld 1986, 366 days | all 14 published variables to 6e-8 |
| BASGRA | BASGRA Fortran, Saerheim 2000, 254 days | 42 of 48 quantities to 1e-5, the rest at the single-precision floor |
| GrasSim | the original GrasSim workbook, 2019 | all 23 daily quantities within 0.42%, most within 0.01% |
| ModVege | none shipped | behavioural tests only |
| GrasProg | none shipped | behavioural tests only |

The comparison is against the reference output of each original
implementation -- PCSE, the BASGRA Fortran, the GrasSim workbook -- not
against GraMinR's R, which could not be run here.

Two discrepancies in GraMinR's own reference material were found while
doing this and are documented in the tests:

- the LINGRA reference parameter file carries `ANGSTA`/`ANGSTB` =
  0.22/0.47, while the PCSE run that produced its reference output used
  the PCSE and WOFOST defaults 0.18/0.55. The pair enters only the Penman
  net long-wave term, so it moves `E0`, `ES0` and through them `EVS`, and
  leaves the Penman-Monteith `ET0` untouched;
- the BASGRA Fortran computes the day length from the *previous* day of
  year. Its reference output's `DAYL` column equals the same formula at
  `doy - 1` on all 253 days after the first. This port follows GraMinR in
  using the current day; the BASGRA parameter `dayl_doy_offset` can be
  set to `-1` to reproduce the Fortran.

## EXAMPLES

### ModVege over a season, provisioning everything itself

```sh
t.graminr model=modvege start=2024-03-01 end=2024-10-31 \
    output=modvege_2024 cut_dates=135,175,215,255 cut_height=0.05
```

The grassland mask comes from ESA WorldCover and the climate from
ERA5-Land. The outputs are `modvege_2024_gvb`, `modvege_2024_hb` and so
on.

### LINGRA on prepared datasets, with a per-pixel soil property

```sh
t.graminr model=lingra start=2023-01-01 end=2023-12-31 \
    output=lingra_2023 climate=strds \
    tmin=era5_tmin tmax=era5_tmax rad_glob=era5_rad rain=era5_rain \
    vp=era5_vp wind=era5_wind \
    grassland=landcover grassland_classes=30 \
    parameter_maps=RDMSOL=soil_rooting_depth \
    cut_dates=135,165,200,240,280
```

### BASGRA with a nitrogen calendar and MODIS leaf area index

```sh
t.graminr -r model=basgra start=2024-01-01 end=2024-12-31 \
    output=basgra_2024 gap_fill=8 \
    cut_dates=2024-05-30,2024-08-04 \
    fert_dates=2024-04-24/14,2024-05-30/8 \
    ndep_dates=2024-01-01/0.0055,2024-12-31/0.0055 \
    variables=lai,dry_matter,yield_last,nmin
```

### Listing what a model can write

```sh
t.graminr -l model=grassim
```

### Aggregating a result to monthly means

```sh
t.rast.aggregate input=modvege_2024_gvb output=modvege_2024_gvb_monthly \
    basename=gvb_month granularity="1 months" method=average
```

## REFERENCES

- Jouven, M., Carrère, P., Baumont, R. (2006). *Model predicting dynamics
  of biomass, structure and digestibility of herbage in managed permanent
  pastures. 1. Model description.* Grass and Forage Science 61(2),
  112-124.
- Schapendonk, A.H.C.M., Stol, W., van Kraalingen, D.W.G., Bouman, B.A.M.
  (1998). *LINGRA, a sink/source model to simulate grassland
  productivity in Europe.* European Journal of Agronomy 9(2-3), 87-100.
- Höglind, M., Van Oijen, M., Cameron, D., Persson, T. (2016).
  *Process-based simulation of growth and overwintering of grassland
  using the BASGRA model.* Ecological Modelling 335, 1-15.
- Hargreaves, G.H., Samani, Z.A. (1985). *Reference crop
  evapotranspiration from temperature.* Applied Engineering in
  Agriculture 1(2), 96-99.
- Allen, R.G., Pereira, L.S., Raes, D., Smith, M. (1998). *Crop
  evapotranspiration: guidelines for computing crop water requirements.*
  FAO Irrigation and Drainage Paper 56.

## SEE ALSO

*[g.region](g.region.md),
[r.in.landcover](r.in.landcover.md),
[r.in.modis](r.in.modis.md),
[t.create](t.create.md),
[t.in.ecmwf](t.in.ecmwf.md),
[t.in.era5](t.in.era5.md),
[t.rast.aggregate](t.rast.aggregate.md),
[t.rast.list](t.rast.list.md),
[t.register](t.register.md),
[t.support](t.support.md)*

## AUTHORS

Yann Chemin

Ported from the [GraMinR](https://github.com/YannChemin/graminr) R
package.
