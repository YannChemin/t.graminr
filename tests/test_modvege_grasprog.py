"""
MODULE:     t.graminr (tests)

AUTHOR(S):  Yann Chemin

PURPOSE:    Behavioural tests for the ModVege and GrasProg ports.

            GraMinR ships no reference output for either model -- only
            inputs (inst/extdata/modvege_weather.csv and params, and
            inst/extdata/original_outputs/modvege/modvege_ts_marcenat1994.csv
            for Marcenat 1994; inst/extdata/grasprog_* and the Kiel weather
            series) and, for GrasProg, a development summary describing how
            model_grasprog() was aligned with the original FOPROQ. There is
            therefore nothing to validate against numerically, so these
            tests check the behaviour the two models are specified to have.

SPDX-License-Identifier: Unlicense
"""

import csv
import os

import numpy as np
import pytest

from graminr_lib import grasprog, modvege
from graminr_lib.params import cast_parameters

GRAMINR = os.environ.get("GRAMINR_DIR", os.path.expanduser("~/dev/graminr"))
EXTDATA = os.path.join(GRAMINR, "inst", "extdata")

pytestmark = pytest.mark.skipif(
    not os.path.isdir(EXTDATA),
    reason="the GraMinR example data is not available; set GRAMINR_DIR",
)


def _params(name, model):
    raw = {}
    with open(os.path.join(EXTDATA, name), newline="") as handle:
        for row in csv.reader(handle):
            if len(row) >= 2 and row[0].strip():
                raw[row[0].strip()] = row[1].strip()
    return cast_parameters(raw, model)


def _forcing(name, columns, cut_doys=(), cut_value=1.0):
    with open(os.path.join(EXTDATA, name), newline="") as handle:
        rows = [
            {key.strip(): value for key, value in row.items()}
            for row in csv.DictReader(handle)
        ]
    ndays = len(rows)
    forcing = {
        "year": np.array([int(row["year"]) for row in rows]),
        "doy": np.array([int(row["doy"]) for row in rows]),
    }
    for column in columns:
        forcing[column] = np.array([[float(row[column])] for row in rows])
    forcing["gcut"] = np.array(
        [[cut_value if int(row["doy"]) in set(cut_doys) else 0.0] for row in rows]
    )
    return forcing, ndays


MODVEGE_COLUMNS = ("tavg", "rad_par", "rain", "et_pot", "et_act_rs", "lai_rs")
GRASPROG_COLUMNS = ("tavg", "rad_par", "rain", "et_pot")

#: GraMinR's gcut does not mean the same thing in every model: ModVege and
#: GrasSim read it as the cutting height in metres, while LINGRA, GrasProg
#: and BASGRA read it as a flag. t.graminr feeds each model the right one
#: (see cut_value in t.graminr.py), and so do these tests.
MODVEGE_CUT = 0.05
GRASPROG_CUT = 1.0


def test_modvege_runs_on_its_example_series():
    """Every published output is finite and the compartments stay positive."""
    params = _params("modvege_params.csv", "modvege")
    forcing, ndays = _forcing("modvege_weather.csv", MODVEGE_COLUMNS)
    out = modvege.run(params, forcing, 1)
    assert set(modvege.OUTPUTS) <= set(out)
    for name in modvege.OUTPUTS:
        values = out[name]
        assert values.shape == (ndays, 1), name
        assert np.all(np.isfinite(values)), name
    for name in ("gvb", "grb", "dvb", "drb"):
        assert np.all(out[name] >= -1e-6), name


def test_modvege_cutting_removes_biomass_and_is_cumulative():
    """A cut lowers standing biomass and raises the harvested total."""
    params = _params("modvege_params.csv", "modvege")
    uncut, _ = _forcing("modvege_weather.csv", MODVEGE_COLUMNS)
    cut, _ = _forcing(
        "modvege_weather.csv",
        MODVEGE_COLUMNS,
        cut_doys=(150, 200, 250),
        cut_value=MODVEGE_CUT,
    )

    out_uncut = modvege.run(params, uncut, 1)
    out_cut = modvege.run(params, cut, 1)
    harvested = out_cut["hb"][:, 0]
    assert np.all(np.diff(harvested) >= -1e-6), "hb must be cumulative"
    assert harvested[-1] > 0.0
    assert float(out_uncut["hb"][-1, 0]) == 0.0

    # The harvest grows only on the cut days themselves.
    grew = np.nonzero(np.diff(harvested) > 1e-6)[0] + 1
    cut_days = np.nonzero(cut["gcut"][:, 0] > 0)[0]
    assert set(grew.tolist()) <= set(cut_days.tolist())

    # Standing biomass is NOT lower after a cut: removing the canopy resets
    # ModVege's age classes, so senescence collapses and the sward regrows
    # faster than the undisturbed one, which is declining at this point.
    # That is GraMinR's behaviour and is kept.
    assert out_cut["gvb"][151, 0] > out_uncut["gvb"][151, 0]

    # Cutting more often removes more over the season.
    more, _ = _forcing(
        "modvege_weather.csv",
        MODVEGE_COLUMNS,
        cut_doys=(150, 175, 200, 225, 250),
        cut_value=MODVEGE_CUT,
    )
    assert float(modvege.run(params, more, 1)["hb"][-1, 0]) > harvested[-1]


def test_modvege_temperature_sum_resets_each_calendar_year():
    """The temperature sum is a within-year accumulation."""
    params = _params("modvege_params.csv", "modvege")
    forcing, _ = _forcing("modvege_weather.csv", MODVEGE_COLUMNS)
    out = modvege.run(params, forcing, 1)
    sum_t = out["stp"][:, 0]
    year = forcing["year"]
    for value in np.unique(year):
        days = np.nonzero(year == value)[0]
        within = sum_t[days]
        assert np.all(np.diff(within) >= -1e-6)
        assert within[0] < within[-1] or within[-1] == 0.0


def test_grasprog_runs_and_waits_for_the_temperature_sum():
    """Growth is pinned at the floor until begsum degree-days accumulate."""
    params = _params("grasprog_params.csv", "grasprog")
    forcing, ndays = _forcing("grasprog_weather.csv", GRASPROG_COLUMNS)
    out = grasprog.run(params, forcing, 1)
    for name in grasprog.OUTPUTS:
        assert out[name].shape == (ndays, 1), name
        assert np.all(np.isfinite(out[name])), name
    # begsum defaults to 250 degree-days on base 0, so the first days of
    # January cannot have grown past the floor.
    weight = out["w"][:, 0]
    assert weight[0] == pytest.approx(weight[5], rel=1e-9)
    assert weight.max() > weight[0]


def test_grasprog_waterstress_off_tops_the_soil_up():
    """waterstress=no irrigates to field capacity and lifts the yield."""
    stressed = _params("grasprog_params.csv", "grasprog")
    watered = dict(stressed)
    watered["waterstress"] = "no"
    forcing, _ = _forcing("grasprog_weather.csv", GRASPROG_COLUMNS)
    out_stressed = grasprog.run(stressed, forcing, 1)
    out_watered = grasprog.run(watered, forcing, 1)
    assert out_watered["irrig"].sum() > 0.0
    assert out_stressed["irrig"].sum() == 0.0
    assert out_watered["w"].max() >= out_stressed["w"].max()


def test_grasprog_cut_lowers_the_standing_weight():
    """A cut resets the growth compartment from the next day on."""
    params = _params("grasprog_params.csv", "grasprog")
    uncut, _ = _forcing("grasprog_weather.csv", GRASPROG_COLUMNS)
    cut, _ = _forcing(
        "grasprog_weather.csv",
        GRASPROG_COLUMNS,
        cut_doys=(150,),
        cut_value=GRASPROG_CUT,
    )
    out_uncut = grasprog.run(params, uncut, 1)
    out_cut = grasprog.run(params, cut, 1)
    # The cut is applied on the day itself: the stored weight drops to the
    # compartment floor and stays below the undisturbed run afterwards.
    assert out_cut["w"][149, 0] == pytest.approx(out_uncut["w"][149, 0], rel=1e-9)
    assert out_cut["w"][150, 0] < out_uncut["w"][150, 0]
    assert out_cut["w"][151, 0] < out_uncut["w"][151, 0]
    # The cut is booked as yield. GrasProg's yield is per cut, not
    # cumulative, so it is the peak that carries it.
    assert float(out_cut["yield"].max()) > float(out_uncut["yield"].max())


def test_gcut_means_a_height_in_modvege_and_a_flag_in_grasprog():
    """The two readings of gcut, pinned so a future change cannot pass.

    Feeding ModVege a flag of 1.0 is read as a 1 m cutting height, which on
    this sward removes next to nothing; feeding it 0.05 m harvests properly.
    GrasProg, by contrast, tests gcut == 1 and ignores the magnitude.
    """
    modvege_params = _params("modvege_params.csv", "modvege")
    as_height, _ = _forcing(
        "modvege_weather.csv",
        MODVEGE_COLUMNS,
        cut_doys=(150, 200, 250),
        cut_value=0.05,
    )
    as_flag, _ = _forcing(
        "modvege_weather.csv",
        MODVEGE_COLUMNS,
        cut_doys=(150, 200, 250),
        cut_value=1.0,
    )
    harvest_height = float(modvege.run(modvege_params, as_height, 1)["hb"][-1, 0])
    harvest_flag = float(modvege.run(modvege_params, as_flag, 1)["hb"][-1, 0])
    assert harvest_height > 10.0 * max(harvest_flag, 1e-9)

    grasprog_params = _params("grasprog_params.csv", "grasprog")
    flag, _ = _forcing(
        "grasprog_weather.csv", GRASPROG_COLUMNS, cut_doys=(150,), cut_value=1.0
    )
    height, _ = _forcing(
        "grasprog_weather.csv", GRASPROG_COLUMNS, cut_doys=(150,), cut_value=0.05
    )
    # Only the flag triggers a GrasProg cut.
    assert grasprog.run(grasprog_params, flag, 1)["yield"].max() > 0.0
    assert grasprog.run(grasprog_params, height, 1)["yield"].max() == 0.0
