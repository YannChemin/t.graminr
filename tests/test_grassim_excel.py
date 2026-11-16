"""
MODULE:     t.graminr (tests)

AUTHOR(S):  Yann Chemin

PURPOSE:    Numerical validation of the GrasSim port against the original
            GrasSim workbook shipped with GraMinR
            (inst/extdata/original_outputs/grassim/), the fertilised and cut
            scenario of 2019.

            All 23 daily quantities the workbook publishes are compared over
            its 365 days. The workbook decides cuts and fertilisation
            adaptively (cut once standing biomass passes 2500 kg/ha, then
            fertilise five days later, at most twice); the port takes the
            dates it is given, so the schedule it actually realised is read
            out of the workbook's own Management_events sheet: cuts on days
            of year 112, 162, 191, 228 and 265, and mineral nitrogen of 80,
            60 and 40 kg/ha on days 60, 117 and 167. That is the 180 kg
            N/ha and five cuts its Scenario_summary reports.

            Note that GrasSim's gcut is a cutting height in metres, not a
            flag: the workbook cuts to 0.05 m.

            The agreement is within 0.42% on every quantity and within 0.01%
            on most, over the whole year, including a total harvest of
            12 735 kg/ha. The residual is not float32 noise -- it is about
            600 times the single-precision floor -- and is attributed to the
            workbook being the original spreadsheet implementation rather
            than GraMinR's R.

SPDX-License-Identifier: Unlicense
"""

import csv
import os

import numpy as np
import pytest

from graminr_lib import grassim
from graminr_lib.params import cast_parameters
from xlsx_reader import read_sheet

GRAMINR = os.environ.get("GRAMINR_DIR", os.path.expanduser("~/dev/graminr"))
EXTDATA = os.path.join(GRAMINR, "inst", "extdata")
WORKBOOK = os.path.join(
    EXTDATA, "original_outputs", "grassim", "grassim_fertilized_cut.xlsx"
)

CUT_DOYS = {112, 162, 191, 228, 265}
CUT_HEIGHT = 0.05
FERTILISER = {60: 80.0, 117: 60.0, 167: 40.0}

#: Workbook column -> ported output name.
COLUMNS = {
    "gvb": "gvb",
    "grb": "grb",
    "dvb": "dvb",
    "drb": "drb",
    "gva": "gva",
    "gra": "gra",
    "dva": "dva",
    "dra": "dra",
    "g": "g",
    "pgr": "pgr",
    "env": "env",
    "sea": "sea",
    "ftm": "ftm",
    "lai": "lai",
    "hb": "hb",
    "stp": "stp",
    "tsb": "tsb",
    "fW": "fw",
    "fN": "fn",
    "RNC": "rnc",
    "Nmin": "nmin",
    "Norg": "norg",
    "swardHeight": "sward_height",
}

#: The largest relative difference accepted on any quantity, over the year.
TOLERANCE = 5e-3

pytestmark = pytest.mark.skipif(
    not os.path.exists(WORKBOOK),
    reason="the GraMinR reference outputs are not available; set GRAMINR_DIR",
)


def _inputs():
    raw = {}
    with open(os.path.join(EXTDATA, "grassim_params.csv"), newline="") as handle:
        for row in csv.reader(handle):
            if len(row) >= 2 and row[0].strip():
                raw[row[0].strip()] = row[1].strip()
    params = cast_parameters(raw, "grassim")

    with open(os.path.join(EXTDATA, "grassim_weather.csv"), newline="") as handle:
        weather = list(csv.DictReader(handle))
    ndays = len(weather)

    def column(key):
        return np.array([[float(row[key])] for row in weather])

    forcing = {
        "year": np.array([int(row["year"]) for row in weather]),
        "doy": np.array([int(row["doy"]) for row in weather]),
        "tavg": column("tavg"),
        "rad_par": column("rad_par"),
        "rain": column("rain"),
        "et_pot": column("et_pot"),
        "st_grass": column("st_grass"),
        # gcut carries the cutting height, in metres.
        "gcut": np.array(
            [[CUT_HEIGHT if int(row["doy"]) in CUT_DOYS else 0.0] for row in weather]
        ),
        "fert_min": np.array(
            [[FERTILISER.get(int(row["doy"]), 0.0)] for row in weather]
        ),
        "fert_org": np.zeros((ndays, 1)),
    }
    return params, forcing


def _reference():
    return read_sheet(WORKBOOK, "Daily_outputs")


def test_the_reference_schedule_is_the_one_used():
    """The cuts and fertiliser applications come from the workbook itself."""
    events = read_sheet(WORKBOOK, "Management_events")
    cuts = {int(row["Julian_day"]) for row in events if row["Event"] == "Cut"}
    fertiliser = {
        int(row["Julian_day"]): float(row["Mineral_N_kg_ha"])
        for row in events
        if row["Event"] == "Mineral fertilization"
    }
    assert cuts == CUT_DOYS
    assert fertiliser == FERTILISER
    assert all(float(row["Cut_height_m"]) in (0.0, CUT_HEIGHT) for row in events)

    summary = read_sheet(WORKBOOK, "Scenario_summary")[0]
    assert int(summary["Number_of_cuts"]) == len(CUT_DOYS)
    assert float(summary["Total_mineral_N_kg_ha"]) == sum(FERTILISER.values())


def test_matches_the_grassim_workbook():
    """Every daily quantity the workbook publishes agrees within 0.5%."""
    params, forcing = _inputs()
    reference = _reference()
    assert len(reference) == len(forcing["doy"]) == 365
    out = grassim.run(params, forcing, 1)
    failures = []
    for ref_name, out_name in COLUMNS.items():
        got = out[out_name][:, 0].astype(np.float64)
        want = np.array([float(row[ref_name]) for row in reference])
        relative = np.abs(got - want) / np.maximum(np.abs(want), 1e-9)
        relative = np.where(np.abs(want) < 1e-8, 0.0, relative)
        if relative.max() >= TOLERANCE:
            failures.append("{} {:.3e}".format(ref_name, relative.max()))
    assert not failures, "worst relative differences: " + ", ".join(failures)


def test_total_harvest_matches_the_scenario_summary():
    """The five cuts remove the biomass the workbook reports.

    hb is the cumulative harvest, so its last value is the season total.
    """
    params, forcing = _inputs()
    out = grassim.run(params, forcing, 1)
    harvested = out["hb"][:, 0].astype(np.float64)
    assert np.all(np.diff(harvested) >= -1e-6), "hb must be cumulative"
    summary = read_sheet(WORKBOOK, "Scenario_summary")[0]
    want = float(summary["Total_harvested_biomass_kg_ha"])
    got = float(harvested[-1])
    assert abs(got - want) / want < 1e-3, (got, want)


def test_cut_height_is_a_length_not_a_flag():
    """gcut is a height: cutting to 1 m takes almost nothing off this sward.

    The 2019 sward never reaches 1 m, so a 1 m cut leaves it standing and
    the season's cumulative harvest collapses from 12.7 t/ha to under
    0.5 t/ha. Passing 1.0 as if it were a flag is therefore silently wrong,
    which is why t.graminr takes the height through its cut_height option.
    """
    params, forcing = _inputs()
    out_low = grassim.run(params, forcing, 1)
    high = np.where(forcing["gcut"] > 0, 1.0, 0.0)
    out_high = grassim.run(params, dict(forcing, gcut=high), 1)
    assert float(out_low["hb"][-1, 0]) > 12000.0
    assert float(out_high["hb"][-1, 0]) < 500.0
