"""
MODULE:     t.graminr (tests)

AUTHOR(S):  Yann Chemin

PURPOSE:    Numerical validation of the LINGRA port against the PCSE
            reference run shipped with GraMinR
            (inst/extdata/original_outputs/lingra/), Zegveld 1986.

            The whole 366-day run is compared, state variable by state
            variable, against pcse.crop.lingra. Two findings are pinned
            here, both established by this comparison:

            1. GraMinR's reference parameter file carries ANGSTA/ANGSTB =
               0.22/0.47, while the PCSE run used the PCSE/WOFOST defaults
               0.18/0.55. The pair enters only the Penman net long-wave
               term, so it moves E0, ES0 and through them EVS, and leaves
               the Penman-Monteith ET0 untouched -- which is exactly what
               is observed. ET0 matches to 6e-8 with the file's own values.

            2. GraMinR suspends photosynthesis on every scheduled cutting
               day, whereas PCSE suspends it only when biomass is actually
               removed. The two agree until a scheduled cut lands on a day
               when green leaf weight is already at or below cwlvg: nothing
               is harvested, yet GraMinR still loses that day's
               assimilation. In this run that happens on day of year 240
               and costs about 11% of green leaf weight for the rest of the
               season. The port's HarvestGate parameter selects between the
               two, and defaults to PCSE's rule.

SPDX-License-Identifier: Unlicense
"""

import csv
import os

import numpy as np
import pytest

from graminr_lib import lingra
from graminr_lib.params import cast_parameters

GRAMINR = os.environ.get("GRAMINR_DIR", os.path.expanduser("~/dev/graminr"))
REFDIR = os.path.join(GRAMINR, "inst", "extdata", "original_outputs", "lingra")

#: The PCSE run's mowing calendar (the agromanagement of
#: test_LINGRA_Netherlands-Zegveld-1986_WLP.yml; the same days as GraMinR's
#: own inst/extdata/lingra_cut_dates.csv).
CUT_DOYS = {135, 165, 200, 240, 280}

#: Reference column -> ported output name.
COLUMNS = {
    "LAI": "lai",
    "TillerNumber": "tiller_number",
    "WeightLVgreen": "weight_lv_green",
    "WeightRE": "weight_re",
    "WeightRT": "weight_rt",
    "WeightHARV": "weight_harv",
    "WeightABG": "weight_abg",
    "TRA": "tra",
    "TRAMX": "tramx",
    "LeafLength": "leaf_length",
    "SLAINT": "slaint",
    "SM": "sm",
    "EVS": "evs",
    "RFTRA": "rftra",
}

pytestmark = pytest.mark.skipif(
    not os.path.isdir(REFDIR),
    reason="the GraMinR reference outputs are not available; set GRAMINR_DIR",
)


def _inputs(harvest_gate):
    raw = {}
    with open(os.path.join(REFDIR, "lingra_zegveld_params.csv"), newline="") as handle:
        for row in csv.reader(handle):
            if len(row) >= 2 and row[0].strip():
                raw[row[0].strip()] = row[1].strip()
    raw["ANGSTA"], raw["ANGSTB"] = "0.18", "0.55"
    raw["HarvestGate"] = harvest_gate
    params = cast_parameters(raw, "lingra")

    with open(
        os.path.join(REFDIR, "lingra_zegveld_weather_1986.csv"), newline=""
    ) as handle:
        weather = list(csv.DictReader(handle))
    ndays = len(weather)

    def column(key):
        return np.array([[float(row[key])] for row in weather])

    forcing = {
        "year": np.array([int(row["year"]) for row in weather]),
        "doy": np.array([int(row["doy"]) for row in weather]),
        "rad_glob": column("rad_glob"),
        "tmin": column("tmin"),
        "tmax": column("tmax"),
        "vp": column("vp"),
        "wind": column("wind"),
        "rain": column("rain"),
        "gcut": np.array(
            [[1.0 if int(row["doy"]) in CUT_DOYS else 0.0] for row in weather]
        ),
        "rirr": np.zeros((ndays, 1)),
    }
    return params, forcing


def _reference():
    with open(
        os.path.join(REFDIR, "lingra_zegveld_pcse_output.csv"), newline=""
    ) as handle:
        return list(csv.DictReader(handle))


def _worst_relative(out, reference):
    worst, where = 0.0, None
    for ref_name, out_name in COLUMNS.items():
        got = out[out_name][:, 0].astype(np.float64)
        want = np.array([float(row[ref_name]) for row in reference])
        relative = np.abs(got - want) / np.maximum(np.abs(want), 1e-9)
        relative = np.where(np.abs(want) < 1e-8, 0.0, relative)
        if relative.max() > worst:
            worst, where = float(relative.max()), ref_name
    return worst, where


def test_matches_pcse_to_single_precision():
    """With PCSE's harvest gating, every state variable matches PCSE."""
    params, forcing = _inputs("harvest")
    out = lingra.run(params, forcing, 1)
    worst, where = _worst_relative(out, _reference())
    # The outputs are stored as float32, so about 1e-7 is the floor.
    assert worst < 1e-6, "worst relative difference {:.3e} on {}".format(worst, where)


def test_reference_et0_is_exact_with_the_files_own_angstrom_pair():
    """Penman-Monteith ET0 does not depend on the Angstrom coefficients."""
    with open(
        os.path.join(REFDIR, "lingra_zegveld_weather_1986_all_variables.csv"),
        newline="",
    ) as handle:
        weather = list(csv.DictReader(handle))
    got = []
    for row in weather:
        _e0, _es0, et0 = lingra.reference_et(
            float(row["doy"]),
            52.12,
            5.0,
            float(row["tmin"]),
            float(row["tmax"]),
            float(row["rad_glob"]) * 1e6,
            float(row["vp"]) * 10.0,
            float(row["wind"]),
            0.22,
            0.47,
            "PM",
        )
        got.append(float(et0))
    want = np.array([float(row["et0_pcse"]) for row in weather])
    assert np.max(np.abs(np.array(got) - want)) < 1e-7


def test_mowflag_gating_reproduces_the_graminr_deviation():
    """GraMinR's own gating differs from PCSE, and only after day 240."""
    params, forcing = _inputs("mowflag")
    out = lingra.run(params, forcing, 1)
    reference = _reference()
    got = out["weight_lv_green"][:, 0].astype(np.float64)
    want = np.array([float(row["WeightLVgreen"]) for row in reference])
    relative = np.abs(got - want) / np.maximum(np.abs(want), 1e-9)
    # Identical up to and including the cut of day of year 240 (index 239),
    # then divergent: that day's assimilation is lost although the cut
    # removed nothing. The deficit is immediate and large, and it never
    # closes over the following weeks, though its sign can flip later on as
    # the regrowth cycles drift apart.
    assert relative[:240].max() < 1e-6
    assert relative[240] > 0.1
    assert np.all(got[241:261] < want[241:261])
