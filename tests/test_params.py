"""
MODULE:     t.graminr (tests)

AUTHOR(S):  Yann Chemin

PURPOSE:    Tests for the shipped default parameter sets and for the
            parameter casting, which turns the CSV's text into the scalars,
            vectors and tables each model expects.

            These import no GRASS module, so they run outside a session.

SPDX-License-Identifier: Unlicense
"""

import numpy as np
import pytest

from graminr_lib import basgra, grasprog, grassim, lingra, modvege
from graminr_lib import params as P
from graminr_lib.msg import GraminrError


def _raw(model):
    """The shipped default set as the raw text mapping the CSV would give."""
    return {
        name: ",".join(str(item) for item in value)
        if isinstance(value, tuple)
        else str(value)
        for name, value in P.read_parameter_csv(P.default_parameter_file(model)).items()
    }


MODELS = {
    "lingra": lingra,
    "modvege": modvege,
    "grasprog": grasprog,
    "grassim": grassim,
    "basgra": basgra,
}


@pytest.mark.parametrize("model", sorted(MODELS))
def test_every_model_ships_a_loadable_default_set(model):
    """params/<model>_params.csv must load and cover the model's defaults."""
    values = P.load(model)
    assert values, "no parameter was loaded for {}".format(model)
    missing = sorted(set(MODELS[model].DEFAULTS) - set(values))
    assert not missing, "{} left out of the default set: {}".format(model, missing)


@pytest.mark.parametrize("model", sorted(MODELS))
def test_the_default_set_runs_the_model_for_a_few_days(model):
    """A default set must be complete enough to integrate, not just parse."""
    ndays, npixels = 10, 3
    values = P.load(model)
    doy = np.arange(100, 100 + ndays)
    ones = np.ones((ndays, npixels))
    weather = {
        "year": np.full(ndays, 2024),
        "doy": doy,
        "tavg": ones * 12.0,
        "tmin": ones * 6.0,
        "tmax": ones * 18.0,
        "rain": ones * 2.0,
        "rad_glob": ones * 15.0,
        "rad_par": ones * 7.5,
        "et_pot": ones * 3.0,
        "vp": ones * 1.0,
        "wind": ones * 2.5,
        "lat": np.full(npixels, 45.0),
        "elev": np.full(npixels, 100.0),
        "gcut": np.zeros((ndays, npixels)),
        "rirr": np.zeros((ndays, npixels)),
    }
    out = MODELS[model].run(values, weather, npixels)
    assert out, "{} returned no output".format(model)
    for name, series in out.items():
        assert series.shape == (ndays, npixels), "{}/{} has shape {}".format(
            model, name, series.shape
        )
        assert np.all(np.isfinite(series)), "{}/{} went non-finite".format(model, name)


def test_a_scalar_is_cast_to_a_float():
    """The CSV holds text; the model needs numbers."""
    raw = dict(_raw("lingra"), LAT="52.12")
    values = P.cast_parameters(raw, "lingra")
    assert isinstance(values["LAT"], float)
    assert values["LAT"] == pytest.approx(52.12)


def test_a_comma_separated_value_is_cast_to_a_tuple():
    """The per-layer and per-class vectors are written as a comma list."""
    raw = dict(_raw("basgra"), FGAS="0.1, 0.2, 0.3, 0.4")
    values = P.cast_parameters(raw, "basgra")
    assert np.allclose(np.asarray(values["FGAS"], dtype=float), [0.1, 0.2, 0.3, 0.4])


def test_a_one_element_list_collapses_to_a_scalar():
    raw = dict(_raw("basgra"), FGAS="0.25")
    assert P.cast_parameters(raw, "basgra")["FGAS"] == pytest.approx(0.25)


def test_a_non_numeric_value_is_kept_as_text_not_coerced():
    """ETMODEL and the like are names, not numbers."""
    raw = dict(_raw("lingra"), ETMODEL="PM")
    assert P.cast_parameters(raw, "lingra")["ETMODEL"] == "PM"


def test_a_missing_required_parameter_is_refused():
    """An incomplete set must fail loudly, not run on silent defaults."""
    with pytest.raises(GraminrError):
        P.cast_parameters({"LAT": "52.12"}, "lingra")


def test_an_unknown_model_has_no_default_file():
    with pytest.raises((GraminrError, KeyError)):
        P.load("pasim")


def test_every_required_forcing_can_actually_be_supplied():
    """A model must not require a forcing t.graminr has no option for.

    preflight() rejects a climate=strds request whose REQUIRED forcings are
    not all given. If a model required something absent from
    FORCING_OPTIONS, that check could never be satisfied and the model
    would be unusable on prepared datasets. The module is read as text
    rather than imported, because importing it needs a GRASS session.
    """
    import os
    import re

    source = os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "t.graminr.py"
    )
    with open(source) as handle:
        text = handle.read()
    block = re.search(r"FORCING_OPTIONS = \((.*?)\)", text, re.DOTALL)
    assert block, "FORCING_OPTIONS is no longer a literal tuple in t.graminr.py"
    options = set(re.findall(r'"([a-z_0-9]+)"', block.group(1)))
    derivable = set(re.findall(r'"([a-z_0-9]+)": "([a-z_0-9]+)"', text))

    for name, module in sorted(MODELS.items()):
        for forcing in module.REQUIRED:
            sources = {forcing} | {
                target for source, target in derivable if source == forcing
            }
            assert sources & options, (
                "the '{}' model requires '{}', which no t.graminr option "
                "supplies".format(name, forcing)
            )
