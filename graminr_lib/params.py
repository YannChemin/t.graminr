"""
MODULE:     t.graminr (graminr_lib.params)

AUTHOR(S):  Yann Chemin

PURPOSE:    Parameter handling for t.graminr: read the two-column CSV that
            GraMinR's read_config() reads, apply the same renames, derived
            quantities and defaults, and let any single parameter be given
            as a raster so that it varies across the region.

SPDX-License-Identifier: Unlicense
"""

import csv
import os

import numpy as np

from . import basgra, grassim, grasprog, lingra, modvege
from . import msg

#: The models this addon ports. PaSim is deliberately absent: it wraps an
#: external executable that is not distributed with GraMinR.
MODELS = {
    "lingra": lingra,
    "modvege": modvege,
    "basgra": basgra,
    "grasprog": grasprog,
    "grassim": grassim,
}

#: Directory holding the default parameter sets, beside this package.
PARAM_DIR = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "params"
)


def default_parameter_file(model):
    """Path of the parameter set shipped for `model`."""
    path = os.path.join(PARAM_DIR, "{}_params.csv".format(model))
    if not os.path.exists(path):
        # Installed layout: $GISBASE/etc/t.graminr/params/
        etc = os.environ.get("GISBASE", "")
        alternative = os.path.join(
            etc, "etc", "t.graminr", "params", "{}_params.csv".format(model)
        )
        if os.path.exists(alternative):
            return alternative
        msg.fatal(
            msg.translate(
                "No parameter set is installed for the '{model}' model; pass "
                "one with parameters=."
            ).format(model=model)
        )
    return path


def read_parameter_csv(path):
    """Read a headerless two-column name,value parameter CSV.

    Values are kept as text here; the per-model casting below turns them
    into floats, into tables (a quoted, comma-separated cell such as
    LUEreductionSoilTempTB) or leaves them as strings (ETMODEL, ET_METHOD).
    """
    values = {}
    with open(path, newline="") as handle:
        for row in csv.reader(handle):
            cells = [cell.strip() for cell in row if cell.strip() != ""]
            if not cells or cells[0].startswith("#"):
                continue
            if len(cells) < 2:
                msg.fatal(
                    msg.translate(
                        "Parameter file <{path}>: the row '{row}' has no value."
                    ).format(path=path, row=",".join(row))
                )
            name = cells[0]
            if len(cells) == 2:
                values[name] = cells[1]
            else:
                # An unquoted multi-value table row, e.g. KDIFTB,0,0.6,2,0.6
                values[name] = ",".join(cells[1:])
    if not values:
        msg.fatal(msg.translate("Parameter file <{}> holds no parameter.").format(path))
    return values


def _to_number(text):
    try:
        return float(text)
    except (TypeError, ValueError):
        return None


def _to_table(text):
    parts = [part.strip() for part in str(text).split(",")]
    numbers = [_to_number(part) for part in parts]
    if any(number is None for number in numbers):
        return None
    return tuple(numbers)


def cast_parameters(raw, model):
    """Cast a raw name->text mapping into the model's parameter types."""
    module = MODELS[model]
    values = {}
    table_names = set(getattr(module, "TABLE_PARAMS", ()))
    text_names = {"ETMODEL", "ET_METHOD", "droughtindex", "waterstress"}

    for name, text in raw.items():
        if name in text_names:
            values[name] = str(text)
            continue
        if name in table_names or "," in str(text):
            table = _to_table(text)
            if table is None:
                values[name] = str(text)
            else:
                values[name] = table[0] if len(table) == 1 else table
            continue
        number = _to_number(text)
        values[name] = str(text) if number is None else number

    # Per-model renames (read_config()'s names(plist)[...] <- ... block).
    for source, target in getattr(module, "ALIASES", {}).items():
        if source in values and target not in values:
            values[target] = values.pop(source)
        elif source in values:
            values.pop(source)

    if model in ("lingra",):
        # SLA and sla are the same parameter under two names; supplying
        # either alone must work.
        if "SLA" in values and "sla" not in values:
            values["sla"] = values["SLA"]
        if "sla" in values and "SLA" not in values:
            values["SLA"] = values["sla"]

    # Defaults, then the derived quantities that depend on them.
    for name, value in getattr(module, "DEFAULTS", {}).items():
        values.setdefault(name, value)

    if model == "grasprog":
        _collect_compartments(values)
    if model == "grassim":
        _grassim_extras(values)
    if model == "lingra":
        # A parameter set may be given either in the legacy LINGRA names
        # (rootdi, sla, LAIi, ...) or directly in the PCSE ones (RDI, SLA,
        # LAIinit, ...). A derivation whose inputs are absent is simply
        # skipped, as in read_config(), where a missing name yields an empty
        # value rather than an error; what the model truly needs is checked
        # by REQUIRED_PARAMS below.
        for name, rule in lingra.DERIVED.items():
            if name in values:
                continue
            try:
                values[name] = rule(values)
            except (KeyError, TypeError):
                continue
    if model == "modvege":
        values["cellSurfaceMeter"] = 10000.0 * values["cellSurface"]
        method = str(values.get("ET_METHOD", "split"))
        if method not in ("split", "linear"):
            msg.fatal(
                msg.translate(
                    "ET_METHOD must be 'split' or 'linear', not '{}'."
                ).format(method)
            )

    _check_required(values, model)
    return values


def _collect_compartments(values):
    """Gather rs_1.., ak_1.., w0_1.. into per-compartment vectors."""
    for prefix in grasprog.COMPARTMENT_PREFIXES:
        found = {}
        for name in list(values):
            parts = name.split("_")
            if len(parts) == 2 and parts[0] == prefix and parts[1].isdigit():
                found[int(parts[1])] = float(values.pop(name))
        if found:
            values["{}_compartments".format(prefix)] = tuple(
                found[key] for key in sorted(found)
            )


def _grassim_extras(values):
    """Monthly Kc and the soil water limits of the GrasSim parameter set."""
    kc_names = ["Kc_{}".format(month) for month in range(1, 13)]
    present = [name for name in kc_names if name in values]
    if present:
        if len(present) != 12:
            msg.fatal(
                msg.translate(
                    "GrasSim monthly Kc must be given as the 12 rows "
                    "Kc_1..Kc_12; found {count}."
                ).format(count=len(present))
            )
        values["kc_monthly"] = tuple(float(values.pop(name)) for name in kc_names)

    hydro = ("waterCapacity", "waterSaturation", "Wiltingpoint")
    if not all(name in values for name in hydro):
        if not all(name in values for name in ("sand", "clay", "org")):
            msg.fatal(
                msg.translate(
                    "GrasSim needs either waterCapacity, waterSaturation and "
                    "Wiltingpoint, or the soil texture sand, clay and org."
                )
            )
        derived = grassim.soil_hydro(values["sand"], values["clay"], values["org"])
        values.update(derived)


def _check_required(values, model):
    """Fail loudly on a parameter the model cannot run without."""
    module = MODELS[model]
    required = getattr(module, "REQUIRED_PARAMS", None)
    if required is None:
        return
    missing = [name for name in required if name not in values]
    if missing:
        msg.fatal(
            msg.translate("The '{model}' parameter set is missing: {names}").format(
                model=model, names=", ".join(sorted(missing))
            )
        )


def apply_parameter_rasters(values, assignments, mask_index):
    """Replace named parameters with per-pixel values read from rasters.

    `assignments` is a sequence of "name=raster" strings, so that a soil
    property or a sward parameter can vary over the region while the rest
    of the set stays scalar.
    """
    from .gio import read_raster

    for assignment in assignments:
        if "=" not in assignment:
            msg.fatal(
                msg.translate(
                    "parameter_maps= expects name=raster pairs, not '{}'."
                ).format(assignment)
            )
        name, raster = assignment.split("=", 1)
        name, raster = name.strip(), raster.strip()
        if name not in values:
            msg.warning(
                msg.translate(
                    "Parameter '{name}' is not in the parameter set; its raster "
                    "<{raster}> is used all the same."
                ).format(name=name, raster=raster)
            )
        data = read_raster(raster).ravel()[mask_index]
        if np.any(~np.isfinite(data)):
            msg.warning(
                msg.translate(
                    "Raster <{raster}> is NULL on {count} grassland pixel(s); "
                    "the scalar value of '{name}' is used there."
                ).format(
                    raster=raster,
                    count=int(np.count_nonzero(~np.isfinite(data))),
                    name=name,
                )
            )
            scalar = values.get(name)
            if isinstance(scalar, (int, float)):
                data = np.where(np.isfinite(data), data, float(scalar))
        values[name] = data
        msg.verbose(
            msg.translate(
                "Parameter '{name}' now varies per pixel, from <{raster}>."
            ).format(name=name, raster=raster)
        )
    return values


def load(model, path=None, parameter_maps=(), mask_index=None, overrides=()):
    """Full parameter set of `model`, ready for its run() function."""
    path = path or default_parameter_file(model)
    values = cast_parameters(read_parameter_csv(path), model)
    for override in overrides:
        if "=" not in override:
            msg.fatal(
                msg.translate(
                    "parameter_values= expects name=value pairs, not '{}'."
                ).format(override)
            )
        name, text = override.split("=", 1)
        number = _to_number(text.strip())
        values[name.strip()] = text.strip() if number is None else number
    if parameter_maps and mask_index is not None:
        values = apply_parameter_rasters(values, parameter_maps, mask_index)
    return values
