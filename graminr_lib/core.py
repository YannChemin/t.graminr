"""
MODULE:     t.graminr (graminr_lib.core)

AUTHOR(S):  Yann Chemin

PURPOSE:    Shared numeric helpers for the vectorised GraMinR grassland
            models. Every model in this package works on NumPy arrays whose
            single axis is the pixel, so each scalar branch of the original
            R code becomes a numpy.where() over that axis.

SPDX-License-Identifier: Unlicense
"""

import numpy as np


def afgen(table, x):
    """Piecewise-linear interpolation of an AFGEN-style table.

    :param table: flat sequence x1,y1,x2,y2,... with strictly increasing x
    :param x: scalar or array of abscissae
    :return: interpolated ordinate(s), clamped to the table's end values
    """
    tbl = np.asarray(table, dtype=np.float64).ravel()
    if tbl.size < 4 or tbl.size % 2 != 0:
        raise ValueError(
            "afgen: table must hold an even number of values (at least 4), "
            "as x1,y1,x2,y2,..."
        )
    xs = tbl[0::2]
    ys = tbl[1::2]
    if np.any(np.diff(xs) <= 0):
        raise ValueError(
            "afgen: the x values of the table must be strictly increasing: "
            + ", ".join(str(v) for v in xs)
        )
    # numpy.interp clamps outside the table, which is what AFGEN does.
    return np.interp(x, xs, ys)


def clamp(vmin, vmax, value):
    """Clamp value into [vmin, vmax], in the argument order used by LINGRA."""
    lo = np.asarray(vmin, dtype=np.float64)
    hi = np.asarray(vmax, dtype=np.float64)
    if np.any(lo > hi):
        raise ValueError("clamp: vmin is greater than vmax")
    return np.minimum(hi, np.maximum(lo, value))


def limit01(x):
    """Clip to [0, 1]."""
    return np.clip(x, 0.0, 1.0)


def safe_div(num, den, fill=0.0):
    """Divide, substituting fill where the denominator is zero or not finite.

    The original R code divides by state variables that can legitimately
    reach zero (bare soil, an empty carbon pool), where R yields NaN or Inf
    and carries it silently through every later day. Over a whole raster
    that would poison entire pixels, so the quotient is defined explicitly
    here instead.
    """
    num = np.asarray(num, dtype=np.float64)
    den = np.asarray(den, dtype=np.float64)
    out = np.full(np.broadcast(num, den).shape, float(fill), dtype=np.float64)
    good = np.isfinite(den) & (den != 0)
    np.divide(num, den, out=out, where=good)
    out[~np.isfinite(out)] = float(fill)
    return out


def as_array(value, n, dtype=np.float64):
    """Broadcast a scalar or sequence to a length-n array of pixels."""
    arr = np.asarray(value, dtype=dtype)
    if arr.ndim == 0:
        return np.full(n, arr, dtype=dtype)
    if arr.shape == (n,):
        return arr.astype(dtype, copy=True)
    if arr.size == 1:
        return np.full(n, arr.reshape(-1)[0], dtype=dtype)
    raise ValueError(
        "as_array: cannot broadcast a value of shape {} to {} pixels".format(
            arr.shape, n
        )
    )


def year_cumsum(values, year):
    """Cumulative sum along the day axis, restarting on each calendar year.

    :param values: (ndays, npixels) array of daily contributions
    :param year: (ndays,) array of calendar years
    :return: (ndays, npixels) running sum, reset at every change of year
    """
    values = np.asarray(values, dtype=np.float64)
    year = np.asarray(year)
    out = np.empty_like(values)
    acc = np.zeros(values.shape[1], dtype=np.float64)
    previous = None
    for i in range(values.shape[0]):
        if previous is not None and year[i] != previous:
            acc = np.zeros_like(acc)
        acc = acc + values[i]
        out[i] = acc
        previous = year[i]
    return out


def trailing_mean(values, window):
    """Causal trailing mean over the day axis, NaN until the window is full.

    Matches the ModVege spreadsheet's AVERAGE(day-9:day): it never looks at
    a future day, and is undefined before `window` days of history exist.
    """
    values = np.asarray(values, dtype=np.float64)
    ndays = values.shape[0]
    out = np.full_like(values, np.nan)
    if ndays < window:
        return out
    csum = np.cumsum(values, axis=0)
    out[window - 1] = csum[window - 1] / window
    if ndays > window:
        out[window:] = (csum[window:] - csum[:-window]) / window
    return out


def day_of_year_length(year):
    """Number of days in a calendar year."""
    year = int(year)
    leap = (year % 4 == 0 and year % 100 != 0) or year % 400 == 0
    return 366 if leap else 365
