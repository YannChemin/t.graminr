"""
MODULE:     t.graminr (graminr_lib.cutdates)

AUTHOR(S):  Yann Chemin

PURPOSE:    Cut and harvest calendars for t.graminr, ported from GraMinR's
            resolve_cut_dates()/resolve_cut_positions() (R/graminr_cutdates.R)
            and extended so that cut days can also come from a space-time
            raster dataset, which is how a remotely sensed mowing detection
            enters the model with a different date on every field.

SPDX-License-Identifier: Unlicense
"""

import datetime as dt

import numpy as np

from . import msg


def _read_raster(*args, **kwargs):
    """Deferred graminr_lib.gio.read_raster, which needs a GRASS session."""
    from .gio import read_raster

    return read_raster(*args, **kwargs)


def parse_cut_dates(text):
    """Split a cut_dates= string into day-of-year numbers and calendar dates.

    A bare number is a day of year and applies to every simulated year; a
    YYYY-MM-DD value applies to that one date only. The two can be mixed.
    """
    doys = []
    dates = []
    if not text:
        return doys, dates
    for item in str(text).split(","):
        item = item.strip()
        if not item:
            continue
        try:
            doys.append(int(item))
            continue
        except ValueError:
            pass
        try:
            dates.append(dt.datetime.strptime(item, "%Y-%m-%d").date())
        except ValueError:
            msg.fatal(
                msg.translate(
                    "cut_dates=: '{item}' is neither a day of year (1-366) nor "
                    "a date (YYYY-MM-DD)."
                ).format(item=item)
            )
    for doy in doys:
        if not 1 <= doy <= 366:
            msg.fatal(
                msg.translate("cut_dates=: the day of year {} is out of range.").format(
                    doy
                )
            )
    return doys, dates


class CutCalendar:
    """Per-day, per-pixel cut height (or cut flag) for a simulation.

    Three mutually exclusive sources, as in GraMinR:

      * an STRDS of cut flags or cut heights, one value per pixel per day;
      * a cut_dates= list, the same for every pixel;
      * nothing, leaving the model to its own internal cutting rule (LINGRA
        with imopt=1 cuts on a biomass threshold, for instance).
    """

    def __init__(
        self, days, mask_index, strds=None, cut_dates=None, cut_value=1.0, gap_fill=0
    ):
        self.days = days
        self.mask_index = mask_index
        self.cut_value = float(cut_value)
        self.index = None
        self._cache = {}
        self._schedule = None

        if strds and cut_dates:
            msg.fatal(
                msg.translate(
                    "Both a cut dataset (cut_strds=) and a cut calendar "
                    "(cut_dates=) were given; supply only one source of "
                    "cut/harvest information."
                )
            )
        if strds:
            # Only the per-pixel STRDS route touches GRASS, so the import is
            # deferred: the calendars themselves are plain logic and stay
            # usable, and testable, outside a GRASS session.
            from .gio import DayIndex

            self.index = DayIndex(strds, max_gap=gap_fill)
            msg.message(
                msg.translate(
                    "Cut/harvest days are taken from the dataset <{}>."
                ).format(strds)
            )
            return
        if cut_dates:
            doys, dates = parse_cut_dates(cut_dates)
            schedule = np.zeros(len(days))
            for position, day in enumerate(days):
                if day.timetuple().tm_yday in doys or day in dates:
                    schedule[position] = self.cut_value
            self._schedule = schedule
            flagged = int(np.count_nonzero(schedule))
            msg.message(
                msg.translate(
                    "Cut/harvest days are taken from cut_dates=, {count} day(s) "
                    "within the simulated period, at a value of {value}."
                ).format(count=flagged, value=self.cut_value)
            )
            if flagged == 0:
                msg.warning(
                    msg.translate(
                        "None of the cut_dates= values falls inside the "
                        "simulated period."
                    )
                )
            return
        msg.message(
            msg.translate(
                "No cut dataset and no cut calendar: the model will rely on its "
                "own cutting rule, where it has one."
            )
        )

    def values(self, position, day):
        """Cut height or flag of `day`, over the masked pixels."""
        size = self.mask_index.size
        if self.index is not None:
            name = self.index.map_for(day)
            if name is None:
                return np.zeros(size)
            if name not in self._cache:
                self._cache = {name: _read_raster(name).ravel()[self.mask_index]}
            data = self._cache[name]
            return np.where(np.isfinite(data), data, 0.0)
        if self._schedule is not None:
            return np.full(size, self._schedule[position])
        return np.zeros(size)


class EventCalendar:
    """A per-day amount from a list of year/doy/amount events.

    Used for the nitrogen fertilisation of BASGRA and GrasSim and for
    BASGRA's background nitrogen deposition.
    """

    def __init__(self, days, text, label="events", interpolate=False):
        self.days = days
        self.label = label
        self.amounts = np.zeros(len(days))
        if not text:
            return
        events = []
        for item in str(text).split(","):
            item = item.strip()
            if not item:
                continue
            parts = item.split("/")
            if len(parts) != 2:
                msg.fatal(
                    msg.translate(
                        "{label}=: '{item}' is not a date/amount pair such as "
                        "2024-05-10/60."
                    ).format(label=label, item=item)
                )
            try:
                when = dt.datetime.strptime(parts[0].strip(), "%Y-%m-%d").date()
            except ValueError:
                msg.fatal(
                    msg.translate(
                        "{label}=: '{date}' is not a date (YYYY-MM-DD)."
                    ).format(label=label, date=parts[0])
                )
            try:
                amount = float(parts[1])
            except ValueError:
                msg.fatal(
                    msg.translate("{label}=: '{amount}' is not a number.").format(
                        label=label, amount=parts[1]
                    )
                )
            events.append((when, amount))

        positions = {day: index for index, day in enumerate(days)}
        if interpolate:
            # A background rate: interpolate linearly between the points.
            events.sort()
            if len(events) == 1:
                self.amounts[:] = events[0][1]
            else:
                xs = [(when - days[0]).days for when, _ in events]
                ys = [amount for _, amount in events]
                self.amounts = np.interp([(day - days[0]).days for day in days], xs, ys)
            return
        for when, amount in events:
            if when not in positions:
                msg.warning(
                    msg.translate(
                        "{label}=: {date} lies outside the simulated period and "
                        "is ignored."
                    ).format(label=label, date=when)
                )
                continue
            self.amounts[positions[when]] += amount

    def values(self, position, size):
        return np.full(size, self.amounts[position])
