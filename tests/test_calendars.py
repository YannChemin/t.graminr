"""
MODULE:     t.graminr (tests)

AUTHOR(S):  Yann Chemin

PURPOSE:    Tests for the cutting and fertilisation calendars, which parse
            the cut_dates=, fert_dates=, fert_org_dates= and ndep_dates=
            options into per-day arrays.

            They import no GRASS module, so they run outside a session;
            graminr_lib.msg raises GraminrError where the module would call
            G_fatal_error().

SPDX-License-Identifier: Unlicense
"""

import datetime as dt

import numpy as np
import pytest

from graminr_lib import cutdates
from graminr_lib.msg import GraminrError

DAYS = [dt.date(2024, 4, 1) + dt.timedelta(days=offset) for offset in range(90)]


def test_parse_cut_dates_separates_days_of_year_from_dates():
    doys, dates = cutdates.parse_cut_dates("135,2024-05-01,200")
    assert doys == [135, 200]
    assert dates == [dt.date(2024, 5, 1)]


@pytest.mark.parametrize("text", ["0", "367", "springtime", "2024-13-01"])
def test_parse_cut_dates_refuses_nonsense(text):
    with pytest.raises(GraminrError):
        cutdates.parse_cut_dates(text)


def test_cut_calendar_repeats_a_day_of_year_every_year():
    days = [dt.date(2023, 1, 1) + dt.timedelta(days=offset) for offset in range(800)]
    calendar = cutdates.CutCalendar(days, np.arange(1), cut_dates="135", cut_value=1.0)
    flagged = [day for day, value in zip(days, calendar._schedule) if value]
    assert [day.year for day in flagged] == [2023, 2024]
    assert {day.timetuple().tm_yday for day in flagged} == {135}


def test_cut_calendar_carries_the_cut_value_not_just_a_flag():
    """ModVege and GrasSim need the height here, not a 1."""
    calendar = cutdates.CutCalendar(DAYS, np.arange(1), cut_dates="135", cut_value=0.05)
    assert float(np.max(calendar._schedule)) == pytest.approx(0.05)


def test_cut_calendar_refuses_two_sources():
    with pytest.raises(GraminrError):
        cutdates.CutCalendar(DAYS, np.arange(1), strds="some_strds", cut_dates="135")


def test_cut_calendar_without_a_source_schedules_nothing():
    calendar = cutdates.CutCalendar(DAYS, np.arange(1))
    assert calendar.index is None
    assert calendar._schedule is None


def test_event_calendar_places_amounts_on_their_days():
    calendar = cutdates.EventCalendar(DAYS, "2024-04-10/60,2024-05-25/40", "fert_dates")
    placed = {
        DAYS[index]: float(calendar.amounts[index])
        for index in np.nonzero(calendar.amounts)[0]
    }
    assert placed == {dt.date(2024, 4, 10): 60.0, dt.date(2024, 5, 25): 40.0}


def test_event_calendar_interpolates_when_asked():
    """Nitrogen deposition is a rate, so it is carried between the points."""
    calendar = cutdates.EventCalendar(
        DAYS, "2024-04-01/0.004,2024-04-11/0.014", "ndep_dates", interpolate=True
    )
    amounts = calendar.amounts
    assert amounts[0] == pytest.approx(0.004)
    assert amounts[10] == pytest.approx(0.014)
    assert amounts[5] == pytest.approx(0.009, abs=1e-9)
    assert np.all(np.diff(amounts[:11]) > 0)


def test_event_calendar_without_interpolation_leaves_zeros_between():
    calendar = cutdates.EventCalendar(DAYS, "2024-04-01/10,2024-04-11/10", "fert_dates")
    assert calendar.amounts[5] == 0.0


@pytest.mark.parametrize(
    "text", ["2024-04-10", "2024-04-10/", "notadate/60", "2024-04-10/lots"]
)
def test_event_calendar_refuses_malformed_pairs(text):
    with pytest.raises(GraminrError):
        cutdates.EventCalendar(DAYS, text, "fert_dates")


def test_event_calendar_empty_text_is_all_zeros():
    calendar = cutdates.EventCalendar(DAYS, "", "fert_dates")
    assert calendar.amounts.shape == (len(DAYS),)
    assert not np.any(calendar.amounts)
