# protocols.py

import numpy as np


def allnighter_protocol(
    baseline_days: int = 30,
    recovery_days: int = 14,
    dt: float = 0.1,
    day_start: float = 7.0,
    day_end: float = 23.0,
    day_lux: float = 250.0,
    night_lux: float = 1.0,
    allnighter_lux: float = 1000.0,
    allnighter_start: float = 22.0,
    allnighter_end: float = 6.0,
):
    """
    Build the common 'baseline + all-nighter + recovery' light schedule.

    Layout:
      - days 0 .. baseline_days-1: normal LD (day_lux / night_lux)
      - night beginning on baseline_days (i.e. baseline_days, 22:00) is an all-nighter
        with allnighter_lux until baseline_days+1 at allnighter_end
      - days baseline_days+1 .. baseline_days+recovery_days: normal LD again

    Returns
    -------
    t : np.ndarray
        Time array in hours.
    lux_all : np.ndarray
        Light intensity (lux) for the entire protocol.
    baseline_days : int
        Number of baseline days (just passed through for convenience).
    recovery_start_day : int
        Day index where recovery begins (first full day after the all-nighter).
    """
    total_days = baseline_days + 1 + recovery_days
    t = np.arange(0.0, 24.0 * total_days, dt)

    # start with all times at night_lux
    lux = np.full_like(t, night_lux, dtype=float)

    for i, ti in enumerate(t):
        day = int(ti // 24.0)
        hour = ti % 24.0

        # ----- baseline days -----
        if day < baseline_days:
            if day_start <= hour < day_end:
                lux[i] = day_lux

        # ----- day with all-nighter -----
        elif day == baseline_days:
            # normal daytime light
            if day_start <= hour < day_end:
                lux[i] = day_lux

            # all-nighter: from allnighter_start on this day
            # until allnighter_end next morning
            if hour >= allnighter_start:
                lux[i] = allnighter_lux

        # ----- day immediately after all-nighter -----
        elif day == baseline_days + 1:
            # allnighter continues until allnighter_end
            if hour < allnighter_end:
                lux[i] = allnighter_lux
            elif day_start <= hour < day_end:
                lux[i] = day_lux
            # otherwise remains at night_lux (early morning / late night)

        # ----- subsequent recovery days (normal LD) -----
        else:
            if day_start <= hour < day_end:
                lux[i] = day_lux
            # otherwise keep night_lux

    recovery_start_day = baseline_days + 1  # first full day after the all-nighter
    return t, lux, baseline_days, recovery_start_day


def branch_ld_after(t, lux_all, recovery_start_day: int):
    """
    LD-recovery branch: just use the original schedule (baseline + all-nighter + LD recovery).
    Provided for symmetry & future flexibility.
    """
    # For now this is a no-op; we return a copy to avoid accidental in-place edits.
    return lux_all.copy()


def branch_dark_after(t, lux_all, recovery_start_day: int, dark_lux: float = 0.0):
    """
    Dark-recovery branch: from recovery_start_day onward, hold light at dark_lux
    for all hours (constant darkness or very dim light).

    Parameters
    ----------
    t : np.ndarray
        Time array in hours.
    lux_all : np.ndarray
        Original lux schedule (baseline + all-nighter + LD recovery).
    recovery_start_day : int
        Day index where constant-darkness recovery begins.
    dark_lux : float
        Light level to use in 'darkness' (0 or 1 lux).

    Returns
    -------
    lux_dark : np.ndarray
        Modified lux schedule for the darkness condition.
    """
    lux_dark = lux_all.copy()
    start_t = recovery_start_day * 24.0
    mask = t >= start_t
    lux_dark[mask] = dark_lux
    return lux_dark
