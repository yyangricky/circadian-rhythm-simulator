

import math
import numpy as np


# -------- Phase-marker helpers --------

def get_phase_markers(model, trajectory, marker="cbt"):
    """
    Extract daily phase-marker times (in hours) from a model trajectory.

    Parameters
    ----------
    model : circadian model instance
        Must expose model.cbt(trajectory) and/or model.dlmos(trajectory).
    trajectory : DynamicalTrajectory
        Output of model(time, initial_state, light_values).
    marker : {"cbt", "dlmo"}
        Which phase marker to use.

    Returns
    -------
    np.ndarray
        1D array of absolute times (in hours) at which the marker occurs.
    """
    marker = marker.lower()
    try:
        if marker == "cbt":
            times = model.cbt(trajectory)
        elif marker == "dlmo":
            times = model.dlmos(trajectory)
        else:
            raise ValueError("marker must be 'cbt' or 'dlmo'")
        return np.array(times) if times is not None else np.array([])
    except Exception as e:
        print(f"[WARNING] Error extracting {marker.upper()} markers: {e}")
        return np.array([])


def circular_mean_hours(hours):
    """
    Circular mean of clock times on a 0–24 h scale.

    Parameters
    ----------
    hours : array-like
        Clock hours in [0, 24).

    Returns
    -------
    float or None
        Mean clock hour in [0, 24), or None if input is empty.
    """
    if len(hours) == 0:
        return None

    ang = np.array(hours) * (2 * np.pi / 24.0)
    mean_sin = np.mean(np.sin(ang))
    mean_cos = np.mean(np.cos(ang))
    mean_ang = math.atan2(mean_sin, mean_cos)

    if mean_ang < 0:
        mean_ang += 2 * np.pi

    return (mean_ang * 24.0) / (2 * np.pi)


def compute_baseline_marker(marker_times, baseline_days=25):
    """
    Baseline phase marker (CBTmin, DLMO, etc.) from last baseline_days.

    Parameters
    ----------
    marker_times : array-like
        Absolute times of phase markers (in hours) over the whole experiment.
    baseline_days : int
        Only markers strictly before baseline_days * 24 h are used.

    Returns
    -------
    float or None
        Baseline clock hour in [0, 24), or None if too few markers.
    """
    if len(marker_times) == 0:
        return None

    base = [t for t in marker_times if t < baseline_days * 24.0]
    if len(base) < 3:
        print(f"[WARNING] Only {len(base)} baseline markers found, need at least 3")
        return None

    hours = np.array(base) % 24.0
    return circular_mean_hours(hours)


def circular_distance(h1, h2):
    """
    Shortest signed distance between two clock hours on a 24 h circle.

    Parameters
    ----------
    h1, h2 : float
        Clock hours in [0, 24).

    Returns
    -------
    float
        Difference in hours in [-12, 12].
    """
    diff = h1 - h2
    while diff > 12:
        diff -= 24
    while diff < -12:
        diff += 24
    return diff


def reentrainment_hours(marker_times, baseline_marker_h,
                        start_day=31, tol_min=15, streak=3):
    """
    Time to re-entrainment, measured from the start of the recovery phase.

    Parameters
    ----------
    marker_times : array-like
        Absolute times (hours) of daily markers after the all-nighter.
    baseline_marker_h : float
         Baseline marker clock time (0–24 h) from compute_baseline_marker.
    start_day : int
        First day that belongs to the recovery phase.
    tol_min : float
        Tolerance in minutes for 'within baseline'.
    streak : int
        Number of consecutive days within tolerance required.

    Returns
    -------
    float or None
        Hours from start_day*24 until first day of the successful streak,
        or None if re-entrainment not reached.
    """
    if len(marker_times) == 0 or baseline_marker_h is None:
        return None

    tol_h = tol_min / 60.0
    start_recovery_time = start_day * 24.0

    # markers after recovery starts
    recovery_markers = [t for t in marker_times if t >= start_recovery_time]
    if len(recovery_markers) == 0:
        return None

    # group by day
    day_markers = {}
    for t in recovery_markers:
        day = int(t // 24)
        day_markers.setdefault(day, []).append(t)

    sorted_days = sorted(day_markers.keys())

    for i, d in enumerate(sorted_days):
        todays = day_markers[d]
        mean_h = circular_mean_hours([ti % 24.0 for ti in todays])
        if mean_h is None:
            continue

        err = circular_distance(mean_h, baseline_marker_h)
        if abs(err) > tol_h:
            continue

        # check streak of consecutive days
        if i + streak > len(sorted_days):
            continue

        good = True
        for k in range(1, streak):
            next_idx = i + k
            if next_idx >= len(sorted_days):
                good = False
                break

            d2 = sorted_days[next_idx]
            if d2 != d + k:
                good = False
                break

            hs2 = [ti % 24.0 for ti in day_markers[d2]]
            mean_h2 = circular_mean_hours(hs2)
            if mean_h2 is None:
                good = False
                break

            err2 = circular_distance(mean_h2, baseline_marker_h)
            if abs(err2) > tol_h:
                good = False
                break

        if good:
            first_marker_time = np.mean(day_markers[d])
            return first_marker_time - start_recovery_time

    return None


# -------- Amplitude helpers --------

def amp_series(model, trajectory):
    """
    Extract an amplitude proxy R(t) from a trajectory.

    Uses:
    - R = state[0] for Hannay-style models, OR
    - R = sqrt(x^2 + y^2) if state[0] looks non-amplitude-like, OR
    - |state[0]| as a fallback.

    Parameters
    ----------
    model : model instance (unused but kept for flexibility)
    trajectory : DynamicalTrajectory

    Returns
    -------
    np.ndarray
        Non-negative amplitude series R(t) of same length as trajectory.time.
    """
    try:
        states = np.asarray(trajectory.states)
        if states.size == 0:
            return np.array([])

        R = states[:, 0]

        # heuristic sanity check
        if np.nanmax(np.abs(R)) > 5.0 or np.nanmin(R) < -2.0:
            if states.shape[1] >= 2:
                R = np.sqrt(states[:, 0] ** 2 + states[:, 1] ** 2)
            else:
                R = np.abs(states[:, 0])

        return np.abs(R)
    except Exception as e:
        print(f"[WARNING] Error extracting amplitude: {e}")
        return np.array([])


def hours_to_90pct(R, t, baseline_days=30, start_day=31):
    """
    Time from start of recovery until amplitude first reaches 90% baseline.

    Parameters
    ----------
    R : array-like
        Amplitude series.
    t : array-like
        Time array corresponding to R (hours).
    baseline_days : int
        Days used to compute baseline amplitude (t < baseline_days*24).
    start_day : int
        First day of recovery.

    Returns
    -------
    float or None
        Hours after start_day*24 until R(t) >= 0.9 * baseline_mean,
        or None if not reached.
    """
    R = np.asarray(R)
    t = np.asarray(t)
    if R.size == 0 or t.size == 0:
        return None

    base_mask = t < baseline_days * 24.0
    if not np.any(base_mask):
        return None

    R_base = np.mean(R[base_mask])
    if R_base <= 0:
        return None

    target = 0.9 * R_base
    start_time = start_day * 24.0

    for Ri, ti in zip(R, t):
        if ti < start_time:
            continue
        if Ri >= target:
            return ti - start_time

    return None
