import numpy as np
'''
__all__ = ["recovery_time"]

def recovery_time(dlmo, baseline, start_day=25, threshold=0.10):
    """
    Compute the first day where the DLMO returns within a threshold of the baseline.

    Parameters
    ----------
    dlmo : array-like
        Array of daily DLMO values (in hours).
    baseline : float
        Baseline DLMO value before perturbation.
    start_day : int
        Day index at which recovery search begins.
    threshold : float
        Fractional threshold (default = 0.10 for 10%).

    Returns
    -------
    int or None
        The day index where recovery occurs, or None if never recovered.
    """
    dlmo = np.asarray(dlmo)
    deviation = np.abs(dlmo[start_day:] - baseline)
    allowed = threshold * baseline

    recovery_days = np.where(deviation <= allowed)[0]

    if len(recovery_days) == 0:
        return None

    return start_day + recovery_days[0]

'''
from metrics import reentrainment_hours

def recovery_time(marker_times, baseline_marker, start_day=25, tol_min=15, streak=3):
    hours = reentrainment_hours(marker_times, baseline_marker, start_day, tol_min, streak)
    if hours is None:
        return None
    return start_day + hours / 24.0
