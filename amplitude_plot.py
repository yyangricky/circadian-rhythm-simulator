import importlib
import numpy as np
import math
import matplotlib.pyplot as plt

# --- Simulation & Experiment Constants ---
MODEL_NAMES = [
    'Jewett99', 'Forger99', 'Hannay19', 'Hannay19TP', 'Hilaire07', 'Breslow13',
    'Skeldon23'
]

HOURS = 24.0 * (30 + 1 + 14)  # 30 baseline days + 1 disruption night + 14 recovery
DT = 0.1                      # hour step; smaller = more accurate, slower
TOL_MIN = 15                  # tolerance in minutes for re-entrainment
STREAK = 3                    # consecutive days within tolerance

# Baseline light levels (lux)
DAY_LUX = 250.0
NIGHT_LUX = 1.0

# All-nighter light level & timing (local clock hours)
ALL_NIGHTER_LUX = 1000.0
ALL_NIGHTER_START = 22.0
ALL_NIGHTER_END = 6.0         # next morning

# Low-lux sanity pulse (optional, but defined in original)
LOW_LUX = 75.0
LOW_PULSE_START = 21.0
LOW_PULSE_END = 24.0

# --- Light Schedule Functions ---

def make_ld_schedule(total_hours, day_lux, night_lux, day_start=7.0, day_end=23.0):
    """Generates a standard light-dark (LD) schedule."""
    t = np.arange(0, total_hours, DT)
    lux = np.zeros_like(t)

    # daily pattern repeating:
    for i, ti in enumerate(t):
        h = ti % 24.0
        in_day = (day_start <= h < day_end)
        lux[i] = day_lux if in_day else night_lux
    return t, lux

def insert_all_nighter(lux, t, start_hr=22.0, end_hr=6.0, lux_level=1000.0, day_index=30):
    """Replace night of baseline day_index with bright light (the all-nighter)."""
    # map day_index's clock interval
    for i, ti in enumerate(t):
        day = int(ti // 24.0)
        if day == day_index:
            h = ti % 24.0
            # Check for night hours, handling the wrap-around (22:00 to 24:00 and 0:00 to 6:00)
            if (h >= start_hr) or (h < end_hr):
                lux[i] = lux_level
    return lux

def branch_darkness_after(t, lux, start_day=31, night_lux=1.0):
    """Enforce constant dark nights during recovery phase (starting at start_day)."""
    out = lux.copy()
    for i, ti in enumerate(t):
        day = int(ti // 24.0)
        if day >= start_day:
            h = ti % 24.0
            # keep days same; enforce dark nights
            if not (7.0 <= h < 23.0):
                out[i] = night_lux
    return out

# --- Metric Functions ---

def daily_marker_times(model, traj):
    """
    Return times (in hours) for daily phase marker.
    Assumes models.py exposes a 'cbt(trajectory)' method to find markers.
    """
    try:
        # Both variants in your models.py expose CBT minima finders.
        cbt_times = model.cbt(traj)  # expected in hours timeline
        return np.array(cbt_times) if cbt_times else np.array([])
    except Exception as e:
        print(f"[WARNING] Error extracting CBT markers: {e}")
        return np.array([])

def circular_mean_hours(hours):
    """
    Compute circular mean of hours (0-24 scale).
    Returns mean hour in [0, 24).
    """
    if len(hours) == 0:
        return None
    
    # Convert hours to radians
    ang = np.array(hours) * (2 * np.pi / 24.0)
    
    # Compute mean angle
    mean_sin = np.mean(np.sin(ang))
    mean_cos = np.mean(np.cos(ang))
    mean_ang = math.atan2(mean_sin, mean_cos)
    
    # Convert back to hours [0, 24)
    if mean_ang < 0:
        mean_ang += 2 * np.pi
    
    return (mean_ang * 24.0) / (2 * np.pi)

def compute_baseline_marker(cbt_times, baseline_days=25):
    """
    Compute the target baseline phase marker (CBTmin) by averaging the
    last few baseline days (up to baseline_days).
    Returns mean clock hour (0-24).
    """
    if len(cbt_times) == 0:
        return None
    
    # Filter to baseline period
    base = [t for t in cbt_times if t < baseline_days * 24.0]
    if len(base) < 3:
        print(f"[WARNING] Only {len(base)} baseline markers found, need at least 3")
        return None

    # Get clock hours
    hours = np.array(base) % 24.0
    
    # Use circular mean
    return circular_mean_hours(hours)

def circular_distance(h1, h2):
    """
    Compute the shortest angular distance between two clock hours.
    Returns value in [-12, 12].
    """
    diff = h1 - h2
    # Wrap to [-12, 12]
    while diff > 12:
        diff -= 24
    while diff < -12:
        diff += 24
    return diff

def reentrainment_hours(cbt_times, baseline_marker_h, start_day=31, tol_min=15, streak=3):
    """
    Return the time to re-entrainment in HOURS (float), measured from
    the start of the recovery phase (start_day * 24 h).
    Requires streak consecutive DAYS within the tolerance.
    """
    if len(cbt_times) == 0 or baseline_marker_h is None:
        return None
    
    tol_h = tol_min / 60.0
    start_recovery_time = start_day * 24.0

    # Get all days that have markers after recovery started
    recovery_markers = [t for t in cbt_times if t >= start_recovery_time]
    if len(recovery_markers) == 0:
        return None
    
    # Group markers by day
    day_markers = {}
    for t in recovery_markers:
        day = int(t // 24)
        if day not in day_markers:
            day_markers[day] = []
        day_markers[day].append(t)
    
    # Sort days
    sorted_days = sorted(day_markers.keys())
    
    # Check each day for the start of a streak
    for i, d in enumerate(sorted_days):
        # Get mean clock hour for this day
        todays = day_markers[d]
        h = circular_mean_hours([ti % 24.0 for ti in todays])
        if h is None:
            continue
        
        # Check if within tolerance
        err = circular_distance(h, baseline_marker_h)
        
        if abs(err) <= tol_h:
            # Check if we have enough consecutive days for the streak
            if i + streak > len(sorted_days):
                # Not enough days left to complete streak
                continue
            
            # Verify the streak
            good = True
            for k in range(1, streak):
                next_day_idx = i + k
                if next_day_idx >= len(sorted_days):
                    good = False
                    break
                
                d2 = sorted_days[next_day_idx]
                
                # Check if days are actually consecutive
                if d2 != d + k:
                    good = False
                    break
                
                # Check if markers are within tolerance
                h2 = circular_mean_hours([ti % 24.0 for ti in day_markers[d2]])
                if h2 is None:
                    good = False
                    break
                
                err2 = circular_distance(h2, baseline_marker_h)
                if abs(err2) > tol_h:
                    good = False
                    break
            
            if good:
                # Re-entrainment achieved!
                # Use the mean marker time on the first day of the streak
                first_marker = np.mean(day_markers[d])
                return first_marker - start_recovery_time
    
    return None

def amp_series(model, traj):
    """
    Extracts the amplitude series from the trajectory states.
    Uses R (state[0]) for Hannay-like models or sqrt(x^2 + xc^2) for Forger-like models.
    """
    try:
        states = np.asarray(traj.states)
        if states.size == 0:
            return np.array([])
        
        R = states[:, 0]
        
        # Heuristic check: if the max value seems too large for an 'R' state
        if np.nanmax(np.abs(R)) > 5.0 or np.nanmin(R) < -2.0:
            # Fallback to sqrt(x^2+xc^2) for Cartesian coordinates
            if states.shape[1] >= 2:
                R = np.sqrt(states[:, 0]**2 + states[:, 1]**2)
            else:
                # If only one state, use absolute value as amplitude proxy
                R = np.abs(states[:, 0])
        
        # Ensure all values are positive (amplitude should be non-negative)
        R = np.abs(R)
        return R
        
    except Exception as e:
        print(f"[WARNING] Error extracting amplitude: {e}")
        return np.array([])

def hours_to_90pct(R, t, baseline_days=30, start_day=31):
    """
    Return the time (in HOURS) from the start of recovery until the
    amplitude first reaches 90% of its baseline mean.
    """
    if len(R) == 0 or len(t) == 0:
        return None
    
    # baseline amplitude = mean over last few baseline days
    base_mask = t < baseline_days * 24.0
    if not np.any(base_mask):
        return None
    
    R_base = np.mean(R[base_mask])
    if R_base <= 0:
        return None
    
    target = 0.9 * R_base
    start_time = start_day * 24.0
    
    for i, ti in enumerate(t):
        if ti < start_time:
            continue
        if R[i] >= target:
            return ti - start_time  # hours since recovery start
    
    return None

# --- Plotting Function ---

def plot_phase_amp(t, model, traj_dark, traj_ld, title):
    """Generates a two-panel plot showing phase error and daily mean amplitude."""
    # Phase markers
    cbt_dark = daily_marker_times(model, traj_dark)
    cbt_ld = daily_marker_times(model, traj_ld)

    # Compute baseline for error calculation
    baseline_h = compute_baseline_marker(cbt_ld, baseline_days=30) 
    
    if baseline_h is None:
        print(f"[WARNING] Cannot plot {title}: no baseline marker")
        return

    def err_series(cbt_times, baseline_h):
        """Calculates the daily phase error relative to the baseline marker."""
        if len(cbt_times) == 0:
            return np.array([]), np.array([])
        
        days = sorted(set([int(tt // 24.0) for tt in cbt_times]))
        errs = []
        xs = []
        for d in days:
            hs = [tt % 24.0 for tt in cbt_times if d * 24.0 <= tt < (d + 1) * 24.0]
            if not hs:
                continue
            xs.append(d)
            # Use circular distance
            mean_h = circular_mean_hours(hs)
            if mean_h is not None:
                e = circular_distance(mean_h, baseline_h)
                errs.append(e)
            else:
                errs.append(np.nan)
        return np.array(xs), np.array(errs)

    xd, ed = err_series(cbt_dark, baseline_h)
    xl, el = err_series(cbt_ld, baseline_h)

    # Amplitude series
    Rd = amp_series(model, traj_dark)
    Rl = amp_series(model, traj_ld)

    fig, axs = plt.subplots(2, 1, figsize=(8, 6), sharex=True)
    
    # --- Phase Plot ---
    if len(xd) > 0:
        axs[0].plot(xd, ed, 'o-', label='Darkness Recovery', markersize=3)
    if len(xl) > 0:
        axs[0].plot(xl, el, 'x-', label='LD Recovery', markersize=3)
    axs[0].axhline(0, color='gray', linestyle='--', lw=0.8, alpha=0.6)
    axs[0].axhline(TOL_MIN/60.0, color='red', linestyle=':', lw=0.8, alpha=0.6)
    axs[0].axhline(-TOL_MIN/60.0, color='red', linestyle=':', lw=0.8, alpha=0.6)
    axs[0].set_ylabel('Phase error (h)')
    axs[0].set_title(f'Circadian Re-entrainment Simulation: {title}')
    axs[0].legend()
    axs[0].grid(True, alpha=0.3)

    # --- Amplitude Plot ---
    if len(t) > 0:
        days = np.arange(0, int(t[-1] // 24.0) + 1)
        
        def daily_mean(R):
            """Calculates the mean amplitude for each full day."""
            if len(R) == 0:
                return np.full(len(days), np.nan)
            out = []
            for d in days:
                mask = (t >= d * 24.0) & (t < (d + 1) * 24.0)
                if np.any(mask):
                    out.append(np.mean(R[mask]))
                else:
                    out.append(np.nan)
            return np.array(out)
            
        daily_Rd = daily_mean(Rd)
        daily_Rl = daily_mean(Rl)
        
        # Baseline for 90% target line
        baseline_R = np.nanmean(daily_Rl[:30]) if len(daily_Rl) >= 30 else np.nanmean(daily_Rl)
        target_90 = 0.9 * baseline_R if not np.isnan(baseline_R) else None

        if len(daily_Rd) > 0:
            axs[1].plot(days, daily_Rd, 'o-', label='Darkness Recovery', markersize=3)
        if len(daily_Rl) > 0:
            axs[1].plot(days, daily_Rl, 'x-', label='LD Recovery', markersize=3)
        if target_90 is not None:
            axs[1].axhline(target_90, color='green', linestyle=':', lw=0.8, alpha=0.6, label='90% Baseline Amp')
        axs[1].set_xlabel('Day')
        axs[1].set_ylabel('Daily Mean Amplitude')
        axs[1].legend()
        axs[1].grid(True, alpha=0.3)

    plt.tight_layout()
    plt.show()

# --- Main Execution ---

def main():
    try:
        import circadian.models as models  # your models.py
    except ImportError:
        print("[CRITICAL ERROR] Could not import 'circadian.models'. Ensure it is in the path.")
        print("Trying alternative import...")
        try:
            import models
        except ImportError:
            print("[CRITICAL ERROR] Could not import 'models' either. Exiting.")
            return

    # 1. Build baseline schedule
    t, lux = make_ld_schedule(HOURS, day_lux=DAY_LUX, night_lux=NIGHT_LUX)

    # 2. Insert all-nighter on day 30 -> 22:00–06:00 bright
    lux_disrupted = insert_all_nighter(lux.copy(), t,
                                       start_hr=ALL_NIGHTER_START,
                                       end_hr=ALL_NIGHTER_END,
                                       lux_level=ALL_NIGHTER_LUX,
                                       day_index=30)

    # 3. Branches: darkness vs LD during recovery (starting day 31)
    lux_dark = branch_darkness_after(t, lux_disrupted, start_day=31, night_lux=NIGHT_LUX)
    lux_ld = lux_disrupted  # already back to normal LD after day 31

    for name in MODEL_NAMES:
        if not hasattr(models, name):
            print(f"\n[!] Could not find model class '{name}' in models.py. Skipping.")
            continue
        
        ModelClass = getattr(models, name)
        model = ModelClass()
        
        print(f"\n=== Running {name} ===")
        
        # Integrate trajectories
        try:
            traj_dark = model.integrate(t, input=lux_dark)
            traj_ld = model.integrate(t, input=lux_ld)
        except Exception as e:
            print(f"[ERROR] Integration failed for {name}: {e}. Skipping metrics.")
            continue

        # Phase markers & baseline
        cbt_ld = daily_marker_times(model, traj_ld)
        baseline_h = compute_baseline_marker(cbt_ld, baseline_days=30)
        
        if baseline_h is None:
            print(f"[ERROR] Could not establish a reliable baseline marker for {name}. Skipping metrics.")
            continue

        print(f"Baseline CBTmin: {baseline_h:.2f} h")

        # Metrics
        rdark = reentrainment_hours(daily_marker_times(model, traj_dark), baseline_h,
                                    start_day=31, tol_min=TOL_MIN, streak=STREAK)
        rld = reentrainment_hours(cbt_ld, baseline_h,
                                  start_day=31, tol_min=TOL_MIN, streak=STREAK)

        Rd = amp_series(model, traj_dark)
        Rl = amp_series(model, traj_ld)

        dA_d = hours_to_90pct(Rd, t, baseline_days=30, start_day=31)
        dA_l = hours_to_90pct(Rl, t, baseline_days=30, start_day=31)
        
        recovery_ratio = None
        if dA_d is not None and dA_l is not None and dA_l != 0:
            recovery_ratio = dA_d / dA_l
        
        print(f"Re-entrainment (Dark): {rdark:.2f} h ({rdark/24:.1f} days)" if rdark is not None else "Re-entrainment (Dark): not reached")
        print(f"Re-entrainment (LD):   {rld:.2f} h ({rld/24:.1f} days)"    if rld is not None else "Re-entrainment (LD):   not reached")
        print(f"90% amplitude (Dark):  {dA_d:.2f} h ({dA_d/24:.1f} days)" if dA_d is not None else "90% amplitude (Dark):  not reached")
        print(f"90% amplitude (LD):    {dA_l:.2f} h ({dA_l/24:.1f} days)" if dA_l is not None else "90% amplitude (LD):    not reached")
        print(f"Recovery Ratio (Dark/LD): {recovery_ratio:.2f}" if recovery_ratio is not None else "Recovery Ratio (Dark/LD): N/A")

        plot_phase_amp(t, model, traj_dark, traj_ld, title=name)

if __name__ == "__main__":
    main()