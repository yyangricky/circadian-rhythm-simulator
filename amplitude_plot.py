
import importlib
import numpy as np
import math
import matplotlib.pyplot as plt
from metrics import get_phase_markers, compute_baseline_marker, reentrainment_hours, amp_series, hours_to_90pct, circular_mean_hours, circular_distance 
from protocols import allnighter_protocol, branch_ld_after, branch_dark_after
from models import Jewett99, Forger99, Hannay19, Hannay19TP

# Which phase marker to use for the analysis ("cbt" or "dlmo")
PHASE_MARKER = "cbt" #"cbt"    "dlmo"
MODEL_CLASSES = [
    ("Jewett99",   Jewett99),
    ("Forger99",   Forger99),
    ("Hannay19",   Hannay19),
    ("Hannay19TP", Hannay19TP),]

# ==== protocol and light schedule ====
HOURS = 24.00* (30 + 1 + 14)           # 30 baseline days + 1 disruption night + 14 recovery
DT = 0.10
BASELINE_DAYS = 30
RECOVERY_DAYS = 14
# Re-entrainment definition
TOL_MIN = 15     # phase tolerance in minutes
STREAK = 3       # consecutive days within tolerance
DAY_START = 7.0
DAY_END = 23.0
# Baseline light levels (lux)
DAY_LUX = 1000.0
NIGHT_LUX = 1.0
DARK_LUX = 100.0
# All-nighter 
ALLNIGHTER_LUX = 1000.0
ALL_NIGHTER_START = 22.0
ALL_NIGHTER_END = 6.0   # next morning

t, lux_all, baseline_days, recovery_start_day = allnighter_protocol(
    baseline_days=BASELINE_DAYS,
    recovery_days=RECOVERY_DAYS,
    dt=DT,
    day_start=7.0,
    day_end=23.0,
    day_lux=DAY_LUX,
    night_lux=NIGHT_LUX,
    allnighter_lux=ALLNIGHTER_LUX,
)
'''

lux_ld   = branch_ld_after(t, lux_all, recovery_start_day)
lux_dark = branch_dark_after(t, lux_all, recovery_start_day, dark_lux=NIGHT_LUX)



def make_ld_schedule(total_hours, day_lux, night_lux, day_start=7.0, day_end=23.0):
    t = np.arange(0, total_hours, DT)
    lux = np.zeros_like(t)
    # daily pattern repeating:
    for i, ti in enumerate(t):
        h = ti % 24.00
        in_day = (day_start <= h < day_end)
        lux[i] = day_lux if in_day else night_lux
    return t, lux

def insert_all_nighter(lux, t, start_hr=22.0, end_hr=6.0, lux_level=1000.0, day_index=30):
    """Replace night of baseline day_index with bright light 22:00–06:00."""
    # map day_index's clock interval
    for i, ti in enumerate(t):
        day = int(ti // 24.00)
        if day == day_index:
            h = ti % 24.00
            if (h >= start_hr) or (h < end_hr):
                lux[i] = lux_level
    return lux

def branch_darkness_after(t, lux, start_day=31, night_lux=1.0):
    out = lux.copy()
    for i, ti in enumerate(t):
        day = int(ti // 24.00)
        if day >= start_day: 
            out[i] = night_lux
    return out

def daily_marker_times(model, traj):
    """Return times (in hours) for daily phase marker.
       For Forger-like models we’ll use CBTmin via model.cbt(); for Hannay-like via model.cbt()."""
    # Both variants in your models.py expose CBT minima finders. We’ll try model.cbt()
    cbt_times = model.cbt(traj)  # expected in hours timeline
    return np.array(cbt_times)
  '''
def daily_mean_marker(marker_times):
    """
    Given absolute marker times (hours), return (days, daily_mean_clock_hour).
    Uses circular mean within each day.
    """
    marker_times = np.asarray(marker_times)
    if marker_times.size == 0:
        return np.array([]), np.array([])

    day_dict = {}
    for t in marker_times:
        d = int(t // 24.0)
        day_dict.setdefault(d, []).append(t)

    days = sorted(day_dict.keys())
    means = []
    for d in days:
        hs = [ti % 24.0 for ti in day_dict[d]]
        means.append(circular_mean_hours(hs))

    return np.array(days), np.array(means)


def reentrainment_hours(cbt_times, baseline_marker_h, start_day=31, tol_min=15, streak=3):
    """
    Return the time to re-entrainment in HOURS (float), measured from
    the start of the recovery phase (start_day*24 h).

    We still require `streak` consecutive DAYS within the tolerance,
    but we place the re-entrainment time at the *mean marker time*
    on the first good day.
    """
    tol_h = tol_min / 60.0
    # consider only markers during the recovery phase
    days = sorted(set(int(t // 24) for t in cbt_times if t >= start_day * 24))

    for d in days:
        # markers on this day
        todays = [t for t in cbt_times if d * 24 <= t < (d + 1) * 24]
        if not todays:
            continue

        # mean clock hour for that day
        h = np.mean([ti % 24.0 for ti in todays])
        err = (((h - baseline_marker_h + 12) % 24) - 12)   # wrapped error in hours
        if abs(err) <= tol_h:
            # check streak on following days
            good = True
            for k in range(1, streak):
                d2 = d + k
                toms = [t for t in cbt_times if d2 * 24 <= t < (d2 + 1) * 24]
                if not toms:
                    good = False
                    break
                h2 = np.mean([ti % 24.0 for ti in toms])
                err2 = (((h2 - baseline_marker_h + 12) % 24) - 12)
                if abs(err2) > tol_h:
                    good = False
                    break

            if good:
                # absolute re-entrainment time in hours (since t=0)
                t_abs = d * 24.0 + h
                # convert to hours since start of recovery
                return t_abs - start_day * 24.0

    return None

'''
def plot_phase_amp(t, model, traj_dark, traj_ld, title):
    # phase markers
    marker_dark = get_phase_markers(model, traj_dark, marker=PHASE_MARKER)
    marker_ld   = get_phase_markers(model, traj_ld, marker=PHASE_MARKER)

    # phase error lines
    # compute baseline
    baseline_h = compute_baseline_marker(marker_ld)  # both had baseline; use ld to get a baseline
    def err_series(cbt_times):
        days = sorted(set([int(tt//24.00) for tt in cbt_times]))
        errs = []
        xs = []
        for d in days:
            hs = [tt%24.00 for tt in cbt_times if d*24.00 <= tt < (d+1)*24.00]
            if not hs: continue
            xs.append(d)
            # circular wrap error
            e = (((np.mean(hs) - baseline_h + 12) % 24.00) - 12)
            errs.append(e)
        return np.array(xs), np.array(errs)
    xd, ed = err_series(marker_dark)
    xl, el = err_series(marker_ld)

    # amplitude
    Rd = amp_series(model, traj_dark)
    Rl = amp_series(model, traj_ld)

    fig, axs = plt.subplots(2, 1, figsize=(8,6), sharex=True)
    axs[0].plot(xd, ed, label='Darkness')
    axs[0].plot(xl, el, label='LD')
    axs[0].axhline(0, lw=0.8, alpha=0.6)
    axs[0].set_ylabel('Phase error (h)')
    axs[0].set_title(title)
    axs[0].legend()

    # downsample amplitude daily mean
    days = np.arange(0, int(t[-1]//24.00)+1)
    def daily_mean(R):
        out=[]
        for d in days:
            mask = (t>=d*24.00)&(t<(d+1)*24.00)
            if np.any(mask): out.append(np.mean(R[mask]))
            else: out.append(np.nan)
        return np.array(out)
    axs[1].plot(days, daily_mean(Rd), label='Darkness')
    axs[1].plot(days, daily_mean(Rl), label='LD')
    axs[1].set_xlabel('Day')
    axs[1].set_ylabel('Daily mean amplitude')
    axs[1].legend()
    plt.tight_layout()
    plt.show()
'''
def plot_phase_amp(
    t,
    model,
    traj_dark,
    traj_ld,
    baseline_h,
    recovery_start_day,
    title,
    phase_marker_name="CBTmin",
):
    """
    Two-panel plot:
      top: daily phase error vs baseline
      bottom: daily mean amplitude vs day (dark vs LD), with 90% baseline line.
    """

    # --- Phase error series ---

    markers_dark = get_phase_markers(model, traj_dark, marker=PHASE_MARKER)
    markers_ld   = get_phase_markers(model, traj_ld,   marker=PHASE_MARKER)

    days_d, mean_d = daily_mean_marker(markers_dark)
    days_l, mean_l = daily_mean_marker(markers_ld)

    err_d = np.array(
        [circular_distance(h, baseline_h) if h is not None else np.nan for h in mean_d]
    )
    err_l = np.array(
        [circular_distance(h, baseline_h) if h is not None else np.nan for h in mean_l]
    )

    # --- Amplitude series (daily means) ---

    Rd = amp_series(model, traj_dark)
    Rl = amp_series(model, traj_ld)

    def daily_mean_amp(R, t):
        if len(R) == 0:
            return np.array([]), np.array([])
        total_days = int(t[-1] // 24) + 1
        days = np.arange(total_days)
        out = []
        for d in days:
            mask = (t >= d * 24.0) & (t < (d + 1) * 24.0)
            if np.any(mask):
                out.append(np.mean(R[mask]))
            else:
                out.append(np.nan)
        return days, np.array(out)

    days_amp_d, daily_Rd = daily_mean_amp(Rd, t)
    days_amp_l, daily_Rl = daily_mean_amp(Rl, t)

    # Baseline amplitude (from LD branch only)
    baseline_mask = days_amp_l < BASELINE_DAYS
    if np.any(baseline_mask):
        baseline_R = np.nanmean(daily_Rl[baseline_mask])
    else:
        baseline_R = np.nanmean(daily_Rl)
    target_90 = 0.9 * baseline_R if not np.isnan(baseline_R) else None

    ##plot
    fig, axs = plt.subplots(2, 1, figsize=(8, 6), sharex=True)

    # Top: phase error
    if days_d.size > 0:
        axs[0].plot(days_d, err_d, "o-", label="Darkness Recovery", markersize=3)
    if days_l.size > 0:
        axs[0].plot(days_l, err_l, "x-", label="LD Recovery", markersize=3)
    axs[0].axhline(0.0, color="k", linestyle="--", linewidth=0.8)
    axs[0].axvline(recovery_start_day, color="grey", linestyle=":", linewidth=0.8)
    axs[0].set_ylabel(f"Phase error ({phase_marker_name}) [h]")
    axs[0].legend()
    axs[0].grid(True, alpha=0.3)

    # Bottom: amplitude
    if days_amp_d.size > 0:
        axs[1].plot(days_amp_d, daily_Rd, "o-", label="Darkness Recovery", markersize=3)
    if days_amp_l.size > 0:
        axs[1].plot(days_amp_l, daily_Rl, "x-", label="LD Recovery", markersize=3)
    if target_90 is not None:
        axs[1].axhline(
            target_90,
            color="green",
            linestyle=":",
            linewidth=0.8,
            alpha=0.6,
            label="90% baseline amplitude",
        )
    axs[1].axvline(recovery_start_day, color="grey", linestyle=":", linewidth=0.8)
    axs[1].set_xlabel("Day")
    axs[1].set_ylabel("Daily mean amplitude")
    axs[1].legend()
    axs[1].grid(True, alpha=0.3)

    fig.suptitle(f"{title} – {phase_marker_name}-based analysis")
    plt.tight_layout()
    plt.show()
'''
def main():
    models = importlib.import_module("models")  # your models.py
    # Build baseline schedule
    t, lux = make_ld_schedule(HOURS, day_lux=DAY_LUX, night_lux=NIGHT_LUX)
    # Insert all-nighter on day 30 -> 22:00–06:00 bright
    lux_disrupted = insert_all_nighter(lux.copy(), t,
                                       start_hr=ALL_NIGHTER_START,
                                       end_hr=ALL_NIGHTER_END,
                                       lux_level=ALLNIGHTER_LUX,
                                       day_index=30)
    # Branches: darkness vs LD during recovery
    lux_dark = branch_darkness_after(t, lux_disrupted, start_day=31, night_lux=0.0)
    lux_ld   = lux_disrupted  # already back to normal LD after day 31

    for name in MODEL_NAMES:
        if not hasattr(models, name):
            print(f"[!] Could not find model class '{name}' in models.py. Skipping.")
            continue
        ModelClass = getattr(models, name)
    

        print(f"\n=== Running {name} ===")

        try:
            model = ModelClass()
            traj_dark = model.integrate(t, input=lux_dark)
            traj_ld   = model.integrate(t, input=lux_ld)

        # continue as normal...
        except Exception as e:
            print(f"[ERROR running {name}]: {e}")
            continue

        # Phase markers & baseline
        marker_ld   = get_phase_markers(model, traj_ld, marker=PHASE_MARKER)
        baseline_h = compute_baseline_marker(marker_ld, baseline_days=30)
        # Metrics
        rdark = reentrainment_hours(daily_marker_times(model, traj_dark), baseline_h,
                                    start_day=31, tol_min=TOL_MIN, streak=STREAK)
        rld   = reentrainment_hours(marker_ld, baseline_h,
                                    start_day=31, tol_min=TOL_MIN, streak=STREAK)

        Rd = amp_series(model, traj_dark)
        Rl = amp_series(model, traj_ld)
        dA_d = hours_to_90pct(Rd, t, baseline_days=30, start_day=31)
        dA_l = hours_to_90pct(Rl, t, baseline_days=30, start_day=31)

        recovery_ratio = None
        if dA_d and dA_l:
            recovery_ratio = dA_d / dA_l if dA_l != 0 else None

        print(f"Re-entrainment (Dark): {rdark:.2f} h" if rdark is not None else "Re-entrainment (Dark): not reached")
        print(f"Re-entrainment (LD):   {rld:.2f} h"   if rld   is not None else "Re-entrainment (LD):   not reached")
        print(f"90% amplitude (Dark):  {dA_d:.2f} h" if dA_d is not None else "90% amplitude (Dark):  not reached")
        print(f"90% amplitude (LD):    {dA_l:.2f} h" if dA_l is not None else "90% amplitude (LD):    not reached")
        print(f"Recovery Ratio (Dark/LD): {recovery_ratio}")


        plot_phase_amp(t, model, traj_dark, traj_ld, title=name)
'''
def main():
    # Build common protocol (baseline + all-nighter + recovery)
    t, lux_all, baseline_days, recovery_start_day = allnighter_protocol(
        baseline_days=BASELINE_DAYS,
        recovery_days=RECOVERY_DAYS,
        dt=DT,
        day_start=DAY_START,
        day_end=DAY_END,
        day_lux=DAY_LUX,
        night_lux=NIGHT_LUX,
        allnighter_lux=ALLNIGHTER_LUX,
        allnighter_start=ALL_NIGHTER_START,
        allnighter_end=ALL_NIGHTER_END,
    )

    lux_ld   = branch_ld_after(t, lux_all, recovery_start_day)
    lux_dark = branch_dark_after(t, lux_all, recovery_start_day, dark_lux=DARK_LUX)

    print(f"Using phase marker: {PHASE_MARKER.upper()}")
    print(f"Baseline days: {BASELINE_DAYS}, recovery days: {RECOVERY_DAYS}")
    print(f"Recovery starts on day {recovery_start_day}")

    for name, cls in MODEL_CLASSES:
        print("\n==============================")
        print(f"Model: {name}")

        try:
            model = cls()
        except Exception as e:
            print(f"[ERROR] Could not construct model {name}: {e}")
            continue

        # Integrate under dark and LD recovery branches
        try:
            traj_dark = model.integrate(t, input=lux_dark)
            traj_ld   = model.integrate(t, input=lux_ld)
        except Exception as e:
            print(f"[ERROR] Integration failed for {name}: {e}")
            continue

        # Phase markers and baseline
        markers_ld = get_phase_markers(model, traj_ld, marker=PHASE_MARKER)
        baseline_h = compute_baseline_marker(markers_ld, baseline_days=BASELINE_DAYS)

        if baseline_h is None:
            print(f"[ERROR] No reliable baseline {PHASE_MARKER.upper()} for {name}, skipping metrics.")
            continue

        print(f"Baseline {PHASE_MARKER.upper()}: {baseline_h:.2f} h")

        # Re-entrainment times (hours from start of recovery)
        markers_dark = get_phase_markers(model, traj_dark, marker=PHASE_MARKER)

        rdark = reentrainment_hours(
            markers_dark,
            baseline_h,
            start_day=recovery_start_day,
            tol_min=TOL_MIN,
            streak=STREAK,
        )
        rld = reentrainment_hours(
            markers_ld,
            baseline_h,
            start_day=recovery_start_day,
            tol_min=TOL_MIN,
            streak=STREAK,
        )

        # Amplitude recovery (hours from start of recovery)
        Rd = amp_series(model, traj_dark)
        Rl = amp_series(model, traj_ld)

        dA_d = hours_to_90pct(
            Rd,
            t,
            baseline_days=BASELINE_DAYS,
            start_day=recovery_start_day,
        )
        dA_l = hours_to_90pct(
            Rl,
            t,
            baseline_days=BASELINE_DAYS,
            start_day=recovery_start_day,
        )

        # Recovery ratio (dark/LD)
        recovery_ratio = None
        if dA_d is not None and dA_l is not None and dA_l != 0:
            recovery_ratio = dA_d / dA_l

        # Print metrics in a consistent format
        if rdark is not None:
            print(f"Re-entrainment (Dark): {rdark:.2f} h ({rdark/24.0:.2f} days)")
        else:
            print("Re-entrainment (Dark): not reached")

        if rld is not None:
            print(f"Re-entrainment (LD):   {rld:.2f} h ({rld/24.0:.2f} days)")
        else:
            print("Re-entrainment (LD):   not reached")

        if dA_d is not None:
            print(f"90% amplitude (Dark):  {dA_d:.2f} h ({dA_d/24.0:.2f} days)")
        else:
            print("90% amplitude (Dark):  not reached")

        if dA_l is not None:
            print(f"90% amplitude (LD):    {dA_l:.2f} h ({dA_l/24.0:.2f} days)")
        else:
            print("90% amplitude (LD):    not reached")

        if recovery_ratio is not None:
            print(f"Recovery ratio (Dark/LD): {recovery_ratio:.2f}")
        else:
            print("Recovery ratio (Dark/LD): N/A")

        # Plot for this model
        phase_label = "CBTmin" if PHASE_MARKER.lower() == "cbt" else "DLMO"
        plot_phase_amp(
            t,
            model,
            traj_dark,
            traj_ld,
            baseline_h,
            recovery_start_day,
            title=name,
            phase_marker_name=phase_label,
        )


if __name__ == "__main__":
    main()

