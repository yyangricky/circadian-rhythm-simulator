
import importlib
import numpy as np
import math
import matplotlib.pyplot as plt
from metrics import get_phase_markers, compute_baseline_marker, reentrainment_hours, amp_series, hours_to_90pct, circular_mean_hours, circular_distance 
from protocols import allnighter_protocol, branch_ld_after, branch_dark_after, branch_shortnight_after
from models import Jewett99, Forger99, Hannay19, Hannay19TP, Breslow13,Skeldon23
'''
ld, LD -> a regular 'light-dark' cycle (16 hour day, 8 hour nights)
sn, SN-> a 'short night' cycle that implements a poor sleep schedule (5 hour nights, 19 hour days)  
d, dark -> simulates a constant darkness condition (reflects how the a external zeitgebers, re-entrainment is impossible)
'''
# Which phase marker to use for the analysis ("cbt" or "dlmo")
PHASE_MARKER = "dlmo" #"cbt"    "dlmo"
MODEL_CLASSES = [
    ("Jewett99",   Jewett99),
    ("Forger99",   Forger99),
    ("Hannay19",   Hannay19),
    ("Hannay19TP", Hannay19TP),
    ("Breslow13", Breslow13),
    ('Skeldon23', Skeldon23)
]

# ==== protocol and light schedule ====
DT = 0.1
BASELINE_DAYS = 25
RECOVERY_DAYS = 15
HOURS = 24* (BASELINE_DAYS + 1 + RECOVERY_DAYS) 
# Re-entrainment definition
TOL_MIN = 15     # tolerance in minutes = 0.25 hours
STREAK = 3       # consecutive days within tolerance
DAY_START = 6.0
DAY_END = 22.0
PCT = 0.9999 
# Baseline light levels (lux)
DAY_LUX = 3000.0
NIGHT_LUX = 0.0
DARK_LUX = 1.0
# All-nighter 
ALLNIGHTER_LUX = 1000.0
ALL_NIGHTER_START = 22.0 
ALL_NIGHTER_END = 6.0   # next morning
# short-night light delay (in hours) - we are assuming a lux level equal to the allnighter_lux
PAST_BEDTIME = 3
t, lux_all, baseline_days, recovery_start_day = allnighter_protocol(
    baseline_days=BASELINE_DAYS,
    recovery_days=RECOVERY_DAYS,
    dt=DT,
    day_start=DAY_START,
    day_end=DAY_END,
    day_lux=DAY_LUX,
    night_lux=NIGHT_LUX,
    allnighter_lux=ALLNIGHTER_LUX,
)

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
                # absolute re-entrainment time in hours 
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
    traj_sn,
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
    markers_sn   = get_phase_markers(model, traj_sn,   marker=PHASE_MARKER)

    days_d, mean_d = daily_mean_marker(markers_dark)
    days_l, mean_l = daily_mean_marker(markers_ld)
    days_s, mean_s = daily_mean_marker(markers_sn)

    err_d = np.array(
        [circular_distance(h, baseline_h) if h is not None else np.nan for h in mean_d]
    )
    err_l = np.array(
        [circular_distance(h, baseline_h) if h is not None else np.nan for h in mean_l]
    )
    err_s = np.array(
        [circular_distance(h, baseline_h) if h is not None else np.nan for h in mean_s])
    # --- Amplitude series (daily means) ---
    Rd = amp_series(model, traj_dark)
    Rl = amp_series(model, traj_ld)
    Rsn = amp_series(model, traj_sn)
    
   
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
    days_amp_s, daily_Rsn = daily_mean_amp(Rsn, t)

    # Baseline amplitude (from LD branch only)
    start_baseline = (recovery_start_day - 10) * 24.0
    end_baseline = recovery_start_day * 24.0
    mask = (t >= start_baseline) & (t < end_baseline)

    if np.any(mask):
        baseline_R = np.nanmean(Rl[mask])
    else:
        baseline_R = np.nan  # or handle as you prefer

    target_90 = PCT * baseline_R if not np.isnan(baseline_R) else None
    '''CHANGED TO 99% FOR NOW'''
    ##plot
    fig, axs = plt.subplots(2, 1, figsize=(8, 6), sharex=True)

    # Top: phase error
    if days_d.size > 0:
        axs[0].plot(days_d, err_d, "o-", label="Darkness", markersize=3)
    if days_l.size > 0:
        axs[0].plot(days_l, err_l, "x-", label="LD", markersize=3)
    if days_s.size > 0:
        axs[0].plot(days_s, err_s, "s-", label="Short-night", markersize=3)

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
    if days_amp_s.size > 0:
        axs[1].plot(days_amp_s, daily_Rsn, "s-", label="Short-night", markersize=3)
    if target_90 is not None:
        axs[1].axhline(
            target_90,
            color="green",
            linestyle=":",
            linewidth=0.8,
            alpha=0.6,
            label="AVG baseline amplitude",
        )
    axs[1].axvline(recovery_start_day, color="grey", linestyle=":", linewidth=0.8)
    axs[1].set_xlabel("Day")
    axs[1].set_ylabel("Daily mean amplitude")
    axs[1].legend()
    axs[1].grid(True, alpha=0.3)

    fig.suptitle(f"{title} – {phase_marker_name}-based analysis")
    plt.tight_layout()
    plt.show()

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
        allnighter_end=ALL_NIGHTER_END)

    lux_ld   = branch_ld_after(t, lux_all, recovery_start_day)
    lux_dark = branch_dark_after(t, lux_all, recovery_start_day, dark_lux=DARK_LUX)
    lux_sn   = branch_shortnight_after(t,lux_all,recovery_start_day,
        day_start=DAY_START,  
        day_lux=DAY_LUX,
        night_lux=NIGHT_LUX,
        allnighter_lux=ALLNIGHTER_LUX,
        past_bedtime = PAST_BEDTIME)

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

        # Integrate under each branch
        try:
            traj_dark = model.integrate(t, input=lux_dark)
            traj_ld   = model.integrate(t, input=lux_ld)
            traj_sn   = model.integrate(t, input=lux_sn)
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
        markers_sn   = get_phase_markers(model, traj_sn,   marker=PHASE_MARKER)

        # For DLMO we anchor re-entrainment to the last full baseline day
        reentrain_start_day = baseline_days if PHASE_MARKER.lower() == "dlmo" else recovery_start_day

        rdark = reentrainment_hours(
            markers_dark,
            baseline_h,
            start_day=reentrain_start_day,
            tol_min=TOL_MIN,
            streak=STREAK,
        )
        rld = reentrainment_hours(
            markers_ld,
            baseline_h,
            start_day=reentrain_start_day,
            tol_min=TOL_MIN,
            streak=STREAK,
        )
        rsn = reentrainment_hours(
            markers_sn,
            baseline_h,
            start_day=reentrain_start_day,
            tol_min=TOL_MIN,
            streak=STREAK,
        )
        # Amplitude recovery (hours from start of recovery)
        Rd = amp_series(model, traj_dark)
        Rl = amp_series(model, traj_ld)
        Rsn = amp_series(model, traj_sn)

        start_baseline = (recovery_start_day - 10) * 24.0
        end_baseline   = recovery_start_day * 24.0
        baseline_mask  = (t >= start_baseline) & (t < end_baseline)
        baseline_R_ld  = np.nanmean(Rl[baseline_mask]) if np.any(baseline_mask) else np.nan

        dA_d  = hours_to_90pct(Rd,  t, pct=PCT, start_day=recovery_start_day, sustain_days=STREAK, baseline_R=baseline_R_ld)
        dA_l  = hours_to_90pct(Rl,  t, pct=PCT, start_day=recovery_start_day, sustain_days=STREAK, baseline_R=baseline_R_ld)
        dA_sn = hours_to_90pct(Rsn, t, pct=PCT, start_day=recovery_start_day, sustain_days=STREAK, baseline_R=baseline_R_ld,)

                
        if dA_d is not None and dA_l is not None and dA_l != 0:
            recovery_ratio = dA_d / dA_l
        # Recovery ratio (short night/LD)
        recovery_ratio_dark = None
        recovery_ratio_sn   = None

        if dA_d is not None and dA_l is not None and dA_l != 0:
            recovery_ratio_dark = dA_d / dA_l

        if dA_sn is not None and dA_l is not None and dA_l != 0:
            recovery_ratio_sn = dA_sn / dA_l

        # Print metrics 
        if rdark is not None:
            print(f"Re-entrainment (Dark): {rdark:.2f} h ({rdark/24.0:.2f} days)")
       # else:
     #       print("Re-entrainment (Dark): N/A (In constant darkness, re-entrainment to an external cycle is impossible as the clock is free-running!")

        if rld is not None:
            print(f"Re-entrainment (LD): {rld:.2f} h ({rld/24.0:.2f} days)")
        else:
            print("Re-entrainment (LD): not reached")

        if rsn is not None:
            print(f"Re-entrainment (SN): {rsn:.2f} h ({rsn/24.0:.2f} days)")
        else:
            print("Re-entrainment (SN): not reached")

        if dA_d is not None:
            print(f"90% amplitude (Dark): {dA_d:.2f} h ({dA_d/24.0:.2f} days)")
        else:
            print("90% amplitude (Dark): not reached")

        if dA_l is not None:
            print(f"90% amplitude (LD): {dA_l:.2f} h ({dA_l/24.0:.2f} days)")
        else:
            print("90% amplitude (LD): not reached")

        if dA_sn is not None:
            print(f"90% amplitude (SN): {dA_sn:.2f} h ({dA_sn/24.0:.2f} days)")
        else:
            print("90% amplitude (SN): not reached")

        if recovery_ratio_dark is not None:
            print(f"Recovery ratio Dark/LD: {recovery_ratio_dark:.2f}")
        else:
            print("Recovery ratio Dark/LD: N/A")

        if recovery_ratio_sn is not None:
            print(f"Recovery ratio SN/LD: {recovery_ratio_sn:.2f}")
        else:
            print("Recovery ratio SN/LD: N/A")

        phase_label = "CBTmin" if PHASE_MARKER.lower() == "cbt" else "DLMO"
        plot_phase_amp(
            t,
            model,
            traj_dark,
            traj_ld,
            traj_sn, 
            baseline_h,
            recovery_start_day,
            title=name,
            phase_marker_name=phase_label,
        )

if __name__ == "__main__":
    main()
