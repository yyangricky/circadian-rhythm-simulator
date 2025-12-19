"""
figure 4 generation - phase deviation analysis
makes the phase deviation figure for social jet lag
"""

import numpy as np
import matplotlib.pyplot as plt
from models import Forger99, Jewett99, Hannay19, Hannay19TP
from metrics import get_phase_markers, circular_distance

# plot settings
plt.rcParams['font.size'] = 11
plt.rcParams['axes.labelsize'] = 12
plt.rcParams['axes.titlesize'] = 13
plt.rcParams['xtick.labelsize'] = 10
plt.rcParams['ytick.labelsize'] = 10
plt.rcParams['legend.fontsize'] = 10
plt.rcParams['figure.titlesize'] = 14
plt.rcParams['figure.dpi'] = 100
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['savefig.bbox'] = 'tight'

# protocol parameters
DT = 0.10
BASELINE_WEEKS = 4
DISRUPTION_WEEKENDS = 3
RECOVERY_WEEKS = 3
PHASE_MARKER = "cbt"

# models to test
MODELS = [
    ("Forger99", Forger99),
    ("Jewett99", Jewett99),
    ("Hannay19", Hannay19),
    ("Hannay19TP", Hannay19TP),
]

# colors for plotting
MODEL_COLORS = {
    "Forger99": "#2E86AB",
    "Jewett99": "#A23B72",
    "Hannay19": "#F18F01",
    "Hannay19TP": "#C73E1D",
}

# marker shapes for plotting
MODEL_MARKERS = {
    "Forger99": "o",
    "Jewett99": "s",
    "Hannay19": "^",
    "Hannay19TP": "D",
}


def social_jetlag_protocol(baseline_weeks, disruption_weekends, recovery_weeks, dt):
    """create the social jet lag light schedule"""
    baseline_days = baseline_weeks * 7
    disruption_days = disruption_weekends * 7
    recovery_days = recovery_weeks * 7
    total_days = baseline_days + disruption_days + recovery_days
    
    t = np.arange(0.0, 24.0 * total_days, dt)
    lux = np.zeros_like(t)
    
    disruption_start = baseline_days
    disruption_end = baseline_days + disruption_days
    
    # go through each timepoint and set light level
    for i, ti in enumerate(t):
        day = int(ti // 24.0)
        hour = ti % 24.0
        day_of_week = day % 7
        is_weekend = (day_of_week == 5 or day_of_week == 6)  # sat/sun
        
        if day < disruption_start:
            # baseline: normal 7am-11pm schedule
            lux[i] = 1000.0 if (7.0 <= hour < 23.0) else 1.0
        elif day < disruption_end:
            # disruption: late weekends (10am-2am)
            if is_weekend:
                lux[i] = 1000.0 if (10.0 <= hour < 24.0) or (0.0 <= hour < 2.0) else 1.0
            else:
                lux[i] = 1000.0 if (7.0 <= hour < 23.0) else 1.0
        else:
            # recovery: back to normal
            lux[i] = 1000.0 if (7.0 <= hour < 23.0) else 1.0
    
    return t, lux, baseline_days, disruption_end


def allnighter_protocol(baseline_days, recovery_days, dt):
    """create one all-nighter light schedule"""
    allnighter_day = baseline_days
    total_days = baseline_days + 1 + recovery_days
    
    t = np.arange(0.0, 24.0 * total_days, dt)
    lux = np.zeros_like(t)
    
    for i, ti in enumerate(t):
        day = int(ti // 24.0)
        hour = ti % 24.0
        
        if day < allnighter_day:
            lux[i] = 1000.0 if (7.0 <= hour < 23.0) else 1.0
        elif day == allnighter_day:
            lux[i] = 1000.0  # all-nighter = 24h light
        else:
            lux[i] = 1000.0 if (7.0 <= hour < 23.0) else 1.0
    
    return t, lux, baseline_days, allnighter_day + 1


def analyze_model(model_name, model_class, t, lux_disrupt, baseline_days, recovery_start_day):
    """run the model for both LD and darkness conditions and calculate all metrics"""
    
    # set up two conditions: LD recovery and darkness recovery
    lux_ld = lux_disrupt.copy()
    lux_dark = lux_disrupt.copy()
    start_t = recovery_start_day * 24.0
    lux_dark[t >= start_t] = 0.0  # turn off light after recovery starts
    
    # run model for LD condition
    model_ld = model_class()
    traj_ld = model_ld.integrate(t, input=lux_ld)
    
    # run model for darkness condition
    model_dark = model_class()
    traj_dark = model_dark.integrate(t, input=lux_dark)
    
    # get phase markers (CBTmin or DLMO)
    markers_ld = get_phase_markers(model_ld, traj_ld, marker=PHASE_MARKER)
    markers_dark = get_phase_markers(model_dark, traj_dark, marker=PHASE_MARKER)
    
    # convert markers to daily phases
    def daily_phases(markers):
        phases = {}
        for m in markers:
            d = int(m // 24)
            h = m % 24
            if d not in phases:
                phases[d] = []
            phases[d].append(h)
        return {d: np.mean(hs) for d, hs in phases.items()}
    
    phases_ld = daily_phases(markers_ld)
    phases_dark = daily_phases(markers_dark)
    
    # calculate baseline phase (average over baseline period)
    baseline_phases = [phases_ld[d] for d in phases_ld.keys() if d < baseline_days]
    baseline_mean = np.mean(baseline_phases) if baseline_phases else None
    
    # calculate deviations from baseline for recovery period
    days_recovery = range(recovery_start_day, recovery_start_day + 21)
    deviations_ld = []
    deviations_dark = []
    days_plot = []
    
    for d in days_recovery:
        if d in phases_ld and baseline_mean is not None:
            dev = circular_distance(phases_ld[d], baseline_mean) * 60
            deviations_ld.append(dev)
            days_plot.append(d)
        if d in phases_dark and baseline_mean is not None:
            dev = circular_distance(phases_dark[d], baseline_mean) * 60
            deviations_dark.append(dev)
    
    # get amplitude for each day
    # hannay models: R is literally the amplitude variable
    # forger/jewett: need to calculate from oscillation range
    def get_daily_amplitudes(model_obj, traj):
        amplitudes = {}
        
        # loop through each day
        for day in range(int(t[-1] // 24) + 1):
            day_start = day * 24.0
            day_end = (day + 1) * 24.0
            
            # collect all values for this day
            day_values = []
            for i in range(len(t)):
                if day_start <= t[i] < day_end:
                    _, state = traj[i]  # traj returns (time, state) so need both
                    
                    # get the right variable depending on model
                    if model_name in ["Hannay19", "Hannay19TP"]:
                        # hannay uses R directly as amplitude
                        if model_name == "Hannay19":
                            val = float(state[0])  # R
                        else:  # hannay19TP has two populations
                            val = float((state[0] + state[1]) / 2.0)  # avg R_v and R_d
                    else:
                        # forger/jewett: x is the oscillating variable
                        val = float(state[0])
                    
                    day_values.append(val)
            
            # calculate amplitude for this day
            if day_values:
                if model_name in ["Hannay19", "Hannay19TP"]:
                    # R is already amplitude so just avg over the day
                    amplitudes[day] = np.mean(day_values)
                else:
                    # forger/jewett: amplitude = half the peak-to-peak distance
                    amplitudes[day] = (max(day_values) - min(day_values)) / 2.0
        
        return amplitudes
    
    amplitudes_ld = get_daily_amplitudes(model_ld, traj_ld)
    amplitudes_dark = get_daily_amplitudes(model_dark, traj_dark)
    
    # get baseline amplitude (avg over baseline period)
    baseline_amps = [amplitudes_ld[d] for d in range(baseline_days) if d in amplitudes_ld]
    baseline_amplitude = np.mean(baseline_amps) if baseline_amps else None
    
    # calculate peak deviation
    peak_deviation = max(abs(d) for d in deviations_ld) if deviations_ld else 0
    
    # find when re-entrainment happens (3 days within +/- 20 min)
    # using linear interpolation to get hour-level precision
    def find_reentrainment_time(deviations, days, threshold=20, streak=3):
        print(f"    {model_name} LD deviations: ", end="")
        for i in range(min(7, len(deviations))):
            print(f"Day {days[i]}: {deviations[i]:.1f}min, ", end="")
        print()
        
        # find first day where it stays within threshold for 3 consecutive days
        for i in range(len(deviations) - streak + 1):
            if all(abs(deviations[i+k]) < threshold for k in range(streak)):
                first_day_idx = i
                
                # interpolate between previous and current day to get exact hour
                if first_day_idx > 0:
                    dev_before = abs(deviations[first_day_idx - 1])
                    dev_after = abs(deviations[first_day_idx])
                    
                    # linear interp: find fraction of day when crossing threshold
                    if dev_before > threshold and dev_after < threshold:
                        fraction = (dev_before - threshold) / (dev_before - dev_after)
                        day_with_fraction = days[first_day_idx - 1] + fraction
                        hours = (day_with_fraction - recovery_start_day) * 24.0
                        print(f"    {model_name}: Re-entrained at day {day_with_fraction:.2f} ({hours:.1f} hours)")
                        return hours
                
                # fallback if no interpolation possible
                hours = (days[first_day_idx] - recovery_start_day) * 24.0
                print(f"    {model_name}: Re-entrained at day {days[first_day_idx]} ({hours:.1f} hours)")
                return hours
        
        print(f"    {model_name}: Never re-entrained")
        return None
    
    reentrainment_ld_hours = find_reentrainment_time(deviations_ld, days_plot)
    reentrainment_dark_hours = find_reentrainment_time(deviations_dark, days_plot)
    
    # check when amplitude recovers to 90% of baseline
    def find_90pct_amplitude_recovery(amplitudes, baseline_amp, recovery_start, threshold=0.90):
        if baseline_amp is None or baseline_amp == 0:
            return None
        target_amp = baseline_amp * threshold
        print(f"    {model_name}: Baseline amp = {baseline_amp:.3f}, 90% target = {target_amp:.3f}")
        
        # look through recovery period to find when it hits 90%
        for day in range(recovery_start, min(recovery_start + 21, max(amplitudes.keys()) + 1)):
            if day in amplitudes:
                day_amp = amplitudes[day]
                if day_amp >= target_amp:
                    print(f"    {model_name}: Reached 90% on day {day} (amp = {day_amp:.3f})")
                    return (day - recovery_start) * 24.0
        print(f"    {model_name}: Did not reach 90% within recovery period")
        return None
    
    amp_90_ld_hours = find_90pct_amplitude_recovery(amplitudes_ld, baseline_amplitude, recovery_start_day)
    amp_90_dark_hours = find_90pct_amplitude_recovery(amplitudes_dark, baseline_amplitude, recovery_start_day)
    
    # calculate drift rate in darkness (avg change per day)
    dark_changes = []
    for i in range(len(deviations_dark) - 1):
        dark_changes.append(abs(deviations_dark[i+1] - deviations_dark[i]))
    avg_drift = np.mean(dark_changes) if dark_changes else 0
    
    # return all the metrics we calculated
    return {
        'model': model_name,
        'baseline_days': baseline_days,
        'recovery_start': recovery_start_day,
        'deviations_ld': deviations_ld,
        'deviations_dark': deviations_dark,
        'days_plot': days_plot,
        'reentrainment_ld_hours': reentrainment_ld_hours,
        'reentrainment_dark_hours': reentrainment_dark_hours,
        'amp_90_ld_hours': amp_90_ld_hours,
        'amp_90_dark_hours': amp_90_dark_hours,
        'peak_deviation': peak_deviation,
        'avg_drift': avg_drift,
        'baseline_mean': baseline_mean,
        'baseline_amplitude': baseline_amplitude,
    }


def generate_figure_4(protocol_type="social_jetlag"):
    """make figure 4 showing phase deviation analysis"""
    print(f"\n=== Generating Figure 4: {protocol_type.upper()} ===")
    
    # set up the protocol
    if protocol_type == "social_jetlag":
        t, lux, baseline_days, recovery_start = social_jetlag_protocol(
            BASELINE_WEEKS, DISRUPTION_WEEKENDS, RECOVERY_WEEKS, DT
        )
        title = 'Circadian Phase Recovery: Social Jet Lag'
        filename = 'figure_4_social_jetlag.png'
    else:  # allnighter
        t, lux, baseline_days, recovery_start = allnighter_protocol(
            baseline_days=30, recovery_days=14, dt=DT
        )
        title = 'Circadian Phase Recovery: One All-Nighter'
        filename = 'figure_4_allnighter.png'
    
    fig, axes = plt.subplots(1, 2, figsize=(16, 7))
    
    # run all models and collect results
    results = []
    for name, cls in MODELS:
        print(f"  Analyzing {name}...")
        result = analyze_model(name, cls, t, lux, baseline_days, recovery_start)
        results.append(result)
    
    # plot left panel: LD recovery
    ax_ld = axes[0]
    for result in results:
        color = MODEL_COLORS[result['model']]
        marker = MODEL_MARKERS[result['model']]
        ax_ld.plot(result['days_plot'], result['deviations_ld'],
                  marker=marker, linestyle='-', label=result['model'],
                  color=color, markersize=6, linewidth=2.5, alpha=0.85)
    
    # add reference lines
    ax_ld.axhline(0, color='green', linestyle='--', linewidth=2, alpha=0.6, label='Baseline')
    ax_ld.axhline(20, color='gray', linestyle=':', linewidth=1.5, alpha=0.5)
    ax_ld.axhline(-20, color='gray', linestyle=':', linewidth=1.5, alpha=0.5, label='±20 min threshold')
    ax_ld.axvline(results[0]['recovery_start'], color='purple', linestyle='--', 
                 linewidth=2, alpha=0.5, label='Recovery Start')
    
    ax_ld.set_xlabel('Day', fontsize=13, fontweight='bold')
    ax_ld.set_ylabel('Phase Deviation (minutes)', fontsize=13, fontweight='bold')
    ax_ld.set_title('LD Recovery: Re-entrainment with Zeitgeber', fontsize=14, fontweight='bold')
    ax_ld.legend(fontsize=10, loc='upper right', framealpha=0.95)
    ax_ld.grid(True, alpha=0.3)
    ax_ld.set_xlim([results[0]['recovery_start'] - 1, results[0]['recovery_start'] + 12])
    
    # plot right panel: darkness (free-running)
    ax_dark = axes[1]
    for result in results:
        color = MODEL_COLORS[result['model']]
        marker = MODEL_MARKERS[result['model']]
        ax_dark.plot(result['days_plot'], result['deviations_dark'],
                    marker=marker, linestyle='-', label=result['model'],
                    color=color, markersize=6, linewidth=2.5, alpha=0.85)
    
    ax_dark.axhline(0, color='green', linestyle='--', linewidth=2, alpha=0.6, label='Baseline')
    ax_dark.axvline(results[0]['recovery_start'], color='purple', linestyle='--',
                   linewidth=2, alpha=0.5, label='Recovery Start')
    
    ax_dark.set_xlabel('Day', fontsize=13, fontweight='bold')
    ax_dark.set_ylabel('Phase Deviation (minutes)', fontsize=13, fontweight='bold')
    ax_dark.set_title('Darkness: Free-Running (No Re-entrainment)', fontsize=14, fontweight='bold')
    ax_dark.legend(fontsize=10, loc='upper left', framealpha=0.95)
    ax_dark.grid(True, alpha=0.3)
    ax_dark.set_xlim([results[0]['recovery_start'] - 1, results[0]['recovery_start'] + 12])
    
    # save the figure
    plt.suptitle(title, fontsize=16, fontweight='bold', y=0.98)
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.savefig(filename)
    print(f"✓ Saved: {filename}")
    plt.close()
    
    # print summary table
    print(f"\n{'Model':<12} {'Peak Dev':<12} {'Reentrain LD':<15} {'Reentrain Dark':<15} {'90% Amp LD':<12} {'90% Amp Dark':<12}")
    print(f"{'':12} {'(min)':<12} {'(hours)':<15} {'(hours)':<15} {'(hours)':<12} {'(hours)':<12}")
    print("-" * 78)
    for r in results:
        peak = f"{r['peak_deviation']:.1f}"
        reent_ld = f"{r['reentrainment_ld_hours']:.1f}" if r['reentrainment_ld_hours'] else "N/A"
        reent_dark = f"{r['reentrainment_dark_hours']:.1f}" if r['reentrainment_dark_hours'] else "N/A"
        amp_ld = f"{r['amp_90_ld_hours']:.1f}" if r['amp_90_ld_hours'] else "N/A"
        amp_dark = f"{r['amp_90_dark_hours']:.1f}" if r['amp_90_dark_hours'] else "N/A"
        print(f"{r['model']:<12} {peak:<12} {reent_ld:<15} {reent_dark:<15} {amp_ld:<12} {amp_dark:<12}")
    print()


def main():
    """run figure 4 generation"""
    print("\n" + "="*80)
    print("FIGURE 4 GENERATION - PHASE DEVIATION ANALYSIS")
    print("="*80)
    print(f"\nPhase Marker: {PHASE_MARKER.upper()}")
    print(f"Protocol: {BASELINE_WEEKS}w baseline → {DISRUPTION_WEEKENDS}w disruption → {RECOVERY_WEEKS}w recovery")
    print("="*80)
    
    generate_figure_4(protocol_type="social_jetlag")
    
    print("\n" + "="*80)
    print("✓ FIGURE GENERATED")
    print("="*80)
    print("\nOutput:")
    print("  figure_4_social_jetlag.png - Social jet lag phase deviations")
    print("\nMetrics calculated:")
    print("  - Re-entrainment time (±20 min for 3 days)")
    print("  - 90% Amplitude recovery time")
    print("  - Peak phase deviation")
    print("="*80 + "\n")


if __name__ == "__main__":
    main()