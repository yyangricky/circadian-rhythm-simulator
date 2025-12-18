"""
Figure 4 Generation - Phase Deviation Analysis
Generates phase deviation figures for social jet lag and one all-nighter protocols
"""

import numpy as np
import matplotlib.pyplot as plt
from models import Forger99, Jewett99, Hannay19, Hannay19TP
from metrics import get_phase_markers, circular_distance

# Publication-quality settings
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

# Parameters
DT = 0.10
BASELINE_WEEKS = 4
DISRUPTION_WEEKENDS = 3
RECOVERY_WEEKS = 3
PHASE_MARKER = "cbt"

MODELS = [
    ("Forger99", Forger99),
    ("Jewett99", Jewett99),
    ("Hannay19", Hannay19),
    ("Hannay19TP", Hannay19TP),
]

MODEL_COLORS = {
    "Forger99": "#2E86AB",
    "Jewett99": "#A23B72",
    "Hannay19": "#F18F01",
    "Hannay19TP": "#C73E1D",
}

MODEL_MARKERS = {
    "Forger99": "o",
    "Jewett99": "s",
    "Hannay19": "^",
    "Hannay19TP": "D",
}


def social_jetlag_protocol(baseline_weeks, disruption_weekends, recovery_weeks, dt):
    """Social Jet Lag Protocol"""
    baseline_days = baseline_weeks * 7
    disruption_days = disruption_weekends * 7
    recovery_days = recovery_weeks * 7
    total_days = baseline_days + disruption_days + recovery_days
    
    t = np.arange(0.0, 24.0 * total_days, dt)
    lux = np.zeros_like(t)
    
    disruption_start = baseline_days
    disruption_end = baseline_days + disruption_days
    
    for i, ti in enumerate(t):
        day = int(ti // 24.0)
        hour = ti % 24.0
        day_of_week = day % 7
        is_weekend = (day_of_week == 5 or day_of_week == 6)
        
        if day < disruption_start:
            lux[i] = 1000.0 if (7.0 <= hour < 23.0) else 1.0
        elif day < disruption_end:
            if is_weekend:
                lux[i] = 1000.0 if (10.0 <= hour < 24.0) or (0.0 <= hour < 2.0) else 1.0
            else:
                lux[i] = 1000.0 if (7.0 <= hour < 23.0) else 1.0
        else:
            lux[i] = 1000.0 if (7.0 <= hour < 23.0) else 1.0
    
    return t, lux, baseline_days, disruption_end


def allnighter_protocol(baseline_days, recovery_days, dt):
    """One all-nighter protocol"""
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
            lux[i] = 1000.0  # 24h light
        else:
            lux[i] = 1000.0 if (7.0 <= hour < 23.0) else 1.0
    
    return t, lux, baseline_days, allnighter_day + 1


def analyze_model(model_name, model_class, t, lux_disrupt, baseline_days, recovery_start_day):
    """Complete analysis for one model"""
    
    lux_ld = lux_disrupt.copy()
    lux_dark = lux_disrupt.copy()
    start_t = recovery_start_day * 24.0
    lux_dark[t >= start_t] = 0.0
    
    model_ld = model_class()
    traj_ld = model_ld.integrate(t, input=lux_ld)
    
    model_dark = model_class()
    traj_dark = model_dark.integrate(t, input=lux_dark)
    
    markers_ld = get_phase_markers(model_ld, traj_ld, marker=PHASE_MARKER)
    markers_dark = get_phase_markers(model_dark, traj_dark, marker=PHASE_MARKER)
    
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
    
    baseline_phases = [phases_ld[d] for d in phases_ld.keys() if d < baseline_days]
    baseline_mean = np.mean(baseline_phases) if baseline_phases else None
    
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
    
    # Calculate metrics
    peak_deviation = max(abs(d) for d in deviations_ld) if deviations_ld else 0
    
    def find_recovery_hours(deviations, days, threshold=20, streak=3):
        for i in range(len(deviations) - streak + 1):
            if all(abs(deviations[i+k]) < threshold for k in range(streak)):
                return (days[i] - recovery_start_day) * 24.0
        return None
    
    recovery_ld_hours = find_recovery_hours(deviations_ld, days_plot)
    
    dark_changes = []
    for i in range(len(deviations_dark) - 1):
        dark_changes.append(abs(deviations_dark[i+1] - deviations_dark[i]))
    avg_drift = np.mean(dark_changes) if dark_changes else 0
    
    return {
        'model': model_name,
        'baseline_days': baseline_days,
        'recovery_start': recovery_start_day,
        'deviations_ld': deviations_ld,
        'deviations_dark': deviations_dark,
        'days_plot': days_plot,
        'recovery_ld_hours': recovery_ld_hours,
        'peak_deviation': peak_deviation,
        'avg_drift': avg_drift,
        'baseline_mean': baseline_mean,
    }


def generate_figure_4(protocol_type="social_jetlag"):
    """Generate Figure 4: Phase Deviation Analysis
    
    Args:
        protocol_type: Either "social_jetlag" or "allnighter"
    """
    print(f"\n=== Generating Figure 4: {protocol_type.upper()} ===")
    
    # Create protocol
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
    
    # Analyze all models
    results = []
    for name, cls in MODELS:
        print(f"  Analyzing {name}...")
        result = analyze_model(name, cls, t, lux, baseline_days, recovery_start)
        results.append(result)
    
    # Left panel: LD Recovery
    ax_ld = axes[0]
    for result in results:
        color = MODEL_COLORS[result['model']]
        marker = MODEL_MARKERS[result['model']]
        ax_ld.plot(result['days_plot'], result['deviations_ld'],
                  marker=marker, linestyle='-', label=result['model'],
                  color=color, markersize=6, linewidth=2.5, alpha=0.85)
    
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
    
    # Right panel: Darkness (Free-Running)
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
    
    plt.suptitle(title, fontsize=16, fontweight='bold', y=0.98)
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.savefig(filename)
    print(f"✓ Saved: {filename}")
    plt.close()
    
    # Print statistics
    print(f"\n{'Model':<12} {'Peak Dev (min)':<16} {'Reentrain (days)':<18} {'Drift Rate (min/d)':<18}")
    print("-" * 64)
    for r in results:
        peak = f"{r['peak_deviation']:.1f}"
        reentrain = f"{r['recovery_ld_hours']/24:.2f}" if r['recovery_ld_hours'] else "Not achieved"
        drift = f"{r['avg_drift']:.2f}"
        print(f"{r['model']:<12} {peak:<16} {reentrain:<18} {drift:<18}")
    print()


def main():
    """Generate both figure 4 versions"""
    print("\n" + "="*80)
    print("FIGURE 4 GENERATION - PHASE DEVIATION ANALYSIS")
    print("="*80)
    print(f"\nPhase Marker: {PHASE_MARKER.upper()}")
    print("="*80)
    
    # Generate social jet lag figure
    print("\n[1/2] Social Jet Lag Protocol")
    print(f"Protocol: {BASELINE_WEEKS}w baseline → {DISRUPTION_WEEKENDS}w disruption → {RECOVERY_WEEKS}w recovery")
    generate_figure_4(protocol_type="social_jetlag")
    
    # Generate all-nighter figure
    print("\n[2/2] One All-Nighter Protocol")
    print("Protocol: 30d baseline → 1 all-nighter → 14d recovery")
    generate_figure_4(protocol_type="allnighter")
    
    print("\n" + "="*80)
    print("✓ ALL FIGURES GENERATED")
    print("="*80)
    print("\nOutputs:")
    print("  1. figure_4_social_jetlag.png - Social jet lag phase deviations")
    print("  2. figure_4_allnighter.png - All-nighter phase deviations")
    print("="*80 + "\n")


if __name__ == "__main__":
    main()