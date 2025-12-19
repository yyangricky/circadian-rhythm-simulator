"""
actogram visualization
shows phase markers with light/dark schedule bars
"""

import numpy as np
import matplotlib.pyplot as plt
from models import Forger99, Jewett99, Hannay19, Hannay19TP
from metrics import get_phase_markers

# plot settings
plt.rcParams['font.size'] = 11
plt.rcParams['savefig.dpi'] = 300

# protocol parameters
DT = 0.10
BASELINE_WEEKS = 4
DISRUPTION_WEEKENDS = 3
RECOVERY_WEEKS = 3
PHASE_MARKER = "dlmo"  # can change to "cbt" if needed

# models with colors
MODELS = [
    ("Forger99", Forger99, "#4472C4"),      # blue
    ("Jewett99", Jewett99, "#A23B72"),      # purple 
    ("Hannay19", Hannay19, "#70AD47"),      # green
    ("Hannay19TP", Hannay19TP, "#C00000"),  # red
]


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
    
    for i, ti in enumerate(t):
        day = int(ti // 24.0)
        hour = ti % 24.0
        day_of_week = day % 7
        is_weekend = (day_of_week == 5 or day_of_week == 6)  # sat/sun
        
        if day < disruption_start:
            # baseline: normal schedule
            lux[i] = 1000.0 if (7.0 <= hour < 23.0) else 1.0
        elif day < disruption_end:
            # disruption: late weekends
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
            # baseline
            lux[i] = 1000.0 if (7.0 <= hour < 23.0) else 1.0
        elif day == allnighter_day:
            # all-nighter = 24h light
            lux[i] = 1000.0
        else:
            # recovery
            lux[i] = 1000.0 if (7.0 <= hour < 23.0) else 1.0
    
    return t, lux, baseline_days, allnighter_day + 1


def create_actogram(protocol_type="social_jetlag"):
    """make the actogram figure"""
    
    print("\n" + "="*70)
    print(f"GENERATING ACTOGRAM - {protocol_type.upper()}")
    print("="*70)
    
    # set up protocol
    if protocol_type == "social_jetlag":
        t, lux, baseline_days, recovery_start = social_jetlag_protocol(
            BASELINE_WEEKS, DISRUPTION_WEEKENDS, RECOVERY_WEEKS, DT
        )
        title = 'Actogram for Social Jet Lag Protocol'
        filename = 'actogram_social_jetlag.png'
    else:  # allnighter
        t, lux, baseline_days, recovery_start = allnighter_protocol(
            baseline_days=30, recovery_days=14, dt=DT
        )
        title = 'Actogram for One All-Nighter Protocol'
        filename = 'actogram_allnighter.png'
    
    total_days = int(t[-1] / 24.0)
    
    # set up figure with two panels
    fig = plt.figure(figsize=(10, 12))
    
    # left panel: phase markers
    ax_main = plt.subplot2grid((1, 20), (0, 0), colspan=14)
    
    # right panel: light schedule
    ax_light = plt.subplot2grid((1, 20), (0, 15), colspan=5, sharey=ax_main)
    
    # run all models
    print("\nAnalyzing models...")
    all_markers = {}
    for model_name, model_class, color in MODELS:
        print(f"  - {model_name}")
        
        model = model_class()
        traj = model.integrate(t, input=lux)
        markers = get_phase_markers(model, traj, marker=PHASE_MARKER)
        
        all_markers[model_name] = (markers, color)
    
    # plot phase markers
    for model_name, (markers, color) in all_markers.items():
        marker_days = [int(m // 24) for m in markers]
        marker_times = [m % 24 for m in markers]
        
        ax_main.plot(marker_times, marker_days, '-', color=color,
                    label=model_name, linewidth=2.0, alpha=0.85)
    
    # add phase boundary lines
    ax_main.axhline(baseline_days, color='red', linestyle='--', 
                   linewidth=1.5, alpha=0.4)
    ax_main.axhline(recovery_start, color='green', linestyle='--',
                   linewidth=1.5, alpha=0.4)
    
    # format main plot
    ax_main.set_xlim([18, 6])  # show 18:00 to 06:00 (wrapping around midnight)
    ax_main.set_ylim([total_days, 0])  # day 0 at top
    ax_main.set_xlabel('ZT', fontsize=13, fontweight='bold')
    ax_main.set_ylabel('Days', fontsize=13, fontweight='bold')
    
    # x-axis ticks
    ax_main.set_xticks([18, 21, 24, 3, 6])
    ax_main.set_xticklabels(['18', '21', '24', '3', '6'])
    
    # legend at top
    ax_main.legend(loc='upper center', bbox_to_anchor=(0.5, 1.08), 
                  ncol=4, fontsize=11, frameon=True)
    
    ax_main.grid(True, alpha=0.2, axis='y')
    ax_main.set_axisbelow(True)
    
    # draw light/dark schedule bars
    for day in range(total_days):
        # figure out light schedule for this day
        if protocol_type == "social_jetlag":
            day_of_week = day % 7
            is_weekend = (day_of_week == 5 or day_of_week == 6)
            
            if day < baseline_days:
                # normal schedule
                light_periods = [(7, 23)]
            elif day < recovery_start:
                # disruption phase
                if is_weekend:
                    # late weekend
                    light_periods = [(10, 24), (0, 2)]
                else:
                    # normal weekday
                    light_periods = [(7, 23)]
            else:
                # recovery
                light_periods = [(7, 23)]
        else:  # allnighter
            if day < baseline_days:
                light_periods = [(7, 23)]
            elif day == baseline_days:
                # all-nighter day
                light_periods = [(0, 24)]
            else:
                # recovery
                light_periods = [(7, 23)]
        
        # draw black background (darkness)
        ax_light.barh(day, 24, left=0, height=1.0, 
                     color='black', edgecolor='black', linewidth=0.3)
        
        # draw white bars for light
        for start, end in light_periods:
            ax_light.barh(day, end - start, left=start, height=1.0,
                        color='white', edgecolor='black', linewidth=0.3)
    
    # format light schedule panel
    ax_light.set_xlim([0, 24])
    ax_light.set_ylim([total_days, 0])
    ax_light.set_xticks([0, 6, 12, 18, 24])
    ax_light.set_xticklabels(['0', '6', '12', '18', '24'], fontsize=9)
    ax_light.yaxis.set_visible(False)
    ax_light.spines['left'].set_visible(False)
    ax_light.spines['top'].set_visible(False)
    ax_light.spines['right'].set_visible(False)
    ax_light.set_xlabel('Light\nSchedule', fontsize=10, fontweight='bold')
    
    # save figure
    fig.suptitle(title, fontsize=14, fontweight='bold', y=0.99)
    
    plt.tight_layout(rect=[0, 0, 1, 0.97])
    plt.savefig(filename, bbox_inches='tight', dpi=300)
    print(f"\n✓ Saved: {filename}")
    
    # print baseline stats
    print("\n" + "="*70)
    print(f"BASELINE {PHASE_MARKER.upper()} TIMES (average)")
    print("="*70)
    
    for model_name, (markers, _) in all_markers.items():
        baseline_markers = [m % 24 for m in markers if m < baseline_days * 24]
        if baseline_markers:
            avg_time = np.mean(baseline_markers)
            hours = int(avg_time)
            minutes = int((avg_time - hours) * 60)
            print(f"{model_name:12} - {hours:02d}:{minutes:02d}")
    
    print("="*70)
    print()
    
    plt.close()


if __name__ == "__main__":
    print("\n" + "="*70)
    print("ACTOGRAM GENERATION - CIRCADIAN RHYTHM PROJECT")
    print("="*70)
    print(f"\nPhase Marker: {PHASE_MARKER.upper()}")
    print(f"Protocol: {BASELINE_WEEKS}w baseline → {DISRUPTION_WEEKENDS}w disruption → {RECOVERY_WEEKS}w recovery")
    print("="*70)
    
    # run actogram generation
    create_actogram(protocol_type="social_jetlag")
    
    print("\n" + "="*70)
    print("ACTOGRAM COMPLETE!")
    print("="*70)
    print("\nOutput:")
    print("  actogram_social_jetlag.png - Social jet lag (late weekends)")
    print("="*70)