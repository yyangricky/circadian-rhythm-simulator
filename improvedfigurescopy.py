"""
IMPROVED Figure Generation for Circadian Rhythm Project
Single all-nighter protocol (Day 49) with publication-quality figures
"""

import numpy as np
import matplotlib.pyplot as plt
from models import Forger99, Jewett99, Hannay19, Hannay19TP
from metrics import get_phase_markers, circular_distance

# =========================
# Publication-quality settings
# =========================
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

# =========================
# Global parameters
# =========================
DT = 0.10
TOTAL_DAYS = 70
ALL_NIGHTER_DAY = 48
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

# =========================
# All-nighter protocol
# =========================
def all_nighter_protocol(dt, total_days, all_nighter_day):
    """
    Single all-nighter protocol.

    Normal schedule:
        Light ON 07:00–23:00

    All-nighter:
        Light ON continuously from 07:00 (day 49)
        until 23:00 (day 50)
    """
    t = np.arange(0.0, 24.0 * total_days, dt)
    lux = np.zeros_like(t)

    for i, ti in enumerate(t):
        day = int(ti // 24.0)
        hour = ti % 24.0

        if day == all_nighter_day and hour >= 7.0:
            lux[i] = 1000.0
        elif day == all_nighter_day + 1 and hour < 23.0:
            lux[i] = 1000.0
        else:
            lux[i] = 1000.0 if (7.0 <= hour < 23.0) else 1.0

    return t, lux

# =========================
# Analysis function
# =========================
def analyze_model_all_nighter(model_name, model_class):
    print(f"Analyzing {model_name}...")

    t, lux_ld = all_nighter_protocol(DT, TOTAL_DAYS, ALL_NIGHTER_DAY)

    # Darkness recovery condition
    lux_dark = lux_ld.copy()
    recovery_start_t = (ALL_NIGHTER_DAY + 1) * 24.0
    lux_dark[t >= recovery_start_t] = 0.0

    # Integrate models
    model_ld = model_class()
    traj_ld = model_ld.integrate(t, input=lux_ld)

    model_dark = model_class()
    traj_dark = model_dark.integrate(t, input=lux_dark)

    # Phase markers
    markers_ld = get_phase_markers(model_ld, traj_ld, marker=PHASE_MARKER)
    markers_dark = get_phase_markers(model_dark, traj_dark, marker=PHASE_MARKER)

    # Daily phase means
    def daily_phases(markers):
        phases = {}
        for m in markers:
            d = int(m // 24)
            phases.setdefault(d, []).append(m % 24)
        return {d: np.mean(hs) for d, hs in phases.items()}

    phases_ld = daily_phases(markers_ld)
    phases_dark = daily_phases(markers_dark)

    # Baseline phase
    baseline_phases = [phases_ld[d] for d in phases_ld if d < ALL_NIGHTER_DAY]
    baseline_mean = np.mean(baseline_phases) if baseline_phases else None

    # Recovery window
    days_recovery = range(ALL_NIGHTER_DAY + 1, ALL_NIGHTER_DAY + 22)

    deviations_ld = []
    deviations_dark = []
    days_plot = []

    for d in days_recovery:
        if d in phases_ld and baseline_mean is not None:
            deviations_ld.append(
                circular_distance(phases_ld[d], baseline_mean) * 60
            )
            days_plot.append(d)

        if d in phases_dark and baseline_mean is not None:
            deviations_dark.append(
                circular_distance(phases_dark[d], baseline_mean) * 60
            )

    # Metrics
    peak_deviation = max(abs(d) for d in deviations_ld) if deviations_ld else 0

    # Time to 50% recovery
    target_50 = peak_deviation * 0.5
    time_to_50pct = None
    for i, dev in enumerate(deviations_ld):
        if abs(dev) <= target_50:
            time_to_50pct = (days_plot[i] - (ALL_NIGHTER_DAY + 1)) * 24.0
            break

    # Re-entrainment time (±20 min for 3 days)
    def find_recovery_hours(deviations, days, threshold=20, streak=3):
        for i in range(len(deviations) - streak + 1):
            if all(abs(deviations[i + k]) < threshold for k in range(streak)):
                return (days[i] - (ALL_NIGHTER_DAY + 1)) * 24.0
        return None

    recovery_ld_hours = find_recovery_hours(deviations_ld, days_plot)

    # Drift in darkness
    dark_changes = [
        abs(deviations_dark[i + 1] - deviations_dark[i])
        for i in range(len(deviations_dark) - 1)
    ]
    avg_drift = np.mean(dark_changes) if dark_changes else 0

    return {
        "model": model_name,
        "t": t,
        "lux_ld": lux_ld,
        "lux_dark": lux_dark,
        "days_plot": days_plot,
        "deviations_ld": deviations_ld,
        "deviations_dark": deviations_dark,
        "peak_deviation": peak_deviation,
        "time_to_50pct": time_to_50pct,
        "recovery_ld_hours": recovery_ld_hours,
        "avg_drift": avg_drift,
    }

# =========================
# Figure 1: Protocol overview
# =========================
def generate_figure_1():
    result = analyze_model_all_nighter("Forger99", Forger99)
    t_days = result["t"] / 24.0

    fig, axes = plt.subplots(2, 1, figsize=(16, 9))

    for idx, (lux, title, color) in enumerate([
        (result["lux_ld"], "LD Recovery", "orange"),
        (result["lux_dark"], "Darkness Recovery", "navy"),
    ]):
        ax = axes[idx]
        ax.plot(t_days, lux, color=color, linewidth=1.0)

        ax.axvline(ALL_NIGHTER_DAY, color="red", linestyle="--",
                   linewidth=2.5, label="All-Nighter")
        ax.axvline(ALL_NIGHTER_DAY + 1, color="green", linestyle="--",
                   linewidth=2.5, label="Recovery Start")

        ax.set_ylabel("Light Intensity (lux)")
        ax.set_title(title, fontweight="bold")
        ax.set_xlim([45, 55])
        ax.set_ylim([-50, 1150])
        ax.grid(True, alpha=0.3)
        ax.legend()

        if idx == 1:
            ax.set_xlabel("Day")

    plt.tight_layout()
    plt.savefig("figure_1_protocol_overview.png")
    plt.close()

# =========================
# Figures 4 & 5 + table
# =========================
def generate_all():
    results = [analyze_model_all_nighter(name, cls) for name, cls in MODELS]

    # Trajectory comparison
    fig, axes = plt.subplots(1, 2, figsize=(16, 7))

    for r in results:
        axes[0].plot(
            r["days_plot"], r["deviations_ld"],
            marker=MODEL_MARKERS[r["model"]],
            color=MODEL_COLORS[r["model"]],
            label=r["model"],
            linewidth=2.5,
        )

    axes[0].axhline(0, color="green", linestyle="--")
    axes[0].set_title("LD Recovery")
    axes[0].set_xlabel("Day")
    axes[0].set_ylabel("Phase Deviation (min)")
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)

    for r in results:
        axes[1].plot(
            r["days_plot"], r["deviations_dark"],
            marker=MODEL_MARKERS[r["model"]],
            color=MODEL_COLORS[r["model"]],
            label=r["model"],
            linewidth=2.5,
        )

    axes[1].set_title("Darkness (Free-Running)")
    axes[1].set_xlabel("Day")
    axes[1].set_ylabel("Phase Deviation (min)")
    axes[1].legend()
    axes[1].grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig("figure_4_trajectory_comparison.png")
    plt.close()

    # Print table
    print("\nMODEL SUMMARY")
    print("-" * 70)
    print(f"{'Model':<12} {'Peak Dev':<10} {'50% (d)':<10} {'Reent (d)':<10} {'Drift':<10}")
    for r in results:
        print(
            f"{r['model']:<12} "
            f"{r['peak_deviation']:<10.1f} "
            f"{(r['time_to_50pct'] or 0)/24:<10.2f} "
            f"{(r['recovery_ld_hours'] or 0)/24:<10.2f} "
            f"{r['avg_drift']:<10.2f}"
        )

# =========================
# Main
# =========================
def main():
    generate_figure_1()
    generate_all()
    print("\n✓ All-nighter analysis complete (Day 49)")

if __name__ == "__main__":
    main()
