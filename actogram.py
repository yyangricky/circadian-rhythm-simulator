import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as lines
'''
from circadian.plots import Actogram
from circadian.lights import LightSchedule
from circadian.models import Forger99, Jewett99, Hannay19, Hannay19TP
'''
from recovery_time import recovery_time

from models import Forger99, Jewett99, Hannay19, Hannay19TP

from circadian.plots import Actogram
from metrics import get_phase_markers, compute_baseline_marker, reentrainment_hours
from protocols import allnighter_protocol

DT = 0.10
BASELINE_DAYS = 30
RECOVERY_DAYS = 14

time, lux_all, baseline_days, recovery_start_day = allnighter_protocol(
    baseline_days=BASELINE_DAYS,
    recovery_days=RECOVERY_DAYS,
    dt=DT,
    day_start=7.0,
    day_end=23.0,
    day_lux=250.0,
    night_lux=1.0,
    allnighter_lux=1000.0,
)
light_values = lux_all  # actogram sees the same schedule as amplitude_plot
'''
days_night = 3
days_day = 2
#slam_shift = LightSchedule.ShiftWork(lux=300.0, days_on=days_night, days_off=days_day)

total_days = 45
time = np.arange(0, 24*total_days, 0.10)
#light_values = slam_shift(time)

jet_lag_shift = LightSchedule.SocialJetlag(lux=1000.0, num_regular_days=20, num_jetlag_days=5, hours_delayed=6)
light_values = jet_lag_shift(time)
'''

f_model = Forger99()
kj_model = Jewett99()
spm_model = Hannay19()
tpm_model = Hannay19TP()

equilibration_reps = 2
initial_conditions_forger = f_model.equilibrate(time, light_values, equilibration_reps)
initial_conditions_kj = kj_model.equilibrate(time, light_values, equilibration_reps)
initial_conditions_spm = spm_model.equilibrate(time, light_values, equilibration_reps)
initial_conditions_tpm = tpm_model.equilibrate(time, light_values, equilibration_reps)

# from f_model(time, initial_conditions_forger, light_values) to  f_model.integrate(time, initial_condition=initial_conditions_forger, input=light_values)
# consistant 
trajectory_f = f_model.integrate(time, initial_condition=initial_conditions_forger, input = light_values)
trajectory_kj = kj_model.integrate(time, initial_condition=initial_conditions_kj, input = light_values)
trajectory_spm = spm_model.integrate(time, initial_condition=initial_conditions_spm, input = light_values)
trajectory_tpm = tpm_model.integrate(time, initial_condition=initial_conditions_tpm, input = light_values)
'''
dlmo_f = f_model.dlmos()
dlmo_kj = kj_model.dlmos()
dlmo_spm = spm_model.dlmos()
dlmo_tpm = tpm_model.dlmos()
'''

dlmo_f  = get_phase_markers(f_model, trajectory_f,  marker="dlmo")
cbt_f   = get_phase_markers(f_model, trajectory_f,  marker="cbt")

dlmo_kj  = get_phase_markers(kj_model, trajectory_kj,  marker="dlmo")
cbt_kj   = get_phase_markers(kj_model, trajectory_kj,  marker="cbt")
dlmo_spm  = get_phase_markers(spm_model, trajectory_spm,  marker="dlmo")
cbt_spm  = get_phase_markers(spm_model, trajectory_spm,  marker="cbt")
dlmo_tpm  = get_phase_markers(tpm_model, trajectory_tpm,  marker="dlmo")
cbt_tpm   = get_phase_markers(tpm_model, trajectory_tpm,  marker="cbt")


acto = Actogram(time, light_vals=light_values, opacity=1.0, smooth=False)
acto.plot_phasemarker(dlmo_f, color='blue')
acto.plot_phasemarker(dlmo_spm, color='darkgreen')
acto.plot_phasemarker(dlmo_tpm, color='red')
acto.plot_phasemarker(dlmo_kj, color='purple')
# legend
blue_line = lines.Line2D([], [], color='blue', label='Forger99')
green_line = lines.Line2D([], [], color='darkgreen', label='Hannay19')
red_line = lines.Line2D([], [], color='red', label='Hannay19TP')
purple_line = lines.Line2D([], [], color='purple', label='Jewett99')

#this is testing
baseline_day = 20
baseline_dlmo_f = compute_baseline_marker(dlmo_f, baseline_days=20) #dlmo_f[baseline_day]      # Forger99
baseline_dlmo_kj = compute_baseline_marker(dlmo_kj, baseline_days=20) #dlmo_kj[baseline_day]    # Jewett99
baseline_dlmo_spm = compute_baseline_marker(dlmo_spm, baseline_days=20) #dlmo_spm[baseline_day]  # Hannay19
baseline_dlmo_tpm = compute_baseline_marker(dlmo_tpm, baseline_days=20) #dlmo_tpm[baseline_day]  # Hannay19TP

# added the % 24, changed from axhline to axvline so the baselines are on the right axis (vertical lines)
acto.ax.axvline(baseline_dlmo_f % 24, color='blue', linestyle='--', alpha=0.6)
acto.ax.axvline(baseline_dlmo_spm% 24, color='darkgreen', linestyle='--', alpha=0.6)
acto.ax.axvline(baseline_dlmo_tpm% 24, color='red', linestyle='--', alpha=0.6)
acto.ax.axvline(baseline_dlmo_kj% 24, color='purple', linestyle='--', alpha=0.6)

# Compute recovery days for each model
rec_f  = recovery_time(dlmo_f, baseline_dlmo_f)
rec_kj = recovery_time(dlmo_kj, baseline_dlmo_kj)
rec_spm = recovery_time(dlmo_spm, baseline_dlmo_spm)
rec_tpm = recovery_time(dlmo_tpm, baseline_dlmo_tpm)

# Print results
print("Days until DLMO returns to within 10% of baseline:")
print(f"Forger99:   Day {rec_f}")
print(f"Jewett99:   Day {rec_kj}")
print(f"Hannay19:   Day {rec_spm}")
print(f"Hannay19TP: Day {rec_tpm}")


plt.legend(handles=[blue_line, purple_line, green_line, red_line], 
           loc='upper center', bbox_to_anchor=(0.5, 1.12), ncol=4)
plt.title("Actogram for DLMO recovery", pad=35)
plt.tight_layout()
plt.show()

