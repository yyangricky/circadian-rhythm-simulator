import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as lines
from circadian.plots import Actogram
from circadian.lights import LightSchedule
from circadian.models import Forger99, Jewett99, Hannay19, Hannay19TP
from recovery_time import recovery_time


days_night = 3
days_day = 2
#slam_shift = LightSchedule.ShiftWork(lux=300.0, days_on=days_night, days_off=days_day)

total_days = 45
time = np.arange(0, 24*total_days, 0.10)
#light_values = slam_shift(time)

jet_lag_shift = LightSchedule.SocialJetlag(lux=150.0, num_regular_days=20, num_jetlag_days=5, hours_delayed=6)
light_values = jet_lag_shift(time)


f_model = Forger99()
kj_model = Jewett99()
spm_model = Hannay19()
tpm_model = Hannay19TP()

equilibration_reps = 2
initial_conditions_forger = f_model.equilibrate(time, light_values, equilibration_reps)
initial_conditions_kj = kj_model.equilibrate(time, light_values, equilibration_reps)
initial_conditions_spm = spm_model.equilibrate(time, light_values, equilibration_reps)
initial_conditions_tpm = tpm_model.equilibrate(time, light_values, equilibration_reps)

trajectory_f = f_model(time, initial_conditions_forger, light_values)
trajectory_kj = kj_model(time, initial_conditions_kj, light_values)
trajectory_spm = spm_model(time, initial_conditions_spm, light_values)
trajectory_tpm = tpm_model(time, initial_conditions_tpm, light_values)

dlmo_f = f_model.dlmos()
dlmo_kj = kj_model.dlmos()
dlmo_spm = spm_model.dlmos()
dlmo_tpm = tpm_model.dlmos()

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
baseline_dlmo_f = dlmo_f[baseline_day]      # Forger99
baseline_dlmo_kj = dlmo_kj[baseline_day]    # Jewett99
baseline_dlmo_spm = dlmo_spm[baseline_day]  # Hannay19
baseline_dlmo_tpm = dlmo_tpm[baseline_day]  # Hannay19TP

acto.ax.axhline(baseline_dlmo_f, color='blue', linestyle='--', alpha=0.6)
acto.ax.axhline(baseline_dlmo_spm, color='darkgreen', linestyle='--', alpha=0.6)
acto.ax.axhline(baseline_dlmo_tpm, color='red', linestyle='--', alpha=0.6)
acto.ax.axhline(baseline_dlmo_kj, color='purple', linestyle='--', alpha=0.6)

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
plt.title("Actogram for a Simulated Shift Worker", pad=35)
plt.tight_layout()
plt.show()

