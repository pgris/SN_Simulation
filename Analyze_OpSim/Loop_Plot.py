import os

for val in ['feature_baseline_10yrs','feature_rolling_half_mask_10yrs','feature_rolling_twoThird_10yrs','alt_sched','alt_sched_rolling']:
    cmd='python Plots.py --simu_name '+val
    os.system(cmd)
