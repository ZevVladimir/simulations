import pandas as pd
import matplotlib.pyplot as plt
import json

snapshot = "9980"
data_path = "/home/zvladimi/ASTR415-Fall24/TERM-PROJECTs/BL_JA_ZM_ZV_project/data_BH_N1000/"
snapshot_path = data_path + "galaxy_step_" + snapshot + ".csv"
config_path = data_path + "sim_params.json"
data = pd.read_csv(snapshot_path)

with open(config_path, "r") as param_file:
    parameters = json.load(param_file)

g1_num_particles = parameters["g1_num_particles"]
g2_num_particles = parameters["g2_num_particles"]

if len(data) < g1_num_particles + g2_num_particles:
    raise ValueError(f"Snapshot {snapshot} has fewer particles than expected.")

first_half = data.iloc[1:g1_num_particles+1]
second_half = data.iloc[g1_num_particles+1:]

fig, ax = plt.subplots(1)
ax.clear()
ax.scatter(first_half['x'], first_half['y'], s=1, alpha=0.7, color='blue', label='Galaxy 1')
ax.scatter(second_half['x'], second_half['y'], s=1, alpha=0.7, color='red', label='Galaxy 2')
ax.scatter(data.iloc[0]['x'], data.iloc[0]['y'], s=30, alpha=1, color='cyan', marker="+", label='BH 1')
ax.scatter(data.iloc[g1_num_particles+1]['x'], data.iloc[g1_num_particles+1]['y'], s=30, alpha=1, color='orange', marker="+", label='BH 2')

ax.set_title(f"Snapshot: " + snapshot)
ax.set_xlabel("X")
ax.set_ylabel("Y")
ax.set_xlim(-200, 200)
ax.set_ylim(-200, 200)
ax.legend(loc="upper right")
fig.savefig(data_path + "static_snap_"+snapshot+".png")