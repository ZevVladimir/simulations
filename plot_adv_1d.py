import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import glob
import os
import imageio.v2 as imageio  # For saving GIFs
from matplotlib.animation import FuncAnimation

# Directory containing the files
data_dir = "/home/zvladimi/simulations/data_adv_1d/"
file_pattern = os.path.join(data_dir, "adv1d_step_*.csv")

# Find and sort all files numerically
files = sorted(glob.glob(file_pattern), key=lambda x: int(x.split("_step_")[-1].split(".csv")[0]))

# Define colormap
cmap = plt.cm.viridis
num_files = len(files)

# Create figure
fig, ax = plt.subplots(figsize=(8, 5))
ax.set_xlabel("x")
ax.set_ylabel("u")
ax.set_title("1D Advection Over Time")
ax.set_xlim(-0.1, 10.1)  # Adjust based on your domain
ax.set_ylim(-1.1, 1.1)   # Adjust based on expected u values

# Initialize scatter plot
scatter = ax.scatter([], [], c=[], cmap=cmap, edgecolor='k')

# Update function for animation
def update(i):
    df = pd.read_csv(files[i], delimiter=",")
    x = df["x"].to_numpy()
    u = df["u"].to_numpy()
    
    scatter.set_offsets(np.column_stack((x, u)))  # Update positions
    scatter.set_array(np.full_like(x, i / num_files))  # Update colors
    
    ax.set_title(f"1D Advection Step {i*5}")
    return scatter,

# Create animation
ani = FuncAnimation(fig, update, frames=len(files), interval=100, blit=True)

# Save as GIF
gif_path = os.path.join(data_dir, "advection_sim.gif")
ani.save(gif_path, writer="pillow", fps=10)

print(f"GIF saved at: {gif_path}")
