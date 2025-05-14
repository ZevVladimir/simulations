import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import glob  # To read multiple files

# Load all files
def load_files(pattern="hydro_data/hydro_output_*.txt"):
    files = sorted(glob.glob(pattern))  # Get all files matching the pattern
    data = []
    for file in files:
        timestep_data = np.loadtxt(file)  # Load data from each file
        data.append(timestep_data)
    return data

# Load data from all files
data = load_files()

# Initialize the figure
fig, ax = plt.subplots(2, 2, figsize=(10, 8))

lines = []
titles = ["Energy", "Density", "Velocity", "Pressure"]
colors = ["blue", "orange", "green", "magenta"]

ax = ax.flat

# Set up each subplot
for i in range(len(titles)):
    ax[i].set_xlim(0, data[0][:, 0].max())  # X range based on position
    ax[i].set_ylim(0, np.max([step[:, i + 1].max() for step in data]))  # Y range based on max value
    ax[i].set_ylabel(titles[i])
    line, = ax[i].plot([], [], color=colors[i], label=titles[i])
    lines.append(line)
    ax[i].legend()

ax[-1].set_xlabel("Position")

# Update function for animation
def update(frame):
    for i, line in enumerate(lines):
        line.set_data(data[frame][:, 0], data[frame][:, i + 1])  # Update position and variable
    return lines

# Create the animation
ani = FuncAnimation(fig, update, frames=len(data), interval=100, blit=True)

# Save or show the animation
ani.save("simulation.gif", fps=10, dpi=200)  # Save as MP4 (requires ffmpeg)
plt.show()  # Or display the animation
