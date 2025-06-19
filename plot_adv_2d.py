import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import glob
import os

# Grid size (must match your simulation)
NX, NY = 100, 100  # example values — replace with yours

# Load all file paths and sort by step number
files = sorted(glob.glob("data_adv_2d/adv2d_step_*.csv"),
               key=lambda f: int(os.path.splitext(f)[0].split('_')[-1]))

# Load first frame to get shape
first_data = np.loadtxt(files[0]).reshape((NX, NY))

# Set up the plot
fig, ax = plt.subplots()
im = ax.imshow(first_data, origin='lower', cmap='viridis', extent=[0, NY, 0, NX])
cbar = fig.colorbar(im)
cbar.set_label("u")

# Update function for animation
def update(frame):
    data = np.loadtxt(files[frame]).reshape((NX, NY))
    im.set_array(data)
    ax.set_title(f"Step {frame}")
    return [im]

# Create animation
ani = animation.FuncAnimation(fig, update, frames=len(files), interval=100, blit=True)

# Show animation
plt.show()
