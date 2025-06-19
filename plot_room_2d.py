import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import glob
import os

# Grid size (must match your simulation)
NX, NY = 11, 11  # example values — replace with yours

# Load all file paths and sort by step number
temp_files = sorted(glob.glob("data_room_2d/room_2d_step_*.csv"),
                    key=lambda f: int(os.path.splitext(f)[0].split('_')[-1]))
vel_files = sorted(glob.glob("data_room_2d/velocity_step_*.csv"),
                   key=lambda f: int(f.split('_')[-1].split('.')[0]))

# Load first frame
first_temp = np.loadtxt(temp_files[0], delimiter=",")
first_vel = np.loadtxt(vel_files[0], delimiter=",")
vx = first_vel[:, ::2]
vy = first_vel[:, 1::2]

# Set up the figure with 3 panels
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15, 5))

# --- Temperature panel ---
im1 = ax1.imshow(first_temp, origin='lower', cmap='viridis', extent=[0, NY, 0, NX])
cbar1 = fig.colorbar(im1, ax=ax1)
cbar1.set_label("Temperature (u)")
ax1.set_title("Temperature Field")

# Velocity arrows on top of temperature
X, Y = np.meshgrid(np.arange(NY), np.arange(NX))
quiv = ax1.quiver(X, Y, vx, vy, color='white', scale=50)

# --- vx panel ---
im2 = ax2.imshow(vx, origin='lower', cmap='coolwarm', extent=[0, NY, 0, NX], vmin=-1, vmax=1)
cbar2 = fig.colorbar(im2, ax=ax2)
cbar2.set_label("vₓ")
ax2.set_title("Velocity X Component")

# --- vy panel ---
im3 = ax3.imshow(vy, origin='lower', cmap='coolwarm', extent=[0, NY, 0, NX], vmin=-1, vmax=1)
cbar3 = fig.colorbar(im3, ax=ax3)
cbar3.set_label("vᵧ")
ax3.set_title("Velocity Y Component")

# --- Update function for animation ---
def update(frame):
    temp = np.loadtxt(temp_files[frame], delimiter=",")
    im1.set_array(temp)
    
    vel = np.loadtxt(vel_files[frame], delimiter=",")
    vx = vel[:, ::2]
    vy = vel[:, 1::2]

    quiv.set_UVC(vx, vy)
    im2.set_array(vx)
    im3.set_array(vy)

    ax1.set_title(f"Temperature Field (Step {frame})")
    ax2.set_title(f"Velocity X (vₓ) (Step {frame})")
    ax3.set_title(f"Velocity Y (vᵧ) (Step {frame})")
    
    return [im1, im2, im3, quiv]

# Create animation
ani = animation.FuncAnimation(fig, update, frames=len(temp_files), interval=500)

# Show
plt.tight_layout()
plt.show()
