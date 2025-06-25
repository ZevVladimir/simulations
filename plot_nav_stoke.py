import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import glob
import os

# Grid size (must match your simulation)
NX, NY = 64,64  # example values — replace with yours
step_size = 10

# Load all file paths and sort by step number
vx_files = sorted(glob.glob("nav_stoke_output/u_step_*.csv"),
                    key=lambda f: int(os.path.splitext(f)[0].split('_')[-1]))
vy_files = sorted(glob.glob("nav_stoke_output/v_step_*.csv"),
                    key=lambda f: int(os.path.splitext(f)[0].split('_')[-1]))
w_files = sorted(glob.glob("nav_stoke_output/omega_step_*.csv"),
                    key=lambda f: int(os.path.splitext(f)[0].split('_')[-1]))


# Load first frame
first_vx = np.loadtxt(vx_files[0], delimiter=",")
first_vy = np.loadtxt(vy_files[0], delimiter=",")
first_w = np.loadtxt(w_files[0], delimiter=",")


Lx = 1
Ly = 1

dx = Lx / NX
dy = Ly / NY

x = np.linspace(0, Lx, NX, endpoint=False)
y = np.linspace(0, Ly, NY, endpoint=False)
X, Y = np.meshgrid(x, y)



# Set up the figure with 3 panels
fig, axs = plt.subplots(2, 2, figsize=(15, 15),tight_layout=True)
ax1, ax2, ax3, ax4 = axs[0,0], axs[0,1], axs[1,0], axs[1,1]

im1 = ax1.imshow(first_w, origin='upper', cmap='coolwarm', extent=[0, NY, 0, NX])
cbar1 = fig.colorbar(im1, ax=ax1)
cbar1.set_label("w")
ax1.set_title("Vorticity")

quiv = ax2.quiver(X, Y, first_vx, first_vy, scale=50, color='blue')
ax2.invert_yaxis()
ax2.set_title("Velocity Field")
ax2.set_xlabel("vₓ")
ax2.set_ylabel("vᵧ")
ax2.set_aspect('equal')

im3 = ax3.imshow(first_vx, origin='upper', cmap='coolwarm', extent=[0, NY, 0, NX])
cbar3 = fig.colorbar(im3, ax=ax3)
cbar3.set_label("vₓ")
ax3.set_title("Velocity X Component")

im4 = ax4.imshow(first_vy, origin='upper', cmap='coolwarm', extent=[0, NY, 0, NX])
cbar4 = fig.colorbar(im4, ax=ax4)
cbar4.set_label("vᵧ")
ax4.set_title("Velocity Y Component")


# --- Update function for animation ---
def update(frame):   
    w = np.loadtxt(w_files[frame], delimiter=",")
    vx = np.loadtxt(vx_files[frame], delimiter=",")
    vy = np.loadtxt(vy_files[frame], delimiter=",")

    im1.set_array(w)
    im3.set_array(vx)
    im4.set_array(vy)
    
    quiv.set_UVC(vx, vy)

    ax1.set_title(f"Vorticity (Step {frame*step_size})")
    ax2.set_title(f"Velocity Field (Step {frame*step_size})")
    ax3.set_title(f"Velocity X (vₓ) (Step {frame*step_size})")
    ax4.set_title(f"Velocity Y (vᵧ) (Step {frame*step_size})")
    
    return [im1, im3, im4, quiv]

# Create animation
ani = animation.FuncAnimation(fig, update, frames=len(w_files), interval=0)

# Show
plt.tight_layout()
plt.show()
