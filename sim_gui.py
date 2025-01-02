import tkinter as tk
from tkinter import ttk, messagebox, simpledialog
import subprocess
import os
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import pandas as pd
import time
import sys
import numpy as np
import json
global data_dir
data_dir = "./data/"

def remove_outliers(data, columns, threshold=3):
    """
    Removes outliers based on the IQR method.
    
    Args:
        data (pd.DataFrame): The dataset containing particle coordinates.
        columns (list): List of column names to check for outliers (e.g., ['x', 'y']).
        threshold (float): The IQR multiplier to define outliers (default: 1.5).
    
    Returns:
        pd.DataFrame: DataFrame with outliers removed.
    """
    mask = np.ones(len(data), dtype=bool)  # Start with all True (no outliers)
    for col in columns:
        Q1 = data[col].quantile(0.25)  # First quartile (25th percentile)
        Q3 = data[col].quantile(0.75)  # Third quartile (75th percentile)
        IQR = Q3 - Q1  # Interquartile range
        lower_bound = Q1 - threshold * IQR
        upper_bound = Q3 + threshold * IQR
        mask &= (data[col] >= lower_bound) & (data[col] <= upper_bound)  # Update mask
        
    return data[mask]

# Function to create tooltip
class ToolTip:
    def __init__(self, widget, text):
        self.widget = widget
        self.text = text
        self.tooltip = None
        self.widget.bind("<Enter>", self.show_tooltip)
        self.widget.bind("<Leave>", self.hide_tooltip)

    def show_tooltip(self, event):
        x, y, _, _ = self.widget.bbox("insert")
        x += self.widget.winfo_rootx() + 25
        y += self.widget.winfo_rooty() + 25
        self.tooltip = tk.Toplevel(self.widget)
        self.tooltip.wm_overrideredirect(1)  # Remove window decorations
        self.tooltip.wm_geometry(f"+{x}+{y}")
        label = tk.Label(self.tooltip, text=self.text, background="yellow", relief="solid", borderwidth=1)
        label.pack()

    def hide_tooltip(self, event):
        if self.tooltip:
            self.tooltip.destroy()


def save_to_gif():
    if not os.path.exists(data_dir):
        messagebox.showerror("Error", f"Data directory '{data_dir}' not found!")
        return

    snapshot_files = sorted(
        [f for f in os.listdir(data_dir) if f.startswith("galaxy_step_")],
        key=lambda x: int(x.split('_')[2].split('.')[0])  # Extract numeric part and sort
    )

    if not snapshot_files:
        messagebox.showinfo("Info", "No snapshots found in the data directory.")
        return

    g1_num_particles = int(g1_num_particles_var.get())
    g2_num_particles = int(g2_num_particles_var.get())
    
    fig, ax = plt.subplots(figsize=(8, 8))

    def update(frame_index):
        snapshot = snapshot_files[frame_index]
        snapshot_path = os.path.join(data_dir, snapshot)
        data = pd.read_csv(snapshot_path)

        if len(data) < g1_num_particles + g2_num_particles:
            raise ValueError(f"Snapshot {snapshot} has fewer particles than expected.")

        first_half = data.iloc[1:g1_num_particles+1]
        second_half = data.iloc[g1_num_particles+1:]

        ax.clear()
        ax.scatter(first_half['x'], first_half['y'], s=1, alpha=0.7, color='blue', label='Galaxy 1')
        ax.scatter(second_half['x'], second_half['y'], s=1, alpha=0.7, color='red', label='Galaxy 2')
        ax.scatter(data.iloc[0]['x'], data.iloc[0]['y'], s=30, alpha=1, color='cyan', marker="+", label='BH 1')
        ax.scatter(data.iloc[g1_num_particles+1]['x'], data.iloc[g1_num_particles+1]['y'], s=30, alpha=1, color='orange', marker="+", label='BH 2')
        snap_num = int(snapshot.split('_')[2].split('.')[0])
        ax.set_title(f"Snapshot: {snap_num}")
        ax.set_xlabel("X")
        ax.set_ylabel("Y")
        ax.set_xlim(-500, 500)
        ax.set_ylim(-500, 500)
        ax.legend(loc="upper right")

    try:
        ani = animation.FuncAnimation(fig, update, frames=len(snapshot_files), repeat=False)
        gif_path = os.path.join(data_dir, "simulation.gif")
        ani.save(gif_path, writer='pillow', fps=30)  # Save the animation as a GIF
        messagebox.showinfo("Success", f"Animation saved as GIF at {gif_path}")
    except Exception as e:
        messagebox.showerror("Error", f"An error occurred while saving GIF: {e}")

# Function to run the simulation
def run_simulation():
    global simulation_process
    try:
        # Clear the data directory
        if os.path.exists(data_dir):
            for file in os.listdir(data_dir):
                file_path = os.path.join(data_dir, file)
                if os.path.isfile(file_path) or os.path.islink(file_path):
                    os.unlink(file_path)  # Delete the file
                elif os.path.isdir(file_path):
                    os.rmdir(file_path)  # Delete the folder (if necessary)
        else:
            os.makedirs(data_dir)  # Create the directory if it doesn't exist
        # Get values from the GUI inputs
        g1_num_particles = int(g1_num_particles_var.get())
        g2_num_particles = int(g2_num_particles_var.get())
        num_timesteps = int(num_timesteps_var.get())
        timestep_size = float(timestep_size_var.get())
        step_output = int(step_output_var.get())
        g1_x = float(g1_x_var.get())
        g1_y = float(g1_y_var.get())
        g1_vx = float(g1_vx_var.get())
        g1_vy = float(g1_vy_var.get())
        g2_x = float(g2_x_var.get())
        g2_y = float(g2_y_var.get())
        g2_vx = float(g2_vx_var.get())
        g2_vy = float(g2_vy_var.get())

        # Save parameters to a file in the ./data folder
        parameters = {
            "g1_num_particles": g1_num_particles,
            "g2_num_particles": g2_num_particles,
            "num_timesteps": num_timesteps,
            "timestep_size": timestep_size,
            "step_output": step_output,
            "g1_x": g1_x,
            "g1_y": g1_y,
            "g1_vx": g1_vx,
            "g1_vy": g1_vy,
            "g2_x": g2_x,
            "g2_y": g2_y,
            "g2_vx": g2_vx,
            "g2_vy": g2_vy
        }

        # Save the parameters to a JSON file (simulation_parameters.json)
        with open(os.path.join(data_dir, "sim_params.json"), "w") as param_file:
            json.dump(parameters, param_file, indent=4)

        # Construct the command
        command = [
            "./sim.o",
            str(g1_num_particles),
            str(g2_num_particles),
            str(num_timesteps),
            str(timestep_size), 
            str(step_output),
            str(g1_x),
            str(g1_y),
            str(g1_vx),
            str(g1_vy),
            str(g2_x),
            str(g2_y),
            str(g2_vx),
            str(g2_vy),
        ]

        # Create a popup window to indicate the simulation is running
        popup = tk.Toplevel(root)
        popup.title("Simulation in Progress")
        label = ttk.Label(popup, text="Simulation is running, please wait...")
        label.pack(padx=20, pady=20)
        popup.geometry("400x150")
        popup.grab_set()  # Make it modal

        # List to track observed snapshots
        observed_snapshots = set()

        # Function to monitor the data directory
        def monitor_snapshots():
            nonlocal observed_snapshots
            if os.path.exists(data_dir):
                snapshot_files = sorted(
                    [f for f in os.listdir(data_dir) if f.startswith("galaxy_step_")],
                    key=lambda x: int(x.split('_')[2].split('.')[0])  # Extract numeric part and sort
                )
                new_snapshots = set(snapshot_files) - observed_snapshots
                if new_snapshots:
                    observed_snapshots.update(new_snapshots)
                    current_snapshot = max(new_snapshots, key=lambda x: int(x.split('_')[2].split('.')[0]))
                    label.config(text=f"Processing Snapshot: {current_snapshot}")
            
            if simulation_process.poll() is None:  # Process is still running
                root.after(100, monitor_snapshots)
            else:
                # Destroy the popup once the process finishes
                popup.destroy()
                play_simulation()  # Automatically start playback

        # Run the simulation process
        simulation_process = subprocess.Popen(command)

        # Start monitoring snapshots
        monitor_snapshots()
        
    except Exception as e:
        messagebox.showerror("Error", f"An error occurred: {e}")

def continue_simulation():
    try:
        # Load the previously saved parameters from the JSON file
        param_file_path = os.path.join(data_dir, "sim_params.json")
        
        if not os.path.exists(param_file_path):
            messagebox.showerror("Error", "No saved simulation parameters found.")
            return

        with open(param_file_path, "r") as param_file:
            parameters = json.load(param_file)

        # Prompt for the additional number of timesteps
        additional_timesteps = simpledialog.askinteger("Additional Timesteps", 
                                                      "How many more timesteps to run for?", 
                                                      parent=root, 
                                                      minvalue=1)
        if additional_timesteps is None:
            # User canceled the input
            return
        
        # Update the num_timesteps with the additional timesteps
        updated_num_timesteps = parameters["num_timesteps"] + additional_timesteps
        parameters["num_timesteps"] = updated_num_timesteps

        # Save the updated parameters back to the JSON file
        with open(param_file_path, "w") as param_file:
            json.dump(parameters, param_file, indent=4)

        # Construct the command with the updated parameters
        command = [
            "./sim.o", 
            str(parameters["g1_num_particles"]),
            str(parameters["g2_num_particles"]),
            str(additional_timesteps),
            str(parameters["timestep_size"]),
            str(parameters["step_output"]),
            str(parameters["g1_x"]),
            str(parameters["g1_y"]),
            str(parameters["g1_vx"]),
            str(parameters["g1_vy"]),
            str(parameters["g2_x"]),
            str(parameters["g2_y"]),
            str(parameters["g2_vx"]),
            str(parameters["g2_vy"]),
            "--resume",
        ]

        # Create a popup window to indicate the simulation is running
        popup = tk.Toplevel(root)
        popup.title("Simulation in Progress")
        label = ttk.Label(popup, text="Simulation is running, please wait...")
        label.pack(padx=20, pady=20)
        popup.geometry("400x150")
        popup.grab_set()  # Make it modal

        # Function to monitor the data directory (same as in run_simulation)
        observed_snapshots = set()

        def monitor_snapshots():
            nonlocal observed_snapshots
            if os.path.exists(data_dir):
                snapshot_files = sorted(
                    [f for f in os.listdir(data_dir) if f.startswith("galaxy_step_")],
                    key=lambda x: int(x.split('_')[2].split('.')[0])  # Extract numeric part and sort
                )
                new_snapshots = set(snapshot_files) - observed_snapshots
                if new_snapshots:
                    observed_snapshots.update(new_snapshots)
                    current_snapshot = max(new_snapshots, key=lambda x: int(x.split('_')[2].split('.')[0]))
                    label.config(text=f"Processing Snapshot: {current_snapshot}")
            
            if simulation_process.poll() is None:  # Process is still running
                root.after(100, monitor_snapshots)
            else:
                # Destroy the popup once the process finishes
                popup.destroy()
                play_simulation()  # Automatically start playback

        # Run the simulation process with the updated parameters
        simulation_process = subprocess.Popen(command)

        # Start monitoring snapshots
        monitor_snapshots()
        
    except Exception as e:
        messagebox.showerror("Error", f"An error occurred: {e}")

def play_simulation():
    if not os.path.exists(data_dir):
        messagebox.showerror("Error", f"Data directory '{data_dir}' not found!")
        return

    snapshot_files = sorted(
        [f for f in os.listdir(data_dir) if f.startswith("galaxy_step_")],
        key=lambda x: int(x.split('_')[2].split('.')[0])  # Extract numeric part and sort
    )

    if not snapshot_files:
        messagebox.showinfo("Info", "No snapshots found in the data directory.")
        return

    g1_num_particles = int(g1_num_particles_var.get())  # Get the number of particles input from the user
    g2_num_particles = int(g2_num_particles_var.get())
    
    def update_plot(index=0):
        if index >= len(snapshot_files):
                return  # Stop when all snapshots are displayed

        snapshot = snapshot_files[index]
        snapshot_path = os.path.join(data_dir, snapshot)

        data = pd.read_csv(snapshot_path)

        if len(data) < g1_num_particles + g2_num_particles:
            messagebox.showerror("Error", f"Snapshot {snapshot} has fewer particles ({len(data)}) than expected ({g1_num_particles + g2_num_particles}).")
            return

        # Split the data
        first_half = data.iloc[1:g1_num_particles+1]
        second_half = data.iloc[g1_num_particles+1:]
    
        # first_half = remove_outliers(first_half, ['x', 'y'])
        # second_half = remove_outliers(second_half, ['x', 'y'])

        # Clear and update the plot
        ax.clear()
        ax.scatter(first_half['x'], first_half['y'], s=1, alpha=0.7, color='blue', label='Galaxy 1')
        ax.scatter(second_half['x'], second_half['y'], s=1, alpha=0.7, color='red', label='Galaxy 2')
        ax.scatter(data.iloc[0]['x'], data.iloc[0]['y'], s=30, alpha=1, color='cyan', marker="+", label='BH 1')
        ax.scatter(data.iloc[g1_num_particles+1]['x'], data.iloc[g1_num_particles+1]['y'], s=30, alpha=1, color='orange', marker="+", label='BH 2')
        snap_num = int(snapshot.split('_')[2].split('.')[0])
        ax.set_title(f"Snapshot: {snap_num}")
        ax.set_xlabel("X")
        ax.set_ylabel("Y")
        ax.set_xlim(-500,500)
        ax.set_ylim(-500,500)
        ax.legend(loc="upper right")
        canvas.draw()
        
        delay = int(playback_speed_var.get() * 1000)  # Convert seconds to milliseconds
        root.after(delay, update_plot, index + 1)
        
    update_plot()

# Handle window close event
def on_close():
    """Clean up resources and exit."""
    if simulation_process:
        simulation_process.terminate()
        simulation_process.wait()  # Ensure the process is cleaned up

    if messagebox.askokcancel("Quit", "Do you want to exit the program?"):
        root.destroy()
        sys.exit(0)  # Ensure Python exits completely

def exit_fullscreen(event=None):
    root.attributes('-fullscreen', False)

# Create the GUI window
root = tk.Tk()
root.title("N-Body Simulation")
root.attributes('-fullscreen', True)

root.bind("<Escape>", exit_fullscreen)

# Attach the close event handler
root.protocol("WM_DELETE_WINDOW", on_close)

playback_speed_frame = ttk.LabelFrame(root, text="Playback Speed")
playback_speed_frame.pack(fill="x", padx=10, pady=5)

playback_speed_var = tk.DoubleVar(value=0.3)  # Default playback speed (seconds per frame)
playback_speed_slider = ttk.Scale(
    playback_speed_frame,
    from_=0.0001,  # Minimum delay (fastest playback)
    to=2.0,  # Maximum delay (slowest playback)
    orient="horizontal",
    variable=playback_speed_var
)
playback_speed_slider.pack(fill="x", padx=10, pady=5)


# Parameters frame
params_frame = ttk.LabelFrame(root, text="Simulation Parameters")
params_frame.pack(fill="x", padx=10, pady=10)

# Create input fields
fields = [
    ("Galaxy 1 X", "0.0"),
    ("Galaxy 1 Y", "0.0"),
    ("Galaxy 1 VX", "0.0"),
    ("Galaxy 1 VY", "0.0"),
    ("Galaxy 2 X", "-50.0"),
    ("Galaxy 2 Y", "-100.0"),
    ("Galaxy 2 VX", "10.0"),
    ("Galaxy 2 VY", "20.0"),
    ("Num of Particles Galaxy 1", "1000"),
    ("Num of Particles Galaxy 2", "1000"),
    ("Number of Time Steps", "5000"),
    ("Time Step Size", "0.0005"),
    ("Step Output Rate", "10"),
]

# Variables to store input data
variables = []

# Number of columns for the grid layout
num_columns = 3

for i, (label, default) in enumerate(fields):
    # Calculate row and column positions
    row, col = divmod(i, num_columns)
    
    # Create label and entry widgets
    lbl = ttk.Label(params_frame, text=label, width=20)
    lbl.grid(row=row, column=col * 2, padx=5, pady=2, sticky="w")
    
    var = tk.StringVar(value=default)
    entry = ttk.Entry(params_frame, textvariable=var)
    entry.grid(row=row, column=col * 2 + 1, padx=5, pady=2, sticky="ew")
    
    variables.append(var)

# Adjust column weights for resizing
for col in range(num_columns * 2):
    params_frame.columnconfigure(col, weight=1)

# Map variables to meaningful names
(
    g1_x_var, g1_y_var, g1_vx_var, g1_vy_var,
    g2_x_var, g2_y_var, g2_vx_var, g2_vy_var,
    g1_num_particles_var, g2_num_particles_var, num_timesteps_var, 
    timestep_size_var, step_output_var,
) = variables


button_frame = ttk.Frame(root)
button_frame.pack(pady=10)

# Run simulation button
run_button = ttk.Button(button_frame, text="Run C Simulation", command=run_simulation)
run_button.pack(side=tk.LEFT, padx=5)  # Align buttons horizontally
run_tooltip = ToolTip(run_button, "Click to start the C simulation from the parameters inputted.")

# Continue simulation button
continue_button = ttk.Button(button_frame, text="Continue C Simulation", command=continue_simulation)
continue_button.pack(side=tk.LEFT, padx=5)  # Align buttons horizontally
continue_tooltip = ToolTip(continue_button, "Click to continue the C simulation from the\nparameters and data in the ./data folder.")

# Play existing data button
play_button = ttk.Button(button_frame, text="Play Existing Data", command=play_simulation)
play_button.pack(side=tk.LEFT, padx=5)  # Align buttons horizontally
play_tooltip = ToolTip(play_button, "Click to play the output simulation data in the ./data folder.")

# Add the "Save to GIF" button to the GUI
save_gif_button = ttk.Button(root, text="Save to GIF", command=save_to_gif)
save_gif_button.pack(pady=10)  # Use pack instead of grid


# Plot frame
plot_frame = ttk.LabelFrame(root, text="Simulation Output")
plot_frame.pack(fill="both", expand=True, padx=10, pady=10)

fig, ax = plt.subplots(figsize=(6, 6))
ax.set_aspect('equal', adjustable='box')  # Ensure the plot is square
canvas = FigureCanvasTkAgg(fig, master=plot_frame)
canvas.get_tk_widget().pack(fill="both", expand=True)

# Global variable for tracking the simulation process
simulation_process = None

# Start the GUI loop
root.mainloop()
