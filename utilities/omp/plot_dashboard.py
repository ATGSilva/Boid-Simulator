import numpy as np
from matplotlib import animation
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from matplotlib import gridspec
import pandas as pd
import configparser
plt.rcParams.update({'font.size':14})


# USER SETTINGS IMPORT ==================================================================
fps = 30
skip_factor = 1 # Only change this if you use small DT, skip_factor 1 does not skip any frames

config = configparser.ConfigParser()
config.read("../../boid_data/settings.cfg")

close_alert = int(config["Neighbour Finding"]["too_close_radius"])
duration = float(config["Time"]["duration"])
dt = float(config["Time"]["timestep"])
wall_bounds = [int(config["Bounds"]["wall_bound_lower"]), int(config["Bounds"]["wall_bound_upper"])]
boids = int(config["Program"]["number_boids"])
threads = int(config["Program"]["number_procs"])
# =======================================================================================


def PlotBox(ax, wall_bounds):
    # Plot the wiremesh of bounding box
    line_ll = [wall_bounds[0], wall_bounds[0]]
    line_lu = [wall_bounds[0], wall_bounds[1]]
    line_uu = [wall_bounds[1], wall_bounds[1]]
    boxcol = (0,0,0,0.4)
    # Plot X axis lines
    ax.plot3D(line_lu , line_ll , line_ll , color=boxcol)
    ax.plot3D(line_lu, line_uu, line_ll, color=boxcol)
    ax.plot3D(line_lu, line_ll, line_uu, color=boxcol)
    ax.plot3D(line_lu, line_uu, line_uu, color=boxcol)
    # Plot Y axis lines
    ax.plot3D(line_ll , line_lu , line_ll , color=boxcol)
    ax.plot3D(line_uu, line_lu, line_ll, color=boxcol)
    ax.plot3D(line_ll, line_lu, line_uu, color=boxcol)
    ax.plot3D(line_uu, line_lu, line_uu, color=boxcol)
    # Plot Z axis lines
    ax.plot3D(line_ll , line_ll , line_lu , color=boxcol)
    ax.plot3D(line_ll, line_uu, line_lu, color=boxcol)
    ax.plot3D(line_uu, line_ll, line_lu, color=boxcol)
    ax.plot3D(line_uu, line_uu, line_lu, color=boxcol)
    

def animate(frame):
    global frame_id
    frame_id += skip_factor
    ax1.view_init(20, (0.2*frame_id % 360))
    if not frame_id % (skip_factor*fps):
        print("Processing frames. Current Frame: ", frame_id, "of", int(duration/dt))
    scatter._offsets3d = (data.loc[data['Frame']==frame_id, 'X'],
                          data.loc[data['Frame']==frame_id, 'Y'],
                          data.loc[data['Frame']==frame_id, 'Z'])
    momentum.set_data(times[:frame_id], mave[:frame_id])
    pos_plot.set_data(pave[0:2, :frame_id])
    pos_plot.set_3d_properties(pave[2, :frame_id])
    fc.set_data(times[:frame_id], flockcoeffs["Coeff"][:frame_id])
    

filename = "B" + str(boids) + "Thr" + str(threads)
frame_count = int(duration//(dt*skip_factor))
data = pd.read_csv("../../boid_data/"+filename+".csv", sep=',')
bound_width = wall_bounds[1] - wall_bounds[0]
plot_bounds = [wall_bounds[0]-(bound_width*0.25), wall_bounds[1]+(bound_width*0.25)]

flockcoeffs = pd.read_csv("../../boid_data/FC_B{0}Thr{1}.csv".format(boids, threads), sep=",", names=["Frame", "Coeff"])

frames_m = [i for i in range(int(duration/dt))]
times = [f/max(frames_m) * duration for f in frames_m]
m0= []
mave = []
pavex = []
pavey = []
pavez = []
for f in frames_m:
    vel_mag = (data["VelX"][f*boids]**2 + data["VelY"][f*boids]**2 + data["VelZ"][f*boids]**2)**0.5
    mom = data["Mass"][f] * vel_mag
    m0.append(mom)
    vel_mag = (data["VelX"][f*boids:f*boids+boids].mean()**2 + 
               data["VelY"][f*boids:f*boids+boids].mean()**2 + 
               data["VelZ"][f*boids:f*boids+boids].mean()**2)**0.5
    mom = data["Mass"][f] * vel_mag
    mave.append(mom)
    pavex.append(data["X"][f*boids:f*boids+boids].mean())
    pavey.append(data["Y"][f*boids:f*boids+boids].mean())
    pavez.append(data["Z"][f*boids:f*boids+boids].mean())
pave = np.array([pavex, pavey, pavez])

# Setup plot for animation
plt.ion()
frame_id = 0
fig = plt.figure()
gs = fig.add_gridspec(7,7)

# 3D Plot Setup -------------------------------------------------------------------------
ax1 = fig.add_subplot(gs[:6,:5], projection="3d")
ax1.set_xlim3d(plot_bounds[0], plot_bounds[1])
ax1.set_ylim3d(plot_bounds[0], plot_bounds[1])
ax1.set_zlim3d(plot_bounds[0], plot_bounds[1])
ax1.set_title("Boid Simulation: {0} Boids".format(boids))
ax1.set_xlabel("x")
ax1.set_ylabel("y")
ax1.set_zlabel("z")
DPI = fig.get_dpi()
fig.set_size_inches(1500.0/float(DPI),1080.0/float(DPI))

PlotBox(ax1, wall_bounds)

# assign initial positions - ie create frame 1
scatter = ax1.scatter(data.loc[data['Frame']==frame_id, 'X'],
                     data.loc[data['Frame']==frame_id, 'Y'],
                     data.loc[data['Frame']==frame_id, 'Z'],
                     marker='o', edgecolor='k', lw=0.4, s=2.5, c="dodgerblue")
pos_plot = ax1.plot(pave[0, :frame_id], pave[1, :frame_id], pave[2, :frame_id], c="crimson", alpha=0.25)[0]

# Lower Momentum Plot Setup -------------------------------------------------------------
ax2 = fig.add_subplot(gs[6:,1:4])
ax2.set_xlim(right=max(times))
ax2.set_ylim(top=60)
ax2.set_ylabel("Average\nMomentum\n(kg m/s)")
ax2.set_xlabel("Time (s)")
momentum, = ax2.plot(times[:frame_id], mave[:frame_id], lw=2, c="dodgerblue")

# Right Flocking Coeff Plot setup -------------------------------------------------------
ax3 = fig.add_subplot(gs[:6, 5:])
ax3.set_xlim(right=max(times))
ax3.set_ylim(bottom = flockcoeffs["Coeff"].min()-0.05, top=flockcoeffs["Coeff"].max()+0.05)
ax3.set_ylabel("Flocking Coefficient")
ax3.set_xlabel("Time (s)")
fc, = ax3.plot(times[:frame_id], flockcoeffs["Coeff"][:frame_id], lw=2, c="dodgerblue")
ax3.yaxis.tick_right()
ax3.yaxis.set_label_position("right")

# Animation -----------------------------------------------------------------------------
# This must be performed in main body to work!
anim = animation.FuncAnimation(fig, animate, frames=frame_count, interval=100)

Writer = animation.writers['ffmpeg']
writer = Writer(fps=fps, metadata=dict(artist='Me'), bitrate=5000)
anim.save("../../animations/DB"+filename+".mp4", writer=writer)
print("Render complete.")