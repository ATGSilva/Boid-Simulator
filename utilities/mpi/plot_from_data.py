from matplotlib import animation
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import axes3d
import pandas as pd
import configparser


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
threads = int(config["Program"]["number_procs_per_node"])
nodes = int(config["Program"]["number_nodes"])
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
    if not frame_id % (skip_factor*fps):
        print("Processing frames. Current Frame: ", frame_id, "of", int(duration/dt))
    ax.view_init(20, (0.2*frame_id % 360))
    scatter._offsets3d = (position_data.loc[position_data['Frame']==frame_id, 'X'],
                          position_data.loc[position_data['Frame']==frame_id, 'Y'],
                          position_data.loc[position_data['Frame']==frame_id, 'Z'])



filename = "B" + str(boids) + "Nodes" + str(nodes) + "Thr" + str(threads)
print("Generating animation...")

frame_count = int(duration//(dt*skip_factor))
position_data = pd.read_csv("../../boid_data/"+filename+".csv", sep=',')
bound_width = wall_bounds[1] - wall_bounds[0]
plot_bounds = [wall_bounds[0]-(bound_width*0.25), wall_bounds[1]+(bound_width*0.25)]

# Setup plot for animation
plt.ion()
fig = plt.figure()
ax = axes3d.Axes3D(fig)
ax.set_xlim3d(plot_bounds[0], plot_bounds[1])
ax.set_ylim3d(plot_bounds[0], plot_bounds[1])
ax.set_zlim3d(plot_bounds[0], plot_bounds[1])
DPI = fig.get_dpi()
fig.set_size_inches(1280.0/float(DPI),960.0/float(DPI))

PlotBox(ax, wall_bounds)

frame_id = 0
# assign initial positions - ie create frame 1
scatter = ax.scatter(position_data.loc[position_data['Frame']==frame_id, 'X'],
                     position_data.loc[position_data['Frame']==frame_id, 'Y'],
                     position_data.loc[position_data['Frame']==frame_id, 'Z'],
                     marker='o', edgecolor='k', lw=0.4, s=2.5)

# This must be performed in main body to work!
anim = animation.FuncAnimation(fig, animate, frames=frame_count, interval=100)

Writer = animation.writers['ffmpeg']
writer = Writer(fps=fps, metadata=dict(artist='Me'), bitrate=5000)
anim.save("../../animations/"+filename+".mp4", writer=writer)
print("Render complete.")