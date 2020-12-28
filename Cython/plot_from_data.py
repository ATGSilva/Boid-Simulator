from matplotlib import animation
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import axes3d
import pandas as pd

def animate(frame):
    global TIME
    TIME += 1
    scatter._offsets3d = (position_data.loc[position_data['Time']==TIME, 'X Position'],
                          position_data.loc[position_data['Time']==TIME, 'Y Position'],
                          position_data.loc[position_data['Time']==TIME, 'Z Position'])


position_data = pd.read_csv("positions.csv", sep=',')

# Setup plot for animation
plt.ion()
fig = plt.figure()
ax = axes3d.Axes3D(fig)
ax.set_xlim3d(-200, 1700)
ax.set_ylim3d(-200, 1700)
ax.set_zlim3d(-200, 1700)

TIME = 0
# assign initial positions - ie create frame 1
scatter = ax.scatter(position_data.loc[position_data['Time']==TIME, 'X Position'],
                     position_data.loc[position_data['Time']==TIME, 'Y Position'],
                     position_data.loc[position_data['Time']==TIME, 'Z Position'],
                     marker='o', edgecolor='k', lw=0.5, s=2)

# This must be performed in main body to work!
anim = animation.FuncAnimation(fig, animate, frames=120, interval=100)

Writer = animation.writers['ffmpeg']
writer = Writer(fps=24, metadata=dict(artist='Me'), bitrate=5000)
anim.save("boids_3d_mpi.mp4", writer=writer)