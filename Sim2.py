import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# Quaternion to DCM
def quat_to_dcm(q):
    q0, q1, q2, q3 = q
    return np.array([
        [1-2*(q2**2+q3**2), 2*(q1*q2 - q0*q3), 2*(q1*q3 + q0*q2)],
        [2*(q1*q2 + q0*q3), 1-2*(q1**2+q3**2), 2*(q2*q3 - q0*q1)],
        [2*(q1*q3 - q0*q2), 2*(q2*q3 + q0*q1), 1-2*(q1**2+q2**2)]
    ])

# Create a cube representing CubeSat
def draw_cubesat(ax, R):
    # Cube vertices (unit cube, centered at origin)
    r = 0.5
    vertices = np.array([[x, y, z] for x in [-r, r] for y in [-r, r] for z in [-r, r]])
    vertices = (R @ vertices.T).T  # rotate
    edges = [(0,1),(0,2),(0,4),(1,3),(1,5),(2,3),(2,6),(3,7),
             (4,5),(4,6),(5,7),(6,7)]
    for e in edges:
        ax.plot(*zip(vertices[e[0]], vertices[e[1]]), color="b")

# Draw reaction wheels as cylinders on each axis
def draw_reaction_wheels(ax, R):
    # Reaction wheels aligned with body X, Y, Z
    axes = np.eye(3)
    colors = ['r','g','k']
    for axis, c in zip(axes, colors):
        # Cylinder parameters
        theta = np.linspace(0, 2*np.pi, 30)
        z = np.linspace(-0.1, 0.1, 2)
        theta, z = np.meshgrid(theta, z)
        x = 0.05*np.cos(theta)
        y = 0.05*np.sin(theta)
        z = z
        cyl = np.vstack((x.flatten(), y.flatten(), z.flatten()))
        # Align cylinder to axis
        from scipy.spatial.transform import Rotation as Rot
        R_align = Rot.align_vectors([axis],[np.array([0,0,1])])[0].as_matrix()
        cyl = (R_align @ cyl).reshape(3, *x.shape)
        # Rotate with satellite
        cyl = (R @ cyl.reshape(3,-1)).reshape(3, *x.shape)
        ax.plot_surface(cyl[0], cyl[1], cyl[2], color=c, alpha=0.7)

# Example quaternion trajectory (rotating about Z)
timesteps = 100
q_traj = []
for t in range(timesteps):
    angle = 0.02 * t
    q = np.array([np.cos(angle/2), 0, 0, np.sin(angle/2)])  # rotation about Z
    q_traj.append(q/np.linalg.norm(q))

# Animate
fig = plt.figure(figsize=(6,6))
ax = fig.add_subplot(111, projection='3d')
ax.set_xlim([-1,1]); ax.set_ylim([-1,1]); ax.set_zlim([-1,1])
ax.set_xlabel("X"); ax.set_ylabel("Y"); ax.set_zlabel("Z")

def update(frame):
    ax.cla()
    ax.set_xlim([-1,1]); ax.set_ylim([-1,1]); ax.set_zlim([-1,1])
    q = q_traj[frame]
    R = quat_to_dcm(q)
    draw_cubesat(ax, R)
    draw_reaction_wheels(ax, R)
    ax.set_title(f"CubeSat + Reaction Wheels (frame {frame})")

ani = FuncAnimation(fig, update, frames=timesteps, interval=100)
ani.save("cubesat_reaction_wheels.gif", writer="pillow")
plt.show()
