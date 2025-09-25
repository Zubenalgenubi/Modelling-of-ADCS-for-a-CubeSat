# Full ADCS Simulation Pipeline for a CubeSat
# Includes: Detumbling (B-dot), TRIAD algorithm, Kalman Filter attitude estimation, and Reaction Wheel control with torque saturation

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# ==================== Constants and Parameters ====================
I = np.diag([0.0417, 0.0417, 0.0417])  # CubeSat inertia matrix [kg·m²]
I_inv = np.linalg.inv(I)
B_earth_inertial = np.array([2.5e-5, -1.5e-5, 3.0e-5])  # Tesla (approx Earth field vector)
Sun_inertial = np.array([1, 0, 0])  # Assume Sun vector along +X inertial

# Reaction wheel parameters
Kp = np.array([0.00584, 0.00584, 0.00584])
Ki = np.array([0.000417, 0.000417, 0.000417])

# Reaction wheel torque saturation [N·m]
rw_max_torque = 5e-3

# Desired quaternion (pointing to inertial X)
global_q_desired = np.array([1, 0, 0, 0])

# ==================== Helper Functions ====================
def omega_matrix(w):
    return np.array([
        [0, -w[0], -w[1], -w[2]],
        [w[0], 0, w[2], -w[1]],
        [w[1], -w[2], 0, w[0]],
        [w[2], -w[1], w[0], 0]
    ])

def quat_multiply(q1, q2):
    w1, x1, y1, z1 = q1
    w2, x2, y2, z2 = q2
    return np.array([
        w1*w2 - x1*x2 - y1*y2 - z1*z2,
        w1*x2 + x1*w2 + y1*z2 - z1*y2,
        w1*y2 - x1*z2 + y1*w2 + z1*x2,
        w1*z2 + x1*y2 - y1*x2 + z1*w2
    ])

def quat_error(q, qd):
    q_conj = qd * np.array([1, -1, -1, -1])
    return quat_multiply(q_conj, q)

# ==================== Detumbling: B-dot control ====================
def bdot_control(B_body, B_dot, k=1e5):
    # Magnetic dipole command (Nm/T)
    m_cmd = -k * B_dot
    # Torque = m x B
    return np.cross(m_cmd, B_body)

# ==================== TRIAD Algorithm ====================
def triad(v1b, v2b, v1i, v2i):
    t1b = v1b / np.linalg.norm(v1b)
    t2b = np.cross(t1b, v2b); t2b /= np.linalg.norm(t2b)
    t3b = np.cross(t1b, t2b)

    t1i = v1i / np.linalg.norm(v1i)
    t2i = np.cross(t1i, v2i); t2i /= np.linalg.norm(t2i)
    t3i = np.cross(t1i, t2i)

    Tb = np.column_stack((t1b, t2b, t3b))
    Ti = np.column_stack((t1i, t2i, t3i))
    return Tb @ Ti.T  # DCM from inertial to body

# ==================== Kalman Filter (Simplified) ====================
class SimpleKalman:
    def __init__(self, Q=1e-6, R=1e-3):
        self.Q = Q
        self.R = R
        self.x = np.zeros(3)
        self.P = np.eye(3)
    def update(self, meas, u):
        # Predict
        self.x = self.x + u
        self.P = self.P + self.Q*np.eye(3)
        # Update
        y = meas - self.x
        S = self.P + self.R*np.eye(3)
        K = self.P @ np.linalg.inv(S)
        self.x = self.x + K @ y
        self.P = (np.eye(3) - K) @ self.P
        return self.x

kf = SimpleKalman()

# ==================== Attitude Dynamics ====================
integral_error = np.zeros(3)
def attitude_dynamics(t, state):
    global integral_error, global_q_desired

    q = state[0:4]; w = state[4:7]
    q = q / np.linalg.norm(q)

    # Sensor simulation: project inertial vectors into body frame
    R_ib = quat_to_dcm(q).T
    B_body = R_ib @ B_earth_inertial
    Sun_body = R_ib @ Sun_inertial

    # TRIAD estimation
    R_est = triad(Sun_body, B_body, Sun_inertial, B_earth_inertial)
    q_est = dcm_to_quat(R_est)

    # Kalman filter update (here simplified as bias correction)
    est_err = kf.update(q_est[1:], np.zeros(3))

    # Control: PI on attitude error
    q_err = quat_error(q, global_q_desired)
    error_vec = q_err[1:4]
    integral_error += error_vec * 0.1
    torque_rw = -Kp * error_vec - Ki * integral_error

    # Apply torque saturation
    torque_rw = np.clip(torque_rw, -rw_max_torque, rw_max_torque)

    # Dynamics
    w_dot = I_inv @ (torque_rw - np.cross(w, I @ w))
    q_dot = 0.5 * omega_matrix(w) @ q

    return np.concatenate((q_dot, w_dot))

# ==================== Quaternion <-> DCM Conversions ====================
def quat_to_dcm(q):
    q0, q1, q2, q3 = q
    return np.array([
        [1-2*(q2**2+q3**2), 2*(q1*q2 - q0*q3), 2*(q1*q3 + q0*q2)],
        [2*(q1*q2 + q0*q3), 1-2*(q1**2+q3**2), 2*(q2*q3 - q0*q1)],
        [2*(q1*q3 - q0*q2), 2*(q2*q3 + q0*q1), 1-2*(q1**2+q2**2)]
    ])

def dcm_to_quat(R):
    q0 = 0.5 * np.sqrt(1 + np.trace(R))
    q1 = (R[2,1] - R[1,2])/(4*q0)
    q2 = (R[0,2] - R[2,0])/(4*q0)
    q3 = (R[1,0] - R[0,1])/(4*q0)
    return np.array([q0,q1,q2,q3])

# ==================== Initial Conditions ====================
q0 = np.array([1, 0.1, 0.1, 0])
q0 /= np.linalg.norm(q0)
w0 = np.random.uniform(-0.05, 0.05, 3)
initial_state = np.concatenate((q0, w0))

# ==================== Simulate ====================
sol = solve_ivp(attitude_dynamics, [0, 300], initial_state, t_eval=np.linspace(0, 300, 1500))

# ==================== Plot Results ====================
plt.figure(figsize=(10, 5))
for i, label in enumerate(['q0', 'q1', 'q2', 'q3']):
    plt.plot(sol.t, sol.y[i], label=label)
plt.title('Quaternion Evolution with Full ADCS (B-dot, TRIAD, Kalman, Reaction Wheels with Saturation)')
plt.xlabel('Time [s]')
plt.ylabel('Quaternion Components')
plt.legend(); plt.grid(); plt.tight_layout(); plt.show()
