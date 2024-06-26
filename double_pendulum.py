import numpy as np
import matplotlib.pyplot as plt

g = 9.81 # m/s^2
m1 = 1.0 # mass rod 1 in kg
m2 = 1.0 # mass rod 2 in kg
l1 = 1.0 # length rod 1 in m
l2 = 1.0 # length rod 2 in m

# 4th order Runge-Kutta
def RK4(y, t, dt):
    k1 = dYdt(y) * dt
    k2 = dYdt(y + k1 / 2) * dt
    k3 = dYdt(y + k2 / 2) * dt
    k4 = dYdt(y + k3) * dt

    return y + (k1 + 2 * k2 + 2 * k3 + k4) / 6

# System of differential equations 
def dYdt(y):
    theta1, theta2, omega1, omega2 = y
    dTheta1_dt = omega1
    dTheta2_dt = omega2

    dOmega1_dt = (-g * (2 * m1 + m2) * np.sin(theta1) - m2 * g * np.sin(theta1 - 2 * theta2) - 2 * np.sin(theta1 - theta2) * m2 * 
                    (omega2**2 * l2 + omega1**2 * l1 * np.cos(theta1 - theta2))) / (l1 * (2 * m1 + m2 - m2 * np.cos(2 * theta1 - 2 * theta2)))
    dOmega2_dt = (2 * np.sin(theta1 - theta2) * (omega1**2 * l1 * (m1 + m2) + g * (m1 + m2) * np.cos(theta1) + 
                    omega2**2 * l2 * m2 * np.cos(theta1 - theta2))) / (l2 * (2 * m1 + m2 - m2 * np.cos(2 * theta1 - 2 * theta2)))
    
    return np.array([dTheta1_dt, dTheta2_dt, dOmega1_dt, dOmega2_dt])

theta1_0 = np.pi / 4 # rod 1 initial angle, radians
theta2_0 = np.pi / 6 # rod 2 initial angle, radians
omega1_0 = 0.0 # rod 1 initial angular velocity, rad/s
omega2_0 = 0.0 # rod 2 initial angular velocity, rad/s

y0 = np.array([theta1_0, theta2_0, omega1_0, omega2_0])

t_end = 10.0
dt = 0.01
n = int(t_end / dt)
t = np.linspace(0, t_end, n)

y = np.zeros((n, 4))
y[0] = y0

for i in range(1, n):
    y[i] = RK4(y[i - 1], t[i - 1], dt)

theta1, theta2, omega1, omega2 = y.T

plt.plot(t, theta1, label='Theta1')
plt.plot(t, theta2, label='Theta2')
plt.xlabel('t')
plt.ylabel('angle')
plt.legend()
plt.show()
