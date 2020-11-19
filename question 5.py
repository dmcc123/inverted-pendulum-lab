import sympy as sym
import numpy as np
import control as ctrl
import matplotlib.pyplot as plt

m, ell, x3, x4, M, g, F, m = sym.symbols('m, ell, x3, x4, M, g, F, m')

# φ(F, x3, x4)
phi = 4*m*ell*x4**2*sym.sin(x3) + 4*F - 3*m*g*sym.sin(x3)*sym.cos(x3)
phi /= 4*(M+m) - 3*m*sym.cos(x3)**2
dphi_x3 = phi.diff(x3)
dphi_x4 = phi.diff(x4)
dphi_F = phi.diff(F)

# Equilibrium point
Feq = 0
x3eq = 0
x4eq = 0
dphi_F_eq = dphi_F.subs([(F, Feq), (x3, x3eq), (x4, x4eq)])
dphi_x3_eq = dphi_x3.subs([(F, Feq), (x3, x3eq), (x4, x4eq)])
dphi_x4_eq = dphi_x4.subs([(F, Feq), (x3, x3eq), (x4, x4eq)])

a = dphi_F_eq
b = -dphi_x3_eq
c = 3/(ell*(4*M + m))
d = 3*(M+m)*g/(ell*(4*M + m))

# GIVEN VALUES!
M_value = 0.3
m_value = 0.1
g_value = 9.81
ell_value = 0.35

def evaluate_at_given_parameters(z):
    """
    :param z:
    :return:
    """
    return float(z.subs([(M, M_value), (m, m_value), (ell, ell_value), (g, g_value)]))


a_value = evaluate_at_given_parameters(a)
b_value = evaluate_at_given_parameters(b)
c_value = evaluate_at_given_parameters(c)
d_value = evaluate_at_given_parameters(d)


transfer_function_G_θ = ctrl.TransferFunction([-c_value], [1, 0, -d_value])


def pid (kp, ki, kd):

    diff    = ctrl.TransferFunction ([1, 0], 1)
    intgr   = ctrl.TransferFunction (1, [1, 0])
    pid_t_f = kp + kd * diff + ki * intgr

    return pid_t_f

kp = 400        # for no oscillations try 400 5 20, to exhibit oscillations try 350 1 5
ki = 5
kd = 20

controller = -pid(kp, ki, kd)

n_points = 500
t_final = 1

closed_loop_t_f = ctrl.feedback(transfer_function_G_θ, controller)
t, y = ctrl.impulse_response(closed_loop_t_f, T = np.linspace(0, t_final, n_points))

y_deg = y * 180 / sym.pi

plt.plot(t, y_deg)
plt.suptitle('Trajectory of pendulum using PID controller', fontsize=14)
plt.xlabel('time(s)')
plt.ylabel('angle of pendulum,θ (degrees)')
plt.grid()
plt.show()
