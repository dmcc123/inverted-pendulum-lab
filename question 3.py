import sympy as sym
import control as ctrl

m, ell, M, g = sym.symbols('m, ell, M, g', real=True, positive=True)
x1, x2, x3, x4, F = sym.symbols('x1, x2, x3, x4, F')

# φ(F, x3, x4)
phi = 4 * m * ell * x4**2 * sym.sin(x3) + 4 * F - 3 * m * g * sym.sin(x3) * sym.cos(x3)
phi /= 4 * (M + m) - 3 * m * sym.cos(x3)**2

#partial derviatives of phi wrt f, x3 and x4
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

a =   dphi_F_eq
b = - dphi_x3_eq
c = 3 / (ell * (4 * M + m))
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


s, t = sym.symbols('s, t')
c, d = sym.symbols('c, d', real=True, positive=True)

# G theta

G_θ = - c / (s**2 - d)

# impulse response

F_imp_θ = 1
x3_imp_θ = G_θ * F_imp_θ
x3_t_imp_θ = sym.inverse_laplace_transform(x3_imp_θ, s, t)

print ("")
print ("Impulse response for G_θ")
print ("")
sym.pprint(x3_t_imp_θ.simplify())
print ("")

# step response

F_step_θ = 1 / s
x3_step_θ = G_θ * F_step_θ
x3_t_step_θ = sym.inverse_laplace_transform(x3_step_θ, s, t)

print ("")
print ("Step response for G_θ")
print ("")
sym.pprint(x3_t_step_θ.simplify())
print ("")

# frequency response

# sin(ω * t) where ω is pi radians

F_freq_θ = sym.sin(sym.pi * t)
x3_freq_θ = G_θ * F_freq_θ
x3_t_freq_θ = sym.inverse_laplace_transform(x3_freq_θ, s, t)

print ("")
print ("Frequency response for G_θ")
print ("")
sym.pprint(x3_t_freq_θ.simplify())
print ("")

# G x

G_x = a + (b * c / (s**4 - d * s**2))

# impulse response

F_imp_x = 1
x3_imp_x = G_x * F_imp_x
x3_t_imp_x = sym.inverse_laplace_transform(x3_imp_x, s, t)

print ("")
print ("Impulse response for G_x")
print ("")
sym.pprint(x3_t_imp_x.simplify())
print ("")

# step response

F_step_x = 1 / s
x3_step_x = G_x * F_step_x
x3_t_step_x = sym.inverse_laplace_transform(x3_step_x, s, t)

print ("")
print ("Step response for G_x")
print ("")
sym.pprint(x3_t_step_x.simplify())
print ("")

# frequency response

F_freq_x = sym.sin(sym.pi * t)
x3_freq_x = G_x * F_freq_x
x3_t_freq_x = sym.inverse_laplace_transform(x3_freq_x, s, t)

print ("")
print ("Frequency response for G_x")
print ("")
sym.pprint(x3_t_freq_x.simplify())
print ("")