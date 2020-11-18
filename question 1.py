import sympy as sym

m, ell, x3, x4, M, g, F, m = sym.symbols('m, ell, x3, x4, M, g, F, m')

phi  = 4*m*ell*x4**2*sym.sin(x3) + 4*F - 3*m*g*sym.sin(x3)*sym.cos(x3)
phi /= 4*(M+m) - 3*m*sym.cos(x3)**2

psi  = -3*(m*ell*x4**2*sym.sin(x3)*sym.cos(x3) + F*sym.cos(x3) - (M + m)*g*sym.sin(x3))
psi /= ell*(4*(M + m) - 3*m*sym.cos(x3)**2)

dphi_x3 = phi.diff(x3)
dphi_x4 = phi.diff(x4)
dphi_F  = phi.diff(F)

dpsi_x3 = psi.diff(x3)
dpsi_x4 = psi.diff(x4)
dpsi_F  = psi.diff(F)

Feq  = 0
x3eq = 0
x4eq = 0

dphi_F_eq  = dphi_F.subs ([(F,Feq), (x3, x3eq), (x4, x4eq)])
dphi_x3_eq = dphi_x3.subs([(F,Feq), (x3, x3eq), (x4, x4eq)])
dphi_x4_eq = dphi_x4.subs([(F,Feq), (x3, x3eq), (x4, x4eq)])

dpsi_F_eq  = dpsi_F.subs ([(F,Feq), (x3, x3eq), (x4, x4eq)])
dpsi_x3_eq = dpsi_x3.subs([(F,Feq), (x3, x3eq), (x4, x4eq)])
dpsi_x4_eq = dpsi_x4.subs([(F,Feq), (x3, x3eq), (x4, x4eq)])

sym.pprint(dphi_F_eq )
print("")
sym.pprint(dphi_x3_eq)
print("")
sym.pprint(dphi_x4_eq)
print("")
sym.pprint(dpsi_F_eq )
print("")
sym.pprint(dpsi_x3_eq)
print("")
sym.pprint(dpsi_x4_eq)