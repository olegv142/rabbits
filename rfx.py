#
# Rabbit-fox population model (Volterra model) with
# logistic growth, covered rabbits and hunting saturation
#
# H  = c * F * (R - R0) * Rc / (R + Rc)
# dR = dt * [ a * R * (1 - R/Rs) - H ]
# dF = dt * [ d * H - b * F ]
#
# a  - rabbit growth
# b  - fox death
# c  - hunting efficiency
# d  - eaten rabbit to fox born coefficient
# R0 - covered rabbits population
# Rc - hunting saturation
# Rs - saturated rabbits population
#

a = 1./20
b = 1./3000
c = 1./100
d = 0.005
R0 = 0
Rc = 200.
Rs = 500.

# Equilibrium
R_ = b / (c * d)
R_e = Rc * (R0 + R_) / (Rc - R_)
F_e = (a / c) * R_e * (1 - R_e/Rs) * (R_e + Rc) / ((R_e - R0) * Rc)

# Time quantum
dt = .1

def population_up(R, F):
	catched = c * F * (R - R0) * Rc / (R + Rc)
	born = d * catched
	vR = a * R * (1 - R/Rs) - catched
	vF = born - b * F
	return R + vR * dt, F + vF * dt

import matplotlib.pyplot as pl

R = R_e + .1
F = F_e
t = 0.
R_ = []
F_ = []
t_ = []
while t < 50000:
	R, F = population_up(R, F)
	R_.append(R)
	F_.append(F)
	t_.append(t)
	t += dt

pl.plot(t_, R_, label='rabbits')
pl.plot(t_, F_, label='foxes')
pl.legend()
pl.show()



