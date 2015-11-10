#
# Rabbit-fox population model (Volterra model)
#
# dR = dt * [ R * a - R * F * c ]
# dF = dt * [ R * F * c * d - F * b ]
#
# a - rabbit growth
# b - fox death
# c - hunting efficiency
# d - eaten rabbit to fox born coefficient
#

a = 1./20
b = 1./3000
c = 1./100
d = 0.005

# Equilibrium
F_e = a / c 
R_e = b / (c * d)

# Time quantum
dt = .1

def population_up(R, F):
	catched = R * F * c
	born = catched * d
	vR = R * a - catched
	vF = born - F * b
	return R + vR * dt, F + vF * dt

import matplotlib.pyplot as pl

R = R_e + 1
F = F_e
t = 0.
R_ = []
F_ = []
t_ = []
while t < 10000:
	R, F = population_up(R, F)
	R_.append(R)
	F_.append(F)
	t_.append(t)
	t += dt

pl.plot(t_, R_, label='rabbits')
pl.plot(t_, F_, label='foxes')
pl.legend()
pl.show()



