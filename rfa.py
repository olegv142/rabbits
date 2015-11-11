#
# Rabbit-fox population model with ageing
#
# only adult foxes are hunting
#
# dR  = dt * [ R * a - R * Fa * c ]
# dFc = dt * [ R * Fa * c * d - Fc * b - Fc / A ]
# dFa = dt * [ Fc / A - Fa * b ]
#
# a - rabbit growth
# b - fox death
# c - hunting efficiency
# d - eaten rabbit to fox born coefficient
# A - adult age
#
# Considering environmental saturation:
# dR  = dt * [ R * a * (1 - R/Rs) - R * Fa * c ]
#
# Rs - saturated population

a = 1./20
b = 1./3000
c = 1./100
d = 0.005

A = 100.
Rs = 1000.

# Equilibrium
R_e = b * (1 + A * b) / (c * d)
Fa_e = (a / c) * (1 - R_e / Rs)
Fc_e = A * b * Fa_e

# Time quantum
dt = .1

def population_up(R, Fc, Fa):
	catched = R * Fa * c
	born = catched * d
	vR = R * a * (1 - R / Rs) - catched
	vFc = born - Fc * b - Fc / A
	vFa = Fc / A - Fa * b
	return R + vR * dt, Fc + vFc * dt, Fa + vFa * dt

import matplotlib.pyplot as pl

R = R_e + 0.01
Fc = Fc_e
Fa = Fa_e
t = 0.
R_  = []
Fc_ = []
Fa_ = []
t_  = []
while t < 100000:
	R, Fc, Fa = population_up(R, Fc, Fa)
	R_.append(R)
	Fc_.append(Fc)
	Fa_.append(Fa)
	t_.append(t)
	t += dt

pl.plot(t_, R_, label='rabbits')
pl.plot(t_, Fc_, label='pups')
pl.plot(t_, Fa_, label='foxes')
pl.legend()
pl.show()



