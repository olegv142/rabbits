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

a = 1./20
b = 1./3000
c = 1./100
d = 0.005

A = 20.

Fa_e = a / c
Fc_e = A * a * b / c
R_e = b * (1 + A * b) / (c * d)

def population_drift(R, Fc, Fa):
	catched = R * Fa * c
	born = catched * d
	dR = R * a - catched
	dFc = born - Fc * b - Fc / A
	dFa = Fc / A - Fa * b
	return dR, dFc, dFa

class Population:
	"""
	Approximated solution with ageing
	"""
	def __init__(self, R, Fc, Fa, dt):
		self.R = R
		self.Fc = Fc
		self.Fa = Fa
		self.t = 0.
		self.dt = dt

	def up(self):
		dR, dFc, dFa = population_drift(self.R, self.Fc, self.Fa)
		self.R  += dR  * self.dt
		self.Fc += dFc * self.dt
		self.Fa += dFa * self.dt
		self.t  += self.dt

def solve(R, Fc, Fa, T, dt):
	p = Population(R, Fc, Fa, dt)
	t, R, F, Fc, Fa = [], [], [], [], []
	while p.t < T:
		p.up()
		t.append(p.t)
		R.append(p.R)
		F.append(p.Fc + p.Fa)
		Fc.append(p.Fc)
		Fa.append(p.Fa)
	return t, R, F, Fc, Fa

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure()
ax = fig.add_subplot(1, 2, 1)
t, R, F, Fc, Fa = solve(R_e + 0.01, Fc_e, Fa_e, 25000, 0.1)
ax.plot(t, R, label='rabbits')
ax.plot(t, F, label='foxes')
ax.legend()
ax = fig.add_subplot(1, 2, 2, projection='3d')
ax.plot(R, Fc, Fa)
ax.set_xlabel('Rabbits')
ax.set_ylabel('Child foxes')
ax.set_zlabel('Adult foxes')
plt.show()



