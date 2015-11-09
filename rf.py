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

Fe = a / c 
Re = b / (c * d)

class Population:
	def __init__(self, R, F, dt):
		self.R = R
		self.F = F
		self.t = 0.
		self.dt = dt
		
	def up(self):
		H  = self.R * self.F * c
		dR = self.dt * ( self.R * a - H )
		dF = self.dt * ( H * d - self.F * b )
		self.R += dR
		self.F += dF
		self.t += self.dt

import matplotlib.pylab as pl

p = Population(Re+1, Fe, 0.1)
t, R, F = [], [], []
while p.t < 10000:
	p.up()
	t.append(p.t)
	R.append(p.R)
	F.append(p.F)

pl.plot(t, R, label='rabbits')
pl.plot(t, F, label='foxes')
#pl.plot(R, F)
pl.legend()
pl.show()



