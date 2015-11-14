#
# Rabbit-fox population model (Volterra model) with
# logistic growth, covered rabbits and hunting saturation
#
# H  = c * F * (R - R0) * Rc / (R + Rc)
# dR = dt * [ a * (R * R / (R + R1)) * (1 - R/Rs) - H ]
# dF = dt * [ d * H - b * F ]
#
# a  - rabbit growth
# b  - fox death
# c  - hunting efficiency
# d  - eaten rabbit to fox born coefficient
# R0 - covered rabbits population
# R1 - lonely rabbits threshold
# Rc - hunting saturation
# Rs - saturated rabbits population
#

import numpy as np

a = 1./20
b = 1./3000
c = 1./100
d = 0.005
R0 = .1
R1 = 1.
Rc = 50.
Rs = 100.

# Equilibrium
R_ = b / (c * d)
R_e = Rc * (R0 + R_) / (Rc - R_)
F_e = (a / c) * R_e * R_e * (1 - R_e/Rs) * (R_e + Rc) / ((R_e + R1) * (R_e - R0) * Rc)

# Time quantum
dt = .1

N = 1000
I = range(N)
Left  = np.array(I[-1:]+I[:-1])
Right = np.array(I[1:]+I[:1])

def migration_rate(P):
	return .5 * P[Left] + .5 * P[Right] - P

def population_up(R, F):
	catched = c * F * (R - R0) * Rc / (R + Rc)
	born = d * catched
	vR = a * R * R * (1 - R/Rs) / (R + R1) - catched + .01 * migration_rate(R)
	vF = born - b * F + 1. * migration_rate(F)
	return R + vR * dt, F + vF * dt

import matplotlib.pyplot as pl
import matplotlib.animation as animation

R = np.array([2*R_e if N/2-5 <= i < N/2+5 else R_e for i in range(N)])
F = np.array([F_e for i in range(N)])
t = 0.
n = 0

fig, ax = pl.subplots()
lR, = ax.plot(I, R, label='rabbits')
lF, = ax.plot(I, F, label='foxes')

y_max = max(R)

def init():
	ax.set_xlim(0, N)
	ax.set_ylim(0, y_max)
	ax.legend()

def run(i):
	global R, F, t, n, y_max
	for i in range(100):
		R, F = population_up(R, F)
		t += dt
		n += 1
	Rmax = max(R)
	if y_max < Rmax:
		y_max = Rmax
		ax.set_ylim(0, y_max)

	lR.set_data(I, R)
	lF.set_data(I, F)
	ax.set_title(str(int(n*dt)))

ani = animation.FuncAnimation(fig, run, init_func=init, interval=0)
pl.show()



