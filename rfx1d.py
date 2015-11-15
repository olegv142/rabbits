# -*- coding: utf-8 -*-
#
# Лисы - кролики
# R - популяция кроликов
# F - популяция лис
# Скорость рождения кроликов:
#  B(R) = a * (R / (1 + R1/R)) * (1 - R/Rs)
#    a  - коэффициент пропорциональности
#    R1 - порог одиночества
#    Rs - порог насыщения популяции (емкость среды)
# Скорость отлова кроликов лисами:
#  H(R,F) = c * F / (1 + Rh / (R - R0))
#    c  - коэффициент пропорциональности 
#    Rh - порог насыщения хищников
#    R0 - порог неуловимости
# Скорости изменения популяций:
#  R' = B - H
#  F' = d * H - b * F
#    d - коэффициент конверсии съеденных кроликов в новых лис
#    b - смертность лис
#

import sys
import numpy as np
import matplotlib.pyplot as pl
import matplotlib.animation as animation

a = 1./20    # один кролик рождается в среднем раз в 20 дней
b = 1./3000  # лиса живет в среднем 3000 дней
c = 1./2     # лисе хватает одного кролика на 2 дня
d = 1./200   # чтобы выросла новая лиса, нужно съесть 200 кроликов
R0 = .1      # этих кроликов невозможно поймать
R1 = 1.      # при меньшем количестве крольчихе трудно найти партнера
Rh = 50.     # при таком количестве кроликов происходит насыщение хищников
Rs = 100.    # столько кроликов выедают всю траву, и популяция перестает расти

# Equilibrium
R_e = R0 + Rh / (d*c/b - 1)
F_e = (a/c) * R_e * (1 - R_e/Rs) * (1 + Rh/(R_e - R0)) / (1 + R1/R_e)

# Time quantum
dt = .1

# Simulation array parameters 
N = 1000
I = range(N)
Left  = np.array(I[-1:]+I[:-1])
Right = np.array(I[1:]+I[:1])

# Migration speed parameters
Vr = .02
Vf = 1.

def M(P):
	return .5 * P[Left] + .5 * P[Right] - P

def population_up(R, F):
	B = a * (R / (1 + R1/R)) * (1 - R/Rs)
	H = c * F / (1 + Rh / (R - R0))
	vR = B - H + Vr * M(R)
	vF = d * H - b * F + Vf * M(F)
	return R + vR * dt, F + vF * dt

presets = {
	# name: (Vr, Vf, random_seed)
	'soliton'      : (.1,  1., 0),
	'waves'        : (.05, .05, None),
	'struct-flash' : (.02, 1., 1),
	'struct'       : (.01, 1., None)
}

def run_simulation():
	class context:
		R, F = None, None
		y_max = None
		t = 0.
		n = 0

	seed = None
	if len(sys.argv) > 1:
		global Vr, Vf
		preset = sys.argv[1]
		print 'Using preset', preset
		Vr, Vf, seed = presets[preset]

	context.F = np.array([F_e for i in range(N)])
	if seed is None:
		print 'Vr=%f, Vf=%f' % (Vr, Vf)
		context.R = np.array([2*R_e if N/2-5 <= i < N/2+5 else R_e for i in range(N)])
	else:
		print 'Vr=%f, Vf=%f, seed=%s' % (Vr, Vf, seed)
		np.random.seed(seed)
		context.R = R_e * (1 + .1 * np.random.randn(N))

	fig, ax = pl.subplots()
	lR, = ax.plot(I, context.R, label='rabbits')
	lF, = ax.plot(I, context.F, label='foxes')

	context.y_max = max(context.R)

	def init():
		ax.set_xlim(0, N)
		ax.set_ylim(0, context.y_max)
		ax.legend()

	def run(i):
		for i in range(100):
			context.R, context.F = population_up(context.R, context.F)
			context.t += dt
			context.n += 1

		Rmax = max(context.R)
		if context.y_max < Rmax:
			context.y_max = Rmax
			ax.set_ylim(0, context.y_max)

		lR.set_data(I, context.R)
		lF.set_data(I, context.F)
		ax.set_title(str(int(context.n*dt)))

	ani = animation.FuncAnimation(fig, run, init_func=init, interval=0)
	pl.show()

if __name__ == '__main__':
	run_simulation()

