# -*- coding: utf-8 -*-
#
# Лисы - кролики
# R - популяция кроликов
# F - популяция лис
# Скорость рождения кроликов:
#  B(R) = a * ((R/2) / (1 + R1/R)) * (1 - R/Rs)
#    a  - сколько кроликов рождается у крольчихи в единицу времени при благоприятных условиях
#    R1 - порог одиночества
#    Rs - порог насыщения популяции (емкость среды)
# Скорость отлова кроликов лисами:
#  H(R,F) = c * F / (1 + Rh/(R - R0))
#    c  - сколько кроликов ловит одна лиса в единицу времени при благоприятных условиях
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
import matplotlib
import matplotlib.pyplot as pl
import matplotlib.animation as animation

a = 1./10    # один кролик рождается в среднем раз в 10 дней
b = 1./3000  # лиса живет в среднем 3000 дней
c = 1./2     # лисе хватает одного кролика на 2 дня
d = 1./200   # чтобы выросла новая лиса, нужно съесть 200 кроликов
R0 = .1      # этих кроликов невозможно поймать
R1 = 1.      # при меньшем количестве крольчихе трудно найти партнера
Rh = 50.     # при таком количестве кроликов происходит насыщение хищников
Rs = 100.    # столько кроликов выедают всю траву, и популяция перестает расти

# Equilibrium
R_e = R0 + Rh / (d*c/b - 1)
F_e = (a/c) * (R_e/2) * (1 - R_e/Rs) * (1 + Rh/(R_e - R0)) / (1 + R1/R_e)

# Time quantum
dt = .1

# Simulation array parameters 
N = 990
I = range(N)
Left  = np.array(I[-1:]+I[:-1])
Right = np.array(I[1:]+I[:1])

# Migration speed parameters
Vr = .02
Vf = 1.

def M(P):
	return .5 * P[Left] + .5 * P[Right] - P

def population_up(R, F):
	B = a * (R * R / (R + R1)) * (1 - R/Rs) / 2
	H = c * F * (R - R0) / (R - R0 + Rh)
	vR = B - H + Vr * M(R)
	vF = d * H - b * F + Vf * M(F)
	return R + vR * dt, F + vF * dt

def single_peak(w):
	return np.array([2*R_e if N/2-w/2 <= i < N/2+w/2 else R_e for i in range(N)])

def island((population, w)):
	return np.array([population if N/2-w/2 <= i < N/2+w/2 else 0 for i in range(N)])

def peaks((per, l, r)):
	return np.array([2*R_e if l[(i/per) % len(l)] <= (i%per) < r[(i/per) % len(r)] else R_e for i in range(N)])

def random(seed):
	np.random.seed(seed)
	return R_e * (1 + .1 * np.random.randn(N))

def foxes(population):
	return np.array([population for i in range(N)])

presets = {
	# name: (Vr, Vf, rabbits_init, rabbits_init_arg, foxes_init, foxes_init_arg)
	'soliton'      : (.1,  1.,  random, 0,       foxes, F_e),
	'waves'        : (.1,  1.,  single_peak, 10, foxes, F_e),
	'struct-flash' : (.02, 1.,  random, 1,       foxes, F_e),
	'struct'       : (.025, 2., single_peak, 10, foxes, F_e),
	'dyn-struct'   : (.025, 1., peaks, (55, (0, 10), (10, 20)), foxes, F_e),
	'no-foxes'     : (.1, 0.,   island, (.1, 10), foxes, 0.),
}

def run_simulation():
	class context:
		R, F = None, None
		t = 0.
		n = 0

	rinit, rarg = single_peak, 10
	finit, farg = foxes, F_e
	title = ''
	args = {
		'step':10.,
		'video':None,
		'fps':30,
		'frames':10000,
	}

	if len(sys.argv) > 1:
		global Vr, Vf
		preset = sys.argv[1]
		print 'Using preset', preset
		Vr, Vf, rinit, rarg, finit, farg = presets[preset]
		title = '[' + preset + ']'

	if len(sys.argv) > 2:
		for arg in sys.argv[2:]:
			exec(arg, args)

	print 'Vr=%f, Vf=%f, %s(%s), %s(%s), step=%f' % (Vr, Vf, rinit.__name__, rarg, finit.__name__, farg, args['step'])

	context.R = rinit(rarg)
	context.F = finit(farg)

	matplotlib.rc('font', family='Arial')
	fig, ax = pl.subplots()
	lR, = ax.plot(I, context.R, label=u'кролики')
	lF, = ax.plot(I, context.F, label=u'лисы')
	ax.set_xlim(0, N)
	ax.set_ylim(0, Rs * 1.02)
	ax.legend()

	def run(i):
		for i in range(int(args['step']/dt)):
			context.R, context.F = population_up(context.R, context.F)
			context.t += dt
			context.n += 1

		lR.set_data(I, context.R)
		lF.set_data(I, context.F)
		ax.set_title(u'%s  %d дней' % (title, context.n*dt))

	if args['video']:
		print "save video to '%s' at %d fps, %d frames" % (args['video'], args['fps'], args['frames'])
		ani = animation.FuncAnimation(fig, run, interval=0, frames=args['frames'])
		ani.save(args['video'], fps=args['fps'], extra_args=['-vcodec', 'libx264'])
	else:
		ani = animation.FuncAnimation(fig, run, interval=0)
		pl.show()

if __name__ == '__main__':
	run_simulation()

