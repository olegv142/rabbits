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
dt = .5

# Simulation array parameters 
N = 25  # height, first index
M = 201 # width, second index

j = range(M)

I = [[i]*M for i in range(N)]
Dn = I[-1:]+I[:-1]
Up = I[1:]+I[:1]

J = [j for i in range(N)]
Le = [j[-1:]+j[:-1] for i in range(N)]
Rt = [j[1:]+j[:1] for i in range(N)]

Dn, Up, Le, Rt = np.array(Dn), np.array(Up), np.array(Le), np.array(Rt)

# Migration speed parameters
Vr = 0.01
Vf = 1.

def Mig(P):
	return .25 * (P[I,Le] + P[I,Rt] + P[Up,J] + P[Dn,J]) - P

def population_up(R, F):
	B = a * (R * R / (R + R1)) * (1 - R/Rs) / 2
	H = c * F * (R - R0) / (R - R0 + Rh)
	vR = B - H + Vr * Mig(R)
	vF = d * H - b * F + Vf * Mig(F)
	return R + vR * dt, F + vF * dt

dR = 1.
R = R_e * np.ones((N,M))
R[N/2,M/2] = R_e + dR
F = F_e * np.ones((N,M))
t = 0
n = 0

def run_simulation():
	args = {
		'step':50.,
		'video':None,
		'fps':30,
		'frames':10000,
	}

	if len(sys.argv) > 1:
		for arg in sys.argv[1:]:
			exec(arg, args)

	matplotlib.rc('font', family='Arial')
	fig, ax = pl.subplots()
	im = ax.imshow(R, vmin=0, vmax=Rs)
	pl.colorbar(im, orientation='horizontal')

	def run(i):
		global R, F, t, n
		for i in range(int(args['step']/dt)):
			R, F = population_up(R, F)
			t += dt
			n += 1
		im.set_data(R)
		ax.set_title(u'%d дней' % (n*dt))

	if args['video']:
		print "save video to '%s' at %d fps, %d frames" % (args['video'], args['fps'], args['frames'])
		ani = animation.FuncAnimation(fig, run, interval=0, frames=args['frames'])
		ani.save(args['video'], fps=args['fps'], extra_args=['-vcodec', 'libx264'])
	else:
		ani = animation.FuncAnimation(fig, run, interval=0)
		pl.show()

if __name__ == '__main__':
	run_simulation()

