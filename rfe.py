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

# Параметры модели
a = 1./10    # один кролик рождается в среднем раз в 10 дней
b = 1./3000  # лиса живет в среднем 3000 дней
c = 1./2     # лисе хватает одного кролика на 2 дня
d = 1./200   # чтобы выросла новая лиса, нужно съесть 200 кроликов
R0 = .1      # этих кроликов невозможно поймать
R1 = 1.      # при меньшем количестве крольчихе трудно найти партнера
Rh = 50.     # при таком количестве кроликов происходит насыщение хищников
Rs = 100.    # столько кроликов выедают всю траву, и популяция перестает расти

# Равновесная численность популяции
R_e = R0 + Rh / (d*c/b - 1)
F_e = (a/c) * (R_e/2) * (1 - R_e/Rs) * (1 + Rh/(R_e - R0)) / (1 + R1/R_e)

# Скорость рождения кроликов как функция численности их популяции
def B(R):
	return a*((R*R/2.)/(R+R1))*(1.-R/Rs)

# Скорость отлова кроликов лисами как функция численности их популяций
def H(R, F):
	return c*F*(R-R0)/(R-R0+Rh)

def up(R, F):
	vR = B(R) - H(R, F)
	vF = d*H(R, F) - b*F
	return R+vR*dt, F+vF*dt

dt = .1

def evolution(T, R=R_e, F=F_e):
	t  = 0
	aT = []
	aR = []
	aF = []
	while t < T:
		aT.append(t)
		aR.append(R)
		aF.append(F)
		for _ in range(100):
			R, F = up(R, F)
			t += dt
	return aT, aR, aF

def plot_evolution(T):
	import matplotlib
	import matplotlib.pyplot as pl
	matplotlib.rc('font', family='Arial')
	aT, aR, aF = evolution(T, R_e+.001)
	pl.plot(aT, aR, label=u'кролики')
	pl.plot(aT, aF, label=u'лисы')
	pl.legend()
	pl.xlabel(u'дни')
	pl.ylabel(u'численность популяции')
	#pl.title(u'эволюция после малого отклонения от равновесия')
	pl.show()

def plot_equlibrium():
	import matplotlib
	import numpy as np
	import matplotlib.pyplot as pl
	matplotlib.rc('font', family='Arial')
	R = np.linspace(0, Rs, 1001)
	pl.plot(R, B(R), label=u'родившиеся')
	pl.plot(R, H(R, F_e), label=u'пойманные')
	pl.plot([R_e], [B(R_e)], 'o')
	pl.text(R_e + 2, B(R_e), u'равновесие (неустойчивое)')
	pl.xlim(0, Rs)
	pl.ylim(0, H(Rs, F_e))
	pl.legend()
	pl.xlabel(u'популяция кроликов')
	pl.ylabel(u'скорость рождения / отлова кроликов')
	#pl.title(u'неустойчивое равновесие')
	pl.show()

plot_equlibrium()
plot_evolution(50000)


