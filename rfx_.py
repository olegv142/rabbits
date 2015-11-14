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

N = 1000

import numpy as np
import matplotlib.pyplot as pl
R = np.linspace(0, Rs, N)
C = c * F_e * (R - R0) * Rc / (R + Rc)
B = a * R * R * (1 - R/Rs) / (R + R1)
pl.plot(R, C, label='catched')
pl.plot(R, B, label='burn')
pl.xlim((0, Rs))
pl.ylim((0, C[-1]))
pl.legend()
pl.show()



