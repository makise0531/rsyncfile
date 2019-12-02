import sys
sys.path.append('/Users/user/data/nopt/b104_2018A0120/work/function/')

import flambaum_def as fl

mode = "phi"
angle = 45
En = 4.53

Ene_1wave = []
g2gamma_n1 = []
gamma_g1 = []
J1 = []

Ene_2wave = []
g2gamma_n2 = []
gamma_g2 = []
J2 = []

jx = []

ns = 0
np = 0
F = 0

# flambaum_def.calcFlambaum(angle, mode)
#fl.make_a1(En, angle)
fl.a1_cal(En, Ene_1wave, g2gamma_n1, gamma_g1, J1, Ene_2wave, g2gamma_n2, gamma_g2, J2, jx, angle, 3, 1, F)
