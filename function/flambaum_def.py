# -*- coding: utf-8 -*-
"""
This is a calculate module.
"""
import sys, math
import ROOT as rt
import numpy as npy
import scipy as sp
from sympy.physics.wigner import wigner_6j
from ROOT import TCanvas, TGraph, TLegend
from array import array
import code

sys.path.append('/Users/user/data/nopt/b104_2018A0120/work/function/')


def setLegend(leg, entries):
    if isinstance(entries, (list, tuple)):
        for ent in entries:
            leg.AddEntry(ent, ent.GetName(), "l")


range_nE = [0., 10.]
nstep = 1000
step = (range_nE[1] - range_nE[0]) / nstep
Npx = 10000
f1 = rt.TF1()


""" Definition of constant parameters """

NeutronMass = 939.5654133e6
protonMass = 938.27813e6
carbonMass = 6. * protonMass + 6. * NeutronMass
spl = 299792458.  # [m / s]
hbar = 6.582119514e-16  # [eV.s]
pi = npy.pi
# barn = 6.50977023e5
barn = 1.0e28  # m^2 → barn

spflag = 0

""" resonance difinition """
Ene_swave = []
gamma_gs = []
g2gamma_ns = []
Js = []

Ene_pwave = []
gamma_gp = []
g2gamma_np = []
Jp = []

ns = 0
np = 0
Esp = 0
""" 111Cd [eV]"""
targetMass = 48. * protonMass + 63. * NeutronMass
Anum = (targetMass * 12. / carbonMass) * 931.494028
Rrad = 1.35 * Anum ** 1/3

k = 1.
j1 = 0.5
j2 = [0.5, 1.5]  # 1/2 or 3/2
I = 1/2
F = 0
x1 = 0
x2 = 0

"""angle type"""

phi = 0.
theta = 0.


def parameterGet():
    """Give information."""
    for l in open('/Users/user/data/nopt/b104_2018A0120/work/study_convolve/flambaum/resopara_swave_111Cd.txt').readlines():
        data_swave = l[:-1].split(' ')
        Ene_swave.append(float(data_swave[0]))
        g2gamma_ns.append(float(data_swave[1]))
        gamma_gs.append(float(data_swave[2]))
        Js.append(float(data_swave[3]))
        global ns
        ns = len(Ene_swave)

    for l in open('/Users/user/data/nopt/b104_2018A0120/work/study_convolve/flambaum/resopara_pwave_111Cd.txt').readlines():
        data_pwave = l[:-1].split(' ')
        Ene_pwave.append(float(data_pwave[0]))
        g2gamma_np.append(float(data_pwave[1]))
        gamma_gp.append(float(data_pwave[2]))
        Jp.append(float(data_pwave[3]))
        global np
        np = len(Ene_pwave)
        print(ns)
        print(Ene_swave)
        print(g2gamma_ns)
        print(gamma_gs)
        print(Js)
        print(np)
        print(Ene_pwave)
        print(g2gamma_np)
        print(gamma_gp)
        print(Jp)
    return



def wig1(k, j1, j2, I, Jp, Js):
    """Return j2 parameter."""
    wigcal1 = wigner_6j(k, j1, j2, I, Jp, Js)
    return wigcal1


def wig2(k, x1, x2, F, Js, Jp):
    """Return x1 x2 parameter."""
    wigcal2 = wigner_6j(k, 1., 1., F, Js, Jp)
    return wigcal2


def P_comp(Js, Jp, j1, j2, I, F):
    """Return j2 parameter."""
    Pcomp1 = (- 1) ** (Js + Jp + j2 + I + F) * 3/2
    Pcomp2 = math.sqrt((2. * Js + 1.)*(2. * Jp + 1.) * (2. * j1 + 1.) * (2. * j2 + 1.))
    Pcompcal = Pcomp1 * Pcomp2
    return Pcompcal


def P_cal(k, j1, j2, I, Jp, Js, F, x1, x2):
    """Return P parameter."""
    p1 = wigner_6j(k, j1, j2, I, Jp, Js)
    p2 = wigner_6j(k, 1., 1., F, Js, Jp)
    p3 = P_comp(Js, Jp, j1, j2, I, F)
    p_cal = p1 * p2 * p3
    return p_cal


def P12_define():
    """Return."""
    p1_12 = wigner_6j(k, j1, 1/2, I, Jp[0], Js[0])
    p2_12 = wigner_6j(k, 1., 1., F, Js[0], Jp[0])
    p3_12 = P_comp(Js[0], Jp[0], j1, 1/2, I, F)
    Pvalue_12 = p1_12 * p2_12 * p3_12
    print(p1_12)
    print(p2_12)
    print(p3_12)
    return Pvalue_12


def P32_define():
    """Return Js Jp j1 j2 k I F."""
    p1_32 = wig1(k * 2, j1 * 2, 3/2 * 2, I * 2, Jp[0] * 2, Js[0] * 2)
    p2_32 = wig2(k, 1., 1., F, Js[0], Jp[0])
    p3_32 = P_comp(Js[0], Jp[0], j1, 3/2, I, F)
    Pvalue_32 = p1_32 * p2_32 * p3_32
    print(p1_32)
    print(p2_32)
    print(p3_32)
    return Pvalue_32


def P_resu():
    """Pcal_result."""
    P_12 = P12_define()
    print(P_12)
    P_32 = P32_define()
    print(P_32)


""" V1 & V2 caluculation. """


def V1(En, Es, g2gamma_ns, gamma_gs, Js):
    """Return V1 equation."""
    #  ks = np.sqrt(2 * Es * NeutronMass) / (hbar * spl)
    kn = math.sqrt(2 * En * NeutronMass) / (hbar * spl)
    gs = (2. * Js + 1.) / (2. * (2. * I + 1.))
    gamma_s = (g2gamma_ns / (2. * gs)) + gamma_gs
    V1_comp1 = - 1 / (2. * kn)
    V1_comp2 = math.sqrt((g2gamma_ns / 2.) * gamma_gs)
    V1_comp3 = (En - Es) + (gamma_s * 1j) / 2.
    V1_cal = V1_comp1 * V1_comp2 / V1_comp3 * math.sqrt(barn)
    return V1_cal


def V2(En, Ep, g2gamma_np, gamma_gp, Jp, j, phi):
    """Return V2 equation."""
    kn = math.sqrt(2. * En * NeutronMass) / (hbar * spl)
    gp = (2. * Jp + 1.) / (2. * (2. * I + 1.))
    gamma_p = (g2gamma_np / (2 * gp)) + gamma_gp
    V2_comp1 = - 1. / (2. * kn)
    V2_comp2 = math.sqrt((g2gamma_np / 2.) * gamma_gp)
    V2_comp3 = (En - Ep) + gamma_p * 1j / 2.
    V2_comp4 = V2_comp1 * V2_comp2 / V2_comp3 * math.sqrt(barn)
    if j == 1./2.:
        V2_cal = V2_comp4 * npy.cos(phi * pi / 180.)
    elif j == 3./2.:
        V2_cal = V2_comp4 * npy.sin(phi * pi / 180.)
    return V2_cal


def a0_cal(En, Esp, g2gamma_n, gamma_g, J, j, phi, spflag):
    """Calculate a0 term."""
    if spflag == 0:
        V1_cal = V1(En, Esp, g2gamma_n, gamma_g, J)
        V_abs = V1_cal * V1_cal.conjugate()
    else:
        V2_cal = V2(En, Esp, g2gamma_n, gamma_g, J, j, phi)
        V_abs = V2_cal * V2_cal.conjugate()
    return V_abs


def a1_cal(En, Ene_swave, g2gamma_ns, gamma_gs, Js, Ene_pwave, g2gamma_np, gamma_gp, Jp, j, phi, ns, np, F):
    """Calculate a1 term."""
    print("We calculate V1")
    print("ns = ", ns)
    print(Ene_swave)
    for i in range(ns):
        for n in range(np):
            V_cal = 2. * V1(En, Ene_swave[i], g2gamma_ns[i], gamma_gs[i], Js[i]) * V2.conjugate(En, Ene_pwave[n], g2gamma_np[n], gamma_gp[n], Jp[n], j, phi) * P_cal(Js[i], Jp[n], 1/2, j, 1., I, F)
            print(V_cal)
            return V_cal


def make_a0(En, phi):
    """Calculate a0."""
    parameterGet()
    a0 = 0.
    j2 = [1/2, 3/2]
    for i in range(len(Ene_swave)):
        a0 += a0_cal(En, Ene_swave[i], g2gamma_ns[i], gamma_gs[i], Js[i], 0, 0, 0)
    for i in range(len(Ene_pwave)):
        for n in range(len(j2)):
            a0 += a0_cal(En, Ene_pwave[i], g2gamma_np[i], gamma_gp[i], Jp[i], j2[n], phi, 1)
    return a0


def make_a1(En, phi):
    """Calculate a0."""
    parameterGet()
    a1 = 0.
    j2 = [1/2, 3/2]
    print(len(j2))
    print("| + Processing makea1...\n")
    for i in range(len(j2)):
        a1 += a1_cal(En, Ene_swave, g2gamma_ns, gamma_gp, Js, Ene_pwave, g2gamma_np, gamma_gp, Jp, j2[i], phi, ns, np, F)
    print(a1)
    return a1


def acal_all(En, phi, theta):
    """Calculate a term."""
    a_all = 0.
    a0 = 0.
    a1 = 0.
    print("| making a0\n")
    a0 = make_a0(En, phi)
    a1 = make_a1(En, phi)   # a3 = make_a3(En, phi);
    a_all = (1/2) * (a0 + a1 * (npy.cos(theta*pi/180)))  # + a3*(pow(cos(theta*pi/180.), 2.) - 1./3.));
    return a_all


def makeFlambaum(En, phi, theta):
    return acal_all(En, phi, theta)


def makeTGraph(nstep, x, y):
    """Make TGraph."""
    x, y = array('d'), array('d')
    for i in range(nstep):
        x.append(float(range_nE[0] + i * step))
        y.append(float(f1.Eval(x)))
        gr1 = TGraph()
        gr1.SetPoint(i, x, y)
    gr1.SetName('g_{}'.format(f1.GetName()))
    gr1.SetTitle('Graph of {}'.format(f1.GetName()))
    gr1.SetMarkerStyle(6)
    gr1.Draw('ACP')
    return makeTGraph(nstep, x, y)


def calcFlambaum(angle, mode):
    """Calculate flambaum."""
    print("| + Processing flambaum...\n")
    parameterGet()
    flambaum1 = rt.TF1()
    flambaum2 = rt.TF1()
    flambaum3 = rt.TF1()
    if mode == phi:
        flambaum1 = rt.TF1("flambaum", "makeFlambaum(x, 0, {})".format(angle), range_nE[0], range_nE[1])
        flambaum2 = rt.TF1("flambaum", 'makeFlambaum(x, 90, {})'.format(angle), range_nE[0], range_nE[1])
        flambaum3 = rt.TF1("flambaum", 'makeFlambaum(x, 180, {})'.format(angle), range_nE[0], range_nE[1])
    elif mode == theta:
        flambaum1 = rt.TF1("flambaum", 'makeFlambaum(x, {}, 36)'.format(angle), range_nE[0], range_nE[1])
        flambaum1 = rt.TF1("flambaum", 'makeFlambaum(x, {}, 90)'.format(angle), range_nE[0], range_nE[1])
        flambaum1 = rt.TF1("flambaum", 'makeFlambaum(x, {}, 144)'.format(angle), range_nE[0], range_nE[1])

    flambaum1.SetNpx(Npx)
    flambaum2.SetNpx(Npx)
    flambaum3.SetNpx(Npx)

    print("| + Set TCanvas\n")
    c1 = TCanvas('can', 'flambaum')
    c1.SetLogy(0)
    c1.SetGrid()
    rt.gStyle.SetOptTitle(0)

    flambaum1.Draw("l")
    flambaum1.SetLineColor(rt.kBlack)
    flambaum1.SetLineWidth(4)
    flambaum1.SetLineStyle(1)

    flambaum2.Draw("lsame")
    flambaum2.SetLineColor(rt.kRed)
    flambaum2.SetLineWidth(4)
    flambaum2.SetLineStyle(2)

    flambaum3.Draw("lsame")
    flambaum3.SetLineColor(rt.kBlue)
    flambaum3.SetLineWidth(4)
    flambaum3.SetLineStyle(3)

    rt.gStyle.SetLabelSize(0.05, "x")
    rt.gStyle.SetLabelSize(0.05, "y")
    rt.gStyle.SetTitleYOffset(1)
    rt.gStyle.SetTitleSize(0.05, "x")
    rt.gStyle.SetTitleSize(0.05, "y")

    print("| Calling double SetRangeUser \n")
    #flambaum1.SetRangeUser(4.03, 5.03)
    #flambaum1.SetRangeUser(0., 2.0)
    print("| Calling double SetRangeUser → True!!!! \n")

    flambaum1.SetTitle("neutron Energy [eV]")
    flambaum1.SetTitle("cross section [barn]")
    leg = TLegend(0.65, 0.70, 0.90, 0.90)
    if mode == phi:
        setLegend(flambaum1, '#phi = 0 , #theta = {}'.format(angle), "l")
        setLegend(flambaum2, '#phi = 90 , #theta = {}'.format(angle), "l")
        setLegend(flambaum3, '#phi = 180 , #theta = {}'.format(angle), "l")
    elif mode == theta:
        setLegend(flambaum1, '#phi = {} , #theta = 36'.format(angle), "l")
        setLegend(flambaum2, '#phi = {} , #theta = 90'.format(angle), "l")
        setLegend(flambaum3, '#phi = {} , #theta = 144'.format(angle), "l")
    leg.SetFillColor(0)
    leg.Draw()
    print("| + Close Process\n")
    print('makeFlambaum(x, 0, {})'.format(angle))
    print(range_nE)
    return code.InteractiveConsole(globals()).interact()


def eachTermDraw(phi):
    """Calculate each term."""
    print("| + Processing flambaum...\n")
    parameterGet()
    a0 = rt.TF1("a0", 'make_a0(x, {})'.format(phi), range_nE[0], range_nE[1])
    a0.SetNpx(Npx)
    a0.SetLineColor(rt.kRed)
    a1x = rt.TF1("a1", 'make_a1(x, {})'.format(phi), 0., range_nE[0], range_nE[1])
    a1x.SetNpx(Npx)
    a1x.SetLineColor(rt.kBlue)
    a1y = rt.TF1("a1", 'make_a1(x, {})'.format(phi), 90., range_nE[0], range_nE[1])
    a1y.SetNpx(Npx)
    a1y.SetLineColor(rt.kGreen)

    print("| + Set TCanvas\n")
    c1 = TCanvas("c1", "c1")
    c1.SetGrid()
    a0.GetXaxis().SetRangeUser(3.53, 5.53)
    a0.GetYaxis().SetRangeUser(-3., 6.)
    rt.gStyle.SetOptTitle(0)

    a0.Draw("l")
    a0.SetLineColor(rt.kBlack)
    a0.SetLineWidth(4)
    a0.SetLineStyle(1)

    a1x.Draw("lsame")
    a1x.SetLineColor(rt.kRed)
    a1x.SetLineWidth(4)
    a1x.SetLineStyle(2)

    a1y.Draw("lsame")
    a1y.SetLineColor(rt.kBlue)
    a1y.SetLineWidth(4)
    a1y.SetLineStyle(3)

    rt.gStyle.SetLabelSize(0.05, "x")
    rt.gStyle.SetLabelSize(0.05, "y")
    rt.gStyle.SetTitleYOffset(1)
    rt.gStyle.SetTitleSize(0.05, "x")
    rt.gStyle.SetTitleSize(0.05, "y")

    a0.GetXaxis().SetTitle("neutron Energy [eV]")
    a0.GetYaxis().SetTitle("cross section [barn]")

    leg = TLegend(0.65, 0.70, 0.90, 0.90)
    setLegend(a0, "a_{0}", "l")
    setLegend(a1x, "a_{1x}", "l")
    setLegend(a1y, "a_{1y}", "l")
    leg.SetFillColor(0)
    leg.Draw()
    print("| + Close Process\n")
    return code.InteractiveConsole(globals()).interact()
