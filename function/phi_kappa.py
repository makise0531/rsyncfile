import sympy
import sys
import math
import ROOT as rt
import numpy as np
import scipy as sp
from ROOT import TCanvas, TGraph, TLegend
from array import array
import code

a_value = -0.282843
b_value = -0.181357
b_error = 0.460868
x = 0
y = 0

Npx = 36000

I = 1./2.
J = 1.
F = 0.

def phi_result():
    a = a_value
    b = b_value
    berr = b_error

    rt.gStyle.SetOptStat(0)

    cphi = TCanvas("c1", "c1", 800, 750)
    cphi.cd()
    cphi.SetGrid()

    h2 = rt.TH2D("h2", "", 10, -5, 5, 10, -5, 5)
    h2.Draw()
    h2.GetXaxis().SetRangeUser(-1.5, 1.5)
    h2.GetYaxis().SetRangeUser(-1.5, 1.5)

    errorbar = rt.TGraphErrors()
    for i in range(40000):
        errorbar.SetPoint(i, 0.0001 * i - 2.0, a * (0.0001 * i - 2.0) + b)
        errorbar.SetPointError(i, 0, berr)
    errorbar.SetFillColor(rt.kRed)
    errorbar.SetLineColor(rt.kRed)
    errorbar.SetFillStyle(3010)
    errorbar.Draw("same a2")

    cir = rt.TArc(0., 0., 1., 0., 360.)
    cir.SetFillStyle(0)
    cir.Draw("same")

    first1 = rt.TF1("first1", "[0]*x + [1]", -2.0, 2.0)
    first1.SetNpx(Npx)
    first1.FixParameter(0, a)
    first1.FixParameter(1, b)
    # first1.SetParError(1,b_error);
    first1.SetLineColor(rt.kRed)
    first1.Draw("same")

    first2 = rt.TF1("first2", "[0]*x + [1]", -2.0, 2.0)
    first2.SetNpx(Npx)
    first2.FixParameter(0, a)
    first2.FixParameter(1, b - berr)
    # firrt1.SetParError(1,b_error);
    first2.SetLineColor(rt.kRed)
    first2.Draw("same")

    first3 = rt.TF1("first3", "[0]*x + [1]", -2.0, 2.0)
    first3.SetNpx(Npx)
    first3.FixParameter(0, a)
    first3.FixParameter(1, b + berr)
    # firrt1.SetParError(1,b_error);
    first3.SetLineColor(rt.kRed)
    first3.Draw("same")

    # upper_circle = rt.TF1("upper_circle", "TMath::Sqrt(1 - x*x)", -0.99999, 0.99999)
    # down_circle = rt.TF1("down_circle", "-TMath::Sqrt(1 - x*x)", -0.99999, 0.99999)
    # h2.Draw()
    # upper_circle.SetNpx(Npx)
    # upper_circle.Draw("sames")
    # upper_circle.SetLineColor(rt.kBlack)
    # down_circle.SetNpx(Npx)
    # down_circle.Draw("sames")
    # down_circle.SetLineColor(rt.kBlack)

    # phi_gr->SetPointError(i, 0, b_error);
    return code.InteractiveConsole(globals()).interact()


def renritsu(error):
    sympy.var('x y')
    a = a_value
    b = b_value
    berr = error
    # f(x)とg(x)を定義
    fx = - y + a * x + (b + berr)
    gx = x**2 + y**2 - 1
    # 連立方程式 f(x)=0, g(x)=0 の解
    s = sympy.solve([fx, gx], [x, y])
    return s


s_maxerror = renritsu(b_error)
s_value = renritsu(0)
s_minerror = renritsu(-b_error)
print(s_maxerror)
print(s_value)
print(s_minerror)


def phi_value():
    x0 = 0.
    x1 = -0.853591272863153
    x2 = 0.999994454098032
    x3 = -0.994980766942492
    x4 = 0.899989007265006
    x5 = -0.924707703675191
    x6 = 0.588321003085340

    y0 = 0.
    y1 =  0.520943316390433
    y2 =  -0.00333043138044954
    y3 = 0.100066345064315
    y4 = -0.435912590781856
    y5 = -0.380677898969398
    y6 = -0.808627477475667

    max_phi1 = math.atan2(y1-y0, x1-x0)
    max_phi2 = math.atan2(y2-y0, x2-x0)
    value_phi1 = math.atan2(y3-y0, x3-x0)
    value_phi2 = math.atan2(y4-y0, x4-x0)
    min_phi1 = math.atan2(y5-y0, x5-x0)
    min_phi2 = math.atan2(y6-y0, x6-x0)

    phi1 = math.degrees(value_phi1)
    phi2 = math.degrees(value_phi2) + 360
    phi1_min = math.degrees(min_phi1) + 360
    phi2_min = math.degrees(min_phi2) + 360
    phi1_max = math.degrees(max_phi1)
    phi2_max = math.degrees(max_phi2)

    phi1_error = abs(phi1_max - phi1)
    phi2_error = abs(phi2_max + 360 - phi2)
    print(" ~ Calculate φ1. ~ ")
    print("")
    print("誤差最大")
    print("max_error deg", phi1_max)
    print("")
    print("誤差最小")
    print("min_error deg", phi1_min)
    print("φ1 value")
    print("true deg", phi1, " ± ", phi1_error)

    print(" ~ Calculate φ2. ~ ")
    print("誤差最大")
    print("max_error deg", phi2_max)
    print("")
    print("誤差最小")
    print("min_error deg", phi2_min)
    print("φ2 value")
    print("true deg", phi2, " ± ", phi2_error)
    print("")
    return


def tan_cal(x):
    tancal1 = ((1/3) * (1. - np.sqrt(2) * np.tan(x * np.pi / 180.)))
    tancal2 = np.tan(0.001 * x * np.pi / 180.)
    return tancal1


for i in range(360):
    phix = tan_cal(i)
    print(i, phix)


def kappa():
    """Calculate kappa."""

    phi_value()

    c_kappa = TCanvas("c_kappa", "c_kappa")
    c_kappa.SetLogy(1)
    c_kappa.SetGrid()

    # h_d = rt.TH2D("hd", "hd", 360, 0, 360, 100000, 0, 1000)

    #h_d.Draw()

    func = rt.TF1("kappa_phi", "abs((-1.) ** (2. * [0] + 1.) * ([0] / ([0] + 1.)) * (1. - (1. / 2.) * TMath::Sqrt((2. * [0] + 3.) / [0]) * tan(x * TMath::Pi() / 180.)))", 0, 360)

    func.FixParameter(0, I)
    func.SetNpx(Npx)

    func.SetMinimum(1.0e-3)
    func.SetMaximum(1.0e3)
    func.SetLineColor(rt.kBlack)
    func.GetXaxis().SetRangeUser(0, 360)
    func.GetYaxis().SetRangeUser(1.0e-3, 1.0e3)
    func.GetXaxis().SetTitle("#phi [deg]")
    func.GetXaxis().SetTitleSize(0.05)
    func.GetXaxis().SetTitleOffset(0.95)
    func.GetYaxis().SetTitle("|#kappa(J)|")
    func.GetYaxis().SetTitleSize(0.05)
    func.GetYaxis().SetTitleOffset(0.90)
    func.SetTitle("#kappa(J) vs #phi")
    func.Draw()

    kappa1_max = func.Eval(174.26 - 25.65)
    kappa1 = func.Eval(174.26)
    kappa1_min = func.Eval(174.26 + 25.65)
    kappa2_max = func.Eval(334.16 - 25.65)
    kappa2 = func.Eval(334.16)
    kappa2_min = func.Eval(334.16 + 25.65 - 360)

    print("")
    print("1σの精度でκ(J = 1)を評価 → φ = 174.26")
    print("")
    print("kappa1_max = ", kappa1_max)
    print("kappa1 = ", kappa1)
    print("kappa1_min = ", kappa1_min)
    print("")
    print("kappa2_maxerror = ", kappa1_max - kappa1)
    print("kappa2_minerror = ", kappa1_min - kappa1)

    print("")
    print("1σの精度でκ(J = 1)を評価 → φ = 334.16")
    print("")
    print("kappa2_max = ", kappa2_max)
    print("kappa2 = ", kappa2)
    print("kappa2_min = ", kappa2_min)
    print("")
    print("kappa2_maxerror = ", kappa2_max - kappa2)
    print("kappa2_minerror = ", kappa2_min - kappa2)

    line1_error = rt.TBox(174.26 - 25.65, 1.01e-3, 174.26 + 25.65, 9.99e2)
    line1_error.SetLineColor(rt.kRed)
    line1_error.SetFillColor(rt.kRed)
    line1_error.SetFillStyle(3004)
    line1_error.SetLineStyle(3001)  # 1=line,2=broken,3=dotted,4=broken-dot,5=long-broken-dot
    line1_error.Draw("same")

    line2_error = rt.TBox(334.16 - 25.65, 1.01e-3, 334.16 + 25.65, 9.99e2)
    line2_error.SetLineColor(rt.kRed)
    line2_error.SetFillColor(rt.kRed)
    line2_error.SetFillStyle(3004)
    line2_error.SetLineStyle(3001)  # 1=line,2=broken,3=dotted,4=broken-dot,5=long-broken-dot
    line2_error.Draw("same")

    line2_error2 = rt.TBox(0, 1.01e-3, 334.16 + 25.65 - 360, 9.99e2)
    line2_error2.SetLineColor(rt.kRed)
    line2_error2.SetFillColor(rt.kRed)
    line2_error2.SetFillStyle(3004)
    line2_error2.SetLineStyle(3001)  # 1=line,2=broken,3=dotted,4=broken-dot,5=long-broken-dot
    line2_error2.Draw("same")

    line1 = rt.TLine(174.70, 1.0e-3, 174.26, 1.0e3)
    line1.SetLineColor(rt.kRed)
    line1.SetLineWidth(1)
    line1.SetLineStyle(1)  # 1=line,2=broken,3=dotted,4=broken-dot,5=long-broken-dot
    line1.Draw("same")

    line1_min = rt.TLine(174.26 - 25.65, 1.0e-3, 174.26 - 25.65, 1.0e3)
    line1_min.SetLineColor(rt.kRed)
    line1_min.SetLineWidth(1)
    line1_min.SetLineStyle(5)  # 1=line,2=broken,3=dotted,4=broken-dot,5=long-broken-dot
    line1_min.Draw("same")

    line1_max = rt.TLine(174.26 + 25.65, 1.0e-3, 174.26 + 25.65, 1.0e3)
    line1_max.SetLineColor(rt.kRed)
    line1_max.SetLineWidth(1)
    line1_max.SetLineStyle(5)  # 1=line,2=broken,3=dotted,4=broken-dot,5=long-broken-dot
    line1_max.Draw("same")

    line2 = rt.TLine(334.16, 1.0e-3, 334.16, 1.0e3)
    line2.SetLineColor(rt.kRed)
    line2.SetLineWidth(1)
    line2.SetLineStyle(1)  # 1=line,2=broken,3=dotted,4=broken-dot,5=long-broken-dot
    line2.Draw("same")

    line2_min = rt.TLine(334.16 - 25.65, 1.0e-3, 334.16 - 25.65, 1.0e3)
    line2_min.SetLineColor(rt.kRed)
    line2_min.SetLineWidth(1)
    line2_min.SetLineStyle(5)  # 1=line,2=broken,3=dotted,4=broken-dot,5=long-broken-dot
    line2_min.Draw("same")

    line2_max = rt.TLine(334.16 + 25.65 - 360, 1.0e-3, 334.16 + 25.65 - 360, 1.0e3)
    line2_max.SetLineColor(rt.kRed)
    line2_max.SetLineWidth(1)
    line2_max.SetLineStyle(5)  # 1=line,2=broken,3=dotted,4=broken-dot,5=long-broken-dot
    line2_max.Draw("same")


    #func2 = rt.TF1("kappa_phi2", "(1/3) * (1. - TMath::Sqrt(2) * tan(x * TMath::Pi() / 180.))", 0., 360.)
    return code.InteractiveConsole(globals()).interact()


kappa()

phi_result()
