import flambaum_def
from ROOT import TCanvas, TGraph
from array import array
import code
from ROOT import TLegend


def setLegend(leg, entries):
    if isinstance(entries, (list, tuple)):
        for ent in entries:
            leg.AddEntry(ent, ent.GetName(), "l")


flambaum_def.a0_cal(4.53, 4.53, 0.00216e-3, 163.0e-3, 1., 1/2, 0., 0)

range_nE = [0., 10.]
nstep = 1000
step = (range_nE[1] - range_nE[0]) / nstep
f1 = rt.TF1()

def TGraph(nstep, x, y):
    """Make TGraph."""
    x, y = array('d'), array('d')
    for i in range(nstep):
        x.append(float(range_nE[0] + i * step))
        y.append(float(f1.Eval(x)))
        gr1 = TGraph()
        gr1.SetPoint(i, x, y)
    gr1.SetName(Form("g_%s", f1.GetName()));
    gr1.SetTitle(Form("Graph of %s", f1.GetName()));
    gr1.SetMarkerStyle(6);
    gr1.Draw('ACP')
    return TGraph(nstep, x, y)

def calcFlambaum(angle, mode):
    """calculate flambaum."""
    print("| + Processing flambaum...\n")
    parameterGet();
    if mode == phi:
        flambaum1 = rt.TF1("flambaum", 'makeFlambaum(x, 0, {})'.format(angle), range_nE[0], range_nE[1])
        flambaum2 = rt.TF1("flambaum", 'makeFlambaum(x, 90, {})'.format(angle), range_nE[0], range_nE[1])
        flambaum3 = rt.TF1("flambaum", 'makeFlambaum(x, 180, {})'.format(angle), range_nE[0], range_nE[1])
    elif mode== theta:
        flambaum1 = rt.TF1("flambaum", 'makeFlambaum(x, {}, 36)'.format(angle), range_nE[0], range_nE[1])
        flambaum1 = rt.TF1("flambaum", 'makeFlambaum(x, {}, 90)'.format(angle), range_nE[0], range_nE[1])
        flambaum1 = rt.TF1("flambaum", 'makeFlambaum(x, {}, 144)'.format(angle), range_nE[0], range_nE[1])
    flambaum1.SetNpx(Npx)
    flambaum2.SetNpx(Npx)
    flambaum3.SetNpx(Npx)
    print("| + Set TCanvas\n")
    c1 = ROOT.TCanvas('can', 'flambaum')
    c1.SetLogy(0)
    c1.SetGrid()
    gStyle.SetOptTitle(0)
    flambaum1.Draw("l");
    flambaum1.SetLineColor(kBlack);
    flambaum1.SetLineWidth(4);
    flambaum1.SetLineStyle(1);

    flambaum2.Draw("lsame");
    flambaum2.SetLineColor(kRed);
    flambaum2.SetLineWidth(4);
    flambaum2.SetLineStyle(2);

    flambaum3.Draw("lsame");
    flambaum3.SetLineColor(kBlue);
    flambaum3.SetLineWidth(4);
    flambaum3.SetLineStyle(3);

    gStyle.SetLabelSize(0.05,"x");
    gStyle.SetLabelSize(0.05,"y");
    gStyle.SetTitleYOffset(1);
    gStyle.SetTitleSize(0.05,"x");
    gStyle.SetTitleSize(0.05,"y");

    flambaum1.GetXaxis().SetRangeUser(4.03, 5.03);
    flambaum1.GetYaxis().SetRangeUser(0., 2.0);

    flambaum1.GetXaxis().SetTitle("neutron Energy [eV]");
    flambaum1.GetYaxis().SetTitle("cross section [barn]");
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
    return


def eachTermDraw(phi):
    """calculate each term."""
    print("| + Processing flambaum...\n")
    parameterGet()
    a0 = rt.TF1("a0", 'make_a0(x, {})'.format(phi), range_nE[0], range_nE[1])
    a0.SetNpx(Npx)
    a1x = rt.TF1("a1", 'make_a1(x, {})'.format(phi), 0., range_nE[0], range_nE[1])
    a1x.SetNpx(Npx)
    a1x.SetLineColor(kBlue)
    a1y = rt.TF1("a1", 'make_a1(x, {})'.format(phi), 90., range_nE[0], range_nE[1])
    a1y.SetNpx(Npx)
    a1y.SetLineColor(kGreen)

    print("| + Set TCanvas\n")
    c1 = TCanvas("c1", "c1");
    c1.SetGrid();
    a0.GetXaxis().SetRangeUser(3.53, 5.53);
    a0.GetYaxis().SetRangeUser(-3., 6.);
    gStyle.SetOptTitle(0);

    a0.Draw("l");
    a0.SetLineColor(kBlack);
    a0.SetLineWidth(4);
    a0.SetLineStyle(1);

    a1x.Draw("lsame");
    a1x.SetLineColor(kRed);
    a1x.SetLineWidth(4);
    a1x.SetLineStyle(2);

    a1y.Draw("lsame");
    a1y.SetLineColor(kBlue);
    a1y.SetLineWidth(4);
    a1y.SetLineStyle(3);

    gStyle.SetLabelSize(0.05,"x");
    gStyle.SetLabelSize(0.05,"y");
    gStyle.SetTitleYOffset(1);
    gStyle.SetTitleSize(0.05,"x");
    gStyle.SetTitleSize(0.05,"y");

    a0.GetXaxis().SetTitle("neutron Energy [eV]");
    a0.GetYaxis().SetTitle("cross section [barn]");

    leg = TLegend(0.65, 0.70, 0.90, 0.90)
    setLegend(a0, "a_{0}", "l")
    setLegend(a1x, "a_{1x}", "l")
    setLegend(a1y, "a_{1y}", "l")
    leg.SetFillColor(0)
    leg.Draw()


    print("| + Close Process\n")
    return


code.InteractiveConsole(globals()).interact()
