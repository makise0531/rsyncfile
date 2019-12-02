
//header
//for main program

 #include <stdio.h>
 #include <time.h>
 #include <sys/stat.h>
 #include <sys/types.h>
 #include <iostream>
 #include <fstream>
 #include <math.h>
 #include <string>
 #include <cmath>
 #include <complex>

 #include <TCanvas.h>
 #include <TFile.h>
 #include <TTree.h>
 #include <TChain.h>
 #include <TH1.h>
 #include <TH2D.h>
 #include <TF1.h>
 #include <TGraph.h>
 #include <TGraphErrors.h>
 #include <TMultiGraph.h>
 #include <TList.h>
 #include <TDirectory.h>
 #include <TLegend.h>
 #include <TStyle.h>
 #include <TLine.h>
 #include <TChain.h>
 #include <TCut.h>
 #include <TEnv.h>
 #include <TNtuple.h>
 #include <TString.h>
 #include <TROOT.h>
 #include <TLegend.h>
 #include <TSpectrum.h>

#include "include/AnglerMomentum.cpp"

const Double_t NeutronMass = 939.5654133e6; // neutron mass [MeV]
const Double_t protonMass  = 938.27813e6; // proton mass [MeV]
const Double_t targetMass  = 48.*protonMass + 63.*NeutronMass; // 111Cd [eV]
const Double_t I = 1./2.; // spin of target nucleus (111Cd)
const Double_t F = 0.; // spin of final state (ground state of 112Cd)

const Double_t spl = 299792458.; // [m/s]　speed of light(photon)
const Double_t hbar = 6.582119514e-16; // eVs
const Double_t barn = 6.50977023e5; // parameter for chaging unit to (eV)
const Double_t pi = TMath::Pi();

//----------- Doppler ---------------
Double_t m1 = NeutronMass;
Double_t m2 = targetMass;
Double_t alpha = m2/m1;
const Double_t kB = 6.2e18 * 1.3806503e-23;
const int ndoppler = 5000;//2019/11/20
Double_t all_doppler  = 10.;//2019/11/20
Double_t divide_doppler = all_doppler / ndoppler;
Double_t TK = 300;
//------------------------------------

const Int_t Npx = 10000;
Double_t range[2] = {0.1, 10.}; //2019/11/20

//----------------------------------------------------------------

const Int_t npara = 10;
Double_t Es[npara] = {}; // resonance energy
Double_t Ggs[npara]       = {}; // resonance width (exit channel)
Double_t Gns2g[npara]     = {}; // resonacne width (entrance channel)
Double_t Js[npara]        = {}; // total spin of compound state
Double_t ns = 0.; // number of s-wave resonance

Double_t Ep[npara] = {}; // resonance energy
Double_t Ggp[npara]       = {}; // resonance width (exit channel)
Double_t Gnp2g[npara]     = {}; // resonacne width (entrance channel)
Double_t Jp[npara]        = {}; // total spin of compound state
Double_t np = 0; // number of p-wave resonance

//----------------------------------------------------------------

//----------------------------------------------------------------
void sp_paratext(){
        Int_t n = 0;

        //Char_t s_file[100] = "./resopara_s994_111Cd.txt";
	//        Char_t s_file[100] = "./resopara_s275_111Cd.txt";
        Char_t s_file[100] = "./resopara_nega_111Cd.txt";
	//        Char_t s_file[100] = "./resopara_swave_111Cd.txt";
        ifstream s_data(s_file);

        while(!s_data.eof()) {
                s_data >> Es[n] >> Gns2g[n] >> Ggs[n] >> Js[n];
                n++;
        }
        ns = n - 1;

        Char_t p_file[100] = "./resopara_pwave_111Cd.txt";
        ifstream p_data(p_file);
        n = 0;
        while(!p_data.eof()) {
                p_data >> Ep[n] >> Gnp2g[n] >> Ggp[n] >> Jp[n];
                n++;
        }
        np = n - 1;

        return;
}


// Calculation of Pulse function
// (Ikeda-Carpenter function)
Double_t t0cal(Double_t x) {
        Double_t val=0;
        Double_t y;
        y = log10(x);
        Double_t parat0[4] = {1.23667e-01, -4.75351e-01, 2.42091e-03, 1.58608e-03};
        val = parat0[0] + parat0[1]*y + parat0[2]*y*y + parat0[3]*y*y*y;
        val = pow(10,val);
        return val;
}

double alphacal(double x){
        double val=0;
        double y;
        y = log10(x);
        double paraA[4] = {-7.18741e-02, 4.77821e-01, 1.59414e-03, -2.65068e-03};
        double paraB[4] = {6.58025e+00, 8.51179e+00, 3.45309e+00,4.66905e-01};
        if ((y < 8.0)&&(y >= -1.5)) {
                val = paraA[0] + paraA[1]*y + paraA[2]*y*y + paraA[3]*y*y*y;
        }
        if (y < -1.5) {
                val =
                        paraA[0] + paraA[1]*y + paraA[2]*y*y + paraA[3]*y*y*y +
                        paraB[0] + paraB[1]*y + paraB[2]*y*y + paraB[3]*y*y*y;
        }
        val = pow(10,val);
        return val;
}

Double_t betacal(Double_t x){
        Double_t val=0;
        Double_t y;
        y = log10(x);
        double paraA[4] = {-1.04912e+00, 5.80854e-01, -4.06260e-02, 3.56565e-03};
        Double_t paraB[4] = {2.86470e+00,4.75718e+00, 2.31579e+00, 3.24623e-01};
        if ((y < 8.0)&&(y >=-1.5)) {
                val = paraA[0] + paraA[1]*y + paraA[2]*y*y + paraA[3]*y*y*y;
        }
        if (y < -1.5) {
                val =
                        paraA[0] + paraA[1]*y + paraA[2]*y*y + paraA[3]*y*y*y +
                        paraB[0] + paraB[1]*y + paraB[2]*y*y + paraB[3]*y*y*y;
        }
        val = pow(10,val);
        return val;
}

Double_t Rcal(Double_t x) {
        Double_t val=0;
        Double_t y;
        y = log10(x);
        Double_t paraA[4] = {-1.04668e+00, -9.10505e-02, 9.59473e-02, -1.07923e-02};
        Double_t paraB[5]={ -3.35694e+00, -8.16440e+00, -6.70276e+00, -2.37025e+00,-3.21353e-01};
        if ((y < 8.0) && (y >=-0.8)) {
                val = paraA[0] + paraA[1]*y + paraA[2]*y*y + paraA[3]*y*y*y;
        }
        if (y < -0.8) {
                val = paraA[0] + paraA[1]*y + paraA[2]*y*y + paraA[3]*y*y*y + paraB[0] + paraB[1]*y + paraB[2]*y*y + paraB[3]*y*y*y + paraB[4]*y*y*y*y;
        }
        val = pow(10,val);
        return val;
}

Double_t Ikeda_Carpenter(Double_t tof, Double_t y){
        Double_t t0cal(Double_t);
        Double_t alphacal(Double_t);
        Double_t betacal(Double_t);
        Double_t Rcal(Double_t);

        Double_t ft;
        Double_t fa;
        Double_t fb;
        Double_t fr;
        Double_t ff;
        Double_t x;

        Double_t TOF;

        TOF = 72.3*21.5/TMath::Sqrt(y) + 2.50;
        x = tof - TOF;

        ft = t0cal(y);
        fa = alphacal(y);
        fb = betacal(y);
        fr = Rcal(y);

        if(x < ft) {
                ff = 0.;
        }else{
                ff = ((1-fr)*(0.5*fa * pow(fa*(x-ft),2)*exp(-fa*(x-ft))) + fr*0.5*fa * 2*pow(fa,2)*fb/pow(fa-fb,3)* (exp(-fb*(x-ft)) - exp(-fa*(x-ft))*(1+(fa-fb)*(x-ft)+0.5*pow(fa-fb,2)*pow(x-fa,2))));
                if(ff <= 0) {
                        ff = 0.;
                }
        }

        return ff;
}

//----------------------------------------------------------------

complex<Double_t> V1cal(Double_t E, Double_t Es, Double_t Gns2g,  Double_t Ggs, Double_t Js){
        // Gns2g: 2gGamma_n(s-wave) resonance width of neutron (Entrance channel)
        // Ggs: Gamma_gamma(s-wave) resonance width of gamma-ray (Exit channel)
        // Js: total spin of compound state of s-wave resonance
        complex<Double_t> V1;
        complex<Double_t> V1_total;

        Double_t ks = sqrt(2. * fabs(Es) * NeutronMass)/(hbar * spl);
        Double_t k;
        Double_t gs = (2.*Js + 1.)/(2.*(2.*I + 1.)); // g-factor
        //
        Double_t Ex = 0.;
        Double_t delta  = sqrt(4.*E*kB*TK/alpha);
        Double_t coefficient   = 1./(delta * sqrt(pi));
        //
        for(int i = 0; i < ndoppler; i++) {

                Ex = (i+1) * divide_doppler;
                k = sqrt(2. * Ex * NeutronMass)/(hbar * spl);

                V1 = complex<Double_t> (Ex - Es, (0.5*Gns2g/gs + Ggs)/2.);
                V1 = -1.0e14 * sqrt(pi/ks/k) * sqrt(0.5*Gns2g*Ggs) / V1;

                V1_total += exp(-pow((Ex - E)/delta, 2.)) * sqrt(Ex/E) * V1;
        }

        return coefficient * V1_total * divide_doppler;
}

complex<Double_t> V2cal(Double_t E, Double_t Ep, Double_t Gnp2g, Double_t Ggp, Double_t Jp, Double_t j, Double_t phi){
        complex<Double_t> V2;
        complex<Double_t> V2_total;

        Double_t kp = sqrt(2. * fabs(Ep) * NeutronMass)/(hbar * spl);
        Double_t k;
        Double_t gp = (2.*Jp + 1.)/(2.*(2.*I + 1.)); // g-factor
        //
        Double_t Ex = 0.;
        Double_t delta  = sqrt(4.*E*kB*TK/alpha);
        Double_t coefficient   = 1./(delta * sqrt(pi));
        //
        for(int i = 0; i < ndoppler; i++) {

                Ex = (i+1) * divide_doppler;
                k = sqrt(2. * Ex * NeutronMass)/(hbar * spl);

                V2 = complex<Double_t> (Ex - Ep, (0.5*Gnp2g/gp + Ggp)/2.);
                V2 = -1.0e14 * sqrt(pi/kp/k) * sqrt(0.5*Gnp2g*Ggp) / V2;

                if(j == 1./2.) {
                        V2 = V2 * cos(phi*pi/180.);
                }else if(j == 3./2.) {
                        V2 = V2 * sin(phi*pi/180.);
                }

                V2_total += exp(-pow((Ex - E)/delta, 2.)) * sqrt(Ex/E) * V2;
        }

        return coefficient * V2_total * divide_doppler;
}

Double_t Pcal(Double_t J1, Double_t J2, Double_t j1, Double_t j2, Double_t k, Double_t I, Double_t F){
        Double_t P = 0.; // Wiggner's 6j signal
        P = pow(-1,J1+J2+j2+I+F) * 3./2. * sqrt((2*J1+1) * (2*J2+1) * (2*j1+1) * (2*j2+1)) * Wigner6j(k,j1,j2,I,J2,J1) * Wigner6j(k,1.,1.,F,J1,J2);
        return P;
}

Double_t a0cal(Double_t E, Double_t Ene, Double_t gGn2, Double_t Gg, Double_t J, Double_t j, Double_t phi, Int_t Vpattern){
        // Vpattern = 0 (s-wave),  Vpattern =1 (p-wave)
        complex<Double_t> V1cal(Double_t, Double_t, Double_t, Double_t, Double_t);
        complex<Double_t> V2cal(Double_t, Double_t, Double_t, Double_t, Double_t, Double_t, Double_t);

        Double_t a0 = 0.;
        if(Vpattern == 0) {
                a0 = abs(V1cal(E, Ene, gGn2, Gg, J) * conj(V1cal(E, Ene, gGn2, Gg, J)));
        }else{
                a0 = abs(V2cal(E, Ene, gGn2, Gg, J, j, phi) * conj(V2cal(E, Ene, gGn2, Gg, J, j, phi)));
        }
        return a0;
}

Double_t a1cal(Double_t E, Double_t *Es, Double_t *Gns2g, Double_t *Ggs, Double_t *Js, Double_t *Ep, Double_t *Gnp2g, Double_t *Ggp, Double_t *Jp, Double_t j, Double_t phi, Int_t ns, Int_t np, Double_t F){
        complex<Double_t> V1cal(Double_t, Double_t, Double_t, Double_t, Double_t);
        complex<Double_t> V2cal(Double_t, Double_t, Double_t, Double_t, Double_t, Double_t, Double_t);
        Double_t Pcal(Double_t, Double_t, Double_t, Double_t, Double_t, Double_t, Double_t);

        Double_t a1 = 0.;
        for(int i=0; i < ns; i++) {
                for(int n=0; n < np; n++) {
                        a1 += 2.*real(V1cal(E, Es[i], Gns2g[i], Ggs[i], Js[i])
                                      * conj(V2cal(E, Ep[n], Gnp2g[n], Ggp[n], Jp[n], j, phi) ))
                              * Pcal(Js[i], Jp[n], 1./2., j, 1., I, F);
                }
        }
        return a1;
}

Double_t a0func(Double_t E, Double_t phi){
        Double_t a0 = 0.;
        Double_t j[2] = {1./2., 3./2.};
        //s波
        for(int i=0; i < ns; i++) {
                a0 += a0cal(E, Es[i], Gns2g[i], Ggs[i], Js[i], 0., 0., 0);
        }
        //p波
        for(int i=0; i < np; i++) {
                for(int n=0; n < 2; n++) {
                        a0 += a0cal(E, Ep[i], Gnp2g[i], Ggp[i], Jp[i], j[n], phi, 1);
                }
        }
        return a0;
}

Double_t a1func(Double_t E, Double_t phi){
        Double_t a1 = 0.;
        Double_t j[2] = {1./2., 3./2.};

        for(int i=0; i < 2; i++) {
                a1 += a1cal(E, Es, Gns2g, Ggs, Js, Ep, Gnp2g, Ggp, Jp, j[i], phi, ns, np, F);
        }
        return a1;
}

Double_t aterm_function(Double_t E, Double_t phi, Double_t theta){
        Double_t a_all = 0.;
        Double_t a0 = 0.;
        Double_t a1 = 0.;

        a0 = a0func(E, phi);
        a1 = a1func(E, phi);
        a_all = (1./2.)*(a0 + a1*(cos(theta*pi/180.)));

        return a_all;
}

Double_t makeFlambaum(Double_t E, Double_t phi, Double_t theta){
        return aterm_function(E, phi, theta);
}

void calcFlambaum(Int_t angle, Int_t mode){
        cout<<"Fight!!"<<endl;
        sp_paratext();
        ////// Function /////
        Double_t DrawRange[2];
        DrawRange[0] = range[0];
        DrawRange[1] = range[1];

        TF1 *flambaum1;
        TF1 *flambaum2;
        TF1 *flambaum3;

        if(mode==0)
        {
                flambaum1 = new TF1("flambaum", Form("makeFlambaum(x, 0, %d)",angle), DrawRange[0], DrawRange[1]);
                flambaum2 = new TF1("flambaum", Form("makeFlambaum(x, 90, %d)",angle), DrawRange[0], DrawRange[1]);
                flambaum3 = new TF1("flambaum", Form("makeFlambaum(x, 180, %d)",angle), DrawRange[0], DrawRange[1]);
        }else if(mode==1)
        {
                flambaum1 = new TF1("flambaum", Form("makeFlambaum(x, %d,36)",angle), DrawRange[0], DrawRange[1]);
                flambaum2 = new TF1("flambaum", Form("makeFlambaum(x, %d,90)",angle), DrawRange[0], DrawRange[1]);
                flambaum3 = new TF1("flambaum", Form("makeFlambaum(x, %d,144)",angle), DrawRange[0], DrawRange[1]);
        }
        flambaum1->SetNpx(Npx);
        flambaum2->SetNpx(Npx);
        flambaum3->SetNpx(Npx);
        cout<<"Fight!!"<<endl;
        ///// TCanvas /////
        TCanvas *c1 = new TCanvas("c1", "c1");
        c1->SetLogy(0);
        c1->SetGrid();

        gStyle->SetOptTitle(0);

        flambaum1->Draw("l");
        flambaum1->SetLineColor(kBlack);
        flambaum1->SetLineWidth(4);
        flambaum1->SetLineStyle(1);

        flambaum2->Draw("lsame");
        flambaum2->SetLineColor(kRed);
        flambaum2->SetLineWidth(4);
        flambaum2->SetLineStyle(2);

        flambaum3->Draw("lsame");
        flambaum3->SetLineColor(kBlue);
        flambaum3->SetLineWidth(4);
        flambaum3->SetLineStyle(3);

        gStyle->SetLabelSize(0.05,"x");
        gStyle->SetLabelSize(0.05,"y");
        gStyle->SetTitleYOffset(1);
        gStyle->SetTitleSize(0.05,"x");
        gStyle->SetTitleSize(0.05,"y");

        flambaum1->GetXaxis()->SetRangeUser(4.03, 5.03);
        //flambaum1->GetYaxis()->SetRangeUser(0., 2.0);
/*
        flambaum1->GetXaxis()->SetRangeUser(6.5, 7.5);
        flambaum1->GetYaxis()->SetRangeUser(-0.5, 1.4);
 */
        flambaum1->GetXaxis()->SetTitle("neutron Energy [eV]");
        flambaum1->GetYaxis()->SetTitle("cross section [barn]");
        cout<<"Fight!!"<<endl;

        TLegend *legend = new TLegend( 0.65, 0.70, 0.90, 0.90);  //（）の中は位置の指定（左下の x , y 、右上の x , y ）
        if(mode==0) {
                legend->AddEntry( flambaum1, Form("#phi=0 , #theta=%d",angle), "l");// AddEntry( pointer , "interpretation" , "option" )
                legend->AddEntry( flambaum2, Form("#phi=90 , #theta=%d",angle), "l"); // option は　"f"=box, "l"="L"=line, "p"=marker
                legend->AddEntry( flambaum3, Form("#phi=180 , #theta=%d",angle), "l");
        } else if(mode==1) {
                legend->AddEntry( flambaum1, Form("#phi=%d , #theta=36",angle), "l");   // AddEntry( pointer , "interpretation" , "option" )
                legend->AddEntry( flambaum2, Form("#phi=%d , #theta=90",angle), "l"); // option は　"f"=box, "l"="L"=line, "p"=marker
                legend->AddEntry( flambaum3, Form("#phi=%d , #theta=144",angle), "l");
        }
        legend->SetFillColor(0);
        legend->Draw();
        cout<<"Fight!!"<<endl;
        return;
}



void eachTermDraw(Double_t phi){

        Double_t DrawRange[2];
        DrawRange[0] = range[0];
        DrawRange[1] = range[1];


        sp_paratext();
        TF1 *a0 = new TF1("a0", Form("a0func(x, %lf)", phi), DrawRange[0], DrawRange[1]);
        a0->SetNpx(Npx);

        TF1 *a1x = new TF1("a1x", Form("a1func(x, %lf)", 0.), DrawRange[0], DrawRange[1]);
        a1x->SetNpx(Npx);
        a1x->SetLineColor(kBlue);
        TF1 *a1y = new TF1("a1y", Form("a1func(x, %lf)", 90.), DrawRange[0], DrawRange[1]);
        a1y->SetNpx(Npx);
        a1y->SetLineColor(kGreen);

        ///// TCanvas /////
        TCanvas *ccc = new TCanvas("ccc", "ccc");
        cout<<"Fight!!"<<endl;
        ccc->cd();
        ccc->SetGrid();
        //ccc->SetLogy();
        cout<<"| Drawing a0 term..."<<endl;

        gStyle->SetOptTitle(0);

        a0->Draw("l");
        a0->SetLineColor(kBlack);
        a0->SetLineWidth(4);
        a0->SetLineStyle(1);


        cout<<"| Drawing a1x term..."<<endl;
        a1x->Draw("lsame");
        a1x->SetLineColor(kRed);
        a1x->SetLineWidth(4);
        a1x->SetLineStyle(2);
        cout<<"| Drawing a1y term..."<<endl;
        a1y->Draw("lsame");
        a1y->SetLineColor(kBlue);
        a1y->SetLineWidth(4);
        a1y->SetLineStyle(3);

        gStyle->SetLabelSize(0.05,"x");
        gStyle->SetLabelSize(0.05,"y");
        gStyle->SetTitleYOffset(1);
        gStyle->SetTitleSize(0.05,"x");
        gStyle->SetTitleSize(0.05,"y");

        a0->GetXaxis()->SetRangeUser(4.53 - 0.7, 4.53 + 0.7);
        a0->GetYaxis()->SetRangeUser(-0.6, 2.8);

        a0->GetXaxis()->SetTitle("neutron Energy [eV]");
        a0->GetYaxis()->SetTitle("cross section [barn]");

        cout<<"| Drawing Legend..."<<endl;

        TLegend *legend = new TLegend( 0.65, 0.70, 0.90, 0.90);          //（）の中は位置の指定（左下の x , y 、右上の x , y ）

        legend->AddEntry( a0, "a_{0}", "l");         // AddEntry( pointer , "interpretation" , "option" )

        legend->AddEntry( a1x, "a_{1x}", "l");            // option は　"f"=box, "l"="L"=line, "p"=marker
        legend->AddEntry( a1y, "a_{1y}", "l");

        legend->SetFillColor(0);
        legend->Draw();

        return;
}

void eachTermAsy(){
        sp_paratext();
        Double_t DrawRange[2];
        DrawRange[0] = range[0];
        DrawRange[1] = range[1];

        TF1 *a0 = new TF1("a0", "a0func(x, 0)", DrawRange[0], DrawRange[1]);
        a0->SetNpx(Npx);
        TF1 *a1x = new TF1("a1x", Form("a1func(x, %lf)", 0.), DrawRange[0], DrawRange[1]);
        a1x->SetNpx(Npx);
        a1x->SetLineColor(kBlue);
        TF1 *a1y = new TF1("a1y", Form("a1func(x, %lf)", 90.), DrawRange[0], DrawRange[1]);
        a1y->SetNpx(Npx);
        a1y->SetLineColor(kGreen);

        cout<<"Fight"<<endl;
        cout<<"| + Set TCanvas"<<endl;
        TCanvas *c1 = new TCanvas("c1", "c1");
        cout<<"Fight!!"<<endl;
        c1->SetGrid();
        cout<<"| Drawing a0 term..."<<endl;

        gStyle->SetOptTitle(0);

        a0->Draw("l");
        a0->SetLineColor(kBlack);
        a0->SetLineWidth(4);
        a0->SetLineStyle(1);
        cout<<"| Drawing a1x term..."<<endl;
        a1x->Draw("lsame");
        a1x->SetLineColor(kRed);
        a1x->SetLineWidth(4);
        a1x->SetLineStyle(2);
        cout<<"| Drawing a1y term..."<<endl;
        a1y->Draw("lsame");
        a1y->SetLineColor(kBlue);
        a1y->SetLineWidth(4);
        a1y->SetLineStyle(3);

        gStyle->SetLabelSize(0.05,"x");
        gStyle->SetLabelSize(0.05,"y");
        gStyle->SetTitleYOffset(1);
        gStyle->SetTitleSize(0.05,"x");
        gStyle->SetTitleSize(0.05,"y");

        a0->GetXaxis()->SetRangeUser(3.83, 5.23);
        a0->GetYaxis()->SetRangeUser(-0.3, 2.8);

        /*
              a0->GetXaxis()->SetRangeUser(6.5, 7.5);
              a0->GetYaxis()->SetRangeUser(-2.0, 2.0);
         */
        a0->GetXaxis()->SetTitle("neutron Energy [eV]");
        a0->GetYaxis()->SetTitle("cross section [barn]");

        cout<<"| Drawing Legend..."<<endl;

        TLegend *legend = new TLegend( 0.65, 0.70, 0.90, 0.90);  //（）の中は位置の指定（左下の x , y 、右上の x , y ）

        legend->AddEntry( a0, "a_{0}", "l"); // AddEntry( pointer , "interpretation" , "option" )
        legend->AddEntry( a1x, "a_{1x}", "l"); // option は　"f"=box, "l"="L"=line, "p"=marker
        legend->AddEntry( a1y, "a_{1y}", "l");

        legend->SetFillColor(0);
        legend->Draw();

        double nL_a0;
        double nH_a0;
        double Asym_a0;

        nL_a0 = a0->Integral(Ep[0] - 2.*Ggp[0], Ep[0]);
        nH_a0 = a0->Integral(Ep[0], Ep[0] + 2.*Ggp[0]);
        Asym_a0 = (nL_a0 + nH_a0);

        double nL_a1x;
        double nH_a1x;
        double Asym_a1x;

        nL_a1x = a1x->Integral(Ep[0] - 2.*Ggp[0], Ep[0]);
        nH_a1x = a1x->Integral(Ep[0], Ep[0] + 2.*Ggp[0]);
        Asym_a1x = (nL_a1x - nH_a1x)/Asym_a0;

        double nL_a1y;
        double nH_a1y;
        double Asym_a1y;

        nL_a1y = a1y->Integral(Ep[0] - 2.*Ggp[0], Ep[0]);
        nH_a1y = a1y->Integral(Ep[0], Ep[0] + 2.*Ggp[0]);
        Asym_a1y = (nL_a1y - nH_a1y)/Asym_a0;

        cout << "Ene_left = " << Ep[0] - 2.*Ggp[0] << " :  Ene_peak = " << Ep[0] << " :  Ene_right = " << Ep[0] + 2.*Ggp[0] <<endl;

        double experimental_Asyvalue = -0.0148;
        double experimental_Asyerror = 0.03761;
        cout << "a1x left= " << nL_a1x << endl;
        cout << "a1x right= " << nH_a1x << endl;
        cout << "a1x 分子 = " << nL_a1x - nH_a1x << endl;
        cout << "a1x 分母 = " << nL_a1x + nH_a1x << endl;
        cout << "a1x = " << Asym_a1x << endl;


        cout << "a1y left= " << nL_a1y << endl;
        cout << "a1y right= " << nH_a1y << endl;
        cout << "a1x 分子 = " << nL_a1y - nH_a1y << endl;
        cout << "a1x 分母 = " << nL_a1y + nH_a1y << endl;
        cout << "a1y = " << Asym_a1y << endl;
        cout << "" << endl;
        cout << "分母のa0 = " << Asym_a0 << endl;

        double a_value = Asym_a1x / Asym_a1y;
        double b_value = experimental_Asyvalue / Asym_a1y;
        double b_error = experimental_Asyerror / Asym_a1y;

        cout << Asym_a1x << "x + " << Asym_a1y << "y = " << experimental_Asyvalue << " ± " << experimental_Asyerror << endl;
        //一次関数
        cout << "y = " << -a_value << "x + " << b_value << " ± " << b_error << endl;

        return;
}
