#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "RooBreitWigner.h"
#include "RooFFTConvPdf.h"
#include "RooGenericPdf.h"
#include "RooAddPdf.h"
#include "RooHist.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "TTree.h"
#include "TH1D.h"
#include "TRandom.h"
#include <RooUnblindUniform.h>
//#include "RooMinuit.h"
//#include "Roo2JBW.h"
#include "RooMinimizer.h"
#include "RooChi2Var.h"
#include <time.h>

using namespace RooFit;

void scnsp(){

    //time_t time1 = clock()/CLOCKS_PER_SEC;	
    //gROOT->ProcessLine(".x lhcbStyle.C");

    // Observables
    // Here we define the observable we are going to study with the code. In this case the observable is the mass of the D^0 + pion 
    // x1 and x2 are the two extremes of the energy interval we are looking in the analysis
    Double_t x1 = 2004.41, x2 = 2020.0;
    RooRealVar m_D0pi("m_D0pi", "m_D0pi", x1, x2, "Mev/c^2");
    //the bin are set as 520 so that each bin has an excact length of 0.03 MeV/c^2
    m_D0pi.setBins(520);

    // This line is used to define an object that we will use as a tagging variable, since the porpouse of the code is to do a simultaneous fit of both D^++ and D^*- 
    // this RooCategory object is what separate a ' + ' event from a ' - ' event 
    RooCategory * q = new RooCategory("q", "q");
    q->defineType("Dsp", 1);
    q->defineType("Dsm", -1);

    // this line defines a set in which we will input our observalbes and parameters later
    RooArgSet * obs = new RooArgSet();
    obs->add(m_D0pi);
    obs->add(*q);

    // This line defines two observables taht are equal for both mass distributions 
    RooRealVar rho("rho", "rho esponente di f1", 0.05, -5, 5);
    RooRealVar m0("m0", "m0", 2004.4, "MeV/c^2"); // threshold parameter, set as the mass energy of D^0+pi

    // Histograms
    // This lines opens a the file.root in which there are the histogram of events used for the data. At first the file is opened then the histograms are extracted and used to 
    // create the TH1D object. At last, the RooDataHist object is created using both TH1D and the tagging object q, preaviously defined
    TFile * file = new TFile("..path..");
    TH1D * h_dm_plus = (TH1D*)file->Get("h_dm_plus");
    TH1D * h_dm_minus = (TH1D*)file->Get("h_dm_minus");
    RooDataHist * hist = new RooDataHist("hist", "hist", m_D0pi, Index(*q), Import("Dsp", *h_dm_plus), Import("Dsm",*h_dm_minus));

    //blinding
    // this string is used for "blinding", which means that the parameter of interest is modified in order not to have a bias while coding
    TString blindString = "SCNSP";
    /*
    if(file.contains("2017") || file.contains("2018") || file.contains("2017_2018")){
        if(file.contains("KK")) blindString = "newDeltaACPRun3KK";
        if(file.contains("pipi")) blindString = "newDeltaACPRun3pipi";
    }
    */

    // this line defines our parameter of interest and then defines a RooUnblindUniform object which blinds the parameter with an offset. this will only modify the value
    // of the parameter, not ist error
    RooRealVar * A_CP = new RooRealVar("A_CP", "A_CP", -0.01, -1, 1);
    RooUnblindUniform * A_CP_blind = new RooUnblindUniform("A_CP_blind", "A_CP_blind", blindString, 0.5, *A_CP);

    RooRealVar * A_CP_bkg = new RooRealVar("A_CP_bkg", "A_CP_bkg", -0.01, -1, 1);


    // The following lines of code defines the functions used in the model, 3 Johnson functions and 1 Gaussian function. Note that theese definitions are made both for 
    // D^*+ distribution and D^*- distributions
    // Model : 3 Johnson functions D*+
    RooRealVar mu1p("mu1p", "mu1p", 2010.0, 2009, 2011);
    RooRealVar sigj1p("sigj1p", "sigj1p", 0.8, 0, 1);
    RooRealVar dj1p("dj1p", "dj1p", 1.7, 0.01, 10);
    RooRealVar gamj1p("gamj1p", "gamj1p",-0.9, -10, 10);
    RooAbsPdf * j1p = new RooGenericPdf("j1p", "(@0>@5)*pow((@0-@5),1+@6)*@4/(@2*TMath::Sqrt(2*TMath::Pi())*TMath::Sqrt(1+((@0-@1)*(@0-@1)/(@2*@2))))*TMath::Exp(-0.5*(@3+@4*TMath::ASinH((@0-@1)/@2))*(@3+@4*TMath::ASinH((@0-@1)/@2)))", 
                        RooArgSet(m_D0pi, mu1p, sigj1p, gamj1p, dj1p, m0, rho));

    RooRealVar mu2p("mu2p", "mu2p", 2010.0, 2009, 2011);
    RooRealVar sigj2p("sigj2p", "sigj2p", 0.8, 0, 1);
    RooRealVar dj2p("dj2p", "dj2p", 1, 0.001, 10);
    RooRealVar gamj2p("gamj2p", "gamj2p", -0.6, -10, 10);
    RooAbsPdf * j2p = new RooGenericPdf("j2p", "(@0>@5)*pow((@0-@5),1+@6)*@4/(@2*TMath::Sqrt(2*TMath::Pi())*TMath::Sqrt(1+((@0-@1)*(@0-@1)/(@2*@2))))*TMath::Exp(-0.5*(@3+@4*TMath::ASinH((@0-@1)/@2))*(@3+@4*TMath::ASinH((@0-@1)/@2)))", 
                        RooArgSet(m_D0pi, mu2p, sigj2p, gamj2p, dj2p, m0, rho));

    RooRealVar mu3p("mu3p", "mu3p", 2010.0, 2009, 2011);
    RooRealVar sigj3p("sigj3p", "sigj3p", 0.8, 0, 1);
    RooRealVar dj3p("dj3p", "dj3p", 1, 0.01, 10);
    RooRealVar gamj3p("gamj3p", "gamj3p", -0.6, -10, 10);
    RooAbsPdf * j3p = new RooGenericPdf("j3p", "(@0>@5)*pow((@0-@5),1+@6)*@4/(@2*TMath::Sqrt(2*TMath::Pi())*TMath::Sqrt(1+((@0-@1)*(@0-@1)/(@2*@2))))*TMath::Exp(-0.5*(@3+@4*TMath::ASinH((@0-@1)/@2))*(@3+@4*TMath::ASinH((@0-@1)/@2)))", 
                        RooArgSet(m_D0pi, mu3p, sigj3p, gamj3p, dj3p, m0, rho));

    RooRealVar mug_p("mug_p", "mug_p", 2010, 2009, 2011);
    RooRealVar sg_p("sg_p", "sg_p", 0.8, 0.01, 5);
    RooAbsPdf * gausp = new RooGenericPdf("gausp", "1/(@2*TMath::Sqrt(2*TMath::Pi()))*TMath::Exp(-0.5*((@0-@1)/@2)*((@0-@1)/@2))",
                        RooArgSet(m_D0pi, mug_p, sg_p));

    // the RooAddPdf used is made so that it does a recursive addition, meaning that the fraction defined here are not the true fraction in the model. Later ther will be a 
    // par of the code dedicated to the output of the true fractions
    RooRealVar frac1p("frac1p", "frac1p", 0.4, 0.01, 0.9);
    RooRealVar frac2p("frac2p", "frac2p", 0.4, 0.01, 0.9);
    RooRealVar frac3p("frac3p", "frac3p", 0.4, 0.01, 0.9); 
    RooAddPdf jjp("jjp", "jjp", RooArgList(*j1p, *j2p, *j3p, *gausp), RooArgList(frac1p, frac2p, frac3p), true);

    // Model : 3 Johnson functions D*-
    RooRealVar mu1m("mu1m", "mu1m", 2010.0, 2009, 2011);
    RooRealVar sigj1m("sigj1m", "sigj1m", 0.8, 0, 1);
    RooRealVar dj1m("dj1m", "dj1m", 1.7, 0.01, 10);
    RooRealVar gamj1m("gamj1m", "gamj1m",-0.6, -10, 10);
    RooAbsPdf * j1m = new RooGenericPdf("j1m", "(@0>@5)*pow((@0-@5),1+@6)*@4/(@2*TMath::Sqrt(2*TMath::Pi())*TMath::Sqrt(1+((@0-@1)*(@0-@1)/(@2*@2))))*TMath::Exp(-0.5*(@3+@4*TMath::ASinH((@0-@1)/@2))*(@3+@4*TMath::ASinH((@0-@1)/@2)))", 
                        RooArgSet(m_D0pi, mu1m, sigj1p, gamj1p, dj1p, m0, rho));

    RooRealVar mu2m("mu2m", "mu2m", 2010.0, 2009, 2011);
    RooRealVar sigj2m("sigj2m", "sugj2m", 0.8, 0, 1);
    RooRealVar dj2m("dj2m", "dj2m", 1.7, 0.001, 10);
    RooRealVar gamj2m("gamj2m", "gamj2m",-0.6, -10, 10);
    RooAbsPdf * j2m = new RooGenericPdf("j2m", "(@0>@5)*pow((@0-@5),1+@6)*@4/(@2*TMath::Sqrt(2*TMath::Pi())*TMath::Sqrt(1+((@0-@1)*(@0-@1)/(@2*@2))))*TMath::Exp(-0.5*(@3+@4*TMath::ASinH((@0-@1)/@2))*(@3+@4*TMath::ASinH((@0-@1)/@2)))", 
                        RooArgSet(m_D0pi, mu2m, sigj2p, gamj2m, dj2m, m0, rho));

    RooRealVar mu3m("mu3m", "mu3m", 2010.0, 2009, 2011);
    RooRealVar sigj3m("sigj3m", "sigj3m", 0.8, 0, 1);
    RooRealVar dj3m("dj3m", "dj3m", 1.7, 0.01, 10);
    RooRealVar gamj3m("gamj3m", "gamj3m",-0.9, -10, 10);
    RooAbsPdf * j3m = new RooGenericPdf("j3m", "(@0>@5)*pow((@0-@5),1+@6)*@4/(@2*TMath::Sqrt(2*TMath::Pi())*TMath::Sqrt(1+((@0-@1)*(@0-@1)/(@2*@2))))*TMath::Exp(-0.5*(@3+@4*TMath::ASinH((@0-@1)/@2))*(@3+@4*TMath::ASinH((@0-@1)/@2)))", 
                        RooArgSet(m_D0pi, mu3m, sigj3p, gamj3m, dj3p, m0, rho));

    RooRealVar mug_m("mug_m", "mug_m", 2010, 2009, 2011);
    RooRealVar sg_m("sg_m", "sg_m", 0.8, 0.01, 5);
    RooAbsPdf * gausm = new RooGenericPdf("gausm", "1/(@2*TMath::Sqrt(2*TMath::Pi()))*TMath::Exp(-0.5*((@0-@1)/@2)*((@0-@1)/@2))",
                        RooArgSet(m_D0pi, mug_m, sg_m));

    RooRealVar frac1m("frac1m", "frac1m", 0.4, 0.01, 0.9);
    RooRealVar frac2m("frac2m", "frac2m", 0.4, 0.01, 0.9);
    RooRealVar frac3m("frac3m", "frac3m", 0.4, 0.01, 0.9);
    RooAddPdf jjm("jjm", "jjm", RooArgList(*j1m, *j2m, *j3m, *gausm), RooArgList(frac1m, frac2m, frac3m), true);
    

    // Tagging function 
    RooGenericPdf * tag_Dsp = new RooGenericPdf("tag_Dsp", "tag_Dsp", "@0==1", RooArgSet(*q));
    RooGenericPdf * tag_Dsm = new RooGenericPdf("tag_Dsm", "tag_Dsm", "@0==-1", RooArgSet(*q));

    // Pdf for D*+ and D*-
    RooProdPdf * pdf_Dsp = new RooProdPdf("pdf_Dsp", "pdf_Dsp", RooArgSet(*tag_Dsp, jjp));
    RooProdPdf * pdf_Dsm = new RooProdPdf("pdf_Dsm", "pdf_Dsm", RooArgSet(*tag_Dsm, jjm)); // NB

    // A_CP fraction 
    /*
    RooFormulaVar * f_Dsp;
    RooFormulaVar * f_Dsm;
    */

    RooFormulaVar * f_Dsp = new RooFormulaVar("f_Dsp", "f_Dsp", "0.5*(1+@0)", RooArgSet(*A_CP_blind));
    RooFormulaVar * f_Dsm = new RooFormulaVar("f_Dsm", "f_Dsm", "0.5*(1-@0)", RooArgSet(*A_CP_blind));

    // Signal PDF
    RooAddPdf * pdf_sig = new RooAddPdf("pdf_sig", "pdf_sig", RooArgSet(*pdf_Dsp, *pdf_Dsm), RooArgSet(*f_Dsp, *f_Dsm));

    // Background : (x-x0)^b + (x-x0)^c 
    /*
    RooRealVar bp("bp", "bp", 0.55,0,1);
    RooRealVar cp("cp", "cp", 0.02, -0.1, 0.1);
    RooGenericPdf * bkg_Dsp = new RooGenericPdf("bkg_Dsp", "(m_D0pi>m0)*(TMath::Power((m_D0pi-m0), bp)*TMath::Exp(-cp*(m_D0pi-m0)))", RooArgList(m_D0pi, m0, bp, cp));
    
    RooRealVar bm("bm", "bm", 0.55,0,1);
    RooRealVar cm("cm", "cm", 0.02, -0.1, 0.1);
    RooGenericPdf * bkg_Dsm = new RooGenericPdf("bkg_Dsm", "(m_D0pi>m0)*(TMath::Power((m_D0pi-m0), bm)*TMath::Exp(-cm*(m_D0pi-m0)))", RooArgList(m_D0pi, m0, bm, cm));
    */
    RooRealVar ap("ap", "ap", 0, -5, 5);
    RooRealVar bp("bp", "bp", 0, -5, 5);
    RooGenericPdf * bkg_Dsp = new RooGenericPdf("bkg_Dsp", "TMath::Power(@0-@1, 0.5) + @2*TMath::Power(@0-@1, 3/2) + @3*TMath::Power(@0-@1, 5/2)", 
                            RooArgList(m_D0pi, m0, ap, bp));

    RooRealVar am("am", "am", 0, -5, 5);
    RooRealVar bm("bm", "bm", 0, -5, 5);
    RooGenericPdf * bkg_Dsm = new RooGenericPdf("bkg_Dm", "TMath::Power(@0-@1, 0.5) + @2*TMath::Power(@0-@1, 3/2) + @3*TMath::Power(@0-@1, 5/2)", 
                            RooArgList(m_D0pi, m0, am, bm)); // NB

    /*
    RooGenericPdf * bkgtagpdf = new RooGenericPdf("bkgtagpdf", "bkgtagpdf", "1+@0*@1/abs(@1)", RooArgSet(*A_CP_bkg, *q));
    RooProdPdf * pdf_bkg = new RooProdPdf("pdf_bkg", "pdf_bkg", RooArgSet(*bkgtagpdf, *bkg_Dsp));
    */

    RooGenericPdf * tag_bkg_Dsp = new RooGenericPdf("tag_bkg_Dsp", "tag_bkg_Dsp", "@0==1", RooArgSet(*q));
    RooGenericPdf * tag_bkg_Dsm = new RooGenericPdf("tag_bkg_Dsm", "tag_bkg_Dsm", "@0==-1", RooArgSet(*q));

    RooProdPdf * pdf_bkg_Dsp = new RooProdPdf("pdf_bkg_Dsp", "pdf_bkg_Dsp", RooArgList(*tag_bkg_Dsp, *bkg_Dsp));
    RooProdPdf * pdf_bkg_Dsm = new RooProdPdf("pdf_bkg_Dsm", "pdf_bkg_Dsm", RooArgList(*tag_bkg_Dsm, *bkg_Dsm));

    RooFormulaVar * f_bkg_Dsp = new RooFormulaVar("f_bkg_Dsp", "f_bkg_Dsp", "0.5*(1+@0)", RooArgSet(*A_CP_bkg));
    RooFormulaVar * f_bkg_Dsm = new RooFormulaVar("f_bkg_Dsm", "f_bkg_Dsm", "0.5*(1-@0)", RooArgSet(*A_CP_bkg));

    RooAddPdf * pdf_bkg = new RooAddPdf("pdf_bkg", "pdf_bkg", RooArgList(*pdf_bkg_Dsp, *pdf_bkg_Dsm), RooArgList(*f_bkg_Dsp, *f_bkg_Dsm));

    // Total model 
    RooRealVar * Nsig = new RooRealVar("Nsig", "Nsig", 1e5, 0, 1e9);
    RooRealVar * Nbkg = new RooRealVar("Nbkg", "Nbkg", 1e5, 0, 1e9);

    RooAddPdf * pdf_tot = new RooAddPdf("pdf_tot", "pdf_tot", RooArgSet(*pdf_sig, *pdf_bkg), RooArgSet(*Nsig, *Nbkg));


    // Fit
    RooArgSet * params = pdf_tot->getParameters(*obs);
    //params->readFromFile("run3_param.txt");
    //params->readFromFile("pisa_param/pisa_test_mu.txt");
    //params->readFromFile("DACP_Run3_params/test_1.txt");

    bool extended = true;
    //RooAbsTestStatistics::Configuration cfg;
    //cfg.nCPU(16);
    RooChi2Var * chi2 = new RooChi2Var("chi2", "chi2", *pdf_tot, *hist, extended, RooAbsData::SumW2);
    RooMinimizer m1(*chi2);
    m1.setVerbose(kTRUE);
    m1.migrad();
    m1.migrad();
    m1.hesse();

    RooFit::OwningPtr<RooFitResult> result = m1.save();
    result->Print("v");


    const TMatrixDSym &cor = result->correlationMatrix();
    cout << "Correlation Matrix: " << endl;
    cor.Print();

    //params->writeToFile("pisa_param/DACP_params.txt");
    params->writeToFile("/Users/gianlucapovero/Documents/Prove_Codeing/SCNSP/test1.txt");    
    //correlation matrix 
    gStyle->SetOptStat(0);
    TH2* hcorr = result->correlationHist();

    Int_t nbin = h_dm_plus->GetNbinsX() + h_dm_minus->GetNbinsX(); // dovrebbe esere (N1 + N2)/2?
    Int_t npar;
    RooArgSet * paramSet = pdf_tot->getParameters(*obs);
    RooArgList paramList(*paramSet);
    RooArgList * _floatParamList = (RooArgList*) paramList.selectByAttrib("Constant", kFALSE);
    npar = _floatParamList->getSize();
    double chi2_plot=0;
    int ndof = nbin - npar;

    cout << "mu1p: " << setprecision(10) << mu1p << "+-" << mu1p.getError() << endl;
    cout << "mu1m: " << setprecision(10) << mu1m << "+-" << mu1m.getError() << endl;
    cout << "mu2p: " << setprecision(10) << mu2p << "+-" << mu2p.getError() << endl;
    cout << "mu2m: " << setprecision(10) << mu2m << "+-" << mu2m.getError() << endl;
    cout << "mu3p: " << setprecision(10) << mu3p << "+-" << mu3p.getError() << endl;
    cout << "mu3m: " << setprecision(10) << mu3m << "+-" << mu3m.getError() << endl;
    cout << "mug_p: " << setprecision(10) << mug_p << "+-" << mug_p.getError() << endl;
    cout << "mug_m: " << setprecision(10) << mug_m << "+-" << mug_m.getError() << endl;
    
    auto mu_arr = RooArgList(mu1m, mu1p, mu2m, mu2p, mu3m, mu3p, mug_m, mug_p);
    const char* mu_a[8] = {"mu1m", "mu1p", "mu2m", "mu2p", "mu3m", "mu3p", "mug_m", "mug_p"};

    for(int i = 0; i <= 6; i=i+2){
        RooFormulaVar a("a", "@0 - @1", RooArgList(mu_arr[i], mu_arr[i+1]));
        cout << "Compatibilità (" << mu_a[i] << "-" << mu_a[i+1] << "): " << a.evaluate() << "+-" << a.getPropagatedError(*result) << endl;
    }
    

    //RooFormulaVar a("a", "@0 - @1", RooArgList(mu1m, mu1p));
    //cout << "Compatibilità (" << mu1m << "-" << mu1p << "): " << a.evaluate() << "+-" << a.getPropagatedError(*result) << endl;

    /*
    auto f2m = frac2m*(1.0 - frac1m);
    auto f3m = frac3m*(1.0 - frac1m)*(1.0 - frac2m);
    auto f2p = frac2p*(1.0 - frac1p);
    auto f3p = frac3p*(1.0 - frac1p)*(1.0 - frac2p);
    */

    //cout << "frac1m: " << frac1m << endl;
    //cout << "frac2m: " << f2m << endl;
    //cout << "frac3m: " << f3m << endl;
    //cout << "frac1p: " << frac1p << endl;
    //cout << "frac2p: " << f2p << endl;
    //cout << "frac3p: " << f3p << endl;
    
    // --- Plotting

    // ----- 
    TCanvas * roofit_c = new TCanvas("roofit_c", "roofit_c", 800, 800);
    roofit_c->cd();
    RooPlot * plot_mass = m_D0pi.frame(Title("Combined signs"));
    hist->plotOn(plot_mass);
    pdf_tot->plotOn(plot_mass, Components("*pdf_bkg"), LineColor(kBlack), LineStyle(kDashed));
    pdf_tot->plotOn(plot_mass, Components("*pdf_sig"), LineColor(kRed));
    //pdf_tot->plotOn(plot_mass, Components("*jjm"), LineColor(kOrange));
    pdf_tot->plotOn(plot_mass, LineColor(kBlue));
    plot_mass->GetXaxis()->SetTitle("#it{m}(#it{D}^{0}#it{#pi}) [MeV/c^{2}]");
    plot_mass->GetYaxis()->SetTitleOffset(1.65);

    //cout << "Chi^2/ndof: " << plot_mass->chiSquare() << endl;
    //cout << "Chi^2: " << chi2->getVal() << endl;

    RooHist * hpull = plot_mass->pullHist();
    RooPlot * pulls = m_D0pi.frame(Title("Pull Distribution"));
    pulls->addPlotable(hpull, "BX");

    TPad * upPad = new TPad("upPad", "upPad", 0.005, 0.2525, 0.995, 0.995);
    TPad * lowPad = new TPad("lowPad", "lowPad", 0.005, 0.005, 0.995, 0.2475);
    lowPad->Draw();
    upPad->Draw();
    upPad->cd();
    upPad->SetLeftMargin(0.25);
    upPad->SetTopMargin(0.1);
    lowPad->SetLeftMargin(0.25);
    chi2_plot = chi2->getVal()/ndof;
    plot_mass->Draw();

    TPaveText * box_chi2 = new TPaveText(0.7, 0.59, 0.9, 0.7, "NDCBR");
    box_chi2->SetFillColor(kWhite);
    box_chi2->SetTextSize(0.035);
    box_chi2->AddText("#chi^{2}/ndof =");
    box_chi2->AddText(Form("%3.2f/%d = %3.3f", chi2->getVal(), ndof, chi2_plot));
    //box_chi2->AddText(plot_mass_plus->chiSquare());
    box_chi2->Draw("same");

    TLegend * leg = new TLegend(0.6, 0.3, 0.89, 0.59);
    //leg->SetTextFont(132);
    //leg->SetTextSize(0.06);
    //leg->SetBorderSize(0);
    leg->AddEntry(plot_mass->getObject(3), "Total PDF", "l");
    leg->AddEntry(plot_mass->getObject(2), "Signal", "l");
    leg->AddEntry(plot_mass->getObject(1), "Background", "l");
    leg->AddEntry(plot_mass->getObject(0), "data", "lep");
    leg->Draw("same");

    lowPad->cd();
    pulls->SetTitle("");
    pulls->GetXaxis()->SetLabelSize(0);
    pulls->GetXaxis()->SetTitle("");
    pulls->GetYaxis()->SetTitle("Pulls");
    pulls->GetYaxis()->SetTitleSize(0.20);
    pulls->GetYaxis()->SetTitleOffset(0.25);
    pulls->GetYaxis()->SetLabelSize(0.15);
    pulls->GetYaxis()->SetNdivisions(503);
    pulls->GetYaxis()->SetRangeUser(-5, 5);
    pulls->Draw("B");

    // log scale 
    TCanvas * roofit_c_log = new TCanvas("roofit_c_log", "roofit_c_log", 800, 800);
    roofit_c_log->cd();
    
    TPad * upPad_log = new TPad("upPad_log", "upPad_log", 0.005, 0.2525, 0.995, 0.995);
    TPad * lowPad_log = new TPad("lowPad_log", "lowPad_log", 0.005, 0.005, 0.995, 0.2475);
    upPad_log->SetLogy();
    lowPad_log->Draw();
    upPad_log->Draw();
    upPad_log->cd();
    upPad_log->SetLeftMargin(0.25);
    upPad_log->SetTopMargin(0.1);
    lowPad_log->SetLeftMargin(0.25);
    chi2_plot = chi2->getVal()/ndof;
    plot_mass->Draw();

    lowPad_log->cd();
    pulls->SetTitle("");
    pulls->GetXaxis()->SetLabelSize(0);
    pulls->GetXaxis()->SetTitle("");
    pulls->GetYaxis()->SetTitle("Pulls");
    pulls->GetYaxis()->SetTitleSize(0.20);
    pulls->GetYaxis()->SetTitleOffset(0.25);
    pulls->GetYaxis()->SetLabelSize(0.15);
    pulls->GetYaxis()->SetNdivisions(503);
    pulls->GetYaxis()->SetRangeUser(-5, 5);
    pulls->Draw("B");


    // ---- Only D*+ events 
    TCanvas * roofit_plus = new TCanvas("roofit_plus", "roofit_plus", 800, 800);
    roofit_plus->cd();
    RooPlot * plot_mass_plus = m_D0pi.frame(Title("#it{D}^{*+}"));
    hist->plotOn(plot_mass_plus, Cut("q==q::Dsp"));
    pdf_tot->plotOn(plot_mass_plus, Components("pdf_bkg"), Slice(*q, "Dsp"), ProjWData(*q, *hist), LineColor(kBlack), LineStyle(kDashed));
    pdf_tot->plotOn(plot_mass_plus, Components("pdf_sig"), Slice(*q, "Dsp"), ProjWData(*q, *hist), LineColor(kRed));
    pdf_tot->plotOn(plot_mass_plus, Slice(*q, "Dsp"), ProjWData(*q, *hist), LineColor(kBlue));
    plot_mass_plus->GetXaxis()->SetTitle("#it{m}(#it{D}^{0}#it{#pi}) [MeV/c^{2}]");
    plot_mass_plus->GetYaxis()->SetTitleOffset(1.65);

    TText *txt_plus = new TText(2, 100, "#it{D}^{*+}");
    txt_plus->SetTextSize(0.03);
    plot_mass_plus->addObject(txt_plus);

    //cout << "Chi^2/ndof (D*+): " << plot_mass_plus->chiSquare() << endl;

    RooHist * hpull_plus = plot_mass_plus->pullHist();
    RooPlot * pulls_plus = m_D0pi.frame(Title("PullDistribution"));
    pulls_plus->addPlotable(hpull_plus, "BX");

    TPad * upPad_plus = new TPad("upPad_plus", "upPad_plus", 0.005, 0.2525, 0.995, 0.995);
    TPad * lowPad_plus = new TPad("lowPad_plus", "lowPad_plus", 0.005, 0.005, 0.995, 0.2475);
    lowPad_plus->Draw();
    upPad_plus->Draw();
    upPad_plus->cd();
    upPad_plus->SetLeftMargin(0.25);
    upPad_plus->SetTopMargin(0.1);
    lowPad_plus->SetLeftMargin(0.25);
    plot_mass_plus->Draw();

    TLegend * leg_plus = new TLegend(0.6, 0.3, 0.89, 0.59);
    //leg_plus->SetTextFont(132);
    //leg_plus->SetTextSize(0.06);
    //leg_plus->SetBorderSize(0);
    leg_plus->AddEntry(plot_mass_plus->getObject(3), "Total PDF", "l");
    leg_plus->AddEntry(plot_mass_plus->getObject(2), "Signal", "l");
    leg_plus->AddEntry(plot_mass_plus->getObject(1), "Background", "l");
    leg_plus->AddEntry(plot_mass_plus->getObject(0), "data", "lep");
    leg_plus->Draw("same");

    lowPad_plus->cd();
    pulls_plus->SetTitle("");
    pulls_plus->GetXaxis()->SetLabelSize(0);
    pulls_plus->GetXaxis()->SetTitle("");
    pulls_plus->GetYaxis()->SetTitle("Pulls");
    pulls_plus->GetYaxis()->SetTitleSize(0.20);
    pulls_plus->GetYaxis()->SetTitleOffset(0.25);
    pulls_plus->GetYaxis()->SetLabelSize(0.15);
    pulls_plus->GetYaxis()->SetNdivisions(503);
    pulls_plus->GetYaxis()->SetRangeUser(-5, 5);
    pulls_plus->Draw("B");

    // log scale 
    TCanvas * roofit_plus_log = new TCanvas("roofit_plus_log", "roofit_plus_log", 800, 800);
    roofit_plus_log->cd();

    TPad * upPad_plus_log = new TPad("upPad_plus_log", "upPad_plus_log", 0.005, 0.2525, 0.995, 0.995);
    TPad * lowPad_plus_log = new TPad("lowPad_plus_log", "lowPad_plus_log", 0.005, 0.005, 0.995, 0.2475);
    upPad_plus_log->SetLogy();
    lowPad_plus_log->Draw();
    upPad_plus_log->Draw();
    upPad_plus_log->cd();
    upPad_plus_log->SetLeftMargin(0.25);
    upPad_plus_log->SetTopMargin(0.1);
    lowPad_plus_log->SetLeftMargin(0.25);
    plot_mass_plus->Draw();

    lowPad_plus_log->cd();
    pulls_plus->SetTitle("");
    pulls_plus->GetXaxis()->SetLabelSize(0);
    pulls_plus->GetXaxis()->SetTitle("");
    pulls_plus->GetYaxis()->SetTitle("Pulls");
    pulls_plus->GetYaxis()->SetTitleSize(0.20);
    pulls_plus->GetYaxis()->SetTitleOffset(0.25);
    pulls_plus->GetYaxis()->SetLabelSize(0.15);
    pulls_plus->GetYaxis()->SetNdivisions(503);
    pulls_plus->GetYaxis()->SetRangeUser(-5, 5);
    pulls_plus->Draw("B");

    // ------- Only D*- events
    TCanvas * roofit_minus = new TCanvas("roofit_minus", "roofit_minus", 800, 800);
    roofit_minus->cd();
    RooPlot * plot_mass_minus = m_D0pi.frame(Title("#it{D}^{*-}"));
    hist->plotOn(plot_mass_minus, Cut("q==q::Dsm"));
    pdf_tot->plotOn(plot_mass_minus, Components("pdf_bkg"), Slice(*q, "Dsm"), ProjWData(*q, *hist), LineColor(kBlack), LineStyle(kDashed));
    pdf_tot->plotOn(plot_mass_minus, Components("pdf_sig"), Slice(*q, "Dsm"), ProjWData(*q, *hist), LineColor(kRed));
    pdf_tot->plotOn(plot_mass_minus, Slice(*q, "Dsm"), ProjWData(*q, *hist), LineColor(kBlue));
    plot_mass_minus->GetXaxis()->SetTitle("#it{m}(#it{D}^{0}#it{#pi}) [MeV/c^{2}]");
    plot_mass_minus->GetYaxis()->SetTitleOffset(1.65);

    //cout << "Chi^2/ndof: " << plot_mass_minus->chiSquare() << endl;

    RooHist * hpull_minus = plot_mass_minus->pullHist();
    RooPlot * pulls_minus = m_D0pi.frame(Title("Pull Distribution"));
    pulls_minus->addPlotable(hpull_minus, "BX");

    TPad * upPad_minus = new TPad("upPad_minus", "upPad_minus", 0.005, 0.2525, 0.995, 0.995);
    TPad * lowPad_minus = new TPad("lowPad_minus", "lowPad_minus", 0.005, 0.005, 0.995, 0.2475);
    lowPad_minus->Draw();
    upPad_minus->Draw();
    upPad_minus->cd();
    upPad_minus->SetLeftMargin(0.25);
    upPad_minus->SetTopMargin(0.1);
    lowPad_minus->SetLeftMargin(0.25);
    plot_mass_minus->Draw();

    TLegend * leg_minus = new TLegend(0.6, 0.3, 0.89, 0.59);
    //leg_plus->SetTextFont(132);
    //leg_plus->SetTextSize(0.06);
    //leg_plus->SetBorderSize(0);
    leg_minus->AddEntry(plot_mass_minus->getObject(3), "Total PDF", "l");
    leg_minus->AddEntry(plot_mass_minus->getObject(2), "Signal", "l");
    leg_minus->AddEntry(plot_mass_minus->getObject(1), "Background", "l");
    leg_minus->AddEntry(plot_mass_minus->getObject(0), "data", "lep");
    leg_minus->Draw("same");

    lowPad_minus->cd();
    pulls_minus->SetTitle("");
    pulls_minus->GetXaxis()->SetLabelSize(0);
    pulls_minus->GetXaxis()->SetTitle("");
    pulls_minus->GetYaxis()->SetTitle("Pulls");
    pulls_minus->GetYaxis()->SetTitleSize(0.20);
    pulls_minus->GetYaxis()->SetTitleOffset(0.25);
    pulls_minus->GetYaxis()->SetLabelSize(0.15);
    pulls_minus->GetYaxis()->SetNdivisions(503);
    pulls_minus->GetYaxis()->SetRangeUser(-5, 5);
    pulls_minus->Draw("B");

    // log scale 
    TCanvas * roofit_minus_log = new TCanvas("roofit_minus_log", "roofit_minus_log", 800, 800);
    roofit_minus_log->cd();

    TPad * upPad_minus_log = new TPad("upPad_minus_log", "upPad_minus_log", 0.005, 0.2525, 0.995, 0.995);
    TPad * lowPad_minus_log = new TPad("lowPad_minus_log", "lowPad_minus_log", 0.005, 0.005, 0.995, 0.2475);
    upPad_minus_log->SetLogy();
    lowPad_minus_log->Draw();
    upPad_minus_log->Draw();
    upPad_minus_log->cd();
    upPad_minus_log->SetLeftMargin(0.25);
    upPad_minus_log->SetTopMargin(0.1);
    lowPad_minus_log->SetLeftMargin(0.25);
    plot_mass_minus->Draw();

    lowPad_minus_log->cd();
    pulls_minus->SetTitle("");
    pulls_minus->GetXaxis()->SetLabelSize(0);
    pulls_minus->GetXaxis()->SetTitle("");
    pulls_minus->GetYaxis()->SetTitle("Pulls");
    pulls_minus->GetYaxis()->SetTitleSize(0.20);
    pulls_minus->GetYaxis()->SetTitleOffset(0.25);
    pulls_minus->GetYaxis()->SetLabelSize(0.15);
    pulls_minus->GetYaxis()->SetNdivisions(503);
    pulls_minus->GetYaxis()->SetRangeUser(-5, 5);
    pulls_minus->Draw("B");

    // Covarinace matrix 
    
    TCanvas * corr_c = new TCanvas("corr_c", "corr_c");
    gPad->SetRightMargin(0.15);
    hcorr->GetXaxis()->SetLabelSize(0.03);
    hcorr->GetYaxis()->SetLabelSize(0.03);
    hcorr->Draw("colz");
    
   
    //roofit_c->SaveAs("plot_DACP3/moreplot/test_4_m.root");
    //roofit_c_log->SaveAs("plot_DACP3/mass_plot_log_lst.pdf");
    //roofit_plus->SaveAs("plot_DACP3/moreplot/test_4_plus.root");
    //roofit_plus_log->SaveAs("plot_DACP3/moreplot/test_1_plus_log.pdf");
    //roofit_minus->SaveAs("plot_DACP3/moreplot/test_4_minus.root");
    //roofit_minus_log->SaveAs("plot_DACP3/mass_minus_log_lst.pdf");
   

}
