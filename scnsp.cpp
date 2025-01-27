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
#include "RooMinuit.h"
#include "RooMinimizer.h"
#include "RooChi2Var.h"

using namespace RooFit;

void scnsp(){

    // Observables
    // Here we define the observable we are going to study with the code. In this case the observable is the mass of the D^0 + pion 
    // x1 and x2 are the two extremes of the energy interval we are looking in the analysis
    Double_t x1 = 2004.41, x2 = 2020.0;
    RooRealVar m_D0pi("m_D0pi", "m_D0pi", x1, x2, "Mev/c^2");
    //the bin are set as 520 so that each bin has an excact length of 0.03 MeV/c^2
    m_D0pi.setBins(520);

    // This line is used to define an object that we will use as a tagging variable, since the porpouse of the code is to do a simultaneous fit of both D^*+ and D^*- 
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
    // theese lines define the tagging functions using the tagging variable "q". 
    RooGenericPdf * tag_Dsp = new RooGenericPdf("tag_Dsp", "tag_Dsp", "@0==1", RooArgSet(*q));
    RooGenericPdf * tag_Dsm = new RooGenericPdf("tag_Dsm", "tag_Dsm", "@0==-1", RooArgSet(*q));

    // Pdf for D*+ and D*-
    // one defined the tagging functions, the two signal functions for plus and minus events are defined by multipying each tag function with each model pdf.
    RooProdPdf * pdf_Dsp = new RooProdPdf("pdf_Dsp", "pdf_Dsp", RooArgSet(*tag_Dsp, jjp));
    RooProdPdf * pdf_Dsm = new RooProdPdf("pdf_Dsm", "pdf_Dsm", RooArgSet(*tag_Dsm, jjm));

    // A_CP fraction
    // in order to define the whole signal model, plus and minus events must be added together but with the particular fraction defined in theese lines
    RooFormulaVar * f_Dsp = new RooFormulaVar("f_Dsp", "f_Dsp", "0.5*(1+@0)", RooArgSet(*A_CP_blind));
    RooFormulaVar * f_Dsm = new RooFormulaVar("f_Dsm", "f_Dsm", "0.5*(1-@0)", RooArgSet(*A_CP_blind));

    // Signal PDF
    // the final signal pdf is the summ of both plus and minus events 
    RooAddPdf * pdf_sig = new RooAddPdf("pdf_sig", "pdf_sig", RooArgSet(*pdf_Dsp, *pdf_Dsm), RooArgSet(*f_Dsp, *f_Dsm));

    // Background : (x-x0)^(1/2) + a*(x-x0)^(3/2) + c*(x-x0)^(5/2) 
    // the following lines defines the background model fot both D^*+ and D^*- events
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

    // as for the signal events, the tagging function is defined also for the background 
    RooGenericPdf * tag_bkg_Dsp = new RooGenericPdf("tag_bkg_Dsp", "tag_bkg_Dsp", "@0==1", RooArgSet(*q));
    RooGenericPdf * tag_bkg_Dsm = new RooGenericPdf("tag_bkg_Dsm", "tag_bkg_Dsm", "@0==-1", RooArgSet(*q));

    // theese lines define the pdf for the background plus and minus events
    RooProdPdf * pdf_bkg_Dsp = new RooProdPdf("pdf_bkg_Dsp", "pdf_bkg_Dsp", RooArgList(*tag_bkg_Dsp, *bkg_Dsp));
    RooProdPdf * pdf_bkg_Dsm = new RooProdPdf("pdf_bkg_Dsm", "pdf_bkg_Dsm", RooArgList(*tag_bkg_Dsm, *bkg_Dsm));

    // here the fraction for plus and minus events is defined
    RooFormulaVar * f_bkg_Dsp = new RooFormulaVar("f_bkg_Dsp", "f_bkg_Dsp", "0.5*(1+@0)", RooArgSet(*A_CP_bkg));
    RooFormulaVar * f_bkg_Dsm = new RooFormulaVar("f_bkg_Dsm", "f_bkg_Dsm", "0.5*(1-@0)", RooArgSet(*A_CP_bkg));

    // this line defines the final background functon 
    RooAddPdf * pdf_bkg = new RooAddPdf("pdf_bkg", "pdf_bkg", RooArgList(*pdf_bkg_Dsp, *pdf_bkg_Dsm), RooArgList(*f_bkg_Dsp, *f_bkg_Dsm));

    // Total model 
    // onec the signal and the background pdf are written, the total model in defined as the sum of them with the expected number of events used as fractions.
    // in theese lines the number of events is defined with an interval between 0 and 1e9 since we expect Nsig of the order of 1e6
    RooRealVar * Nsig = new RooRealVar("Nsig", "Nsig", 1e6, 0, 1e9);
    RooRealVar * Nbkg = new RooRealVar("Nbkg", "Nbkg", 1e6, 0, 1e9);
    
    RooAddPdf * pdf_tot = new RooAddPdf("pdf_tot", "pdf_tot", RooArgSet(*pdf_sig, *pdf_bkg), RooArgSet(*Nsig, *Nbkg));


    // ------ Fit ------

    // this lines calls 
    // this line calls for the creation of a set which contains the parameters of the total pdf
    RooArgSet * params = pdf_tot->getParameters(*obs);

    // This line tells our fit to read the values of its parameters from a specific file.txt
    // Useful when one has to run the program multiple times 
    //params->readFromFile("..path/to/file.txt..");

    // in theese following lines the RooChi2Var object is crated giving the total pdf and the histogram of data
    bool extended = true;
    RooChi2Var * chi2 = new RooChi2Var("chi2", "chi2", *pdf_tot, *hist, extended, RooAbsData::SumW2);

    // this lines construct the MINUIT interface on the RooChi2Var object
    RooMinimizer m1(*chi2);
    //this line gives instruction to the program to print each computation on the logbook, this will result in a very rich logbook (and also very long)
    m1.setVerbose(kTRUE);
    // theese lines apply the minimum gradient method, computing the asteroid vaues of the parameters, and the hessian matrix method, computing the parameters errors, upon 
    // the RooMinimizer object
    m1.migrad();
    m1.migrad();
    m1.hesse();

    // theese lines defines an object RooFitResult extracted from the RooMinimizer and then print it. The "v" options tells the program to print also 
    // the error matrix and the correlation
    RooFit::OwningPtr<RooFitResult> result = m1.save();
    result->Print("v");

    // in theese lines a matrix object is defined and is filled with the correlation matrix given by the fit, it is then printed 
    const TMatrixDSym &cor = result->correlationMatrix();
    cout << "Correlation Matrix: " << endl;
    cor.Print();

    // this line tells the code to print the RooùargSet containing the fit parameter on a file.txt
    params->writeToFile("..path/to/file.txt..");    
    //correlation matrix 
    // this line create a 2 dimentsion histogram which is then filled with the correlation matrix of the fit exprted from the result
    TH2* hcorr = result->correlationHist();

    // the following lines are used to compute the nuber of degree of freedom
    Int_t nbin = h_dm_plus->GetNbinsX() + h_dm_minus->GetNbinsX();
    Int_t npar;
    RooArgSet * paramSet = pdf_tot->getParameters(*obs);
    RooArgList paramList(*paramSet);
    // this line select the argument that respect a certain condition from an existing list of arguments
    RooArgList * _floatParamList = (RooArgList*) paramList.selectByAttrib("Constant", kFALSE);
    npar = _floatParamList->getSize();
    // this line will be used later
    double chi2_plot=0;
    int ndof = nbin - npar;

    // the following lines print in output the means parameters with a specific precision and theis errors, useful beacus the result printed 
    // from the fit in line 235 don't have  high enough precision 
    cout << "mu1p: " << setprecision(10) << mu1p << "+-" << mu1p.getError() << endl;
    cout << "mu1m: " << setprecision(10) << mu1m << "+-" << mu1m.getError() << endl;
    cout << "mu2p: " << setprecision(10) << mu2p << "+-" << mu2p.getError() << endl;
    cout << "mu2m: " << setprecision(10) << mu2m << "+-" << mu2m.getError() << endl;
    cout << "mu3p: " << setprecision(10) << mu3p << "+-" << mu3p.getError() << endl;
    cout << "mu3m: " << setprecision(10) << mu3m << "+-" << mu3m.getError() << endl;
    cout << "mug_p: " << setprecision(10) << mug_p << "+-" << mug_p.getError() << endl;
    cout << "mug_m: " << setprecision(10) << mug_m << "+-" << mug_m.getError() << endl;

    // the following lines defines two arrays, the first contains the values of the RooVar of the means while the second contains a list of char
    auto mu_arr = RooArgList(mu1m, mu1p, mu2m, mu2p, mu3m, mu3p, mug_m, mug_p);
    const char* mu_a[8] = {"mu1m", "mu1p", "mu2m", "mu2p", "mu3m", "mu3p", "mug_m", "mug_p"};

    // this for loop defines a RooFormulaVar that subtract two elements of the first array and then prints in aoutput also the propagated error of the operation
    for(int i = 0; i <= 6; i=i+2){
        RooFormulaVar a("a", "@0 - @1", RooArgList(mu_arr[i], mu_arr[i+1]));
        cout << "Compatibilità (" << mu_a[i] << "-" << mu_a[i+1] << "): " << a.evaluate() << "+-" << a.getPropagatedError(*result) << endl;
    }

    // the folloquinf lines explicitate the recursive operstion  made in the RooAddPdf in order to compute the real fractions of the pdfs
    auto f2m = frac2m*(1.0 - frac1m);
    auto f3m = frac3m*(1.0 - frac1m)*(1.0 - frac2m);
    auto f2p = frac2p*(1.0 - frac1p);
    auto f3p = frac3p*(1.0 - frac1p)*(1.0 - frac2p);
    
    // theese lines print in output the values of the real fractions of the pdf, useful to get an immediate idea of how the pdfs are distributed
    cout << "frac1m: " << frac1m << endl;
    cout << "frac2m: " << f2m << endl;
    cout << "frac3m: " << f3m << endl;
    cout << "frac1p: " << frac1p << endl;
    cout << "frac2p: " << f2p << endl;
    cout << "frac3p: " << f3p << endl;

    
    // ------ Plotting ------
    // in this section all pdfs will be plotted, meaning that we will have a plot for the combined sign distribution, which is the one dostribution it will be made a fit of, 
    // a plot for the D^*+ distribution and for the D^*- distribution. Each of theese will have its own pull distibution plot drawn underneath. Moreover, for each of theese distributions
    // an additional log scale plot will be drawn since some effects at the extremes tìof the interval are noticible only with such a scale. Finally the covariant matrix 2 dimesional
    // histogram will be drawn.
    // The structure used for drawing a plot is the same for each of the distribution mentioned above, therefore only the first one will be discussed in ditailes 

    // this line calls for the cration of a new TCanvas object, it will be the base the distribution will be drawn upon. The dimension is given now is of little importance 
    // since onece printed it can be made bigger, here "800, 800" is given 
    TCanvas * roofit_c = new TCanvas("roofit_c", "roofit_c", 800, 800);
    // this lines tells our program to change directory to the freshly created TCanvas, meaning that everithing created from now will be automatically drawn on the canvas, until
    // it will be said otherwise
    roofit_c->cd();
    // this line calls for the creation of a RooPlot object called frame, it is a container object in which will put all plottable pdfs and data. It will than be drawn 
    // as a whole on the TCanvas
    RooPlot * plot_mass = m_D0pi.frame(Title("Combined signs"));
    // the following lines give our code the order to put the data (hist) and the pdf, both components and the total, inside the frame 
    // note that the order in which the components are called is importatn for the pullHist() function, see line 319
    hist->plotOn(plot_mass);
    pdf_tot->plotOn(plot_mass, Components("*pdf_bkg"), LineColor(kBlack), LineStyle(kDashed));
    pdf_tot->plotOn(plot_mass, Components("*pdf_sig"), LineColor(kRed));
    //pdf_tot->plotOn(plot_mass, Components("*jjm"), LineColor(kOrange));
    pdf_tot->plotOn(plot_mass, LineColor(kBlue));
    // the following line acces the axis of the frame and set the title as we desire 
    plot_mass->GetXaxis()->SetTitle("#it{m}(#it{D}^{0}#it{#pi}) [MeV/c^{2}]");
    plot_mass->GetYaxis()->SetTitleOffset(1.65);

    //cout << "Chi^2/ndof: " << plot_mass->chiSquare() << endl;
    //cout << "Chi^2: " << chi2->getVal() << endl;
    // this lines call for the cration of the RooHist object that contains the pull distribution of the pdf, since it is not specified which instogram and which pdf to use 
    // in order to compute the pull distribution, the code will automatically use the last defined histogram and the last defined pdf
    RooHist * hpull = plot_mass->pullHist();
    // this line calls for the creation of the frame that will contain the pull distribution
    RooPlot * pulls = m_D0pi.frame(Title("Pull Distribution"));
    // this line adds the RooHist object to the frame object
    pulls->addPlotable(hpull, "BX");

    // in order to have a cleaner output, the pull distribution is drawn under its respective pdf, meaning that in one TCanvaas there has to be two frames
    // to specify to each frame how much space to take inside the TCanvas two TPad object are created. theese will create sections on the TCanvas 
    // inside which the frames are said to be drawn 
    TPad * upPad = new TPad("upPad", "upPad", 0.005, 0.2525, 0.995, 0.995);
    TPad * lowPad = new TPad("lowPad", "lowPad", 0.005, 0.005, 0.995, 0.2475);
    lowPad->Draw();
    upPad->Draw();
    // this line call for the code to change dorectory to the freshly created upPad 
    upPad->cd();
    upPad->SetLeftMargin(0.25);
    upPad->SetTopMargin(0.1);
    lowPad->SetLeftMargin(0.25);
    
    // this line compute the reduced chi^2 "by hand" using the RooChi2Var object value
    chi2_plot = chi2->getVal()/ndof;
    // this line calls for the draw of the frame on the canvas 
    plot_mass->Draw();

    // the following lines define TPaveText which is something that is drawn on top of the TCanvas indipendent from the frame objects. It will contain the value of the reduced chi^2 
    // as well as the value of the chi^2 and the degrees of freedom, so that one wille have all the important informations by looking at the plot 
    TPaveText * box_chi2 = new TPaveText(0.7, 0.59, 0.9, 0.7, "NDCBR");
    box_chi2->SetFillColor(kWhite);
    box_chi2->SetTextSize(0.035);
    box_chi2->AddText("#chi^{2}/ndof =");
    box_chi2->AddText(Form("%3.2f/%d = %3.3f", chi2->getVal(), ndof, chi2_plot));
    //box_chi2->AddText(plot_mass_plus->chiSquare());
    box_chi2->Draw("same");

    // the following lines define a TLegend object containing the information for the pdfs and data drawn on the TCanvas
    TLegend * leg = new TLegend(0.6, 0.3, 0.89, 0.59);
    leg->AddEntry(plot_mass->getObject(3), "Total PDF", "l");
    leg->AddEntry(plot_mass->getObject(2), "Signal", "l");
    leg->AddEntry(plot_mass->getObject(1), "Background", "l");
    leg->AddEntry(plot_mass->getObject(0), "data", "lep");
    leg->Draw("same");

    // this line calls for the change directory in the lowPad 
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
    // the following lines defines the same TCanvas and plot, whitout the TPave text and with the logaritmic scale, see line 381
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
    

    // the following lines call for the save of the plots as files. If one possesses the extention, it is preferable to save only the linear scale plots as file.root
    // so that one can change the scale of the y axis later
    roofit_c->SaveAs("..path/to/file.root");
    //roofit_c_log->SaveAs("plot_DACP3/mass_plot_log_lst.pdf");
    roofit_plus->SaveAs("..path/to/file.root.root");
    //roofit_plus_log->SaveAs("plot_DACP3/moreplot/test_1_plus_log.pdf");
    roofit_minus->SaveAs("..path/to/file.root.root");
    //roofit_minus_log->SaveAs("plot_DACP3/mass_minus_log_lst.pdf");
   

}
