//RooFit example
//Philippe d'Argent

#include <cmath>
#include <iostream>
#include <vector>
#include "TSystem.h"
#include <TChain.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TF1.h>
#include "TDatime.h"
#include "TFile.h"
#include <RooDataSet.h>
#include "RooGaussModel.h"
#include "RooExponential.h"
#include "RooChebychev.h"
#include "RooPolynomial.h"
#include "RooAddModel.h"
#include "RooTruthModel.h"
#include "RooFitResult.h"
#include "RooDecay.h"
#include "RooPlot.h"
#include "RooGaussian.h"
#include "RooDstD0BG.h"
#include "RooAddPdf.h"
#include "RooExtendPdf.h"
#include "RooDataHist.h"
#include "RooCBShape.h"
#include "RooCategory.h"
#include "RooFormulaVar.h"
#include "RooSimultaneous.h"
#include "RooHist.h"
#include "RooStats/SPlot.h"
#include "RooTreeDataStore.h"

#include "TH1D.h"
#include "TH2D.h"

#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooConstVar.h"

#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif

using namespace RooFit ;
using namespace std;

//int main(int argc, char** argv)
void fitB_2Gauss_Cheby()
{   
    //Binned or unbinned fit?
    bool binned=false;
    
	///Load file
	TFile* file= new TFile("../Selection/output.root ");
	TTree* tree = (TTree*) file->Get("DecayTree");
    
    //Disable all but needed branches
    tree->SetBranchStatus("*",0);  
    tree->SetBranchStatus("B_DTF_MM",1);  

    //Define branches
	RooRealVar B_DTF_MM("B_DTF_MM", "m(DK_{s}#pi)", 5000., 5800.,"MeV");

    //Create RooDataSet
	RooArgList list =  RooArgList(B_DTF_MM);
	RooDataSet* data = new RooDataSet("data", "data", tree, list,"");

	RooDataSet* data_small = (RooDataSet*) data->reduce(SelectVars(RooArgSet(B_DTF_MM)));
	RooDataHist* data_binned = data_small->binnedClone();

	
    ///Define fit model
	///----------------

	///Signal model
	///-----------------------

    ///Gaussian
	RooRealVar mean("\\mu", "mean", 5280.,5150.,5350.); 
	RooRealVar sigma1("\\sigma_{1}", "sigma1", 20.,10.,30.);	

	RooGaussian Gauss1("Gauss1", "Gauss1", B_DTF_MM, mean, sigma1);

	RooRealVar mean2("mean2", "mean2", 5080.,4950.,5200.);
	RooRealVar sigma2("sigma2", "sigma2", 50.,30.,70.);
    
	RooGaussian Gauss2("Gauss2", "Gauss2", B_DTF_MM, mean2, sigma2);
	//Add a second Gaussian 
    //...
    
	///Background model
	///-------------------------

    //Chebychev
    RooRealVar xa1 = RooRealVar("xa1","xa1",0,-1.,1.);
    RooRealVar xa2 = RooRealVar("xa2","xa2",0,-0.5,0.5);
    RooRealVar xa3 = RooRealVar("xa3","xa3",0,-0.5,0.5);
    RooRealVar xa4 = RooRealVar("xa4","xa4",0,-0.5,0.5);
    RooChebychev cheby("cheby", "cheby", B_DTF_MM, RooArgList(xa1, xa2, xa3, xa4));
	
    ///Total pdf
	///----------------------
	RooRealVar n_sig1("n_sig1", "n_sig1", data->numEntries()/2., 0., data->numEntries());
	RooRealVar n_bkg("n_bkg", "n_bkg", data->numEntries()/2. , 0., data->numEntries());
	RooRealVar n_sig2("n_sig2", "n_sig2", data->numEntries()/4. ,0., data->numEntries()/2);

	//pdf list
	RooAbsPdf* bkg = new RooAddPdf("bkg", "bkg", cheby, n_bkg);
	RooAbsPdf* sig1 = new RooAddPdf("sig1", "sig1", Gauss1, n_sig1);
	RooAbsPdf* sig2 = new RooAddPdf("sig2", "sig2", Gauss2, n_sig2);

	RooArgList pdf0;
	pdf0.add(*bkg);
	pdf0.add(*sig1);
	pdf0.add(*sig2);

	//N list
	RooArgList num;
	num.add(n_bkg);
	num.add(n_sig1);
	num.add(n_sig2);
	
	RooAddPdf pdf("pdf","pdf",pdf0,num);
	///Do an extended Fit
    ///--------
	RooFitResult *result;
	
    //binned Fit
    if(binned) result = pdf.fitTo(*data_binned,Save(kTRUE),Extended());
	
    //unbinned Fit
    else result = pdf.fitTo(*data,Save(kTRUE),Extended());
	
    cout << "result is --------------- "<<endl;
	result->Print(); 
 
	///Plot 
	///----------
	TCanvas* c1= new TCanvas("");
	c1->SetGrid();

	RooPlot* frame_m= B_DTF_MM.frame();
	
	data->plotOn(frame_m,Name("Data"));
	pdf.plotOn(frame_m,Name("FullModel"));
	pdf.plotOn(frame_m,Components(*bkg),LineColor(1), LineStyle(1), Name("bkg"));
	pdf.plotOn(frame_m,Components(*sig1),LineColor(3), LineStyle(1), Name("sig1"));
	pdf.plotOn(frame_m,Components(*sig2),LineColor(5), LineStyle(1), Name("sig2"));

	TLegend leg(0.7, 0.7, 0.9, 0.9);
	leg.AddEntry(frame_m->findObject("FullModel"), "Full Model", "L");
	leg.AddEntry(frame_m->findObject("bkg"), "bkg", "L");
	leg.AddEntry(frame_m->findObject("sig1"), "sig1", "L");
	leg.AddEntry(frame_m->findObject("sig2"), "sig2", "L");
    
	cout<<"chi2 = "<<frame_m->chiSquare("FullModel","Data")<<endl;
	frame_m->Draw();
	leg.DrawClone();
	c1->Print("fitB_2Gauss_Cheby.eps");
	
    //return 0;
  
}


