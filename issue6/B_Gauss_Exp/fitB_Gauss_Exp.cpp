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
void fitB_Gauss_Exp()
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
	RooRealVar sigma1("\\sigma_{1}", "sigma1", 15.45,0.,30.);	

	RooGaussian Gauss1("Gauss1", "Gauss1", B_DTF_MM, mean, sigma1);

    //Add a second Gaussian 
    //...
    
	///Background model
	///-------------------------

    //Exponential
    RooRealVar exp_par("exp_par", "exp_par", -0.01,-1.,0.); 
    RooExponential exp("exp","exp",B_DTF_MM,exp_par); 
    
	
    ///Total pdf
	///----------------------
	RooRealVar n_sig("n_sig", "n_sig", data->numEntries()/2., 0., data->numEntries());
	RooRealVar n_bkg("n_bkg", "n_bkg", data->numEntries()/2. , 0., data->numEntries());

	RooAbsPdf* pdf=new RooAddPdf("pdf", "pdf", RooArgList(Gauss1, exp), RooArgList(n_sig,n_bkg));

	///Do an extended Fit
    ///--------
	RooFitResult *result;
	
    //binned Fit
    if(binned) result = pdf->fitTo(*data_binned,Save(kTRUE),Extended());
	
    //unbinned Fit
    else result = pdf->fitTo(*data,Save(kTRUE),Extended());
	
    cout << "result is --------------- "<<endl;
	result->Print(); 
	///Plot 
	///----------
	TCanvas* c1= new TCanvas("");
	
	RooPlot* frame_m= B_DTF_MM.frame();

	data->plotOn(frame_m,Name("Data"));
	pdf->plotOn(frame_m,Name("FullModel"));

    frame_m->Draw();
	c1->Print("plot.eps");
	cout << "chi2 = " << frame_m->chiSquare("FullModel","Data")<<endl;	
    //return 0;
  
}


