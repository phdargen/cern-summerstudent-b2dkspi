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
#include <TLegend.h>
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

int main(int argc, char** argv)
{   
    //Binned or unbinned fit?
    bool binned=false;
    
	///Load file
	TFile* file= new TFile("/afs/cern.ch/work/m/mbuhayeu/public/cern-summerstudent-b2dkspi-master/Selection/BDTG/B2DKspi_data.root");
	TTree* tree = (TTree*) file->Get("DecayTree");
    
    //Disable all but needed branches
    tree->SetBranchStatus("*",0);  
    tree->SetBranchStatus("Ks_MM",1);  
    tree->SetBranchStatus("BDTG",1);
    tree->SetBranchStatus("KsCat",1);  
    //Define branches
	RooRealVar Ks_MM("Ks_MM", "m(DK_{s}#pi)", 468., 526.,"MeV");

    //Create RooDataSet
	//RooArgList list =  RooArgList(Ks_MM);
	RooRealVar BDTG("BDTG", "BDTG", 0.);
	RooCategory KsCat("KsCat","KsCat") ;
	KsCat.defineType("LL",0);
	RooArgList list =  RooArgList(Ks_MM,BDTG,KsCat);
	RooDataSet*  data = new RooDataSet("data","data",tree,list, " KsCat == 0 && BDTG > -0.8512" );

	RooDataSet* data_small = (RooDataSet*) data->reduce(SelectVars(RooArgSet(Ks_MM)));
	RooDataHist* data_binned = data_small->binnedClone();

	
    ///Define fit model
	///----------------

	///CB
	RooRealVar mean("mean", "mean", 497.8,465.,530.);
	RooRealVar sigma("sigma", "sigma", 0.5,0.2,6.);
	RooRealVar alpha("alpha","alpha",7.0559,0.,10.);
	RooRealVar n("n","n",1.10,0.,10.);
	RooCBShape CB("CB","CB for signal",Ks_MM, mean, sigma, alpha, n);

    //Exponential
    //RooRealVar exp_par("exp_par", "exp_par", -0.01,-1.,0.); 	
    //RooExponential exp("exp","exp",Ks_MM,exp_par); 

	//Polynomial
	RooRealVar coeff("p1","coeff", 1, 0., 10.);
	RooPolynomial Poly1("Poly1","Poly signal",Ks_MM,RooArgList(coeff),1);

    ///Total pdf
	///----------------------
	RooRealVar n_sig("n_sig", "n_sig", data->numEntries()/2., 0., data->numEntries());
	RooRealVar n_bkg("n_bkg", "n_bkg", data->numEntries()/2. , 0., data->numEntries());

	RooAbsPdf* pdf=new RooAddPdf("pdf", "pdf", RooArgList(CB, Poly1), RooArgList(n_sig,n_bkg));

	RooFitResult *result;
	
    //binned Fit
    if(binned) result = pdf->fitTo(*data_binned,Save(kTRUE),Extended());
	
    //unbinned Fit
    else result = pdf->fitTo(*data,Save(kTRUE),Extended());
	
    cout << "result is --------------- "<<endl;
	result->Print(); 
 
	///Plot 
	///----------

	RooPlot* frame_m= Ks_MM.frame();
	
	data->plotOn(frame_m,Name("Data"));
	pdf->plotOn(frame_m,Name("FullModel"));

	cout<<"chi2 = "<<frame_m->chiSquare("FullModel","Data")<<endl;

	RooHist* hpull = frame_m->pullHist("Data","FullModel");
	RooPlot* frame3 = Ks_MM.frame(Title("Pull Distribution")) ;
	frame3->addPlotable(hpull,"P") ;	

	TCanvas* c1= new TCanvas("");
	
//	RooPlot* frame_m= Ks_MM.frame();

	data->plotOn(frame_m);
	pdf->plotOn(frame_m );

    frame_m->Draw();
	c1->Print("fitKs_MC_LL.eps");
	

	auto leg = new TLegend(0.7, 0.7, 0.9, 0.9);
	leg->AddEntry(frame_m->findObject("FullModel"), "Full Model", "L");
	leg->DrawClone();

	TCanvas* c = new TCanvas("fit_chi2residpull","Ks Mass Fit",800,600) ;
	c->Divide(1,2) ;
	c->cd(1) ; 
	gPad->SetBorderMode(1); 
	gPad->SetTopMargin(1); 
	gPad->SetBottomMargin(0.15); 
	gPad->SetRightMargin(0.03); 
	gPad->SetPad(0.03, 0.27, 0.95, 0.95); 
	frame_m->GetXaxis()->SetTitle("m(DK_{s}#pi) [MeV]");frame_m->Draw() ; leg->DrawClone();
	c->cd(2) ; gPad->SetTopMargin(0); gPad->SetPad(0.03, 0.02, 0.95, 0.27); gPad->SetRightMargin(0.03); gPad->SetBottomMargin(0.2); frame3->Draw();

	c->SaveAs("/afs/cern.ch/work/m/mbuhayeu/public/cern-summerstudent-b2dkspi-master/MassFit/fit_LL.png");	
	

    return 0;
}


