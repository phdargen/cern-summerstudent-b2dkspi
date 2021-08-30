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
//void fitKs_MC_LL()
{   
    //Binned or unbinned fit?
    bool binned=false;
    
	///Load file
	TFile* file= new TFile("/afs/cern.ch/work/m/mbuhayeu/public/cern-summerstudent-b2dkspi-master/Data/MC_B2DKspi_DD_12.root");
	TTree* tree = (TTree*) file->Get("DecayTree");
    
    //Disable all but needed branches
    tree->SetBranchStatus("*",0);  
    tree->SetBranchStatus("Ks_MM",1);  

    //Define branches
	RooRealVar Ks_MM("Ks_MM", "m(DK_{s}#pi)", 465., 530.,"MeV");

    //Create RooDataSet
	RooArgList list =  RooArgList(Ks_MM);
	RooDataSet* data = new RooDataSet("data", "data", tree, list,"");

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
	
	//CB1
	RooRealVar mean1("mean1", "mean1", 507.,480.,510.); 
	RooRealVar sigma1("sigma1", "sigma1", 4,1.,10.);
	RooRealVar alpha1("alpha1","alpha1",-2,-80.,0.);
	RooRealVar n1("n1","n1",5,0.,100.);
	RooCBShape CB1("CB1","CB for signal",Ks_MM, mean, sigma1, alpha1, n1);



	//Y: There's no bkg in MC file. Everything is the "signal"
	// BKG
	//Polynomial
	//RooRealVar coeff("p1","coeff", 1, 0., 10.);
	//RooPolynomial Poly1("Poly1","Poly signal",Ks_MM,RooArgList(coeff),1);
	
	
    ///Total pdf
	///----------------------
	RooRealVar n_sig("n_sig", "n_sig", data->numEntries()/2., 0., data->numEntries());
	RooRealVar n_bkg("n_bkg", "n_bkg", data->numEntries()/2. , 0., data->numEntries());

	RooAbsPdf* pdf=new RooAddPdf("pdf", "pdf", RooArgList(CB, CB1), RooArgList(n_sig,n_bkg));

	//RooAddPdf pdf("pdf","pdf",RooArgList(*sig),RooArgList(n_sig1));
	//RooAddPdf pdf("pdf","pdf",RooArgList(CB,CB2),RooArgList(n_sig1,n_sig2));
	
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
	pdf->plotOn(frame_m,Components(CB),LineColor(3), LineStyle(1), Name("sig1"));
	pdf->plotOn(frame_m,Components(CB1),LineColor(5),LineStyle(1), Name("sig2"));

	cout<<"chi2 = "<<frame_m->chiSquare("FullModel","Data")<<endl;

	RooHist* hpull = frame_m->pullHist("Data","FullModel");
	RooPlot* frame3 = Ks_MM.frame(Title("Pull Distribution")) ;
	frame3->addPlotable(hpull,"P") ;	

	TCanvas* c1= new TCanvas("");
	data->plotOn(frame_m);
	pdf->plotOn(frame_m );

    frame_m->Draw();
	//c1->Print("fitKs_MC_DD.eps");
	

	auto leg = new TLegend(0.7, 0.7, 0.9, 0.9);
	leg->AddEntry(frame_m->findObject("FullModel"), "Full Model", "L");
	leg->AddEntry(frame_m->findObject("sig1"),"sig1","L");
	leg->AddEntry(frame_m->findObject("sig2"),"sig2","L");
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
	//c->cd(2) ; gPad->SetLeftMargin(0.15) ; frame2->GetYaxis()->SetTitleOffset(1.6) ; frame2->Draw() ;
	c->cd(2) ; gPad->SetTopMargin(0); gPad->SetPad(0.03, 0.02, 0.95, 0.27); gPad->SetRightMargin(0.03); gPad->SetBottomMargin(0.2); frame3->Draw();

	c->SaveAs("/afs/cern.ch/work/m/mbuhayeu/public/cern-summerstudent-b2dkspi-master/MassFit/fitKs_MC_DD.png");	
	

    return 0;
}

