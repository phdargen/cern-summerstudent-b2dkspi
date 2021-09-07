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
#include "RooCBShape.h"
//#include "RooCrystalBall.h"
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
void Section6_1_DD()
{   
    //Binned or unbinned fit?
    bool binned=false;
    
	///Load file
	TFile* file= new TFile("/afs/cern.ch/work/y/yuxiao/public/newtree.root");
	TTree* tree = (TTree*) file->Get("DecayTree");

    //Disable all but needed branches
    tree->SetBranchStatus("*",0);  
    tree->SetBranchStatus("B_SUB",1);  

    //Define branches
	RooRealVar B_SUB("B_SUB","B_SUB",-0.2,0.2,"time");

	//Create RooDataSet
	RooArgList list =  RooArgList(B_SUB);
	RooDataSet* data = new RooDataSet("data", "data", tree, list, "");

	RooDataSet* data_small = (RooDataSet*) data->reduce(SelectVars(RooArgSet(B_SUB)));
	RooDataHist* data_binned = data_small->binnedClone();

	
    ///Define fit model
	///----------------
	RooRealVar mean("mean","mean",0.,-0.2,0.2);
	RooRealVar sigma1("sigma1","sigma1",0.1,0.,0.2);
	RooGaussian Gauss1("Gauss1","Gauss1",B_SUB, mean,sigma1);

	RooRealVar mean2("mean2","mean2",0.,-0.2,0.2);
	RooRealVar sigma2("sigma2","sigma2",0.1,0.,0.2);
	RooGaussian Gauss2("Gauss2","Gauss2",B_SUB, mean,sigma2);

	RooRealVar n_sig1("n_sig1", "n_sig1", data->numEntries()/2 ,0.,data->numEntries());
	RooRealVar n_sig2("n_sig2", "n_sig2", data->numEntries()/2 ,0.,data->numEntries());

	//RooAddPdf pdf("pdf","pdf",RooArgList(Gauss1),RooArgList(n_sig1));
	RooAddPdf pdf("pdf","pdf",RooArgList(Gauss1,Gauss2),RooArgList(n_sig1,n_sig2));
	
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

	RooPlot* frame_m= B_SUB.frame();

	data->plotOn(frame_m,Name("Data"));
	pdf.plotOn(frame_m,Name("FullModel"));
	pdf.plotOn(frame_m,Components(Gauss1),LineColor(3), LineStyle(1), Name("sig"));
	pdf.plotOn(frame_m,Components(Gauss2),LineColor(5),LineStyle(1), Name("sig2"));

	TLegend leg(0.7, 0.7, 0.9, 0.9);
	leg.AddEntry(frame_m->findObject("FullModel"), "Full Model", "L");
	leg.AddEntry(frame_m->findObject("sig"),"sig","L");
	leg.AddEntry(frame_m->findObject("sig2"),"sig2","L");


	cout<<"chi2 = "<<frame_m->chiSquare("FullModel","Data")<<endl;
        
	frame_m->Draw();
	leg.DrawClone();
	
    //return 0;
  
}


