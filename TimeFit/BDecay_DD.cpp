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
#include "RooBDecay.h"
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
void BDecay_DD()
{   
    //Binned or unbinned fit?
    bool binned=false;
    
	///Load file
	TFile* file= new TFile("/afs/cern.ch/work/y/yuxiao/public/newtree.root");
	TTree* tree = (TTree*) file->Get("DecayTree");
    
    //Disable all but needed branches
    tree->SetBranchStatus("*",0);  
    tree->SetBranchStatus("B_DTF_TAU",1);  

    //Define branches
	RooRealVar B_DTF_TAU("B_DTF_TAU", "Time", 0.0001, 0.01,"time");

    //Create RooDataSet
	RooArgList list =  RooArgList(B_DTF_TAU);
	RooDataSet* data = new RooDataSet("data", "data", tree, list,"");

	RooDataSet* data_small = (RooDataSet*) data->reduce(SelectVars(RooArgSet(B_DTF_TAU)));
	RooDataHist* data_binned = data_small->binnedClone();

	
    ///Define fit model
	///----------------
	RooRealVar tau("tau","tau",1.519e-03/*,1.515e-03,1.623e-03*/);//ns;1.519+-0.04
	RooRealVar m("m","m",1.1403e-06);//ns
	RooRealVar nsig1("nsig1","nsig1",1.2682e+03,1.2682e+03-3*1.44e+02,1.2682e+03+3*1.44e+02);
	RooRealVar nsig2("nsig2","nsig2",1.1019e+03,1.1019e+03-3*1.44e+02,1.1019e+03+3*1.44e+02);
	RooFormulaVar f("f","f","nsig1/(nsig1+nsig2)",RooArgList(nsig1,nsig2));
	RooRealVar s1("s1","s1",3.0273e-02/*,3.0273e-02-3*2.17e-03,3.0273e-02+3*2.17e-03*/);//e-3ns
	RooRealVar s2("s2","s2",7.1769e-02/*,7.1769e-02-3*3.74e-03,7.1769e-02+3*3.74e-03*/);//e-3ns
	//RooRealVar s("s","s",0.0007);
	RooFormulaVar s("s","s","sqrt(f*pow(s1,2)+(1-f)*pow(s2,2))/1000",RooArgList(f,s1,s2));//ns
	RooRealVar t("t","t",0.03,0.,100.);
	RooResolutionModel* model = new RooGaussModel("model","model",B_DTF_TAU, m, s, t);
	RooRealVar dgamma("dgamma","dgamma",0.001/tau.getVal()/*,-0.009/tau.getVal(),0.011/tau.getVal()*/);//ns^-1
	RooRealVar dm("dm","dm",0.5065e+03/*,0.5065e+03-0.0019e+03,0.5065e+03+0.0019e+03*/);//ns^-1
	RooRealVar f0("f0","f0",1);
	RooRealVar f1("f1","f1",0);
	RooRealVar f2("f2","f2",0);
	RooRealVar f3("f3","f3",0);
	RooBDecay Decay("Decay","Decay",B_DTF_TAU,tau,dgamma,f0,f1,f2,f3,dm,*model,RooBDecay::SingleSided);

	RooRealVar tau2("tau2","tau2",1.,0.,3.);
	RooRealVar m2("m2","m2",5.,-5.,100.);
	//RooRealVar s2("s2","s2",10.,0.,200.);
	RooRealVar t2("t2","t2",0.03,0.,1.);
	RooResolutionModel* model2 = new RooGaussModel("model2","model2",B_DTF_TAU, m2, s2, t2);
	RooDecay Decay2("Decay2","Decay2",B_DTF_TAU,tau2,*model2,RooDecay::SingleSided);


	RooRealVar n_sig1("n_sig1", "n_sig1", data->numEntries()/2,0.,data->numEntries()*2);
	RooRealVar n_sig2("n_sig2", "n_sig2", data->numEntries()/2,0.,data->numEntries());

	//RooRealVar f("f","f",0.5,0.,1.);
	//RooAbsPdf* sig = new RooAddPdf("sig","sig",Decay,n_sig1);

	//RooAddPdf pdf("pdf","pdf",RooArgList(Decay,Decay2),RooArgList(n_sig1,n_sig2));
	RooAddPdf pdf("pdf","pdf",Decay,n_sig1);
	
	RooFitResult *result;
	
    //binned Fit
    if(binned) result = pdf.fitTo(*data_binned,Save(kTRUE),Extended());
	
    //unbinned Fit
    else result = pdf.fitTo(*data,Save(kTRUE),Extended());
	
    cout << "result is --------------- "<<endl;
	result->Print(); 
 
	///Plot 
	///----------

	RooPlot* frame_m= B_DTF_TAU.frame();
	
	data->plotOn(frame_m,Name("Data"));
	pdf.plotOn(frame_m,Name("FullModel"));
	pdf.plotOn(frame_m,Components(Decay),LineColor(3), LineStyle(1), Name("sig"));
	//pdf.plotOn(frame_m,Components(Decay2),LineColor(5),LineStyle(1), Name("sig2"));

	TLegend leg(0.7, 0.7, 0.9, 0.9);
	leg.AddEntry(frame_m->findObject("FullModel"), "Full Model", "L");
	leg.AddEntry(frame_m->findObject("sig"),"sig","L");
	//leg.AddEntry(frame_m->findObject("sig2"),"sig2","L");
	frame_m->Draw();
	leg.DrawClone();

	cout<<"chi2 = "<<frame_m->chiSquare("FullModel","Data")<<endl;
/*
	RooHist* hpull = frame_m->pullHist("Data","FullModel");
	RooPlot* frame3 = B_DTF_TAU.frame(Title("Pull Distribution")) ;
	frame3->addPlotable(hpull,"P") ;	

	TCanvas* c = new TCanvas("fit_chi2residpull","B0 Mass Fit",800,600) ;
	c->Divide(1,2) ;
	c->cd(1) ; gPad->SetBorderMode(1); gPad->SetTopMargin(1); gPad->SetBottomMargin(0.15) ; gPad->SetRightMargin(0.03); gPad->SetPad(0.03, 0.27, 0.95, 0.95); frame_m->GetXaxis()->SetTitle("m(DK_{s}#pi) [MeV]");frame_m->Draw() ; leg.DrawClone();
	//c->cd(2) ; gPad->SetLeftMargin(0.15) ; frame2->GetYaxis()->SetTitleOffset(1.6) ; frame2->Draw() ;
	c->cd(2) ; gPad->SetTopMargin(0); gPad->SetPad(0.03, 0.02, 0.95, 0.27); gPad->SetRightMargin(0.03); gPad->SetBottomMargin(0.2); frame3->Draw();

	c->SaveAs("/afs/cern.ch/work/y/yuxiao/public/FitResults/MC_DD.png");	
*/	
    //return 0;
  
}


