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
#include <fstream>
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
#include "RooFormulaVar.h"
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

int main(int argc, char** argv)
//void fit_LL()
{   
    //Binned or unbinned fit?
    bool binned=false;
    
	///Load file
	
    	TChain* tree = new TChain("DecayTree");
		tree->Add("/afs/cern.ch/work/m/mbuhayeu/public/cern-summerstudent-b2dkspi-master/Selection/BDTG/B2DKspi_data.root");

 
    //Disable all but needed branches
    tree->SetBranchStatus("*",0);  
    tree->SetBranchStatus("Ks_MM",1);  
    tree->SetBranchStatus("BDTG",1);
    tree->SetBranchStatus("KsCat",1);  

    //Define branches
	RooRealVar Ks_MM("Ks_MM", "", 468., 525.,"MeV");

    //Create RooDataSet
	//RooArgList list =  RooArgList(Ks_MM);
	RooRealVar BDTG("BDTG", "BDTG", 0.);
	RooCategory KsCat("KsCat","KsCat") ;
	KsCat.defineType("DD",1);
	RooArgList list =  RooArgList(Ks_MM,BDTG,KsCat);
	RooDataSet*  data = new RooDataSet("data","data",tree,list, " (KsCat == 1) && (BDTG > 0.6422)" );
	//RooDataSet* data = new RooDataSet("data", "data", tree, list,"");

	RooDataSet* data_small = (RooDataSet*) data->reduce(SelectVars(RooArgSet(Ks_MM)));
	RooDataHist* data_binned = data_small->binnedClone();

	
    ///Define fit model
	///----------------

	///Signal model: a single CB function
	///-----------------------
	RooRealVar mean("mean", "mean", 497.67,468.,525.);
	RooRealVar sigma("sigma", "sigma", 3.8379,0.2,6.);
	RooRealVar alpha("alpha","alpha",1.7009,0.,10.);
	RooRealVar n("n","n",1.9225,0.,10.);
	RooCBShape CB("CB","CB for signal",Ks_MM, mean, sigma, alpha, n);

	RooFormulaVar meanb("meanb","meanb","mean",mean);
	RooRealVar sigmab("sigmab","sigmab",5.0,0.,5.);
	RooRealVar alpha2("alpha2","alpha2",0.78303);
	RooRealVar n2("n2","n2",9.7360);
	RooCBShape CB2("CB2","CB2 for signal",Ks_MM,meanb,sigmab,alpha2,n);
        
	RooRealVar f("f","f",0.59567);

	RooAbsPdf* sig = new RooAddPdf("sig","sig",RooArgList(CB,CB2), f);

	///Background model
    ///-------------------------
	RooRealVar mean3("mean3", "mean", 501.18,480.,510.);
	RooRealVar sigma3("sigma3", "sigma", 5,1.,8.);
	RooRealVar alpha3("alpha3","alpha",0.96009,0.,80.);
	RooRealVar n3("n3","n",9.9984,0.,10.);
	RooCBShape CB3("CB3","CB for signal",Ks_MM, mean3, sigma3, alpha3, n3);
	//CB1
	RooRealVar mean4("mean4", "mean4", 497.23,480.,510.); 
	RooRealVar sigma4("sigma4", "sigma4", 2.2839,1.,5.);
	RooRealVar alpha4("alpha4","alpha4",1.0350,0.,8.);
	RooRealVar n4("n4","n4",9.9998,0.,10.);
	RooCBShape CB4("CB4","CB for signal",Ks_MM, mean4, sigma4, alpha4, n4);

    ///Total pdf
	///----------------------
	RooRealVar n_sig("n_sig", "n_sig", data->numEntries()/2., 0., data->numEntries());
	RooRealVar n_bkg("n_bkg", "n_bkg", data->numEntries()/2., 0., data->numEntries());
	RooRealVar n_sig1("n_sig1", "n_sig1", data->numEntries()/2., 0., data->numEntries());
	RooRealVar n_sig2("n_sig2", "n_sig2", data->numEntries()/2., 0., data->numEntries());
	
	//pdf list
	RooAbsPdf* sig1 = new RooAddPdf("sig1", "sig1", RooArgList(CB, CB2), RooArgList(n_sig, n_bkg));
	RooAbsPdf* sig2 = new RooAddPdf("sig2", "sig2", RooArgList(CB3, CB4), RooArgList(n_sig1, n_sig2));

	RooArgList pdf0;
	pdf0.add(*sig);
	pdf0.add(*sig1);
	pdf0.add(*sig2);

	//N list
	RooArgList num;
	num.add(n_sig);
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

	RooPlot* frame_m= Ks_MM.frame();
	
	data->plotOn(frame_m,Name("Data"));
	pdf.plotOn(frame_m,Name("FullModel"));
	pdf.plotOn(frame_m,Components(*sig1),LineColor(3), LineStyle(1), Name("sig1"));
	pdf.plotOn(frame_m,Components(*sig2),LineColor(9), LineStyle(1), Name("sig2"));
	pdf.plotOn(frame_m,Components(*sig),LineColor(11), LineStyle(1), Name("sig"));
	//pdf.plotOn(frame_m,Components(*BsDst),LineColor(5), LineStyle(1), Name("BsDst"));
	//pdf.plotOn(frame_m,Components(*Bsbkg),LineColor(7), LineStyle(1), Name("Bsbkg"));
	//pdf.plotOn(frame_m,Components(CB),LineColor(3), LineStyle(1), Name("sig"));
	//pdf.plotOn(frame_m,Components(cheby),LineColor(11), LineStyle(1), Name("bkg"));
	//pdf.plotOn(frame_m,Components(*sigBs),LineColor(7), LineStyle(1), Name("Bsbkg"));


	auto leg = new TLegend(0.7, 0.7, 0.9, 0.9);
	//leg.SetFillColor(0);
	leg->AddEntry(frame_m->findObject("FullModel"), "Full Model", "f");
	leg->AddEntry(frame_m->findObject("sig"), "sig", "f");
	leg->AddEntry(frame_m->findObject("sig1"), "bkg1", "f");
	leg->AddEntry(frame_m->findObject("sig2"), "bkg2", "f");
	//leg.AddEntry(frame_m->findObject("bkg"), "bkg", "f");
	//leg.AddEntry(frame_m->findObject("Bsbkg"), "Bsbkg", "f");

	cout<<"chi2 = "<<frame_m->chiSquare("FullModel","Data")<<endl;

	RooHist* hresid = frame_m->residHist("Data","FullModel");

	RooHist* hpull = frame_m->pullHist("Data","FullModel");

	RooPlot* frame2 = Ks_MM.frame(Title("Residual Distribution")) ;
	frame2->addPlotable(hresid,"P") ;

	RooPlot* frame3 = Ks_MM.frame(Title("Pull Distribution")) ;
	frame3->addPlotable(hpull,"P") ;

	TCanvas* c = new TCanvas("fit_chi2residpull","Ks Mass Fit",800,600) ;
	c->Divide(1,2) ;
	c->cd(1) ; gPad->SetBorderMode(1); gPad->SetTopMargin(1); gPad->SetBottomMargin(0.15) ; gPad->SetRightMargin(0.03); gPad->SetPad(0.03, 0.27, 0.95, 0.95); frame_m->GetXaxis()->SetTitle("m(DK_{s}#pi) [MeV]");frame_m->Draw() ; leg -> DrawClone();
	//c->cd(2) ; gPad->SetLeftMargin(0.15) ; frame2->GetYaxis()->SetTitleOffset(1.6) ; frame2->Draw() ;
	c->cd(2) ; gPad->SetTopMargin(0); gPad->SetPad(0.03, 0.02, 0.95, 0.27); gPad->SetRightMargin(0.03); gPad->SetBottomMargin(0.2); frame3->Draw();
	 
	c->SaveAs("/afs/cern.ch/work/m/mbuhayeu/public/cern-summerstudent-b2dkspi-master/MassFit/fit_DD.png");

	Ks_MM.setRange("signal", 497.67 - 50 , 497.67 + 50) ;
	double bkg1yield = sig1 -> createIntegral(Ks_MM,NormSet(Ks_MM),Range("signal")) ->getVal();	
	double bkg2yield = sig2 -> createIntegral(Ks_MM,NormSet(Ks_MM),Range("signal")) ->getVal();
	double sigyield = sig -> createIntegral(Ks_MM,NormSet(Ks_MM),Range("signal")) ->getVal();

	ofstream myfile;
	myfile.open("/afs/cern.ch/work/m/mbuhayeu/public/cern-summerstudent-b2dkspi-master/MassFit/KsDD.txt");
	//myfile<<"the BdDst yield: "<<BdDstyield*n_BdDst.getVal()<<"\n";
	//myfile<<"the BsDst yield: "<<BsDstyield*n_BsDst.getVal()<<"\n";
	//myfile<<"the Bs yield: "<<Bsyield*n_Bs.getVal()<<"\n";
	myfile<<"the bkg yield: "<<(bkg1yield+bkg2yield)*n_sig.getVal()<<"<- Bkg"<<"\n";
	myfile<<"sig yield: "<<sigyield*n_sig.getVal()<<"<- Sig"<<"\n";

	myfile.close();
	frame_m->Draw();
	leg-> DrawClone();

    return 0;
  
}
