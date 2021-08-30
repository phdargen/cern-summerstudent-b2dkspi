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
	RooRealVar BDTG("BDTG", "BDTG", 0.);
	RooCategory KsCat("KsCat","KsCat") ;
	KsCat.defineType("LL",0);
	RooArgList list =  RooArgList(Ks_MM,BDTG,KsCat);
	RooDataSet*  data = new RooDataSet("data","data",tree,list, " (KsCat == 0) && (BDTG > -0.8512)" );

	RooDataSet* data_small = (RooDataSet*) data->reduce(SelectVars(RooArgSet(Ks_MM)));
	RooDataHist* data_binned = data_small->binnedClone();

	
    ///Define fit model
	///----------------

	///Signal model: a double CB
	///-----------------------
	RooRealVar mean("mean", "mean", 4.9814e+02); //4.9788e+02);
	RooRealVar sigma("sigma", "sigma", 2.7864e-01);//2.0860e+00);
	RooRealVar alpha("alpha","alpha",3.9923e-02);//2.7019e-01);
	RooRealVar n("n","n",9.9999e+00);//9.9958e+00);
	RooCBShape CB0("CB0","CB0 for signal",Ks_MM, mean, sigma, alpha, n);

    RooRealVar sigma1("sigma1", "sigma1", 3.6738);//3.6868);
    RooRealVar alpha1("alpha1","alpha1",-1.1752e+00);//-1.0718e+00);
    RooRealVar n1("n1","n1",1.1904e+01);//2.0778e+01);
    RooCBShape CB1("CB1","CB1 for signal",Ks_MM, mean, sigma1, alpha1, n1);

	RooRealVar f("f","f",1.7101e+02/(8.2600e+02+1.7101e+02));//1.8369e+02/(8.1322e+02+1.8369e+02));//n_sig/n_total

	RooAbsPdf* sig = new RooAddPdf("sig","sig",RooArgList(CB0, CB1), f);

	///Background model
    ///-------------------------

    RooRealVar mean6("mean6", "mean6", 4.9814e+02, 400.,500.); //4.9788e+02,400.,500.);
    RooCBShape CB6("CB6","CB6 for Bsbkg",Ks_MM, mean6, sigma, alpha, n);
    RooCBShape CB7("CB7","CB7 for Bsbkg",Ks_MM, mean6, sigma1, alpha1, n1);

    RooAbsPdf* Bsbkg = new RooAddPdf("Bsbkg","Bsbkg",RooArgList(CB6, CB7), f);
	
	//Bd part.reco.bkg
	//CB
	RooRealVar mean2("mean2", "mean2", 4.9759e+02); 
	RooRealVar sigma2("sigma2", "sigma2", 5.2479e+00);//5.2499e+00);
	RooRealVar alpha2("alpha2","alpha2",1.5091e+00);//1.5104e+00);
	RooRealVar n2("n2","n2",3.2497e+00);//3.2392e+00);
	RooCBShape CB2("CB2","CB2 for signal",Ks_MM, mean2, sigma2, alpha2, n2);
	
	RooRealVar sigma3("sigma3", "sigma3",2.2582e+00); //2.2569e+00);
	RooRealVar alpha3("alpha3","alpha3",-1.8920e+00);//-1.8925e+00);
	RooRealVar n3("n3","n3",7.7643e-01);//7.7517e-01);
	RooCBShape CB3("CB3","CB3 for signal",Ks_MM, mean2, sigma3, alpha3, n3);
	
	RooRealVar f2("f2","f2",1.1758e+03/(1.1758e+03+9.3864e+02));//1.1759e+03/(1.1759e+03+9.3722e+02));
	RooAbsPdf* BdDst = new RooAddPdf("BdDst","BdDst",RooArgList(CB2,CB3),f2);
	
	//Bs part.reco.bkg
	//CB
	RooRealVar mean4("mean4", "mean4",4.9779e+02);//4.9768e+02); 
	RooRealVar sigma4("sigma4", "sigma4", 2.9540);//5.2179);
	RooRealVar alpha4("alpha4","alpha4",8.1196e-01);//1.4733);
	RooRealVar n4("n4","n4", 6.8592);//4.1241);
	RooCBShape CB4("CB4","CB4 for signal",Ks_MM, mean4, sigma4, alpha4, n4);
	
	RooRealVar sigma5("sigma5", "sigma5", 3.3055);//2.5739);
	RooRealVar alpha5("alpha5","alpha5", -7.2261e-01);//-1.5394);
	RooRealVar n5("n5","n5", 9.9575);//1.3978);
	RooCBShape CB5("CB5","CB5 for signal",Ks_MM, mean4, sigma5, alpha5, n5);
	
	RooRealVar f4("f4","f4", 1.0862e+03/(1.0862e+03+1.0870e+03));//1.0869e+03/(1.0869e+03+1.0870e+03));
	RooAbsPdf* BsDst = new RooAddPdf("BsDst","BsDst",RooArgList(CB4,CB5),f4);
	
	//Polynomial
	RooRealVar coeff3("p3","coeff", 1, 0., 10.);
	RooPolynomial Poly3("Poly3","Poly signal",Ks_MM,RooArgList(coeff3),1);


    ///Total pdf
	///----------------------
	RooRealVar n_sig("n_sig", "n_sig", data->numEntries()/2., 0., data->numEntries());
	RooRealVar n_bkg("n_bkg", "n_bkg", data->numEntries()/2., 0., data->numEntries());
	RooRealVar n_BdDst("n_BdDst", "n_BdDst", data->numEntries()/2., 0., data->numEntries());
	RooRealVar n_BsDst("n_BsDst", "n_BsDst", data->numEntries()/2., 0., data->numEntries());	
	RooRealVar n_Bs("n_Bs", "n_Bs", data->numEntries()/2., 0., data->numEntries());
	
	//pdf list
	RooAbsPdf* bkg = new RooAddPdf("bkg", "bkg", Poly3, n_bkg);

	RooArgList pdf0;
	pdf0.add(*sig);
	pdf0.add(*BsDst);
	pdf0.add(*BdDst);
	pdf0.add(*bkg);
	pdf0.add(*Bsbkg);

	//N list
	RooArgList num;
	num.add(n_sig);
	num.add(n_BsDst);
	num.add(n_BdDst);
	num.add(n_bkg);
	num.add(n_Bs);

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
	pdf.plotOn(frame_m,Components(*sig),LineColor(3), LineStyle(1), Name("sig"));
	pdf.plotOn(frame_m,Components(*bkg),LineColor(11), LineStyle(1), Name("bkg"));
	pdf.plotOn(frame_m,Components(*BdDst),LineColor(9), LineStyle(1), Name("BdDst"));
	pdf.plotOn(frame_m,Components(*BsDst),LineColor(5), LineStyle(1), Name("BsDst"));
	pdf.plotOn(frame_m,Components(*Bsbkg),LineColor(7), LineStyle(1), Name("Bsbkg"));


	auto leg = new TLegend(0.7, 0.7, 0.9, 0.9);
	//leg.SetFillColor(0);
	leg->AddEntry(frame_m->findObject("FullModel"), "Full Model", "f");
	leg->AddEntry(frame_m->findObject("sig"), "sig", "f");
	leg->AddEntry(frame_m->findObject("BsDst"), "BsDst", "f");
	leg->AddEntry(frame_m->findObject("BdDst"), "BdDst", "f");
	leg->AddEntry(frame_m->findObject("bkg"), "bkg", "f");
	leg->AddEntry(frame_m->findObject("Bsbkg"), "Bsbkg", "f");

	cout<<"chi2 = "<<frame_m->chiSquare("FullModel","Data")<<endl;

	RooHist* hresid = frame_m->residHist("Data","FullModel");

	RooHist* hpull = frame_m->pullHist("Data","FullModel");

	RooPlot* frame2 = Ks_MM.frame(Title("Residual Distribution")) ;
	frame2->addPlotable(hresid,"P") ;

	RooPlot* frame3 = Ks_MM.frame(Title("Pull Distribution")) ;
	frame3->addPlotable(hpull,"P") ;

	TCanvas* c = new TCanvas("fit_chi2residpull","Ks Mass Fit",800,600) ;
	c->Divide(1,2) ;
	c->cd(1) ; gPad->SetBorderMode(1); gPad->SetTopMargin(1); gPad->SetBottomMargin(0.15) ; gPad->SetRightMargin(0.03); gPad->SetPad(0.03, 0.27, 0.95, 0.95); frame_m->GetXaxis()->SetTitle("m(DK_{s}#pi) [MeV]");frame_m->Draw() ; leg->DrawClone();
	c->cd(2) ; gPad->SetLeftMargin(0.15) ; frame2->GetYaxis()->SetTitleOffset(1.6) ; frame2->Draw() ;
	c->cd(2) ; gPad->SetTopMargin(0); gPad->SetPad(0.03, 0.02, 0.95, 0.27); gPad->SetRightMargin(0.03); gPad->SetBottomMargin(0.2); frame3->Draw();
	 
	c->SaveAs("/afs/cern.ch/work/m/mbuhayeu/public/cern-summerstudent-b2dkspi-master/MassFit/fit_LL.png");

	Ks_MM.setRange("signal", 497.90 - 50 , 497.90 + 50) ;
	double BdDstyield = BdDst -> createIntegral(Ks_MM,NormSet(Ks_MM),Range("signal")) ->getVal();	
	double BsDstyield = BsDst -> createIntegral(Ks_MM,NormSet(Ks_MM),Range("signal")) ->getVal();
	double Bsyield = Bsbkg -> createIntegral(Ks_MM,NormSet(Ks_MM),Range("signal")) ->getVal();
	double bkgyield = Poly3.createIntegral(Ks_MM,NormSet(Ks_MM),Range("signal")) ->getVal();
	double sigyield = sig -> createIntegral(Ks_MM,NormSet(Ks_MM),Range("signal")) ->getVal();

	ofstream myfile;
	myfile.open("/afs/cern.ch/work/m/mbuhayeu/public/cern-summerstudent-b2dkspi-master/MassFit/Ks_LL.txt");
	myfile<<"the BdDst yield: "<<BdDstyield*n_BdDst.getVal()<<"\n";
	myfile<<"the BsDst yield: "<<BsDstyield*n_BsDst.getVal()<<"\n";
	myfile<<"the Bs yield: "<<Bsyield*n_Bs.getVal()<<"\n";
	myfile<<"the 1-order bkg yield: "<<bkgyield*n_bkg.getVal()<<"\n";
	myfile<<"sig yield: "<<sigyield*n_sig.getVal()<<"\n";

	myfile.close();
	frame_m->Draw();
	leg-> DrawClone();

    return 0;
  
}

