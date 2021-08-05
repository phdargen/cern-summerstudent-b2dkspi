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

//int main(int argc, char** argv)
void fitBsd_MC_data()
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
	RooRealVar B_DTF_MM("B_DTF_MM", "m(DK_{s}#pi)", 5100., 5800.,"MeV");

    //Create RooDataSet
	RooArgList list =  RooArgList(B_DTF_MM);
	RooDataSet* data = new RooDataSet("data", "data", tree, list,"");

	RooDataSet* data_small = (RooDataSet*) data->reduce(SelectVars(RooArgSet(B_DTF_MM)));
	RooDataHist* data_binned = data_small->binnedClone();

	
    ///Define fit model
	///----------------

	///Signal model: a single CB function
	///-----------------------
	RooRealVar mean("mean", "mean", 5280.,5260.,5300.);
	RooRealVar sigma("sigma", "sigma", 20. ,0.,50.);
	RooRealVar alpha("alpha","alpha",5.,0.,10.);
	RooRealVar n("n","n",0.5,0.,10.);
	RooCBShape CB("CB","CB for signal",B_DTF_MM, mean, sigma, alpha, n);

        ///Background model
        ///-------------------------
	
	   //Bs background
	RooRealVar Bs("Bs","Bs",5366.88);//MeV
	RooRealVar Bd("Bd","Bd",5279.65);//MeV
	RooFormulaVar means("means","means","mean+Bs-Bd",RooArgList(mean,Bs,Bd));
	RooCBShape CBBs("CBBs","CBBs",B_DTF_MM,means,sigma, alpha, n);

	//Chebychev
	RooRealVar xa1 = RooRealVar("xa1","xa1",0,-1.,1.);
	RooRealVar xa2 = RooRealVar("xa2","xa2",0,-0.5,0.5);
	RooChebychev cheby("cheby", "cheby", B_DTF_MM, RooArgList(xa1));

	//Bd part.reco.bkg
        RooRealVar mean1("mean1", "mean1", 5085.);
        RooRealVar sigma1("sigma1", "sigma1", 47.179);
        RooGaussian Gauss1("Gauss1", "Gauss1", B_DTF_MM, mean1, sigma1);

        RooRealVar mean2("mean2", "mean2", 5104.5);
        RooRealVar sigma2("sigma2", "sigma2", 28.948);
        RooGaussian Gauss2("Gauss2", "Gauss2", B_DTF_MM, mean2, sigma2);

        RooRealVar mean3("mean3", "mean3", 5260.);
        RooRealVar sigma3("sigma3", "sigma3", 70.);
        RooGaussian Gauss3("Gauss3","Gauss3",B_DTF_MM, mean3, sigma3);

	float dn1 = 560.53;
	float dn2 = 586.00;
	float dn3 = 18.240;
	float df1 = dn1/(dn1+dn2+dn3);
	float df2 = dn2/(dn1+dn2+dn3);

	RooRealVar Bdf1("Bdf1","Bd fraction 1",df1);
	RooRealVar Bdf2("Bdf2","Bd fraction 2",df2);
	
	//Bs part.reco.bkg
	RooRealVar mean4("mean4", "mean4", 5173.1);
	RooRealVar sigma4("sigma4", "sigma4", 120.);
	RooGaussian Gauss4("Gauss4", "Gauss4", B_DTF_MM, mean4, sigma4);

	RooRealVar mean5("mean5", "mean5", 5182.2);
	RooRealVar sigma5("sigma5", "sigma5", 33.99);
	RooGaussian Gauss5("Gauss5", "Gauss5", B_DTF_MM, mean5, sigma5);

	float sn1 = 142.89;
	float sn2 = 942.11;
	float sf = sn1/(sn1+sn2);

	RooRealVar Bsf("Bsf","Bs fraction",sf);

    ///Total pdf
	///----------------------
	RooRealVar n_sig("n_sig", "n_sig", data->numEntries()/2., 0., data->numEntries());
	RooRealVar n_bkg("n_bkg", "n_bkg", data->numEntries()/2., 0., data->numEntries());
	RooRealVar n_BdDst("n_BdDst", "n_BdDst", data->numEntries()/2., 0., data->numEntries());
	RooRealVar n_BsDst("n_BsDst", "n_BsDst", data->numEntries()/10., 0., data->numEntries()/5);	
	RooRealVar n_Bs("n_Bs", "n_Bs", data->numEntries()/10., 0., data->numEntries()/5);
	
	//pdf list
	RooAbsPdf* sig = new RooAddPdf("sig", "sig", CB, n_sig);
	RooAbsPdf* BsDst = new RooAddPdf("BsDst","BsDst",RooArgList(Gauss4,Gauss5),RooArgList(Bsf));
	RooAbsPdf* BdDst = new RooAddPdf("BdDst","BdDst",RooArgList(Gauss1,Gauss2,Gauss3),RooArgList(Bdf1,Bdf2));
	RooAbsPdf* bkg = new RooAddPdf("bkg", "bkg", cheby, n_bkg);
	RooAbsPdf* Bsbkg = new RooAddPdf("Bsbkg", "Bsbkg", CBBs, n_Bs);

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
	TCanvas* c1= new TCanvas("");
	c1->SetGrid();

	RooPlot* frame_m= B_DTF_MM.frame();
	
	data->plotOn(frame_m,Name("Data"));
	pdf.plotOn(frame_m,Name("FullModel"));
	pdf.plotOn(frame_m,Components(*sig),LineColor(3), LineStyle(1), Name("sig"));
	pdf.plotOn(frame_m,Components(*bkg),LineColor(11), LineStyle(1), Name("bkg"));
	pdf.plotOn(frame_m,Components(*BdDst),LineColor(9), LineStyle(1), Name("BdDst"));
	pdf.plotOn(frame_m,Components(*BsDst),LineColor(5), LineStyle(1), Name("BsDst"));
	pdf.plotOn(frame_m,Components(*Bsbkg),LineColor(7), LineStyle(1), Name("Bsbkg"));

	TLegend leg(0.7, 0.7, 0.9, 0.9);
	leg.AddEntry(frame_m->findObject("FullModel"), "Full Model", "L");
	leg.AddEntry(frame_m->findObject("sig"), "sig", "L");
	leg.AddEntry(frame_m->findObject("BsDst"), "BsDst", "L");
	leg.AddEntry(frame_m->findObject("BdDst"), "BdDst", "L");
	leg.AddEntry(frame_m->findObject("bkg"), "bkg", "L");
	leg.AddEntry(frame_m->findObject("Bsbkg"), "Bsbkg", "L");

	cout<<"chi2 = "<<frame_m->chiSquare("FullModel","Data")<<endl;
        
	frame_m->Draw();
	leg.DrawClone();
	c1->Print("fitB_MC_data.png");
	
    //return 0;
  
}


