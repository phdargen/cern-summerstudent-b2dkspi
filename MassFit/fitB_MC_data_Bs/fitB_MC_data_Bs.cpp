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
void fitB_MC_data_Bs()
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
	RooRealVar mean1("mean1", "mean1", 5280.,5150.,5350.); 
	RooRealVar sigma1("\\sigma_{1}", "sigma1", 20.,10.,30.);	

	RooGaussian Gauss1("Gauss1", "Gauss1", B_DTF_MM, mean1, sigma1);

	RooRealVar mean0("mean0", "mean0", 5280.,5150.,5350.);
	RooRealVar sigma0("sigma0", "sigma0", 30., 20., 50.);
	RooGaussian Gauss2("Gauss2", "Gauss2", B_DTF_MM, mean0, sigma0);

	RooRealVar mean2("mean2", "mean2", 5080.,4950.,5200.);
	//RooRealVar m2("m2","m2",100,0.1,100);
	RooRealVar sigma2("sigma2", "sigma2", 100. ,20.,150.);
    	RooRealVar alpha("alpha","alpha",1.,0.,10.);
	RooRealVar n("n","n",0.5,0.,10.);
	RooCBShape Crys2("Crys2","Crys2",B_DTF_MM, mean2, sigma2, alpha, n);
	//Add a second Gaussian 
    //...

	//part. reco. bkg
        RooRealVar mean01("mean01", "mean01", 5085.);
        RooRealVar sigma01("sigma01", "sigma01", 50.);
        RooRealVar width1("width1","width1",5.4e-5,0,0.1);
        //RooGaussian Gauss1("Gauss1", "Gauss1", B_DTF_MM, mean1, sigma1);
        RooBreitWigner Gauss01("Gauss01", "Gauss01", B_DTF_MM, mean01, sigma01);

        RooRealVar mean02("mean02", "mean02",5115.);
        RooRealVar sigma02("sigma02", "sigma02", 21.287);
        RooGaussian Gauss02("Gauss02", "Gauss02", B_DTF_MM, mean02, sigma02);
        //RooBreitWigner Gauss2("Gauss2", "Gauss2", B_DTF_MM, mean2, sigma2);

        RooRealVar mean03("mean03", "mean03", 5100.);
        RooRealVar sigma03("sigma03", "sigma03", 38.869);
        RooGaussian Gauss03("Gauss03","Gauss03",B_DTF_MM, mean03, sigma03);
        //RooBreitWigner Gauss3("Gauss3","Gauss3",B_DTF_MM, mean3, sigma3);
        ///----------------------

	RooRealVar Bs("Bs","Bs",5366.88);//MeV
	RooRealVar Bd("Bd","Bd",5279.65);//MeV
	//RooRealVar mean1("mean1", "mean1", 5280.,5150.,5350.);
	RooFormulaVar means("means","means","mean1+Bs-Bd",RooArgList(mean1,Bs,Bd));
	RooGaussian GaussBs("GaussBs","GaussBs",B_DTF_MM,means,sigma1);
	RooFormulaVar means2("means2","means2","mean0+Bs-Bd",RooArgList(mean0,Bs,Bd));
	RooGaussian GaussBs2("GaussBs2","GaussBs2",B_DTF_MM,means2,sigma0);
	///Background model
	///-------------------------

    //Chebychev
    RooRealVar xa1 = RooRealVar("xa1","xa1",0,-1.,1.);
    RooRealVar xa2 = RooRealVar("xa2","xa2",0,-0.5,0.5);
    RooRealVar xa3 = RooRealVar("xa3","xa3",0,-0.5,0.5);
    RooRealVar xa4 = RooRealVar("xa4","xa4",0,-0.5,0.5);
    RooRealVar xa5 = RooRealVar("xa5","xa5",0,-0.5,0.5);
    RooChebychev cheby("cheby", "cheby", B_DTF_MM, RooArgList(xa1, xa2, xa3, xa4, xa5));
	
    ///Total pdf
	///----------------------
	RooRealVar n_sig1("n_sig1", "n_sig1", data->numEntries()/2., 0., data->numEntries());
	float n1 = 3.5630e+02*2;
	//RooRealVar n1 = RooRealVar("n1","n1",3.5630e+02,0,data->numEntries()/2.);
	float n2 = n1*2.0849e+02/3.5630e+02;
	float n3 = n1*5.7500e+02/3.5630e+02;
	//RooRealVar n2 = RooRealVar("n2","n2",n1*2.0849e+02/3.5630e+02);
	//RooRealVar n3 = RooRealVar("n3","n3",n1*5.7500e+02/3.5630e+02);

	RooRealVar n_bkg("n_bkg", "n_bkg", n1);
	RooRealVar n_sig2("n_sig2", "n_sig2", data->numEntries()/4. ,0., data->numEntries()/2);
	RooRealVar n_bkg2("n_bkg2", "n_bkg2", n2);
	RooRealVar n_bkg3("n_bkg3", "n_bkg3", n3);
	RooRealVar n_bkg4("n_bkg4", "n_bkg4", data->numEntries()/2., 0., data->numEntries());
	
	RooRealVar n_Bs("n_Bs","n_Bs",data->numEntries()/4. ,0., data->numEntries()/2);
	RooFormulaVar n_Bs2("n_Bs2","n_Bs2","n_Bs*n_sig2/n_sig1",RooArgList(n_Bs,n_sig2,n_sig1));
	
	//pdf list
	RooAbsPdf* bkg = new RooAddPdf("bkg", "bkg", Gauss01, n_bkg);
	RooAbsPdf* sig1 = new RooAddPdf("sig1", "sig1", Gauss1, n_sig1);
	RooAbsPdf* sig2 = new RooAddPdf("sig2", "sig2", Gauss2, n_sig2);
	RooAbsPdf* bkg2 = new RooAddPdf("bkg2", "bkg2", Gauss02, n_bkg2);
	RooAbsPdf* bkg3 = new RooAddPdf("bkg3", "bkg3", Gauss03, n_bkg3);
	RooAbsPdf* bkg4 = new RooAddPdf("bkg4", "bkg4", cheby, n_bkg4);
	RooAbsPdf* bkgs = new RooAddPdf("bkgs", "bkgs", GaussBs, n_Bs);
	RooAbsPdf* bkgs2 = new RooAddPdf("bkgs2", "bkgs2", GaussBs2, n_Bs2);

	RooArgList pdf0;
	pdf0.add(*bkg);
	pdf0.add(*bkg2);
	pdf0.add(*sig1);
	pdf0.add(*sig2);
	pdf0.add(*bkg3);
	pdf0.add(*bkg4);
	pdf0.add(*bkgs);
	pdf0.add(*bkgs2);

	//N list
	RooArgList num;
	num.add(n_bkg);
	num.add(n_bkg2);
	num.add(n_sig1);
	num.add(n_sig2);
	num.add(n_bkg3);
	num.add(n_bkg4);
	num.add(n_Bs);
	num.add(n_Bs2);

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
	pdf.plotOn(frame_m,Components(*bkg2),LineColor(7), LineStyle(1), Name("bkg2"));
	pdf.plotOn(frame_m,Components(*sig1),LineColor(3), LineStyle(1), Name("sig1"));
	pdf.plotOn(frame_m,Components(*sig2),LineColor(5), LineStyle(1), Name("sig2"));
	pdf.plotOn(frame_m,Components(*bkg3),LineColor(8), LineStyle(1), Name("bkg3"));
	pdf.plotOn(frame_m,Components(*bkg4),LineColor(9), LineStyle(1), Name("bkg4"));
	pdf.plotOn(frame_m,Components(*bkgs),LineColor(11), LineStyle(1), Name("bkgs"));
	pdf.plotOn(frame_m,Components(*bkgs2),LineColor(13), LineStyle(1), Name("bkgs2"));
	TLegend leg(0.7, 0.7, 0.9, 0.9);
	leg.AddEntry(frame_m->findObject("FullModel"), "Full Model", "L");
	leg.AddEntry(frame_m->findObject("bkg"), "bkg", "L");
	leg.AddEntry(frame_m->findObject("bkg2"), "bkg2", "L");
	leg.AddEntry(frame_m->findObject("sig1"), "sig1", "L");
	leg.AddEntry(frame_m->findObject("sig2"), "sig2", "L");
        leg.AddEntry(frame_m->findObject("bkg3"), "bkg3", "L");
	leg.AddEntry(frame_m->findObject("bkg4"), "bkg4", "L");
	leg.AddEntry(frame_m->findObject("bkgs"), "bkgs", "L");
	leg.AddEntry(frame_m->findObject("bkgs2"), "bkgs2", "L");

	cout<<"chi2 = "<<frame_m->chiSquare("FullModel","Data")<<endl;
        
	frame_m->Draw();
	leg.DrawClone();
	c1->Print("fitB_MC_data.pdf");
	
    //return 0;
  
}


