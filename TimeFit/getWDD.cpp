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
void getWDD()
{   
    //Binned or unbinned fit?
    bool binned=false;
    
	///Load file
        /*
    TChain* tree = new TChain("DecayTree");    
    tree->Add("/afs/cern.ch/work/y/yuxiao/cern-summerstudent-b2dkspi/withvetocuts/output_DD11.root");
        tree->Add("/afs/cern.ch/work/y/yuxiao/cern-summerstudent-b2dkspi/withvetocuts/output_DD12.root");
        tree->Add("/afs/cern.ch/work/y/yuxiao/cern-summerstudent-b2dkspi/withvetocuts/output_DD15.root");
        tree->Add("/afs/cern.ch/work/y/yuxiao/cern-summerstudent-b2dkspi/withvetocuts/output_DD16.root");
        tree->Add("/afs/cern.ch/work/y/yuxiao/cern-summerstudent-b2dkspi/withvetocuts/output_DD17.root");
        tree->Add("/afs/cern.ch/work/y/yuxiao/cern-summerstudent-b2dkspi/withvetocuts/output_DD18.root");	
*/	
	TFile* file= new TFile("/afs/cern.ch/work/m/mbuhayeu/public/cern-summerstudent-b2dkspi-master/Selection/BDTG/B2DKspi_data.root");
	TTree* tree = (TTree*) file->Get("DecayTree");
    
    //Disable all but needed branches
    tree->SetBranchStatus("*",0);  
    tree->SetBranchStatus("B_DTF_MM",1);  
    tree->SetBranchStatus("BDTG",1);
    tree->SetBranchStatus("KsCat",1);  
    tree->SetBranchStatus("B_DTF_TAU",1); 
    //Define branches
	RooRealVar B_DTF_MM("B_DTF_MM", "", 5100., 5800.,"MeV");

    //Create RooDataSet
	//RooArgList list =  RooArgList(B_DTF_MM);
	RooRealVar BDTG("BDTG", "BDTG", 0.);
	RooCategory KsCat("KsCat","KsCat") ;
	KsCat.defineType("DD",1);
	RooArgList list =  RooArgList(B_DTF_MM,BDTG,KsCat);
	RooDataSet*  data = new RooDataSet("data","data",tree,list, " KsCat == 1 && BDTG > 0.6422" );
	//RooDataSet* data = new RooDataSet("data", "data", tree, list,"");

	RooDataSet* data_small = (RooDataSet*) data->reduce(SelectVars(RooArgSet(B_DTF_MM)));
	RooDataHist* data_binned = data_small->binnedClone();

	
    ///Define fit model
	///----------------

	///Signal model: a single CB function
	///-----------------------
	RooRealVar mean("mean", "mean", 5277.3,5260.,5300.);
	RooRealVar sigma("sigma", "sigma", 24.4 ,0.,50.);
	RooRealVar alpha("alpha","alpha",1.9902);
	RooRealVar n("n","n",1.0374);
	RooCBShape CB("CB","CB for signal",B_DTF_MM, mean, sigma, alpha, n);

	RooFormulaVar meanb("meanb","meanb","mean",mean);
	RooRealVar sigmab("sigmab","sigmab",11.1,0.,50.);
	RooRealVar alpha2("alpha2","alpha2",-2.1031);
	RooRealVar n2("n2","n2",0.26483);
	RooCBShape CB2("CB2","CB2 for signal",B_DTF_MM,meanb,sigmab,alpha2,n);
        
	RooRealVar f("f","f",0.75051);

	RooAbsPdf* sig = new RooAddPdf("sig","sig",RooArgList(CB,CB2), f);

	///Background model
        ///-------------------------
	
	   //Bs background
	RooRealVar Bs("Bs","Bs",5366.88);//MeV
	RooRealVar Bd("Bd","Bd",5279.65);//MeV
	RooFormulaVar means("means","means","mean+Bs-Bd",RooArgList(mean,Bs,Bd));
	RooFormulaVar means2("means2","means2","meanb+Bs-Bd",RooArgList(meanb,Bs,Bd));

	RooCBShape CBBs("CBBs","CBBs",B_DTF_MM,means,sigma, alpha, n);
	RooCBShape CBBs2("CBBs2","CBBs2",B_DTF_MM,means2,sigmab,alpha2,n2);

	RooAbsPdf* sigBs = new RooAddPdf("sigBs","sigBs",RooArgList(CBBs,CBBs2),f);
	
	
	//Chebychev
	RooRealVar xa1 = RooRealVar("xa1","xa1",0,-1.,1.);
	RooRealVar xa2 = RooRealVar("xa2","xa2",0,-0.5,0.5);
	RooChebychev cheby("cheby", "cheby", B_DTF_MM, RooArgList(xa1));

	//Bd part.reco.bkg
        RooRealVar mean1("mean1", "mean1", 5078.8);
        RooRealVar sigma1("sigma1", "sigma1", 99.797);
        RooGaussian Gauss1("Gauss1", "Gauss1", B_DTF_MM, mean1, sigma1);

        RooRealVar mean2("mean2", "mean2", 5098.3);
        RooRealVar sigma2("sigma2", "sigma2", 32.549);
        RooGaussian Gauss2("Gauss2", "Gauss2", B_DTF_MM, mean2, sigma2);

        RooRealVar mean3("mean3", "mean3", 5469.5);
        RooRealVar sigma3("sigma3", "sigma3", 399.97);
        RooGaussian Gauss3("Gauss3","Gauss3",B_DTF_MM, mean3, sigma3);

	float dn1 = 308.57;
	float dn2 = 2.6344e+03;
	float dn3 = 1.7201e+02;
	float df1 = dn1/(dn1+dn2+dn3);
	float df2 = dn2/(dn1+dn2+dn3);

	RooRealVar Bdf1("Bdf1","Bd fraction 1",df1);
	RooRealVar Bdf2("Bdf2","Bd fraction 2",df2);
	
	//Bs part.reco.bkg
	RooRealVar mean4("mean4", "mean4", 5207.3, 5207.3-1.5*2.48,5207.3+1.5*2.48);
	RooRealVar sigma4("sigma4", "sigma4", 18.963,18.963-3*1.43,18.963+3*1.43);
	RooGaussian Gauss4("Gauss4", "Gauss4", B_DTF_MM, mean4, sigma4);

	RooRealVar mean5("mean5", "mean5", 5164.2,5164.2-3*3.08,5164.2+3*3.08);
	RooRealVar sigma5("sigma5", "sigma5", 26.953,26.953-2.5*2.21,26.953+2.5*2.21);
	RooGaussian Gauss5("Gauss5", "Gauss5", B_DTF_MM, mean5, sigma5);

	RooRealVar alpha4("alpha4","alpha4",7.2323e-01,7.2323e-01-1*0.16,7.2323e-01+1*0.16);
	RooRealVar alpha5("alpha5","alpha5",-1.5078,-1.5078-3*0.157,-1.5078+3*0.157);

	RooRealVar n4("n4","n4",2.9811,2.9811-3*1.02,2.9811+3*1.02);
	RooRealVar n5("n5","n5",1.1430,1.1430-3*0.159,1.1430+3*0.159);

	RooCBShape CB4("CB4","CB4",B_DTF_MM,mean4,sigma4, alpha4, n4);
	RooCBShape CB5("CB5","CB5",B_DTF_MM,mean5,sigma5, alpha5, n5);

	RooRealVar mean6("mean6", "mean6",5105.6/*,5105.6-1*4.53,5105.6+1*4.53*/);
	RooRealVar sigma6("sigma6", "sigma6", 5.4254,5.4254-2*3.24,5.4254+2*3.24);
	RooGaussian Gauss6("Gauss6", "Gauss6", B_DTF_MM, mean6, sigma6);

	//float sn1 = 1.6683e+03;
	RooRealVar sn1("sn1","sn1",1.6683e+03,1.6683e+03-3*187.,1.6683e+03+3*187.);
	//float sn2 = 1.6567e+03;
	RooRealVar sn2("sn2","sn2",1.6567e+03,1.6567e+03-3*188.,1.6567e+03+3*188.);
	//float sn3 = 1.8110e+01;
	RooRealVar sn3("sn3","sn3",1.8110e+01,1.8110e+01-3*14.2,1.8110e+01+3*14.2);
	//float sf = sn1/(sn1+sn2+sn3);
	RooFormulaVar sf("sf","sf","sn1/(sn1+sn2+sn3)",RooArgList(sn1,sn2,sn3));
	//float sf2 = sn2/(sn1+sn2+sn3);
	RooFormulaVar sf2("sf2","sf2","sn2/(sn1+sn2+sn3)",RooArgList(sn1,sn2,sn3));

	RooRealVar Bsf("Bsf","Bs fraction",sf.getVal());
	RooRealVar Bsf2("Bsf2", "Bs2 fraction",sf2.getVal());

    ///Total pdf
	///----------------------
	RooRealVar n_sig("n_sig", "n_sig", data->numEntries()/2., 0., data->numEntries());
	RooRealVar n_bkg("n_bkg", "n_bkg", data->numEntries()/2., 0., data->numEntries());
	RooRealVar n_BdDst("n_BdDst", "n_BdDst", data->numEntries()/2., 0., data->numEntries());
	RooRealVar n_BsDst("n_BsDst", "n_BsDst", data->numEntries()/10., 0., data->numEntries()/5);	
	RooRealVar n_Bs("n_Bs", "n_Bs", data->numEntries()/10., 0., data->numEntries()/5);
	
	//pdf list
	//RooAbsPdf* sig = new RooAddPdf("sig", "sig", CB, n_sig);
	RooAbsPdf* BsDst = new RooAddPdf("BsDst","BsDst",RooArgList(CB4,CB5,Gauss6),RooArgList(Bsf,Bsf2));
	RooAbsPdf* BdDst = new RooAddPdf("BdDst","BdDst",RooArgList(Gauss1,Gauss2,Gauss3),RooArgList(Bdf1,Bdf2));
	//RooAbsPdf* bkg = new RooAddPdf("bkg", "bkg", cheby, n_bkg);
	//RooAbsPdf* Bsbkg = new RooAddPdf("Bsbkg", "Bsbkg", CBBs, n_Bs);

	RooArgList pdf0;
	pdf0.add(*sig);
	//pdf0.add(CB);
	pdf0.add(*BsDst);
	pdf0.add(*BdDst);
	//pdf0.add(*bkg);
	//pdf0.add(*Bsbkg);
	pdf0.add(cheby);
	pdf0.add(*sigBs);

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
    if(binned) result = pdf.fitTo(*data_binned,Save(kTRUE),Extended(true),SumW2Error(kTRUE));
	
    //unbinned Fit
    else result = pdf.fitTo(*data,Save(kTRUE),Extended(true),SumW2Error(kTRUE));
	
    cout << "result is --------------- "<<endl;
	result->Print(); 
 
	///Plot 
	///----------
/*
	RooPlot* frame_m= B_DTF_MM.frame();
	
	data->plotOn(frame_m,Name("Data"));
	pdf.plotOn(frame_m,Name("FullModel"));
	pdf.plotOn(frame_m,Components(*sig),LineColor(3), LineStyle(1), Name("sig"));
	//pdf.plotOn(frame_m,Components(*bkg),LineColor(11), LineStyle(1), Name("bkg"));
	pdf.plotOn(frame_m,Components(*BdDst),LineColor(9), LineStyle(1), Name("BdDst"));
	pdf.plotOn(frame_m,Components(*BsDst),LineColor(5), LineStyle(1), Name("BsDst"));
	//pdf.plotOn(frame_m,Components(*Bsbkg),LineColor(7), LineStyle(1), Name("Bsbkg"));
	//pdf.plotOn(frame_m,Components(CB),LineColor(3), LineStyle(1), Name("sig"));
	pdf.plotOn(frame_m,Components(cheby),LineColor(11), LineStyle(1), Name("bkg"));
	pdf.plotOn(frame_m,Components(*sigBs),LineColor(7), LineStyle(1), Name("Bsbkg"));
	data->plotOn(frame_m,Name("Data"));

	TLegend leg(0.58,0.55,0.94,0.89);
	//leg.SetFillColor(0);
	leg.AddEntry(frame_m->findObject("FullModel"), "Full Model", "f");
	leg.AddEntry(frame_m->findObject("sig"), "sig", "f");
	leg.AddEntry(frame_m->findObject("BsDst"), "BsDst", "f");
	leg.AddEntry(frame_m->findObject("BdDst"), "BdDst", "f");
	leg.AddEntry(frame_m->findObject("bkg"), "bkg", "f");
	leg.AddEntry(frame_m->findObject("Bsbkg"), "Bsbkg", "f");

	cout<<"chi2 = "<<frame_m->chiSquare("FullModel","Data")<<endl;

	RooHist* hresid = frame_m->residHist("Data","FullModel");

	RooHist* hpull = frame_m->pullHist("Data","FullModel");

	RooPlot* frame2 = B_DTF_MM.frame(Title("Residual Distribution")) ;
	frame2->addPlotable(hresid,"P") ;

	RooPlot* frame3 = B_DTF_MM.frame(Title("Pull Distribution")) ;
	frame3->addPlotable(hpull,"P") ;

	TCanvas* c = new TCanvas("fit_chi2residpull","B0 Mass Fit",800,600) ;
	c->Divide(1,2) ;
	c->cd(1) ; gPad->SetBorderMode(1); gPad->SetTopMargin(1); gPad->SetBottomMargin(0.15) ; gPad->SetRightMargin(0.03); gPad->SetPad(0.03, 0.27, 0.95, 0.95); frame_m->Draw() ; leg.DrawClone();
	//c->cd(2) ; gPad->SetLeftMargin(0.15) ; frame2->GetYaxis()->SetTitleOffset(1.6) ; frame2->Draw() ;
	c->cd(2) ; gPad->SetTopMargin(0); gPad->SetPad(0.03, 0.02, 0.95, 0.27); gPad->SetRightMargin(0.03); gPad->SetBottomMargin(0.2); frame3->Draw();
	 
	c->SaveAs("/afs/cern.ch/work/y/yuxiao/public/FitResults/BDD.png");

	B_DTF_MM.setRange("signal", 5.2793e+03 - 50 , 5.2793e+03 + 50) ;
	double BdDstyield = BdDst -> createIntegral(B_DTF_MM,NormSet(B_DTF_MM),Range("signal")) ->getVal();	
	double BsDstyield = BsDst -> createIntegral(B_DTF_MM,NormSet(B_DTF_MM),Range("signal")) ->getVal();
	double Bsyield = sigBs -> createIntegral(B_DTF_MM,NormSet(B_DTF_MM),Range("signal")) ->getVal();
	double bkgyield = cheby.createIntegral(B_DTF_MM,NormSet(B_DTF_MM),Range("signal")) ->getVal();
	double sigyield = sig -> createIntegral(B_DTF_MM,NormSet(B_DTF_MM),Range("signal")) ->getVal();

	ofstream myfile;
	myfile.open("/afs/cern.ch/work/y/yuxiao/public/FitResults/BDD.txt");
	myfile<<"the BdDst yield: "<<BdDstyield*n_BdDst.getVal()<<"\n";
	myfile<<"the BsDst yield: "<<BsDstyield*n_BsDst.getVal()<<"\n";
	myfile<<"the Bs yield: "<<Bsyield*n_Bs.getVal()<<"\n";
	myfile<<"the 1-order bkg yield: "<<bkgyield*n_bkg.getVal()<<"<- Mary needs this one!"<<"\n";
	myfile<<"sig yield: "<<sigyield*n_sig.getVal()<<"<- Mary needs this one!"<<"\n";

	myfile.close();
	//frame_m->Draw();
	//leg.DrawClone();
	
    //return 0;
*/
        RooStats::SPlot* sPlot = new RooStats::SPlot("sPlot","sPlot",*data,&pdf,RooArgList(n_sig,n_BsDst,n_BdDst,n_bkg,n_Bs));

        TFile* fi = new TFile("/afs/cern.ch/work/y/yuxiao/public/WeightTreeDD.root","RECREATE");
        TTree* newtree = new TTree("DecayTree","DecayTree");

        Double_t BMass0,  BDTG0, BMass, Weight, B_DTF_TAU, B_DTF_TAU2;
        Int_t KsCat0;

        tree->SetBranchAddress("B_DTF_MM",&BMass0);
        tree->SetBranchAddress("KsCat",&KsCat0);
        tree->SetBranchAddress("BDTG",&BDTG0);
        tree->SetBranchAddress("B_DTF_TAU",&B_DTF_TAU);

        newtree->Branch("B_DTF_MM",&BMass,"B_DTF_MM/D");
        newtree->Branch("Weight",&Weight,"Weight/D");
        newtree->Branch("B_DTF_TAU",&B_DTF_TAU2,"B_DTF_TAU/D");

        int nEvents = tree->GetEntries();
        int j = 0;
        for(int i=0; i<nEvents; i++){
                tree->GetEntry(i) ;
                if (0ul == (i % 10000ul)) cout << "Read event " << i << "/" << nEvents << " B_DTF_TAU="<<B_DTF_TAU<<endl;

                if( (BMass0<=5100)||(BMass0>=5800) )continue;
                if( KsCat0!=1 )continue;
                if( BDTG0 <= 0.6422)continue;

                B_DTF_TAU2 = B_DTF_TAU/1000;
                if (0ul == (i % 10000ul)) cout << " B_DTF_TAU="<<B_DTF_TAU2<<endl;
		BMass = BMass0;
                Weight = sPlot->GetSWeight(j,"n_sig");

                newtree->Fill();
                j++;
        }


        fi->Write();
        fi->Close();
  
}


