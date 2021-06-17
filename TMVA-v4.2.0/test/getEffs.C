#include <iostream>
#include <iomanip>
#include "tmvaglob.C"
#include "TH1.h"
#include "TROOT.h"
#include "TList.h"
#include "TIterator.h"
#include "TStyle.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TH2.h"
#include "TFormula.h"
#include "TFile.h"

using namespace std;

void eff_for_BDT(double BDT, int nS= 161943, int nB=51495,  string file = "PIDK_myTMVAs2.root"){

  TFile* f= TFile::Open(file.c_str());
  TH1D* effS=(TH1D*)f->Get("Method_BDT/BDT/MVA_BDT_effS");
  TH1D* effB=(TH1D*)f->Get("Method_BDT/BDT/MVA_BDT_effB");

  cout << "-----------------------------------------" << endl ;
  cout << "BDT > " << BDT << ":" << endl;
 
  double eff_s = effS->GetBinContent(effS->FindBin(BDT));
  double eff_b = effB->GetBinContent(effB->FindBin(BDT));
  double n_s =  eff_s*nS ;
  double n_b =  eff_b*nB ;

  cout << "effS = " << eff_s << endl;
  cout << "effB = " << eff_b << endl;
  cout << "nS = " << n_s << endl;
  cout << "nB = " << n_b << endl;
  cout << "S/B = " << n_s/n_b << endl;
  cout << "S/(S+B) = " << n_s/(n_s+n_b) << endl;
  cout << "Sig = " << n_s/sqrt(n_s+n_b) << endl;
  cout << endl;

}

void eff_BDT(string file = "PIDK_myTMVAs2.root"){

  TFile* f= TFile::Open(file.c_str());
  TH1D* effS=(TH1D*)f->Get("Method_BDT/BDT/MVA_BDT_effS");
  TH1D* effB=(TH1D*)f->Get("Method_BDT/BDT/MVA_BDT_effB");

  cout << "# bins =" << effS->GetNbinsX() << endl;

  for(int i=1; i< effS->GetNbinsX(); i++){
    cout << "-----------------------------------------" << endl ;
    cout << "BDT = " << effS->GetXaxis()->GetBinCenter(i)  << ":" << endl;
    cout << "effS = " << effS->GetBinContent(i) << endl;
    cout << "effB = " << effB->GetBinContent(i) << endl;
    cout << endl;
  }
}

/*
double BDT_for_eff(double eff,bool signal = true ,string file = "myTMVA_ee.root"){

  TFile* f= TFile::Open(file.c_str());
  TH1D* effS=(TH1D*)f->Get("Method_BDT/BDT/MVA_BDT_effS");
  TH1D* effB=(TH1D*)f->Get("Method_BDT/BDT/MVA_BDT_effB");

  cout << "effS = " << eff  << endl;

  double BDT;

  for(int i=1; i< effS->GetNbinsX(); i++){
    if(abs(eff - effS->GetBinContent(effS->FindBin(i)))){
      BDT= effS->GetXaxis().GetBinCenter(i);
      continue;
    }
  }

  cout << "effB = " << effB->GetBinContent(effB->FindBin(BDT)) << endl;


  cout << "BDT > " << BDT << ":" << endl;


  return BDT;

}
*/

void getEffs(){
  
  //eff_BDT();

  eff_for_BDT(-0.067);
  eff_for_BDT(0.0);
  eff_for_BDT(0.1);
  eff_for_BDT(0.16);
  eff_for_BDT(0.2);


}
