//Dalitz plots for B0->DKspi
//Yuxiao Wang
#include <TChain.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TNtuple.h>

#include <iostream>
using namespace std;

void BsDstKspi_DD() {

	gStyle->SetOptStat(0);
	
	TH1::SetDefaultSumw2();
	TH2::SetDefaultSumw2();

 	//Load file	
    TChain* tree=new TChain("DecayTree");
    tree->Add("/eos/lhcb/user/p/phdargen/summerStudents21/Stripped/MC_B2DKspi_DD_12_BsDstKspi.root");

        float massPion = 139.570; //MeV
	float massKaon = 493.677; //MeV
  	float massProton = 938.272; //MeV
	//Needed branches (add some more)
	Double_t D_PX, D_PY, D_PZ, D_MM, D_PE;
	Double_t Ks_PX, Ks_PY, Ks_PZ, Ks_MM, Ks_PE;
	Double_t pi_PX, pi_PY, pi_PZ, pi_PE;
	Double_t pi1_D_PX, pi1_D_PY, pi1_D_PZ;
	Double_t pi2_D_PX, pi2_D_PY, pi2_D_PZ;
	Double_t K_D_PX, K_D_PY, K_D_PZ;
	Double_t pi_ProbNNpi,pi1_D_ProbNNpi,pi2_D_ProbNNpi;

	tree->SetBranchAddress("pi_ProbNNpi",&pi_ProbNNpi);
	tree->SetBranchAddress("pi1_D_ProbNNpi",&pi1_D_ProbNNpi);
	tree->SetBranchAddress("pi2_D_ProbNNpi",&pi2_D_ProbNNpi);
	tree->SetBranchAddress("D_PX",&D_PX);
	tree->SetBranchAddress("D_PY",&D_PY);
	tree->SetBranchAddress("D_PZ",&D_PZ);
	tree->SetBranchAddress("D_MM",&D_MM);
	tree->SetBranchAddress("D_PE",&D_PE);
	tree->SetBranchAddress("Ks_PX",&Ks_PX);
	tree->SetBranchAddress("Ks_PY",&Ks_PY);
	tree->SetBranchAddress("Ks_PZ",&Ks_PZ);
	tree->SetBranchAddress("Ks_MM",&Ks_MM);
	tree->SetBranchAddress("Ks_PE",&Ks_PE);
	tree->SetBranchAddress("pi_PX",&pi_PX);
	tree->SetBranchAddress("pi_PY",&pi_PY);
	tree->SetBranchAddress("pi_PZ",&pi_PZ);
	tree->SetBranchAddress("pi_PE",&pi_PE);
	tree->SetBranchAddress("pi1_D_PX",&pi1_D_PX);
	tree->SetBranchAddress("pi1_D_PY",&pi1_D_PY);
	tree->SetBranchAddress("pi1_D_PZ",&pi1_D_PZ);
	tree->SetBranchAddress("pi2_D_PX",&pi2_D_PX);
	tree->SetBranchAddress("pi2_D_PY",&pi2_D_PY);
	tree->SetBranchAddress("pi2_D_PZ",&pi2_D_PZ);
	tree->SetBranchAddress("K_D_PX",&K_D_PX);
	tree->SetBranchAddress("K_D_PY",&K_D_PY);
	tree->SetBranchAddress("K_D_PZ",&K_D_PZ);
	//Create output file
	TFile* output = new TFile("/afs/cern.ch/work/y/yuxiao/public/Dalitz/dalitz_BsDstKspi_DD.root","RECREATE");
    	//TTree* summary_tree = tree->CloneTree(0);

	TH2F* plot = new TH2F("Dalitz",";m^{2}(D^{-}K_{s})[MeV^{2}/c^{4}];m^{2}(K_{s}#pi)[MeV^{2}/c^{4}]",100,0,30e+06,100,0,13e+06);

	int nEvents = tree->GetEntries();
	for ( Int_t j = 0 ; j < nEvents ; j++ ) {
		tree->GetEntry(j) ;
		if (0ul == (j % 10000ul)) cout << "Read event " << j << "/" << nEvents << endl;
		TLorentzVector p_Ks, p_D, p_pi;
		
		p_Ks.SetXYZM(Ks_PX,Ks_PY,Ks_PZ,Ks_MM);
		p_D.SetXYZM(D_PX,D_PY,D_PZ,D_MM);
		p_pi.SetXYZM(pi_PX, pi_PY, pi_PZ, massPion);

		//p_Ks.SetPxPyPzE(Ks_PX,Ks_PY,Ks_PZ,Ks_PE);
		//p_D.SetPxPyPzE(D_PX,D_PY,D_PZ,D_PE);
		//p_pi.SetPxPyPzE(pi_PX, pi_PY, pi_PZ, pi_PE);

		plot->Fill((p_Ks+p_D).M2(), (p_Ks+p_pi).M2());

	}
	plot->Draw("colz");	
	output->Write();

	TCanvas* c = new TCanvas ("c","");
	
	c->cd();
	plot->Draw("colz");
	//plot->GetXaxis()->SetTitle("m^{2}(D^{-}K_{s})[MeV^{2}/c^{4}]");
	//plot->GetYaxis()->SetTitle("m^{2}(K_{s}#pi)[MeV^{2}/c^{4}]");
	c->SaveAs("Dalitz_BsDstKspi_DD.png");

	output->Close();
}

