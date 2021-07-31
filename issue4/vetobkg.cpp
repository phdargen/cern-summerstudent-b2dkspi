//Select events and fill them into a new tree
//Philippe d'Argent
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

void vetobkg() {

	TH1::SetDefaultSumw2();
	TH2::SetDefaultSumw2();

 	//Load file	
    TChain* tree=new TChain("DecayTree");
    tree->Add("/eos/lhcb/user/p/phdargen/summerStudents21/Stripped/Data_B2DKspi_LL_11.root");

        float massPion = 139.570; //MeV
	float massKaon = 493.677; //MeV
  	float massProton = 938.272; //MeV
	//Needed branches (add some more)
	Double_t D_PX, D_PY, D_PZ, D_MM;
	Double_t pi1_D_PX, pi1_D_PY, pi1_D_PZ;
	Double_t pi2_D_PX, pi2_D_PY, pi2_D_PZ;
	Double_t K_D_PX, K_D_PY, K_D_PZ;

	tree->SetBranchAddress("D_PX",&D_PX);
	tree->SetBranchAddress("D_PY",&D_PY);
	tree->SetBranchAddress("D_PZ",&D_PZ);
	tree->SetBranchAddress("D_MM",&D_MM);
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
	TFile* output = new TFile("vetobkg.root","RECREATE");
    	//TTree* summary_tree = tree->CloneTree(0);

    //Add new branches
    //...
    
	TH1F *B0Mass = new TH1F("B0Mass","B0Mass",100,1500,3000);
	TH1F *Misspi1 = new TH1F("Misspi1","Misspi1",100,1840,3000);
	TH1F *Misspi2 = new TH1F("Misspi2","Misspi2",100,1840,3000);
	TH1F *pitoproton = new TH1F("pitoproton","pitoproton",100,1800,7000);

	int nEvents = tree->GetEntries();
	for ( Int_t j = 0 ; j < nEvents ; j++ ) {
		tree->GetEntry(j) ;
		if (0ul == (j % 10000ul)) cout << "Read event " << j << "/" << nEvents << endl;
		TLorentzVector D_D, pi1_D, pi1_D_asK_MissID, K_D, pi2_D, pi2_D_asK_MissID, pi1_D_asP_MissID;
		D_D.SetXYZM(D_PX,D_PY,D_PZ,D_MM);
		pi1_D.SetXYZM(pi1_D_PX,pi1_D_PY,pi1_D_PZ,massPion);
		pi1_D_asK_MissID.SetXYZM(pi1_D_PX,pi1_D_PY,pi1_D_PZ,massKaon);
		K_D.SetXYZM(K_D_PX,K_D_PY,K_D_PZ,massKaon);
		pi2_D.SetXYZM(pi2_D_PX,pi2_D_PY,pi2_D_PZ,massPion);
		pi2_D_asK_MissID.SetXYZM(pi2_D_PX,pi2_D_PY,pi2_D_PZ,massKaon);
		pi1_D_asP_MissID.SetXYZM(pi1_D_PX,pi1_D_PY,pi1_D_PZ,massProton);	
		B0Mass->Fill((D_D + pi1_D).M());
		Misspi1->Fill((K_D + pi1_D_asK_MissID + pi2_D).M());
		Misspi2->Fill((K_D + pi2_D_asK_MissID + pi1_D).M());
		pitoproton->Fill((K_D + pi1_D_asP_MissID + pi2_D).M());

	}

	output->Write();

	TCanvas* c1 = new TCanvas ("c1","B0Mass");
	TCanvas* c2 = new TCanvas ("c2","Pi1ToKaon");
	TCanvas* c3 = new TCanvas ("c3","Pi2ToKaon");
	TCanvas* c4 = new TCanvas ("c4","PionToProton");
	c1->SetGrid();
	c2->SetGrid();
	c3->SetGrid();
	c4->SetGrid();

	c1->cd();
	B0Mass->Draw("HIST");
	c1->SaveAs("B0Mass.pdf");

	c2->cd();
	Misspi1->Draw("HIST");
	c2->SaveAs("Pi1ToKaon.pdf");

	c3->cd();
	Misspi2->Draw("HIST");
	c3->SaveAs("Pi2ToKaon.pdf");

	c4->cd();
	pitoproton->Draw("HIST");
	c4->SaveAs("PionToProton.pdf");

	output->Close();	
}

