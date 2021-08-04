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

void veto() {

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
	Double_t pi_PX, pi_PY, pi_PZ;
	Double_t pi_ProbNNpi,pi1_D_ProbNNpi,pi2_D_ProbNNpi;

	tree->SetBranchAddress("pi_ProbNNpi",&pi_ProbNNpi);
	tree->SetBranchAddress("pi1_D_ProbNNpi",&pi1_D_ProbNNpi);
	tree->SetBranchAddress("pi2_D_ProbNNpi",&pi2_D_ProbNNpi);
	tree->SetBranchAddress("D_PX",&D_PX);
	tree->SetBranchAddress("D_PY",&D_PY);
	tree->SetBranchAddress("D_PZ",&D_PZ);
	tree->SetBranchAddress("D_MM",&D_MM);
	tree->SetBranchAddress("pi_PX",&pi_PX);
	tree->SetBranchAddress("pi_PY",&pi_PY);
	tree->SetBranchAddress("pi_PZ",&pi_PZ);
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
    
	TH1F *B0Mass = new TH1F("B0Mass","B0Mass",100,5100,5350);
	TH1F *B0Mass1 = new TH1F("B0Mass1","B0Mass1",100,5100,5350);
	TH1F *Misspi1 = new TH1F("Misspi1","Misspi1",100,1769,2169);
	TH1F *Misspi11 = new TH1F("Misspi11","Misspi11",100,1769,2169);
	TH1F *Misspi2 = new TH1F("Misspi2","Misspi2",100,1769,2169);
	TH1F *Misspi21 = new TH1F("Misspi21","Misspi21",100,1769,2169);
	TH1F *pitoproton = new TH1F("pitoproton","pitoproton",100,1800,7000);
	TH1F *pitoproton1 = new TH1F("pitoproton1","pitoproton1",100,1800,7000);

	int nEvents = tree->GetEntries();
	for ( Int_t j = 0 ; j < nEvents ; j++ ) {
		tree->GetEntry(j) ;
		if (0ul == (j % 10000ul)) cout << "Read event " << j << "/" << nEvents << endl;
		TLorentzVector D_D, pi1_D, pi1_D_asK_MissID, K_D, pi2_D, pi2_D_asK_MissID, pi1_D_asP_MissID, pi_D;
		D_D.SetXYZM(D_PX,D_PY,D_PZ,D_MM);
		pi_D.SetXYZM(pi_PX, pi_PY, pi_PZ, massPion);
		pi1_D.SetXYZM(pi1_D_PX,pi1_D_PY,pi1_D_PZ,massPion);
		pi1_D_asK_MissID.SetXYZM(pi1_D_PX,pi1_D_PY,pi1_D_PZ,massKaon);
		K_D.SetXYZM(K_D_PX,K_D_PY,K_D_PZ,massKaon);
		pi2_D.SetXYZM(pi2_D_PX,pi2_D_PY,pi2_D_PZ,massPion);
		pi2_D_asK_MissID.SetXYZM(pi2_D_PX,pi2_D_PY,pi2_D_PZ,massKaon);
		pi1_D_asP_MissID.SetXYZM(pi1_D_PX,pi1_D_PY,pi1_D_PZ,massProton);	
		//after cut
		Double_t missB = (D_D + pi_D).M();
		Double_t misspi1 = (K_D + pi1_D_asK_MissID + pi2_D).M();
		Double_t misspi2 = (K_D + pi2_D_asK_MissID + pi1_D).M();
		Double_t missp = (K_D + pi1_D_asP_MissID + pi2_D).M();

		float mass_Ds = 1968.35;//MeV
		float mass_Lambda_c = 2286.46;//MeV
		float mass_B0 = 5279.65;//MeV

		if((abs(misspi1 - mass_Ds) > 50) || (pi1_D_ProbNNpi > 0.7)){
			Misspi11->Fill(misspi1);
		}
		if((abs(misspi2 - mass_Ds) > 50) || (pi2_D_ProbNNpi > 0.7)){
			Misspi21->Fill(misspi2);
		}
		if((abs(missB - mass_B0) > 50) || (pi_ProbNNpi > 0.7)){
			B0Mass1->Fill(missB);
		}
		if((abs(missp - mass_Lambda_c) > 50) || (pi1_D_ProbNNpi > 0.7)){
			pitoproton1->Fill(missp);
		}

		//before cut
		B0Mass->Fill((D_D + pi_D).M());
		Misspi1->Fill((K_D + pi1_D_asK_MissID + pi2_D).M());
		Misspi2->Fill((K_D + pi2_D_asK_MissID + pi1_D).M());
		pitoproton->Fill((K_D + pi1_D_asP_MissID + pi2_D).M());

	}
	B0Mass->Draw("HIST");	
	//output->Write();

	float norm = 1000;
	B0Mass->Scale(norm/B0Mass->Integral());
	Misspi1->Scale(norm/Misspi1->Integral());
	Misspi2->Scale(norm/Misspi2->Integral());
	pitoproton->Scale(norm/pitoproton->Integral());
	B0Mass1->Scale(norm/B0Mass1->Integral());
	Misspi11->Scale(norm/Misspi11->Integral());
	Misspi21->Scale(norm/Misspi21->Integral());
	pitoproton1->Scale(norm/pitoproton1->Integral());

	TCanvas* c1 = new TCanvas ("c1","B0Mass");
	TCanvas* c2 = new TCanvas ("c2","Pi1ToKaon");
	TCanvas* c3 = new TCanvas ("c3","Pi2ToKaon");
	TCanvas* c4 = new TCanvas ("c4","PionToProton");
	c1->SetGrid();
	c2->SetGrid();
	c3->SetGrid();
	c4->SetGrid();

	c1->cd();
	B0Mass1->SetLineColor(kRed);
	B0Mass1->Draw("HIST");
	B0Mass->Draw("HISTsames");
	TLegend leg1(0.7, 0.7, 0.9, 0.9);
	leg1.AddEntry(B0Mass,"B0 Mass before cut","L");
	leg1.AddEntry(B0Mass1,"B0 Mass after cut","L");
	leg1.DrawClone();
	c1->SaveAs("B0Mass.pdf");

	c2->cd();
	Misspi11->SetLineColor(kRed);
	Misspi11->Draw("HIST");
	Misspi1->Draw("HISTsames");
	TLegend leg2(0.7, 0.7, 0.9, 0.9);
	leg2.AddEntry(Misspi1,"D Mass before cut","L");
	leg2.AddEntry(Misspi11,"D Mass after cut","L");
	leg2.DrawClone();
	c2->SaveAs("Pi1ToKaon.pdf");

	c3->cd();
	Misspi21->SetLineColor(kRed);
	Misspi21->Draw("HIST");
	Misspi2->Draw("HISTsames");
	TLegend leg3(0.7, 0.7, 0.9, 0.9);
	leg3.AddEntry(Misspi2,"D Mass before cut","L");
	leg3.AddEntry(Misspi21,"D Mass after cut","L");
	leg3.DrawClone();
	c3->SaveAs("Pi2ToKaon.pdf");

	c4->cd();
	pitoproton->Draw("HIST");
	pitoproton1->SetLineColor(kRed);
	pitoproton1->Draw("HISTsames");
	TLegend leg4(0.7, 0.7, 0.9, 0.9);
	leg4.AddEntry(pitoproton,"D Mass before cut","L");
	leg4.AddEntry(pitoproton1,"D Mass after cut","L");
	leg4.DrawClone();
	c4->SaveAs("PionToProton.pdf");

	//output->Close();	
}


