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

void selection() {

	TH1::SetDefaultSumw2();
	TH2::SetDefaultSumw2();

 	//Load file	
    TChain* tree=new TChain("DecayTree");
    tree->Add("/eos/lhcb/user/p/phdargen/summerStudents21/Stripped/Data_B2DKspi_LL_11.root");
    //tree->Add("/eos/lhcb/user/p/phdargen/summerStudents21/Stripped/Data_B2DKspi_DD_11.root");

  	//Needed branches (add some more)
  	Double_t B_DTF_MM,B_PT,B_IPCHI2_OWNPV,B_FDCHI2_OWNPV,B_TAU,D_ENDVERTEX_Z,B_ENDVERTEX_Z,D_FDCHI2_ORIVX,D_DIRA_OWNPV,Ks_FDCHI2_ORIVX,Ks_DIRA_OWNPV,Ks_PT,TRACK_GhostProb,B_ENDVERTEX_CHI2;
	Double_t D_MM, Ks_MM;
	Double_t K_D_ProbNNk,K_D_ProbNNpi,pi_ProbNNk,pi_ProbNNpi,pi1_D_ProbNNk,pi1_D_ProbNNpi,pi2_D_ProbNNk,pi2_D_ProbNNpi,pim_Ks_ProbNNk,pim_Ks_ProbNNpi,pip_Ks_ProbNNk,pip_Ks_ProbNNpi;
	Bool_t K_D_hasRich,pi1_D_hasRich,pi2_D_hasRich,pi_hasRich,pim_Ks_hasRich,pip_Ks_hasRich;
	Bool_t B_L0Global_TIS,B_L0HadronDecision_TOS,B_Hlt1TrackAllL0Decision_TOS,B_Hlt1TrackMVADecision_TOS,B_Hlt1TwoTrackMVADecision_TOS,B_Hlt2Topo2BodyBBDTDecision_TOS,B_Hlt2Topo3BodyBBDTDecision_TOS,B_Hlt2Topo4BodyBBDTDecision_TOS,B_Hlt2Topo2BodyDecision_TOS,B_Hlt2Topo3BodyDecision_TOS,B_Hlt2Topo4BodyDecision_TOS;
	Int_t B_ENDVERTEX_NDOF,B_BKGCAT;
	
	tree->SetBranchAddress("B_L0Global_TIS",&B_L0Global_TIS);
	tree->SetBranchAddress("B_L0HadronDecision_TOS",&B_L0HadronDecision_TOS);
	tree->SetBranchAddress("B_Hlt1TrackAllL0Decision_TOS",&B_Hlt1TrackAllL0Decision_TOS);
	tree->SetBranchAddress("B_Hlt1TrackMVADecision_TOS",&B_Hlt1TrackMVADecision_TOS);
	tree->SetBranchAddress("B_Hlt1TwoTrackMVADecision_TOS",&B_Hlt1TwoTrackMVADecision_TOS);
	tree->SetBranchAddress("B_Hlt2Topo2BodyBBDTDecision_TOS",&B_Hlt2Topo2BodyBBDTDecision_TOS);
	tree->SetBranchAddress("B_Hlt2Topo3BodyBBDTDecision_TOS",&B_Hlt2Topo3BodyBBDTDecision_TOS);
	tree->SetBranchAddress("B_Hlt2Topo4BodyBBDTDecision_TOS",&B_Hlt2Topo4BodyBBDTDecision_TOS);
	tree->SetBranchAddress("B_Hlt2Topo2BodyDecision_TOS",&B_Hlt2Topo2BodyDecision_TOS);
	tree->SetBranchAddress("B_Hlt2Topo3BodyDecision_TOS",&B_Hlt2Topo3BodyDecision_TOS);
	tree->SetBranchAddress("B_Hlt2Topo4BodyDecision_TOS",&B_Hlt2Topo4BodyDecision_TOS);
	tree->SetBranchAddress("B_BKGCAT",&B_BKGCAT);
	tree->SetBranchAddress("B_DTF_MM",&B_DTF_MM) ;
    	tree->SetBranchAddress("B_PT",&B_PT) ;
	tree->SetBranchAddress("B_IPCHI2_OWNPV",&B_IPCHI2_OWNPV);
	tree->SetBranchAddress("B_FDCHI2_OWNPV",&B_FDCHI2_OWNPV);
	tree->SetBranchAddress("B_TAU",&B_TAU);
	tree->SetBranchAddress("D_ENDVERTEX_Z",&D_ENDVERTEX_Z);
	tree->SetBranchAddress("B_ENDVERTEX_Z",&B_ENDVERTEX_Z);
	tree->SetBranchAddress("D_FDCHI2_ORIVX",&D_FDCHI2_ORIVX);
	tree->SetBranchAddress("D_DIRA_OWNPV",&D_DIRA_OWNPV);
	tree->SetBranchAddress("Ks_FDCHI2_ORIVX",&Ks_FDCHI2_ORIVX);
	tree->SetBranchAddress("Ks_DIRA_OWNPV",&Ks_DIRA_OWNPV);
	tree->SetBranchAddress("Ks_PT",&Ks_PT);
	tree->SetBranchAddress("B_ENDVERTEX_CHI2",&B_ENDVERTEX_CHI2);
	tree->SetBranchAddress("K_D_ProbNNk",&K_D_ProbNNk);
	tree->SetBranchAddress("K_D_ProbNNpi",&K_D_ProbNNpi);
	tree->SetBranchAddress("pi_ProbNNk",&pi_ProbNNk);
	tree->SetBranchAddress("pi_ProbNNpi",&pi_ProbNNpi);
	tree->SetBranchAddress("pi1_D_ProbNNk",&pi1_D_ProbNNk);
	tree->SetBranchAddress("pi1_D_ProbNNpi",&pi1_D_ProbNNpi);
	tree->SetBranchAddress("pi2_D_ProbNNk",&pi2_D_ProbNNk);
	tree->SetBranchAddress("pi2_D_ProbNNpi",&pi2_D_ProbNNpi);
	tree->SetBranchAddress("pim_Ks_ProbNNk",&pim_Ks_ProbNNk);
	tree->SetBranchAddress("pim_Ks_ProbNNpi",&pim_Ks_ProbNNpi);
	tree->SetBranchAddress("pip_Ks_ProbNNk",&pip_Ks_ProbNNk);
	tree->SetBranchAddress("pip_Ks_ProbNNpi",&pip_Ks_ProbNNpi);
	tree->SetBranchAddress("K_D_hasRich",&K_D_hasRich);
	tree->SetBranchAddress("pi1_D_hasRich",&pi1_D_hasRich);
	tree->SetBranchAddress("pi2_D_hasRich",&pi2_D_hasRich);
	tree->SetBranchAddress("pi_hasRich",&pi_hasRich);
	tree->SetBranchAddress("pim_Ks_hasRich",&pim_Ks_hasRich);
	tree->SetBranchAddress("pip_Ks_hasRich",&pip_Ks_hasRich);
	tree->SetBranchAddress("D_MM",&D_MM);
	tree->SetBranchAddress("Ks_MM",&Ks_MM);
	tree->SetBranchAddress("B_ENDVERTEX_NDOF",&B_ENDVERTEX_NDOF);

	//Create output file
	TFile* output = new TFile("output.root","RECREATE");
    	TTree* summary_tree = tree->CloneTree(0);

    //Add new branches
    //...
	TBranch *BM = summary_tree->Branch("B_DTF_MM",&B_DTF_MM,"B_DTF_MM/F");
	TBranch *DM = summary_tree->Branch("D_MM",&D_MM,"D_MM/F");
	TBranch *KsM = summary_tree->Branch("Ks_MM",&Ks_MM,"Ks_MM/F");
    
    //Define some histograms
    //...
    //before selection
	//particle-mass
	TH1F *BMass0 = new TH1F("BMass0","BMass0",100,4700,6100);
	TH1F *DMass0 = new TH1F("DMass0","DMass0",100,1760,1980);
	TH1F *KsMass0 = new TH1F("KsMass0","KsMass0",100,460,530);
	//PID
	TH1F *KK0 = new TH1F("K_D_ProbNNk0","K_D_ProbNNk",100,0,1);
	TH1F *pipi0 = new TH1F("pi_ProbNNpi0","pi_ProbNNpi",100,0,1);
	TH1F *pi1pi0 = new TH1F("pi1_D_ProbNNpi0","pi1_D_ProbNNpi",100,0,1);
	TH1F *pi2pi0 = new TH1F("pi2_D_ProbNNpi0","pi2_D_ProbNNpi",100,0,1);
	TH1F *pimpi0 = new TH1F("pim_Ks_ProbNNpi0","pim_Ks_ProbNNpi",100,0,1);
	TH1F *pippi0 = new TH1F("pip_Ks_ProbNNpi0","pip_Ks_ProbNNpi",100,0,1);
    //after selection
	//particle-mass
	TH1F *BMass = new TH1F("BMass","BMass",100,4700,6100);
	TH1F *DMass = new TH1F("DMass","DMass",100,1760,1980);
	TH1F *KsMass = new TH1F("KsMass","KsMass",100,460,530);
	//PID
	TH1F *KK = new TH1F("K_D_ProbNNk","K_D_ProbNNk",100,0,1);
	TH1F *pipi = new TH1F("pi_ProbNNpi","pi_ProbNNpi",100,0,1);
	TH1F *pi1pi = new TH1F("pi1_D_ProbNNpi","pi1_D_ProbNNpi",100,0,1);
	TH1F *pi2pi = new TH1F("pi2_D_ProbNNpi","pi2_D_ProbNNpi",100,0,1);
	TH1F *pimpi = new TH1F("pim_Ks_ProbNNpi","pim_Ks_ProbNNpi",100,0,1);
	TH1F *pippi = new TH1F("pip_Ks_ProbNNpi","pip_Ks_ProbNNpi",100,0,1);
    //background
	TH1F *bkg = new TH1F("bkg","bkg",100,4700,6100);

    //Loop over tree
    int nEvents = tree->GetEntries();
  	for ( Int_t j = 0 ; j < nEvents ; j++ ) {
    		tree->GetEntry(j) ;
    		if (0ul == (j % 10000ul)) cout << "Read event " << j << "/" << nEvents << endl;
		
		BMass0->Fill(B_DTF_MM);
		DMass0->Fill(D_MM);
		KsMass0->Fill(Ks_MM);

		KK0->Fill(K_D_ProbNNk);
		pipi0->Fill(pi_ProbNNpi);
		pi1pi0->Fill(pi1_D_ProbNNpi);
		pi2pi0->Fill(pi2_D_ProbNNpi);
		pimpi0->Fill(pim_Ks_ProbNNpi);
		pippi0->Fill(pip_Ks_ProbNNpi);

	    if(B_DTF_MM>5500) 
	    { 
		    bkg->Fill(B_DTF_MM);
		    KK->Fill(K_D_ProbNNk);
		    pipi->Fill(pi_ProbNNpi);
		    pi1pi->Fill(pi1_D_ProbNNpi);
		    pi2pi->Fill(pi2_D_ProbNNpi);
		    pimpi->Fill(pim_Ks_ProbNNpi);
		    pippi->Fill(pip_Ks_ProbNNpi);
	    }

	   
            //Add your cuts
	    if((D_MM - 1869.61 >= 30) || (1869.61-D_MM >= 30))continue;
            if(B_PT <= 2000)continue;
            if(B_IPCHI2_OWNPV >= 20)continue;
            if(B_FDCHI2_OWNPV <= 100)continue;
            if(B_TAU<=0.0001)continue;
            if(D_ENDVERTEX_Z - B_ENDVERTEX_Z <= 0)continue;
            if(D_FDCHI2_ORIVX <= 0)continue;
            if(D_DIRA_OWNPV <= 0)continue;
            if(Ks_FDCHI2_ORIVX <= 0)continue;
            if(Ks_DIRA_OWNPV <= 0)continue;
            if(Ks_PT <= 200)continue;
            if(B_ENDVERTEX_CHI2/B_ENDVERTEX_NDOF >= 10)continue;
            if(K_D_ProbNNk <= 0.3)continue;
            //if (K_D_ProbNNpi <= 0.1) continue;
            //if (pi_ProbNNk <= 0.15) continue;
            if (pi_ProbNNpi <= 0.1) continue;
            //if (pi1_D_ProbNNk <= 0.15)continue;
            if (pi1_D_ProbNNpi <= 0.1) continue;
            //if (pi2_D_ProbNNk <= 0.15) continue;
            if (pi2_D_ProbNNpi <= 0.1) continue;
            //if (pim_Ks_ProbNNk <= 0.15) continue;
            if (pim_Ks_ProbNNpi <= 0.2) continue;
            //if (pip_Ks_ProbNNk <= 0.15) continue;
            if (pip_Ks_ProbNNpi <= 0.2) continue;
            if (K_D_hasRich != true) continue;
            if (pi_hasRich != true) continue;
            if (pi1_D_hasRich != true) continue;
            if (pi2_D_hasRich != true) continue;
            if (pim_Ks_hasRich != true) continue;
            if (pip_Ks_hasRich != true) continue;
 
	    //Trigger Selection
	    if (!(B_L0Global_TIS || B_L0HadronDecision_TOS)) continue;
	    //2011-2012
	    if (!B_Hlt1TrackAllL0Decision_TOS) continue;
	    if (!(B_Hlt2Topo2BodyBBDTDecision_TOS || B_Hlt2Topo3BodyBBDTDecision_TOS || B_Hlt2Topo4BodyBBDTDecision_TOS)) continue;
	    //2015-2018
	    //if (!(B_Hlt1TrackMVADecision_TOS || B_Hlt1TwoTrackMVADecision_TOS)) continue;
	    //if (!(B_Hlt2Topo2BodyDecision_TOS || B_Hlt2Topo3BodyDecision_TOS || B_Hlt2Topo4BodyDecision_TOS)) continue;
	    
	    //MC
	    //if (B_BKGCAT >= 30) continue;
	    
	    //Fill histograms
            //...
	    BMass -> Fill(B_DTF_MM);
	    DMass -> Fill(D_MM);
	    KsMass -> Fill(Ks_MM);
            
            
	    //Fill output tree
            summary_tree->Fill();
	}

            //scale the plots
            float norm = 1000;
            float norm2 = 1;

            BMass0->Scale(norm/BMass0->Integral());
            DMass0->Scale(norm/DMass0->Integral());
            KsMass0->Scale(norm/KsMass0->Integral());

            KK0->Scale(norm2/KK0->Integral());
            pipi0->Scale(norm2/pipi0->Integral());
            pi1pi0->Scale(norm2/pi1pi0->Integral());
            pi2pi0->Scale(norm2/pi2pi0->Integral());
            pimpi0->Scale(norm2/pimpi0->Integral());
            pippi0->Scale(norm2/pippi0->Integral());

            BMass->Scale(norm/BMass->Integral());
            DMass->Scale(norm/DMass->Integral());
            KsMass->Scale(norm/KsMass->Integral());

            KK->Scale(norm2/KK->Integral());
            pipi->Scale(norm2/pipi->Integral());
            pi1pi->Scale(norm2/pi1pi->Integral());
            pi2pi->Scale(norm2/pi2pi->Integral());
            pimpi->Scale(norm2/pimpi->Integral());
            pippi->Scale(norm2/pippi->Integral());

            //color the plots
            BMass->SetLineColor(kRed);
            DMass->SetLineColor(kRed);
            KsMass->SetLineColor(kRed);

            KK->SetLineColor(kRed);
            pipi->SetLineColor(kRed);
            pi1pi->SetLineColor(kRed);
            pi2pi->SetLineColor(kRed);
            pimpi->SetLineColor(kRed);
            pippi->SetLineColor(kRed);

	output->Write();

	TCanvas* c1 = new TCanvas ("c1","B Mass before and after selection");
	TCanvas* c2 = new TCanvas ("c2","D Mass before and after selection");
	TCanvas* c3 = new TCanvas ("c3","Ks Mass before and after selection");
	TCanvas* c4 = new TCanvas ("c4","with combinatorial bkg",1200,600);
	c1->SetGrid();
	c2->SetGrid();
	c3->SetGrid();
	c4->SetGrid();
	
	c1->cd();
	BMass->Draw("HIST");
	BMass0->Draw("HISTsames");
	c1->SaveAs("BMass.pdf");

	c2->cd();
	DMass->Draw("HIST");
	DMass0->Draw("HISTsames");
	c2->SaveAs("DMass.pdf");

	c3->cd();
	KsMass->Draw("HIST");
	KsMass0->Draw("HISTsames");
	c3->SaveAs("KsMass.pdf");
	
	c4->Divide(3,2);
	c4->cd(1);
	KK->Draw("HIST");
	KK0->Draw("HISTsames");
	c4->cd(2);
	pipi->Draw("HIST");
	pipi0->Draw("HISTsames");
	c4->cd(3);
	pi1pi->Draw("HIST");
	pi1pi0->Draw("HISTsames");
	c4->cd(4);
	pi2pi->Draw("HIST");
	pi2pi0->Draw("HISTsames");
	c4->cd(5);
	pimpi->Draw("HIST");
	pimpi0->Draw("HISTsames");
	c4->cd(6);
	pippi->Draw("HIST");
	pippi0->Draw("HISTsames");
	c4->SaveAs("bkg.pdf");
	//three canvas
	output->Close();
}
