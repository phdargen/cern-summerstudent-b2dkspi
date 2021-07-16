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

 	//Load file	
    TChain* tree=new TChain("DecayTree");
    tree->Add("/eos/lhcb/user/p/phdargen/summerStudents21/Stripped/Data_B2DKspi_LL_11.root");
   //tree->Add("/eos/lhcb/user/p/phdargen/summerStudents21/Stripped/Data_B2DKspi_DD_11.root");

  	//Needed branches (add some more)
        //Int_t KsCat;
        //bool hasRich;
  	Double_t B_MM,D_MM,Ks_MM,B_PT,B_IPCHI2_OWNPV,B_FDCHI2_OWNPV,B_TAU,D_ENDVERTEX_Z, B_ENDVERTEX_Z,D_FDCHI2_ORIVX,D_DIRA_OWNPV,Ks_FDCHI2_ORIVX,Ks_DIRA_OWNPV,Ks_PT,B_ENDVERTEX_CHI2,B_ENDVERTEX_NDOF,B_BKGCAT,ProbNNpi,ProbNNk;
	//tree->SetBranchAddress("KsCat",&KsCat);
        //tree->SetBranchAddress("hasRich",&hasRich);
        tree->SetBranchAddress("B_MM",&B_MM) ;
        tree->SetBranchAddress("D_MM",&D_MM) ;
        tree->SetBranchAddress("Ks_MM",&Ks_MM) ;
        tree->SetBranchAddress("B_PT",&B_PT) ;
        tree->SetBranchAddress("B_IPCHI2_OWNPV",&B_IPCHI2_OWNPV);
        tree->SetBranchAddress("B_FDCHI2_OWNPV",&B_FDCHI2_OWNPV);
        tree->SetBranchAddress("B_TAU",&B_TAU);
        tree->SetBranchAddress("D_ENDVERTEX_Z",&D_ENDVERTEX_Z);
        tree->SetBranchAddress("B_ENDVERTEX_Z",&B_ENDVERTEX_Z);
        tree->SetBranchAddress("D_FDCHI2_ORIVX",&D_FDCHI2_ORIVX) ;
        tree->SetBranchAddress("D_DIRA_OWNPV",&D_DIRA_OWNPV) ;
        tree->SetBranchAddress("Ks_FDCHI2_ORIVX",&Ks_FDCHI2_ORIVX) ;
        tree->SetBranchAddress("Ks_DIRA_OWNPV",&Ks_DIRA_OWNPV) ;
        tree->SetBranchAddress("Ks_PT",&Ks_PT) ;
        tree->SetBranchAddress("B_ENDVERTEX_CHI2",&B_ENDVERTEX_CHI2) ;
        tree->SetBranchAddress("B_ENDVERTEX_NDOF",&B_ENDVERTEX_NDOF) ;
        tree->SetBranchAddress("B_BKGCAT",&B_BKGCAT) ;
        //tree->SetBranchAddress("ProbNNpi",&ProbNNpi) ;
        //tree->SetBranchAddress("ProbNNk",&ProbNNk) ;

	//Create output file
	TFile* output = new TFile("output.root","RECREATE");
    TTree* summary_tree = tree->CloneTree(0);

    //Add new branches
    //...
    summary_tree->Branch("B_MM",&B_MM,"B_MM/F");
    summary_tree->Branch("D_MM",&D_MM,"D_MM/F");
    summary_tree->Branch("Ks_MM",&Ks_MM,"Ks_MM/F");
    //summary_tree->Branch("B_PT",&B_PT,"B_PT/F");
    //summary_tree->Branch("B_IPCHI2_OWNPV",&B_IPCHI2_OWNPV,"B_IPCHI2_OWNPV/F");
    //summary_tree->Branch("B_FDCHI2_OWNPV",&B_FDCHI2_OWNPV,"B_FDCHI2_OWNPV/F");
    //summary_tree->Branch("B_TAU",&B_TAU,"B_TAU/F");
    //summary_tree->Branch("D_ENDVERTEX_Z",&D_ENDVERTEX_Z,"D_ENDVERTEX_Z/F");
    //summary_tree->Branch("B_ENDVERTEX_Z",&B_ENDVERTEX_Z,"B_ENDVERTEX_Z/F");
    //summary_tree->Branch("D_FDCHI2_ORIVX",&D_FDCHI2_ORIVX,"D_FDCHI2_ORIVX/F");
    //summary_tree->Branch("D_DIRA_OWNPV",&D_DIRA_OWNPV,"D_DIRA_OWNPV/F");
    //summary_tree->Branch("Ks_FDCHI2_ORIVX",&Ks_FDCHI2_ORIVX,"Ks_FDCHI2_ORIVX/F");
    //summary_tree->Branch("Ks_DIRA_OWNPV",&Ks_DIRA_OWNPV,"Ks_DIRA_OWNPV/F");
    //summary_tree->Branch("Ks_PT",&Ks_PT,"Ks_PT/F");
    //summary_tree->Branch("B_ENDVERTEX_CHI2",&B_ENDVERTEX_CHI2,"B_ENDVERTEX_CHI2/F");
    //summary_tree->Branch("B_ENDVERTEX_NDOF",&B_ENDVERTEX_NDOF,"B_ENDVERTEX_NDOF/F");
    //summary_tree->Branch("ProbNNpi",&ProbNNpi,"ProbNNpi/F");
    //summary_tree->Branch("ProbNNk",&ProbNNk,"ProbNNk/F");

    //Define some histograms
    TH1F *BMass = new TH1F("B_Mass","B_Mass",100,0,7000);
    TH1F *DMass = new TH1F("D_Mass","D_Mass",100,0,7000);
    TH1F *KsMass = new TH1F("Ks_Mass","Ks_Mass",100,0,7000);

    //Loop over tree
    int nEvents = tree->GetEntries();
  	for ( Int_t j = 0 ; j < nEvents ; j++ ) {
    		tree->GetEntry(j) ;
    		if (0ul == (j % 1000ul)) cout << "Read event " << j << "/" << nEvents << endl;
		
            //Add your cuts
            //
            //Kinematic
            if(!B_PT>2000)continue;
            if(!B_IPCHI2_OWNPV<17)continue;
            if(!B_FDCHI2_OWNPV>150)continue;
            if(!B_TAU>0.00025)continue;
            if(!D_ENDVERTEX_Z - B_ENDVERTEX_Z > 0)continue;
            if(!D_FDCHI2_ORIVX>0)continue;
            if(!D_DIRA_OWNPV>0)continue;
            if(!Ks_FDCHI2_ORIVX>0)continue;
            if(!Ks_DIRA_OWNPV>0)continue;
            if(!Ks_PT>400)continue;
            //Quality
            if(!B_ENDVERTEX_CHI2/B_ENDVERTEX_NDOF<7)continue;
            //if(!ProbNNpi>0.15)continue;
            //if(!ProbNNk>0.23)continue;
            //if(!KsCat==0)continue;
            //if(!hasRich==true)continue;

            //Fill histograms
            BMass->Fill(B_MM);
            DMass->Fill(D_MM);
            KsMass->Fill(Ks_MM);

            //Fill output tree
            summary_tree->Fill();
	}
	
	//TCanvas *myC1 = new TCanvas("myC1","Masses",10,10,800,600);
	
	//myC1->Divide(3,1);
	//myC1->cd(1);
	BMass->Draw();
	//myC1->cd(2);
	//DMass->Draw();
	//myC1->cd(3);
	//KsMass->Draw();

	//myC1->SaveAs("Masses.pdf");	
	
	output->Write();
	output->Close();

}
