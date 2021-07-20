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

int main() {

	 	//P: Load file		
    	TChain* tree=new TChain("DecayTree","RECREATE");
		tree->Add("/home/maria/Work/Data/Data_B2DKspi_DD_11.root");
		//tree->Add("/home/maria/Work/Data//Data_B2DKspi_DD_12.root");
		//tree->Add("/home/maria/Work/Data//Data_B2DKspi_DD_15.root");
		//tree->Add("/home/maria/Work/Data//Data_B2DKspi_DD_16.root");
		//tree->Add("/home/maria/Work/Data//Data_B2DKspi_DD_17.root");
		//tree->Add("/home/maria/Work/Data//Data_B2DKspi_DD_18.root");
		//tree->Add("/home/maria/Work/Data//MC_B2DKspi_DD_12_BdDstKspi.root");
		//tree->Add("/home/maria/Work/Data//MC_B2DKspi_DD_12_BsDstKsK.root");
		//tree->Add("/home/maria/Work/Data//MC_B2DKspi_DD_12_BsDstKspi.root");
		//tree->Add("/home/maria/Work/Data//DMC_B2DKspi_DD_12.root");


	  	//P: Needed branches (add some more) 
  		Double_t B_MM,B_PT,B_IPCHI2_OWNPV,B_FDCHI2_OWNPV,B_TAU,D_ENDVERTEX_Z, B_ENDVERTEX_Z;
  		Double_t D_FDCHI2_ORIVX,D_DIRA_OWNPV,Ks_FDCHI2_ORIVX, D_MM, Ks_MM, ProbNNpi, ProbNNk;
  		Double_t Ks_PT, Ks_DIRA_OWNPV,B_ENDVERTEX_CHI2, TRACK_GhostProb;
  		Int_t B_BKGCAT, B_ENDVERTEX_NDOF, KsCat;
  		bool hasRich; 
  		
  		
		tree->SetBranchAddress("B_MM",&B_MM) ;
    	tree->SetBranchAddress("B_PT",&B_PT) ;
    	tree->SetBranchAddress("B_IPCHI2_OWNPV",&B_IPCHI2_OWNPV) ;
	   	tree->SetBranchAddress("B_FDCHI2_OWNPV",&B_FDCHI2_OWNPV) ;
    	tree->SetBranchAddress("B_TAU",&B_TAU) ;
    	tree->SetBranchAddress("D_ENDVERTEX_Z",&D_ENDVERTEX_Z) ;
    	tree->SetBranchAddress("B_ENDVERTEX_Z",&B_ENDVERTEX_Z) ;
    	tree->SetBranchAddress("D_FDCHI2_ORIVX",&D_FDCHI2_ORIVX) ;
 	  	tree->SetBranchAddress("D_DIRA_OWNPV",&D_DIRA_OWNPV) ;
    	tree->SetBranchAddress("Ks_FDCHI2_ORIVX",&Ks_FDCHI2_ORIVX) ;
    	tree->SetBranchAddress("Ks_PT",&Ks_PT) ;
     	tree->SetBranchAddress("Ks_DIRA_OWNPV",&Ks_DIRA_OWNPV) ;
 		tree->SetBranchAddress("B_ENDVERTEX_CHI2",&B_ENDVERTEX_CHI2) ;
 		tree->SetBranchAddress("B_ENDVERTEX_NDOF",&B_ENDVERTEX_NDOF) ;
 		tree->SetBranchAddress("KsCat",&KsCat) ;
 		tree->SetBranchAddress("B_BKGCAT",&B_BKGCAT) ;


		//P: Create output file 
		TFile* output = new TFile("output.root","RECREATE");
  		TTree* summary_tree = tree->CloneTree(0);

    	//P: Add new branches
    	summary_tree->Branch("B_MM",&B_MM);
    	summary_tree->Branch("D_MM",&D_MM);
    	summary_tree->Branch("Ks_MM",&Ks_MM);

    
		//P: Define some histograms 
    	TH1F *hB_MM = new TH1F("B_MM","B_MM",100,4700,6100);
   		TH1F *hD_MM = new TH1F("D_MM","D_MM",100,0,2000);
    	TH1F *hKs_MM = new TH1F("Ks_MM","Ks_MM",100,0,2000);
    	
    	TH1F *hB_MM2 = new TH1F("B_MM","B_MM",100,4700,7000);

    
    	//P: Loop over tree
    	int nEvents = tree->GetEntries();
  		for ( Int_t j = 0 ; j < nEvents ; j++ ) {
    		tree->GetEntry(j) ;
    		if (0ul == (j % 2000ul)) cout << "Read event " << j << "/" << nEvents << endl;
		
		    //P: Add your cuts
		    //M: Ks category
		    if (KsCat != 1) continue;
		    
		    //M: Kinematic		    
            if (B_PT<=2000) continue;
            if (B_IPCHI2_OWNPV>=20) continue;
			if (B_FDCHI2_OWNPV<=200) continue;
		    if (B_TAU<=0.0001) continue;
		   	if (D_ENDVERTEX_Z - B_ENDVERTEX_Z <= 0) continue;
		   	if (D_FDCHI2_ORIVX <= 0) continue;
		   	if (D_DIRA_OWNPV <= 0) continue;
		   	if (Ks_FDCHI2_ORIVX <= 0) continue; 
		   	if (Ks_DIRA_OWNPV <= 0) continue; 
		   	if (Ks_PT <= 200) continue;
		   	
		    //M: Quality
		   	if (TRACK_GhostProb >= 0.5) continue;
		    if ((B_ENDVERTEX_CHI2/B_ENDVERTEX_NDOF) >= 10) continue;
		    //if (B_BKGCAT >= 30) continue; // Only for MC!
		    		    
		    //M: PID
		    //if (ProbNNpi <= 0.2) continue;
		    //if (ProbNNk <= 0.3) continue;
		    //if (hasRich != true) continue;
		    
		    		       
		    //P: Fill histograms
            hB_MM ->Fill(B_MM);
            hB_MM -> SetFillColor(kBlue-7);
            
            hD_MM->Fill(D_MM);
            hB_MM -> SetFillColor(kBlue-7);
            
            hKs_MM->Fill(Ks_MM);
            hB_MM -> SetFillColor(kBlue-7);

        
            //Fill output tree
            summary_tree->Fill();
		}

		cout << "before selection: " << tree->GetEntries() << endl;
		cout << "after selection: " << summary_tree->GetEntries() << endl;

		output->Write();
				
		//M: Canvas
		
		TCanvas *cB_MM = new TCanvas("cB_MM","cB_MM");
		tree->Draw("B_MM");
		summary_tree -> SetLineColor(kRed-5);
		summary_tree->Draw("B_MM");
		
		TCanvas *cD_MM = new TCanvas("cD_MM","cD_MM");
		tree->Draw("D_MM");
		summary_tree -> SetLineColor(kGreen-5);
		summary_tree->Draw("D_MM");
		
		TCanvas *cKs_MM = new TCanvas("cKs_MM","cKs_MM");
		tree->Draw("Ks_MM");
		summary_tree -> SetLineColor(kGreen-5);
		summary_tree->Draw("Ks_MM");

		output -> Close();

return 0;
}
