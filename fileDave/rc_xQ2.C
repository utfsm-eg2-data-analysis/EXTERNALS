#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TVirtualIndex.h>
#include <TFile.h>
#include <TCanvas.h>
#include <iostream>             // std::cout, std::endl
#include <fstream>              // std::ifstream
#include <sstream>
#include <iomanip>
#include <stdlib.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TSystem.h>
#include <TH2F.h>
#include <TF1.h>
#include <TStyle.h>
#include <TPostScript.h>
#include <TPaveText.h>
#include <TGraphErrors.h>
#include <TLatex.h>

using namespace std;


//Float_t Eprime[nqbins][nvbins], Theta[nqbins][nvbins], Xb[nqbins][nvbins];

Float_t sigma_Born_D[308], sigma_Rad_D[308];
Float_t sigma_BornIn_D[308], sigma_BornQE_D[308];
Float_t sigma_RadEl_D[308], sigma_RadQE_D[308], sigma_RadDIS_D[308];
Float_t C_corr_D[308];

Float_t sigma_Born_C[308], sigma_Rad_C[308];
Float_t sigma_BornIn_C[308], sigma_BornQE_C[308];
Float_t sigma_RadEl_C[308], sigma_RadQE_C[308], sigma_RadDIS_C[308];
Float_t C_corr_C[308];

Float_t sigma_Born_Fe[308], sigma_Rad_Fe[308];
Float_t sigma_BornIn_Fe[308], sigma_BornQE_Fe[308];
Float_t sigma_RadEl_Fe[308], sigma_RadQE_Fe[308], sigma_RadDIS_Fe[308];
Float_t C_corr_Fe[308];

Float_t sigma_Born_Pb[308], sigma_Rad_Pb[308];
Float_t sigma_BornIn_Pb[308], sigma_BornQE_Pb[308];
Float_t sigma_RadEl_Pb[308], sigma_RadQE_Pb[308], sigma_RadDIS_Pb[308];
Float_t C_corr_Pb[308];

//Float_t RC_D[nqbins][nvbins],RC_C[nqbins][nvbins],RC_Fe[nqbins][nvbins],RC_Pb[nqbins][nvbins];
//Float_t ratioRC_CD[nqbins][nvbins], ratioRC_FeD[nqbins][nvbins], ratioRC_PbD[nqbins][nvbins];

void rc_xQ2(void){
 gROOT->SetStyle("Plain");
 gStyle->SetPalette(1);
 gStyle->SetOptFit(1);
 gStyle->SetOptStat(0);
 gStyle->SetOptTitle(1);
 gStyle->SetOptDate(0);

 cout << " READING " <<endl;
 ////////////////////// READ RC ////////////////////////
Float_t temp1, Ebeam, x, Q2, Eprime,Theta; 
 ifstream in;
 
 in.open("5.766_40.0_12.out");
 for(int n=0; n<308;n++){
      in >> Ebeam
      >> Eprime
      >> Theta
      >> x 
      >> Q2
      >> sigma_Born_D[n]
      >> sigma_BornIn_D[n] 
      >> sigma_BornQE_D[n] 
      >> temp1
      >> sigma_Rad_D[n]
      >> sigma_RadEl_D[n]
      >> sigma_RadQE_D[n] 
      >> sigma_RadDIS_D[n]
      >> C_corr_D[n];
	 if(sigma_Rad_D[n]>0) {
   cout << " n=" << n << " Eprime=" << Eprime <<" Theta="<< Theta<< " Q2=" << Q2 
        << " x="<< x<<  " nu="<< Q2/(2*0.938*x)  
         <<" RC=" << sigma_Born_D[n] / sigma_Rad_D[n] <<endl;
	 }
	 }
 in.close();
}
