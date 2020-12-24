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

const int nqbins=3;
const int nvbins=3;
const int default_binning = 100;

const float qmin= 1.;
const float qmax= 4.1; 
const float vmin= 2.2;
const float vmax= 4.25;

const float q_ave[3]={1.18, 1.54,2.23};
const float q_values[4]= {1.0, 1.33, 1.76, 4.1 }; 
const float v_values[4]= {2.2, 3.2, 3.73, 4.25 };
const float center_q[3] = {1.165 , 1.55, 2.93 };
const float center_v[3] = {2.7, 3.46, 3.99 };

Float_t Eprime[nqbins][nvbins], Theta[nqbins][nvbins];
Float_t delta_pb[nqbins][nvbins];

Float_t sigma_Born_Pb[nqbins][nvbins], sigma_Rad_Pb[nqbins][nvbins];
Float_t sigma_BornIn_Pb[nqbins][nvbins], sigma_BornQE_Pb[nqbins][nvbins];
Float_t sigma_RadEl_Pb[nqbins][nvbins], sigma_RadQE_Pb[nqbins][nvbins], sigma_RadDIS_Pb[nqbins][nvbins];
Float_t C_corr_Pb[nqbins][nvbins];

Float_t RC_D[nqbins][nvbins],RC_C[nqbins][nvbins],RC_Fe[nqbins][nvbins],RC_Pb[nqbins][nvbins];
Float_t ratioRC_CD[nqbins][nvbins], ratioRC_FeD[nqbins][nvbins], ratioRC_PbD[nqbins][nvbins];

TGraphErrors *g_RC_D[nqbins],*g_RC_C[nqbins],*g_RC_Fe[nqbins],*g_RC_Pb[nqbins];
TGraphErrors *g_ratioRC_CD[nqbins], *g_ratioRC_FeD[nqbins], *g_ratioRC_PbD[nqbins];

void rc(void){
 gROOT->SetStyle("Plain");
 gStyle->SetPalette(1);
 gStyle->SetOptFit(1);
 gStyle->SetOptStat(0);
 gStyle->SetOptTitle(1);
 gStyle->SetOptDate(0);

 ////////////////////// READ RC ////////////////////////
Float_t Ebeam, x, Q2; 
 ifstream in;
 //in.open("clasd2.out");
 in.open("clasPb208.out");
 for(int nq=0; nq<nqbins;nq++){
  for(int nv=0;  nv<nvbins; nv++){  
   in >> Ebeam >> Eprime[nq][nv] >> Theta[nq][nv] >> x >> Q2 >> sigma_Born_Pb[nq][nv] >> sigma_BornIn_Pb[nq][nv] >> sigma_BornQE_Pb[nq][nv]
      >> sigma_Rad_Pb[nq][nv] >>  sigma_RadEl_Pb[nq][nv] >> sigma_RadQE_Pb[nq][nv] >> sigma_RadDIS_Pb[nq][nv] >> C_corr_Pb[nq][nv];
   
 cout << Ebeam << " " << Eprime[nq][nv]<< " " << Theta[nq][nv]<< " "<< x<< " " << Q2 << " " << sigma_Born_Pb[nq][nv]<< " " 
	   << sigma_BornIn_Pb[nq][nv]<< " " << sigma_BornQE_Pb[nq][nv]<< " " << sigma_Rad_Pb[nq][nv]<< " " 
	   << sigma_RadEl_Pb[nq][nv]<< " "<< sigma_RadQE_Pb[nq][nv]
	   << " " << sigma_RadDIS_Pb[nq][nv]<< " " <<C_corr_Pb[nq][nv]<< endl;
   }
  }
 in.close();
 
 // E.Borie (1971) ******************************** delta = pi * alpha * Z * sin(theta/2)/(1+sin(theta/2)) *********************************
 float const alpha = 1/137.036;
 int const z_pb = 82;
 int const z_d = 2;
 for(int nq=0; nq<nqbins;nq++){
  for(int nv=0;  nv<nvbins; nv++){  
   delta_pb[nq][nv] =100* TMath::Pi() * alpha * z_pb * TMath::Sin(Theta[nq][nv]/2 *0.0174533) / (1+TMath::Sin(Theta[nq][nv]/2*0.0174533));
       cout << nq << " " << nv << "theta=" << Theta[nq][nv]<< " " << " delta(Pb)=" << delta_pb[nq][nv]<<endl;
  }
 }

////////////////////////////////// WRITE OUT ////////////////////////////////////// 
/*
 ofstream out;
out.open("eRC_qv_all.txt");
for(int nq=0;  nq<nqbins; nq++){
 for(int nv=0; nv<nvbins;nv++){
  out << nq << " " << nv << " " <<  RC_D[nq][nv] << " " << RC_C[nq][nv] << " " << RC_Fe[nq][nv] << " " << RC_Pb[nq][nv] << endl;
 }
}
out.close();
ofstream out;
out.open("eRC_v_all.txt");
 for(int nv=0; nv<nvbins;nv++){
  out << nv << " " <<  RC_vD[nv] << " " << RC_vC[nv] << " " << RC_vFe[nv] << " " << RC_vPb[nv] << endl;
 }
out.close();
*/

 ////////////////////////// Plot ////////////////////////

 cout << "Plotting" <<endl;  
 gStyle->SetOptTitle(0); 
 
 TCanvas *c = new TCanvas("c","c", 500, 400);
 c->SetGridy();
 TH2F *h = new TH2F("h","h",default_binning, 15, 45, default_binning, 20, 65);
 h->GetXaxis()->SetTitle("#theta, deg");
 //h->GetYaxis()->SetTitle("#delta_{D}, %");
 h->GetYaxis()->SetTitle("#delta_{Pb}, %");
 h->GetYaxis()->SetTitleOffset(0.85);
 h->GetXaxis()->SetTitleOffset(0.9);
 h->GetXaxis()->SetTitleSize(0.055);
 h->GetYaxis()->SetTitleSize(0.055);
 h->GetYaxis()->SetLabelSize(0.052);
 h->GetXaxis()->SetLabelSize(0.052);
 h->Draw();
 for(int nq=0; nq<nqbins;nq++){
  g_RC_Pb[nq] = new TGraphErrors(nvbins, Theta[nq], delta_pb[nq],0, 0);
  if(nq==0)g_RC_Pb[nq]->SetName("Pb0");
  if(nq==1)g_RC_Pb[nq]->SetName("Pb1");
  if(nq==2)g_RC_Pb[nq]->SetName("Pb2");
  
  g_RC_Pb[nq]->SetMarkerColor(nq+1);
  g_RC_Pb[nq]->SetMarkerStyle(20);
  g_RC_Pb[nq]->SetMarkerSize(1.7);
  g_RC_Pb[nq]->Draw("Psame");
 
  TLatex a;
  a.SetNDC();
  a.SetTextSize(0.06);
  //a.DrawLatex(0.2,0.92, Form("%3.2f<Q^{2}<%3.2f and <Q^{2}>=%3.2f", q_values[nq], q_values[nq+1],
//		  q_ave[nq]));
  auto legend = new TLegend(0.7,0.95,0.92,0.75);
  legend->AddEntry("Pb0","1.0<Q^{2}<1.33","p");
  legend->AddEntry("Pb1","1.33<Q^{2}<1.76","p");
  legend->AddEntry("Pb2","1.76<Q^{2}<4.2","p");
	        legend->Draw();
 }
 c->Print("~/secure/ExternalsRC/systemPb_SecondBorn.png");
 c->Clear();
 
}
