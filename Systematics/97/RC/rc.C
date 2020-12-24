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

const int nqbins=9;
const int nvbins=7;

const int default_binning = 100;

const float qmin= 1.;
const float qmax= 4.1; 
const float vmin= 2.2;
const float vmax= 4.25;

const float center_q[9] = {1.055, 1.155, 1.275, 1.4, 1.55, 1.685,  2.03, 2.6, 3.5};
//const float center_v[9] = {2.5, 2.9, 3.1, 3.285,3.46,  3.64, 3.815, 3.985,  4.155};
//const float center_v[7] = {2.5, 2.9, 3.1, 3.285,3.46,  3.64, 3.99};
//const float q_values[10]= {1.0, 1.11, 1.22, 1.33, 1.47, 1.61, 1.76, 2.3, 2.9, 4.1};
//const float v_values[10]= {2.2, 2.8, 3.0, 3.2 , 3.37, 3.55, 3.73, 3.90, 4.07, 4.25 };
const float v_values[8]= {2.2, 2.8, 3.0, 3.2 , 3.37, 3.55, 3.73, 4.25 };

const float q_values[4]= {1.0, 1.33, 1.76, 4.1 }; 
//const float v_values[4]= {2.2, 3.2, 3.73, 4.25 };
//const float center_q[3] = {1.165 , 1.55, 2.93 };
const float center_v[3] = {2.7, 3.46, 3.99 };

////////////// (9,9) ///////////////
Float_t Eprime[nqbins][nvbins], Theta[nqbins][nvbins];

Float_t sigma_Born_D[nqbins][nvbins], sigma_Rad_D[nqbins][nvbins];
Float_t sigma_BornIn_D[nqbins][nvbins], sigma_BornQE_D[nqbins][nvbins];
Float_t sigma_RadEl_D[nqbins][nvbins], sigma_RadQE_D[nqbins][nvbins], sigma_RadDIS_D[nqbins][nvbins];
Float_t C_corr_D[nqbins][nvbins];

Float_t sigma_Born_C[nqbins][nvbins], sigma_Rad_C[nqbins][nvbins];
Float_t sigma_BornIn_C[nqbins][nvbins], sigma_BornQE_C[nqbins][nvbins];
Float_t sigma_RadEl_C[nqbins][nvbins], sigma_RadQE_C[nqbins][nvbins], sigma_RadDIS_C[nqbins][nvbins];
Float_t C_corr_C[nqbins][nvbins];

Float_t sigma_Born_Fe[nqbins][nvbins], sigma_Rad_Fe[nqbins][nvbins];
Float_t sigma_BornIn_Fe[nqbins][nvbins], sigma_BornQE_Fe[nqbins][nvbins];
Float_t sigma_RadEl_Fe[nqbins][nvbins], sigma_RadQE_Fe[nqbins][nvbins], sigma_RadDIS_Fe[nqbins][nvbins];
Float_t C_corr_Fe[nqbins][nvbins];

Float_t sigma_Born_Pb[nqbins][nvbins], sigma_Rad_Pb[nqbins][nvbins];
Float_t sigma_BornIn_Pb[nqbins][nvbins], sigma_BornQE_Pb[nqbins][nvbins];
Float_t sigma_RadEl_Pb[nqbins][nvbins], sigma_RadQE_Pb[nqbins][nvbins], sigma_RadDIS_Pb[nqbins][nvbins];
Float_t C_corr_Pb[nqbins][nvbins];

////////////// (3,3)///////////////////
Float_t sigmaSum_Born_D[3][3], sigmaSum_Rad_D[3][3];
//Float_t C_corr_D[3][3];

Float_t sigmaSum_Born_C[3][3], sigmaSum_Rad_C[3][3];
//Float_t C_corr_C[3][3];

Float_t sigmaSum_Born_Fe[3][3], sigmaSum_Rad_Fe[3][3];
//Float_t C_corr_Fe[3][3];

Float_t sigmaSum_Born_Pb[3][3], sigmaSum_Rad_Pb[3][3];
//Float_t C_corr_Pb[3][3];

Float_t rc_D[3][3],rc_C[3][3],rc_Fe[3][3],rc_Pb[3][3];
Float_t ratioRC_CD[3][3], ratioRC_FeD[3][3], ratioRC_PbD[3][3];

TGraphErrors *g_RC_D[3],*g_RC_C[3],*g_RC_Fe[3],*g_RC_Pb[3];
TGraphErrors *g_ratioRC_CD[3], *g_ratioRC_FeD[3], *g_ratioRC_PbD[3];

TGraphErrors *g_RC_vD,*g_RC_vC,*g_RC_vFe,*g_RC_vPb;
TGraphErrors *g_ratioRC_vCD, *g_ratioRC_vFeD, *g_ratioRC_vPbD;

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

 cout << " **************** D ********************** " << endl; 
 in.open("clasd2.out");
 for(int nq=0; nq<nqbins;nq++){
  for(int nv=0;  nv<nvbins; nv++){  
   in >> Ebeam
      >> Eprime[nq][nv]
      >> Theta[nq][nv]
      >> x 
      >> Q2
      >> sigma_Born_D[nq][nv]
      >> sigma_BornIn_D[nq][nv] 
      >> sigma_BornQE_D[nq][nv] 
      >> sigma_Rad_D[nq][nv]
      >> sigma_RadEl_D[nq][nv]
      >> sigma_RadQE_D[nq][nv] 
      >> sigma_RadDIS_D[nq][nv]
      >> C_corr_D[nq][nv];
   
   cout << Ebeam << " " << Eprime[nq][nv]<< " " << Theta[nq][nv]<< " "<< x<< " " << Q2 << " " << sigma_Born_D[nq][nv]<< " " 
	   << sigma_BornIn_D[nq][nv]<< " " << sigma_BornQE_D[nq][nv]<< " " << sigma_Rad_D[nq][nv]<< " " 
	   << sigma_RadEl_D[nq][nv]<< " "<< sigma_RadQE_D[nq][nv]
	   << " " << sigma_RadDIS_D[nq][nv]<< " " <<C_corr_D[nq][nv]<< endl;
   }
  }
 in.close();

 cout << " **************** C ********************** " << endl; 
 in.open("clasC12.out");
 for(int nq=0; nq<nqbins;nq++){
  for(int nv=0;  nv<nvbins; nv++){  
   in >> Ebeam >> Eprime[nq][nv] >> Theta[nq][nv] >> x >> Q2 >> sigma_Born_C[nq][nv] >> sigma_BornIn_C[nq][nv] >> sigma_BornQE_C[nq][nv]
      >> sigma_Rad_C[nq][nv] >>  sigma_RadEl_C[nq][nv] >> sigma_RadQE_C[nq][nv] >> sigma_RadDIS_C[nq][nv] >> C_corr_C[nq][nv];
   
//   cout << Ebeam << " " << Eprime[nq][nv]<< " " << Theta[nq][nv]<< " "<< x<< " " << Q2 << " " << sigma_Born_C[nq][nv]<< " " 
//	   << sigma_BornIn_C[nq][nv]<< " " << sigma_BornQE_C[nq][nv]<< " " << sigma_Rad_C[nq][nv]<< " " 
//	   << sigma_RadEl_C[nq][nv]<< " "<< sigma_RadQE_C[nq][nv]
//	   << " " << sigma_RadDIS_C[nq][nv]<< " " <<C_corr_C[nq][nv]<< endl;
   }
  }
 in.close();

 cout << " **************** Fe ********************** " << endl; 
 in.open("clasFe56.out");
 for(int nq=0; nq<nqbins;nq++){
  for(int nv=0;  nv<nvbins; nv++){  
   in >> Ebeam >> Eprime[nq][nv] >> Theta[nq][nv] >> x >> Q2 >> sigma_Born_Fe[nq][nv] >> sigma_BornIn_Fe[nq][nv] >> sigma_BornQE_Fe[nq][nv]
      >> sigma_Rad_Fe[nq][nv] >>  sigma_RadEl_Fe[nq][nv] >> sigma_RadQE_Fe[nq][nv] >> sigma_RadDIS_Fe[nq][nv] >> C_corr_Fe[nq][nv];
   
 /*cout << Ebeam << " " << Eprime[nq][nv]<< " " << Theta[nq][nv]<< " "<< x<< " " << Q2 << " " << sigma_Born_Fe[nq][nv]<< " " 
	   << sigma_BornIn_Fe[nq][nv]<< " " << sigma_BornQE_Fe[nq][nv]<< " " << sigma_Rad_Fe[nq][nv]<< " " 
	   << sigma_RadEl_Fe[nq][nv]<< " "<< sigma_RadQE_Fe[nq][nv]
	   << " " << sigma_RadDIS_Fe[nq][nv]<< " " <<C_corr_Fe[nq][nv]<< endl;
   */
   }
  }
 in.close();


 cout << " **************** Pb ********************** " << endl; 
 in.open("clasPb208.out");
 for(int nq=0; nq<nqbins;nq++){
  for(int nv=0;  nv<nvbins; nv++){  
   in >> Ebeam >> Eprime[nq][nv] >> Theta[nq][nv] >> x >> Q2 >> sigma_Born_Pb[nq][nv] >> sigma_BornIn_Pb[nq][nv] >> sigma_BornQE_Pb[nq][nv] >> sigma_Rad_Pb[nq][nv] >>  sigma_RadEl_Pb[nq][nv] >> sigma_RadQE_Pb[nq][nv] >> sigma_RadDIS_Pb[nq][nv] >> C_corr_Pb[nq][nv];
   
 /*cout << Ebeam << " " << Eprime[nq][nv]<< " " << Theta[nq][nv]<< " "<< x<< " " << Q2 << " " << sigma_Born_Pb[nq][nv]<< " " 
	   << sigma_BornIn_Pb[nq][nv]<< " " << sigma_BornQE_Pb[nq][nv]<< " " << sigma_Rad_Pb[nq][nv]<< " " 
	   << sigma_RadEl_Pb[nq][nv]<< " "<< sigma_RadQE_Pb[nq][nv]
	   << " " << sigma_RadDIS_Pb[nq][nv]<< " " <<C_corr_Pb[nq][nv]<< endl;
   */
   }
  }
 in.close();

///////////////////////  Calculate //////////////// 
for(int nv=0; nv<3; nv++){
 for(int nq=0; nq<nqbins/3; nq++){
  sigmaSum_Born_D[nq][nv]=0;sigmaSum_Rad_D[nq][nv]=0; 
  sigmaSum_Born_C[nq][nv]=0;sigmaSum_Rad_C[nq][nv]=0; 
  sigmaSum_Born_Fe[nq][nv]=0;sigmaSum_Rad_Fe[nq][nv]=0; 
  sigmaSum_Born_Pb[nq][nv]=0;sigmaSum_Rad_Pb[nq][nv]=0; 
 }
} 
int i, j;
for(int nq=0; nq<nqbins; nq++){
 for(int nv=0; nv<nvbins; nv++){
  i=nq/3;
  if(nv==6) j=2;
  else j=nv/3; 
  sigmaSum_Born_D[i][j] += sigma_Born_D[nq][nv]; 
  sigmaSum_Rad_D[i][j] += sigma_Rad_D[nq][nv];
  sigmaSum_Born_C[i][j] += sigma_Born_C[nq][nv]; 
  sigmaSum_Rad_C[i][j] += sigma_Rad_C[nq][nv];
  sigmaSum_Born_Fe[i][j] += sigma_Born_Fe[nq][nv]; 
  sigmaSum_Rad_Fe[i][j] += sigma_Rad_Fe[nq][nv];
  sigmaSum_Born_Pb[i][j] += sigma_Born_Pb[nq][nv]; 
  sigmaSum_Rad_Pb[i][j] += sigma_Rad_Pb[nq][nv];
  }
}
 for(int i=0; i<3; i++){
  for(int j=0; j<3; j++){ 
   rc_D[i][j] = sigmaSum_Rad_D[i][j]/sigmaSum_Born_D[i][j];
   rc_C[i][j] = sigmaSum_Rad_C[i][j]/sigmaSum_Born_C[i][j];
   rc_Fe[i][j] = sigmaSum_Rad_Fe[i][j]/sigmaSum_Born_Fe[i][j];
   rc_Pb[i][j] = sigmaSum_Rad_Pb[i][j]/sigmaSum_Born_Pb[i][j];
 
   cout << "nq="<< i << " nv=" << j<<" D = " << rc_D[i][j]<<" C = "<<rc_C[i][j]<<" Fe = "<<rc_Fe[i][j]<<" Pb = "<< rc_Pb[i][j]<<endl; 
  }
 }

 for(int i=0; i<3; i++){
  for(int j=0; j<3; j++){ 
   ratioRC_CD[i][j]=rc_C[i][j]/rc_D[i][j];
   ratioRC_FeD[i][j]=rc_Fe[i][j]/rc_D[i][j];
   ratioRC_PbD[i][j]=rc_Pb[i][j]/rc_D[i][j];
  }
 }  
 cout << " "<< endl;

 cout << " ******************************** RC are calculated **************************"<< endl;
 
 ofstream out;
 out.open("eRC_qv97_D.txt");
 for(int i=0;  i<3; i++){
  for(int j=0; j<3;j++){
   out << i << " " << j << " " <<  rc_D[i][j] <<endl;
  }
 }
 out.close();
 
 out.open("eRC_qv97_C.txt");
 for(int i=0;  i<3; i++){
  for(int j=0; j<3;j++){
   out << i << " " << j << " " <<  rc_C[i][j] << " " << ratioRC_CD[i][j]<< endl;
  }
 }
 out.close();
 out.open("eRC_qv97_Fe.txt");
 for(int i=0;  i<3; i++){
  for(int j=0; j<3;j++){
   out << i << " " << j << " " <<  rc_Fe[i][j] << " " << ratioRC_FeD[i][j] <<endl;
  }
 }
 out.close();
 out.open("eRC_qv97_Pb.txt");
 for(int i=0;  i<3; i++){
  for(int j=0; j<3;j++){
   out << i << " " << j << " " <<  rc_Pb[i][j] << " " << ratioRC_PbD[i][j]<<endl;
  }
 }
 out.close();

 //////////////////////////////////////////////////////////////////////////////////////


 ////////////////////////// Plot ////////////////////////
 cout << "Plotting" <<endl;  
 gStyle->SetOptTitle(0); 
 
 ////////////////////// 2D /////////////////////////////
 TCanvas *c = new TCanvas("c","c", 3*500, 400);
 TH2F *h = new TH2F("h","h",default_binning, 2.2, 4.2, default_binning, 1., 1.5);
 h->GetXaxis()->SetTitle("#nu, GeV");
 h->GetYaxis()->SetTitle("#delta_{RC}");
 h->GetYaxis()->SetTitleOffset(0.9);
 h->GetXaxis()->SetTitleOffset(0.9);
 h->GetXaxis()->SetTitleSize(0.055);
 h->GetYaxis()->SetTitleSize(0.055);
 h->GetYaxis()->SetLabelSize(0.052);
 h->GetXaxis()->SetLabelSize(0.052);
 c->Divide(3,1, 0.00001, 0.00001);
 for(int i=0; i<3;i++){
  g_RC_D[i] = new TGraphErrors(nvbins, center_v, rc_D[i],0, 0);
  g_RC_C[i] = new TGraphErrors(nvbins, center_v, rc_C[i],0, 0);
  g_RC_Fe[i] = new TGraphErrors(nvbins, center_v, rc_Fe[i],0, 0);
  g_RC_Pb[i] = new TGraphErrors(nvbins, center_v, rc_Pb[i],0,0);
  g_RC_D[i]->SetName("D");
  g_RC_C[i]->SetName("C");
  g_RC_Fe[i]->SetName("Fe");
  g_RC_Pb[i]->SetName("Pb");
 
  c->cd(i+1)->SetGridy();
  h->Draw();
  g_RC_D[i]->SetMarkerColor(1);
  g_RC_D[i]->SetMarkerStyle(20);
  g_RC_D[i]->SetMarkerSize(1.7);
  g_RC_D[i]->Draw("Psame");
  g_RC_C[i]->SetMarkerColor(2);
  g_RC_C[i]->SetMarkerStyle(20);
  g_RC_C[i]->SetMarkerSize(1.7);
  g_RC_C[i]->Draw("Psame");
  g_RC_Fe[i]->SetMarkerColor(4);
  g_RC_Fe[i]->SetMarkerStyle(20);
  g_RC_Fe[i]->SetMarkerSize(1.7);
  g_RC_Fe[i]->Draw("Psame");
  g_RC_Pb[i]->SetMarkerColor(8);
  g_RC_Pb[i]->SetMarkerStyle(20);
  g_RC_Pb[i]->SetMarkerSize(1.7);
  g_RC_Pb[i]->Draw("Psame");
 
  TLatex a;
  a.SetNDC();
  a.SetTextSize(0.07);
  a.DrawLatex(0.33,0.92, Form("%3.2f < Q^{2} < %3.2f", q_values[i], q_values[i+1]));
	 
  auto legend = new TLegend(0.8,0.95,0.92,0.75);
  legend->AddEntry("D","D","p");
  legend->AddEntry("C","C","p");
  legend->AddEntry("Fe","Fe","p");
  legend->AddEntry("Pb","Pb","p");
	        legend->Draw();
 }
 c->Print("~/secure/ExternalsRC/rc_qv97.png");
 c->Clear();

 /*   RATIO   */
 TH2F *hr = new TH2F("hr","hr",default_binning, 2.2, 4.2, default_binning, 0.96, 1.03);
 hr->GetXaxis()->SetTitle("#nu, GeV");
 hr->GetYaxis()->SetTitle("#delta_{A}/#delta_{D}");
 hr->GetYaxis()->SetLabelSize(0.031);
 hr->GetYaxis()->SetTitleOffset(0.97);
 hr->GetXaxis()->SetTitleOffset(0.97);
 hr->GetXaxis()->SetTitleSize(0.05);
 hr->GetYaxis()->SetTitleSize(0.05);
 hr->GetYaxis()->SetLabelSize(0.04);//(0.052);
 hr->GetXaxis()->SetLabelSize(0.045);
 c->Divide(3,1, 0.0001, 0.0001);
 for(int i=0; i<3;i++){
  g_ratioRC_CD[i] = new TGraphErrors(nvbins, center_v, ratioRC_CD[i],0, 0);
  g_ratioRC_FeD[i] = new TGraphErrors(nvbins, center_v, ratioRC_FeD[i],0,0);
  g_ratioRC_PbD[i] = new TGraphErrors(nvbins, center_v, ratioRC_PbD[i],0, 0);
  g_ratioRC_CD[i]->SetName("cd");
  g_ratioRC_FeD[i]->SetName("fed");
  g_ratioRC_PbD[i]->SetName("pbd");
  c->cd(i+1)->SetGridy();
  hr->Draw();
  g_ratioRC_CD[i]->SetMarkerColor(2);
  g_ratioRC_CD[i]->SetMarkerStyle(20);
  g_ratioRC_CD[i]->SetMarkerSize(1.7);
  g_ratioRC_CD[i]->Draw("Psame");
  g_ratioRC_FeD[i]->SetMarkerColor(4);
  g_ratioRC_FeD[i]->SetMarkerStyle(20);
  g_ratioRC_FeD[i]->SetMarkerSize(1.7);
  g_ratioRC_FeD[i]->Draw("Psame");
  g_ratioRC_PbD[i]->SetMarkerColor(8);
  g_ratioRC_PbD[i]->SetMarkerStyle(20);
  g_ratioRC_PbD[i]->SetMarkerSize(1.7);
  g_ratioRC_PbD[i]->Draw("Psame");

  TLatex a;
  a.SetNDC(); 
  a.SetTextColor(1);
  a.SetTextSize(0.07);
  a.DrawLatex(0.33,0.92, Form("%3.2f < Q^{2} < %3.2f", q_values[i], q_values[i+1]));
  
  auto legend = new TLegend(0.8,0.95,0.92,0.75);
  legend->AddEntry("cd","C/D","p");
  legend->AddEntry("fed","Fe/D","p");
  legend->AddEntry("pbd","Pb/D","p");
  legend->Draw();
 }
  c->Print("~/secure/ExternalsRC/Ratiorc_qv97.png");
 c->Clear();

 /*
 
//////////////////////// 1D in v //////////////////////////
 TCanvas *c1 = new TCanvas("c1","c1", 500, 400);
 TH2F *h1 = new TH2F("h1","h1",default_binning, 2.2, 4.2, default_binning, 1., 1.4);
 h1->GetXaxis()->SetTitle("#nu, GeV");
 h1->GetYaxis()->SetTitle("#delta_{RC}");
 h1->GetYaxis()->SetTitleOffset(0.9);
 h1->GetXaxis()->SetTitleOffset(0.9);
 h1->GetXaxis()->SetTitleSize(0.055);
 h1->GetYaxis()->SetTitleSize(0.055);
 h1->GetYaxis()->SetLabelSize(0.052);
 h1->GetXaxis()->SetLabelSize(0.052);
 c1->Divide(nqbins,1, 0.00001, 0.00001);
 c1->SetGridy();

  g_RC_vD = new TGraphErrors(nvbins, center_v, RC_vD,0, 0);
  g_RC_vC = new TGraphErrors(nvbins, center_v, RC_vC,0, 0);
  g_RC_vFe = new TGraphErrors(nvbins, center_v, RC_vFe,0, 0);
  g_RC_vPb = new TGraphErrors(nvbins, center_v, RC_vPb,0,0);
  g_RC_vD->SetName("vD");
  g_RC_vC->SetName("vC");
  g_RC_vFe->SetName("vFe");
  g_RC_vPb->SetName("vPb");
 
  h1->Draw();
  g_RC_vD->SetMarkerColor(1);
  g_RC_vD->SetMarkerStyle(20);
  g_RC_vD->SetMarkerSize(1.7);
  g_RC_vD->Draw("Psame");
  g_RC_vC->SetMarkerColor(2);
  g_RC_vC->SetMarkerStyle(20);
  g_RC_vC->SetMarkerSize(1.7);
  g_RC_vC->Draw("Psame");
  g_RC_vFe->SetMarkerColor(4);
  g_RC_vFe->SetMarkerStyle(20);
  g_RC_vFe->SetMarkerSize(1.7);
  g_RC_vFe->Draw("Psame");
  g_RC_vPb->SetMarkerColor(8);
  g_RC_vPb->SetMarkerStyle(20);
  g_RC_vPb->SetMarkerSize(1.7);
  g_RC_vPb->Draw("Psame");
 
  TLatex a;
  a.SetNDC();
  a.SetTextSize(0.06);
  a.DrawLatex(0.27,0.92, "1.0 < Q^{2} < 4.1 GeV");
	 
  auto legend = new TLegend(0.8,0.95,0.92,0.75);
  legend->AddEntry("vD","D","p");
  legend->AddEntry("vC","C","p");
  legend->AddEntry("vFe","Fe","p");
  legend->AddEntry("vPb","Pb","p");
  legend->Draw();
  c1->Print("~/secure/ExternalsRC/rc_v.png");
  c1->Clear();


 TH2F *hr1 = new TH2F("hr1","hr1",default_binning, 2.2, 4.2, default_binning, 0.96, 1.03);
 hr1->GetXaxis()->SetTitle("#nu, GeV");
 hr1->GetYaxis()->SetTitle("#delta_{A}/#delta_{D}");
 hr1->GetYaxis()->SetLabelSize(0.031);
 hr1->GetYaxis()->SetTitleOffset(0.97);
 hr1->GetXaxis()->SetTitleOffset(0.97);
 hr1->GetXaxis()->SetTitleSize(0.05);
 hr1->GetYaxis()->SetTitleSize(0.05);
 hr1->GetYaxis()->SetLabelSize(0.04);//(0.052);
 hr1->GetXaxis()->SetLabelSize(0.045);
  
 g_ratioRC_vCD = new TGraphErrors(nvbins, center_v, ratioRC_vCD,0, 0);
 g_ratioRC_vFeD = new TGraphErrors(nvbins, center_v, ratioRC_vFeD,0,0); 
 g_ratioRC_vPbD = new TGraphErrors(nvbins, center_v, ratioRC_vPbD,0, 0);
 g_ratioRC_vCD->SetName("vcd");
 g_ratioRC_vFeD->SetName("vfed");
 g_ratioRC_vPbD->SetName("vpbd");
 
  c1->SetGridy();
  hr1->Draw();
  g_ratioRC_vCD->SetMarkerColor(2);
  g_ratioRC_vCD->SetMarkerStyle(20);
  g_ratioRC_vCD->SetMarkerSize(1.7);
  g_ratioRC_vCD->Draw("Psame");
  g_ratioRC_vFeD->SetMarkerColor(4);
  g_ratioRC_vFeD->SetMarkerStyle(20);
  g_ratioRC_vFeD->SetMarkerSize(1.7);
  g_ratioRC_vFeD->Draw("Psame");
  g_ratioRC_vPbD->SetMarkerColor(8);
  g_ratioRC_vPbD->SetMarkerStyle(20);
  g_ratioRC_vPbD->SetMarkerSize(1.7);
  g_ratioRC_vPbD->Draw("Psame");

  a.SetTextSize(0.06);
  a.DrawLatex(0.33,0.92, "1.0 < Q^{2} < 4.1 GeV^{2}");
  
  auto legend2 = new TLegend(0.8,0.95,0.92,0.75);
  legend2->AddEntry("vcd","C/D","p");
  legend2->AddEntry("vfed","Fe/D","p");
  legend2->AddEntry("vpbd","Pb/D","p");
  legend2->Draw();
 c1->Print("~/secure/ExternalsRC/Ratiorc_v.png");
 c1->Clear();
*/

}
