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

const float q_values[4]= {1.0, 1.33, 1.76, 4.1 }; 
const float v_values[4]= {2.2, 3.2, 3.73, 4.25 };
const float center_q[3] = {1.165 , 1.55, 2.93 };
const float center_v[3] = {2.7, 3.46, 3.99 };


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

 in.open("clasd2.out");
 for(int nv=0;  nv<nvbins; nv++){  
  for(int nq=0; nq<nqbins;nq++){
   in>>Ebeam>>Eprime[nq][nv]>>Theta[nq][nv]>>x>>Q2>>sigma_Born_D[nq][nv]>>sigma_BornIn_D[nq][nv]>>sigma_BornQE_D[nq][nv]
     >>sigma_Rad_D[nq][nv]>> sigma_RadEl_D[nq][nv]>>sigma_RadQE_D[nq][nv]>>sigma_RadDIS_D[nq][nv]>>C_corr_D[nq][nv];
   
   cout << " **************** D ********************** " << endl; 
   cout << Ebeam << " " << Eprime[nq][nv]<< " " << Theta[nq][nv]<< " "<< x<< " " << Q2 << " " << sigma_Born_D[nq][nv]<< " " 
	   << sigma_BornIn_D[nq][nv]<< " " << sigma_BornQE_D[nq][nv]<< " " << sigma_Rad_D[nq][nv]<< " " 
	   << sigma_RadEl_D[nq][nv]<< " "<< sigma_RadQE_D[nq][nv]
	   << " " << sigma_RadDIS_D[nq][nv]<< " " <<C_corr_D[nq][nv]<< endl;
   }
  }
 in.close();

 in.open("clasC12.out");
 for(int nv=0;  nv<nvbins; nv++){  
  for(int nq=0; nq<nqbins;nq++){
   in>>Ebeam>>Eprime[nq][nv]>>Theta[nq][nv]>>x>>Q2>>sigma_Born_C[nq][nv]>>sigma_BornIn_C[nq][nv]>>sigma_BornQE_C[nq][nv]
     >>sigma_Rad_C[nq][nv]>> sigma_RadEl_C[nq][nv]>>sigma_RadQE_C[nq][nv]>>sigma_RadDIS_C[nq][nv]>>C_corr_C[nq][nv];
   
   cout << " **************** C ********************** " << endl; 
   cout << Ebeam << " " << Eprime[nq][nv]<< " " << Theta[nq][nv]<< " "<< x<< " " << Q2 << " " << sigma_Born_C[nq][nv]<< " " 
	   << sigma_BornIn_C[nq][nv]<< " " << sigma_BornQE_C[nq][nv]<< " " << sigma_Rad_C[nq][nv]<< " " 
	   << sigma_RadEl_C[nq][nv]<< " "<< sigma_RadQE_C[nq][nv]
	   << " " << sigma_RadDIS_C[nq][nv]<< " " <<C_corr_C[nq][nv]<< endl;
   }
  }
 in.close();

 cout << " **************** Fe ********************** " << endl; 
 in.open("clasFe56.out");
 for(int nv=0;  nv<nvbins; nv++){  
  for(int nq=0; nq<nqbins;nq++){
   in>>Ebeam>>Eprime[nq][nv]>>Theta[nq][nv]>>x>>Q2>>sigma_Born_Fe[nq][nv]>>sigma_BornIn_Fe[nq][nv]>>sigma_BornQE_Fe[nq][nv]
     >>sigma_Rad_Fe[nq][nv]>> sigma_RadEl_Fe[nq][nv]>>sigma_RadQE_Fe[nq][nv]>>sigma_RadDIS_Fe[nq][nv]>>C_corr_Fe[nq][nv];
   
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
 for(int nv=0;  nv<nvbins; nv++){  
  for(int nq=0; nq<nqbins;nq++){
   in>>Ebeam>>Eprime[nq][nv]>>Theta[nq][nv]>>x>>Q2>>sigma_Born_Pb[nq][nv]>>sigma_BornIn_Pb[nq][nv]>>sigma_BornQE_Pb[nq][nv]
     >>sigma_Rad_Pb[nq][nv]>> sigma_RadEl_Pb[nq][nv]>>sigma_RadQE_Pb[nq][nv]>>sigma_RadDIS_Pb[nq][nv]>>C_corr_Pb[nq][nv];
   
 /*cout << Ebeam << " " << Eprime[nq][nv]<< " " << Theta[nq][nv]<< " "<< x<< " " << Q2 << " " << sigma_Born_Pb[nq][nv]<< " " 
	   << sigma_BornIn_Pb[nq][nv]<< " " << sigma_BornQE_Pb[nq][nv]<< " " << sigma_Rad_Pb[nq][nv]<< " " 
	   << sigma_RadEl_Pb[nq][nv]<< " "<< sigma_RadQE_Pb[nq][nv]
	   << " " << sigma_RadDIS_Pb[nq][nv]<< " " <<C_corr_Pb[nq][nv]<< endl;
   */
   }
  }
 in.close();

 
 ////////////////////// RC observed/absolute ////////////////// 
 
 cout << " "<< endl;
 cout << " ******************************** RC *********************************"<< endl;
 float deltaq;
 for(int nv=0;  nv<nvbins; nv++){  
  for(int nq=0; nq<nqbins;nq++){
   //deltaq = (q_values[nq+1]-q_values[nq]); 
   //nD_v_rad[nv] += nD_rad[nq][nv]   * deltaq;
   RC_D[nq][nv]= sigma_Rad_D[nq][nv]/sigma_Born_D[nq][nv];
   RC_C[nq][nv]= sigma_Rad_C[nq][nv]/sigma_Born_C[nq][nv];
   RC_Fe[nq][nv]= sigma_Rad_Fe[nq][nv]/sigma_Born_Fe[nq][nv];
   RC_Pb[nq][nv]= sigma_Rad_Pb[nq][nv]/sigma_Born_Pb[nq][nv];

   ratioRC_CD[nq][nv] = RC_C[nq][nv]/RC_D[nq][nv];
   ratioRC_FeD[nq][nv] = RC_Fe[nq][nv]/RC_D[nq][nv];
   ratioRC_PbD[nq][nv] = RC_Pb[nq][nv]/RC_D[nq][nv];

   cout << nq <<" " <<nv<< " RC(D)" << RC_D[nq][nv] << " RC(Pb)" << RC_Pb[nq][nv]<< " "<<
	 "Ratio(Pb/D)" << ratioRC_PbD[nq][nv]<<endl; 
  
  }
 }

 
 ofstream out;
 out.open("eRC_qv_D.txt");
 for(int nq=0;  nq<nqbins; nq++){
  for(int nv=0; nv<nvbins;nv++){
   out << nq << " " << nv << " " <<  RC_D[nq][nv] <<endl;
  }
 }
 out.close();
 
 out.open("eRC_qv_C.txt");
 for(int nq=0;  nq<nqbins; nq++){
  for(int nv=0; nv<nvbins;nv++){
   out << nq << " " << nv << " " <<  RC_C[nq][nv] << " " << ratioRC_CD[nq][nv]<< endl;
  }
 }
 out.close();
 out.open("eRC_qv_Fe.txt");
 for(int nq=0;  nq<nqbins; nq++){
  for(int nv=0; nv<nvbins;nv++){
   out << nq << " " << nv << " " <<  RC_Fe[nq][nv] << " " << ratioRC_FeD[nq][nv] <<endl;
  }
 }
 out.close();
 out.open("eRC_qv_Pb.txt");
 for(int nq=0;  nq<nqbins; nq++){
  for(int nv=0; nv<nvbins;nv++){
   out << nq << " " << nv << " " <<  RC_Pb[nq][nv] << " " << ratioRC_PbD[nq][nv]<<endl;
  }
 }
 out.close();

 ////////////////////////////////////////////////////////
 /*
 out.open("eRC_v_D.txt");
 for(int nv=0; nv<nvbins;nv++){
  out << nv << " " << nD_v_rad[nv]/n_v_norad[nv] <<endl;
  }
 out.close();

 out.open("eRC_v_C.txt");
 for(int nv=0; nv<nvbins;nv++){
  out << nv << " " << nC_v_rad[nv]/n_v_norad[nv] <<endl;
  }
 out.close(); 

 out.open("eRC_v_Fe.txt");
 for(int nv=0; nv<nvbins;nv++){
  out << nv << " " << nFe_v_rad[nv]/n_v_norad[nv] <<endl;
  }
 out.close();

 out.open("eRC_v_Pb.txt");
 for(int nv=0; nv<nvbins;nv++){
  out << nv << " " << nPb_v_rad[nv]/n_v_norad[nv] <<endl;
  }
 out.close();

  */

 ////////////////////////// Plot ////////////////////////
 cout << "Plotting" <<endl;  
 gStyle->SetOptTitle(0); 
 ////////////////////// 2D /////////////////////////////
 TCanvas *c = new TCanvas("c","c", nqbins*500, 400);
 TH2F *h = new TH2F("h","h",default_binning, 2.2, 4.2, default_binning, 1., 1.4);
 h->GetXaxis()->SetTitle("#nu, GeV");
 h->GetYaxis()->SetTitle("#delta_{RC}");
 h->GetYaxis()->SetTitleOffset(1.03);
 h->GetXaxis()->SetTitleOffset(0.92);
 h->GetXaxis()->SetTitleSize(0.055);
 h->GetYaxis()->SetTitleSize(0.045);
 h->GetYaxis()->SetLabelSize(0.052);
 h->GetXaxis()->SetLabelSize(0.052);
 c->Divide(nqbins,1, 0.00001, 0.00001);
 for(int nq=0; nq<nqbins;nq++){
   c->cd(nq+1)->SetGridy();
  g_RC_D[nq] = new TGraphErrors(nvbins, center_v, RC_D[nq],0, 0);
  g_RC_C[nq] = new TGraphErrors(nvbins, center_v, RC_C[nq],0, 0);
  g_RC_Fe[nq] = new TGraphErrors(nvbins, center_v, RC_Fe[nq],0, 0);
  g_RC_Pb[nq] = new TGraphErrors(nvbins, center_v, RC_Pb[nq],0,0);
  g_RC_D[nq]->SetName("D");
  g_RC_C[nq]->SetName("C");
  g_RC_Fe[nq]->SetName("Fe");
  g_RC_Pb[nq]->SetName("Pb");
 
  h->Draw();
  g_RC_D[nq]->SetMarkerColor(1);
  g_RC_D[nq]->SetMarkerStyle(20);
  g_RC_D[nq]->SetMarkerSize(1.7);
  g_RC_D[nq]->Draw("Psame");
  g_RC_C[nq]->SetMarkerColor(2);
  g_RC_C[nq]->SetMarkerStyle(20);
  g_RC_C[nq]->SetMarkerSize(1.7);
  g_RC_C[nq]->Draw("Psame");
  g_RC_Fe[nq]->SetMarkerColor(4);
  g_RC_Fe[nq]->SetMarkerStyle(20);
  g_RC_Fe[nq]->SetMarkerSize(1.7);
  g_RC_Fe[nq]->Draw("Psame");
  g_RC_Pb[nq]->SetMarkerColor(8);
  g_RC_Pb[nq]->SetMarkerStyle(20);
  g_RC_Pb[nq]->SetMarkerSize(1.7);
  g_RC_Pb[nq]->Draw("Psame");
 
  TLatex a;
  a.SetNDC();
  a.SetTextSize(0.07);
  a.DrawLatex(0.33,0.92, Form("%3.2f < Q^{2} < %3.2f", q_values[nq], q_values[nq+1]));

 
  /********
  a.SetTextSize(0.05);
  a.SetTextColor(1);
  a.DrawLatex(0.91,0.87, "D");
  a.SetTextColor(2);
  a.DrawLatex(0.91,0.82, "C");
  a.SetTextColor(4);
  a.DrawLatex(0.91,0.77, "Fe");
  a.SetTextColor(8);
  a.DrawLatex(0.91,0.72, "Pb");
*****/

	 
  auto legend = new TLegend(0.8,0.95,0.92,0.75);
  legend->AddEntry("D","D","p");
  legend->AddEntry("C","C","p");
  legend->AddEntry("Fe","Fe","p");
  legend->AddEntry("Pb","Pb","p");
	        legend->Draw();
 }
 c->Print("~/secure/ExternalsRC/rc_Delta100MeV.png");
 c->Clear();

 TH2F *hr = new TH2F("hr","hr",default_binning, 2.2, 4.2, default_binning, 0.96, 1.03);
 hr->GetXaxis()->SetTitle("#nu, GeV");
 hr->GetYaxis()->SetTitle("#delta_{C}/#delta_{D}");
 hr->GetYaxis()->SetLabelSize(0.031);
 hr->GetYaxis()->SetTitleOffset(1.05);
 hr->GetXaxis()->SetTitleOffset(0.97);
 hr->GetXaxis()->SetTitleSize(0.055);
 hr->GetYaxis()->SetTitleSize(0.045);
 hr->GetYaxis()->SetLabelSize(0.045);//(0.052);
 hr->GetXaxis()->SetLabelSize(0.045);
 c->Divide(nqbins,1, 0.0001, 0.0001);
 for(int nq=0; nq<nqbins;nq++){
  g_ratioRC_CD[nq] = new TGraphErrors(nvbins, center_v, ratioRC_CD[nq],0, 0);
  g_ratioRC_FeD[nq] = new TGraphErrors(nvbins, center_v, ratioRC_FeD[nq],0,0);
  g_ratioRC_PbD[nq] = new TGraphErrors(nvbins, center_v, ratioRC_PbD[nq],0, 0);
  g_ratioRC_CD[nq]->SetName("cd");
  g_ratioRC_FeD[nq]->SetName("fed");
  g_ratioRC_PbD[nq]->SetName("pbd");
 
  c->cd(nq+1)->SetGridy();
  hr->Draw();
  g_ratioRC_CD[nq]->SetMarkerColor(2);
  g_ratioRC_CD[nq]->SetMarkerStyle(20);
  g_ratioRC_CD[nq]->SetMarkerSize(1.7);
  g_ratioRC_CD[nq]->Draw("Psame");
  g_ratioRC_FeD[nq]->SetMarkerColor(4);
  g_ratioRC_FeD[nq]->SetMarkerStyle(20);
  g_ratioRC_FeD[nq]->SetMarkerSize(1.7);
  g_ratioRC_FeD[nq]->Draw("Psame");
  g_ratioRC_PbD[nq]->SetMarkerColor(8);
  g_ratioRC_PbD[nq]->SetMarkerStyle(20);
  g_ratioRC_PbD[nq]->SetMarkerSize(1.7);
  g_ratioRC_PbD[nq]->Draw("Psame");

  TLatex a;
  a.SetNDC(); 
  a.SetTextColor(1);
  a.SetTextSize(0.07);
  a.DrawLatex(0.33,0.92, Form("%3.2f < Q^{2} < %3.2f", q_values[nq], q_values[nq+1]));

     	/*******
  a.SetTextSize(0.05);
  a.SetTextColor(2);
  a.DrawLatex(0.91,0.82, "C/D");
  a.SetTextColor(4);
  a.DrawLatex(0.91,0.77, "Fe/D");
  a.SetTextColor(8);
  a.DrawLatex(0.91,0.72, "Pb/D");
  ******/
  
  auto legend = new TLegend(0.8,0.95,0.92,0.75);
  legend->AddEntry("cd","C/D","p");
  legend->AddEntry("fed","Fe/D","p");
  legend->AddEntry("pbd","Pb/D","p");
  legend->Draw();
 }
  c->Print("~/secure/ExternalsRC/Ratiorc_Delta100MeV.png");

}
