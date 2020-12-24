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

Float_t ratioRC_CD[nqbins][nvbins], ratioRC_FeD[nqbins][nvbins], ratioRC_PbD[nqbins][nvbins];
Float_t ratioRC1_CD[nqbins][nvbins], ratioRC1_FeD[nqbins][nvbins], ratioRC1_PbD[nqbins][nvbins];
Float_t ratioRC100_CD[nqbins][nvbins], ratioRC100_FeD[nqbins][nvbins], ratioRC100_PbD[nqbins][nvbins];

Float_t systRC_CD[nqbins][nvbins], systRC_FeD[nqbins][nvbins], systRC_PbD[nqbins][nvbins];
Float_t syst_percRC_CD[nqbins][nvbins], syst_percRC_FeD[nqbins][nvbins], syst_percRC_PbD[nqbins][nvbins];

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
float rc;

cout << " **************** C ********************** " << endl;
 in.open("../RC/eRC_qv_C.txt");
 for(int nv=0;  nv<nvbins; nv++){  
  for(int nq=0; nq<nqbins;nq++){
   in>>nq>>nv>>rc>>ratioRC_CD[nq][nv]; 
   }
  }
 in.close();

 cout << " **************** Fe ********************** " << endl; 
 in.open("../RC/eRC_qv_Fe.txt");
 for(int nv=0;  nv<nvbins; nv++){  
  for(int nq=0; nq<nqbins;nq++){
   in>>nq>>nv>>rc>>ratioRC_FeD[nq][nv]; 
   }
  }
 in.close(); 

 cout << " **************** Pb ********************** " << endl; 
 in.open("../RC/eRC_qv_Pb.txt");
 for(int nv=0;  nv<nvbins; nv++){  
  for(int nq=0; nq<nqbins;nq++){
   in>>nq>>nv>>rc>>ratioRC_PbD[nq][nv]; 
   }
  }
 in.close();
///////////////////////// read Delta 1 MeV //////////////////

cout << " **************** C 1MeV ********************** " << endl;
 in.open("../externals_Delta1/RC/eRC_qv_C.txt");
 for(int nv=0;  nv<nvbins; nv++){  
  for(int nq=0; nq<nqbins;nq++){
   in>>nq>>nv>>rc>>ratioRC1_CD[nq][nv]; 
   }
  }
 in.close();

 cout << " **************** Fe ********************** " << endl; 
 in.open("../externals_Delta1/RC/eRC_qv_Fe.txt");
 for(int nv=0;  nv<nvbins; nv++){  
  for(int nq=0; nq<nqbins;nq++){
   in>>nq>>nv>>rc>>ratioRC1_FeD[nq][nv]; 
   }
  }
 in.close(); 

 cout << " **************** Pb ********************** " << endl; 
 in.open("../externals_Delta1/RC/eRC_qv_Pb.txt");
 for(int nv=0;  nv<nvbins; nv++){  
  for(int nq=0; nq<nqbins;nq++){
   in>>nq>>nv>>rc>>ratioRC1_PbD[nq][nv]; 
   }
  }
 in.close();

///////////////////////// read Delta 100 MeV //////////////////

cout << " **************** C 100MeV ********************** " << endl;
 in.open("../externals_Delta100/RC/eRC_qv_C.txt");
 for(int nv=0;  nv<nvbins; nv++){  
  for(int nq=0; nq<nqbins;nq++){
   in>>nq>>nv>>rc>>ratioRC100_CD[nq][nv]; 
   }
  }
 in.close();

 cout << " **************** Fe ********************** " << endl; 
 in.open("../externals_Delta100/RC/eRC_qv_Fe.txt");
 for(int nv=0;  nv<nvbins; nv++){  
  for(int nq=0; nq<nqbins;nq++){
   in>>nq>>nv>>rc>>ratioRC100_FeD[nq][nv]; 
   }
  }
 in.close(); 

 cout << " **************** Pb ********************** " << endl; 
 in.open("../externals_Delta100/RC/eRC_qv_Pb.txt");
 for(int nv=0;  nv<nvbins; nv++){  
  for(int nq=0; nq<nqbins;nq++){
   in>>nq>>nv>>rc>>ratioRC100_PbD[nq][nv]; 
   }
  }
 in.close();

 /////////////////////// Calculate systematic error /////////////
 
 for(int nv=0;  nv<nvbins; nv++){
  for(int nq=0; nq<nqbins;nq++){	 
   systRC_CD[nq][nv] = TMath::Sqrt((TMath::Power(ratioRC_CD[nq][nv]-ratioRC1_CD[nq][nv],2)
		                  +TMath::Power(ratioRC_CD[nq][nv]-ratioRC100_CD[nq][nv],2))/2);
   syst_percRC_CD[nq][nv] = systRC_CD[nq][nv] / ratioRC_CD[nq][nv] *100;
   cout << "Carbon syst = " << systRC_CD[nq][nv]<< "  and is "  << syst_percRC_CD[nq][nv] << " %" <<endl;
 
   systRC_FeD[nq][nv] = TMath::Sqrt((TMath::Power(ratioRC_FeD[nq][nv]-ratioRC1_FeD[nq][nv],2)
		                  +TMath::Power(ratioRC_FeD[nq][nv]-ratioRC100_FeD[nq][nv],2))/2);
   syst_percRC_FeD[nq][nv] = systRC_FeD[nq][nv] / ratioRC_FeD[nq][nv] *100;
   cout << "Iron syst = " << systRC_FeD[nq][nv]<< "  and is "  << syst_percRC_FeD[nq][nv] << " %" <<endl;
 
   systRC_PbD[nq][nv] = TMath::Sqrt((TMath::Power(ratioRC_PbD[nq][nv]-ratioRC1_PbD[nq][nv],2)
		                  +TMath::Power(ratioRC_PbD[nq][nv]-ratioRC100_PbD[nq][nv],2))/2);
   syst_percRC_PbD[nq][nv] = systRC_PbD[nq][nv] / ratioRC_PbD[nq][nv] *100;
   cout << "Lead syst = " << systRC_PbD[nq][nv]<< "  and is "  << syst_percRC_PbD[nq][nv] << " %" <<endl;
 
  }
 }

 ////////////////////// RC observed/absolute ////////////////// 
/* 
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
/*

 ////////////////////////// Plot ////////////////////////
 cout << "Plotting" <<endl;  
 gStyle->SetOptTitle(0); 
 ////////////////////// 2D /////////////////////////////
 TCanvas *c = new TCanvas("c","c", nqbins*500, 400);
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
  
  auto legend = new TLegend(0.8,0.95,0.92,0.75);
  legend->AddEntry("cd","C/D","p");
  legend->AddEntry("fed","Fe/D","p");
  legend->AddEntry("pbd","Pb/D","p");
  legend->Draw();
 }
  c->Print("~/secure/ExternalsRC/systemRC_ratio.png");
*/
}
