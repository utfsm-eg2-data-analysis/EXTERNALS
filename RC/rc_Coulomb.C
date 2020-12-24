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
Float_t Eprime_v[nvbins], Theta_v[nvbins];

Float_t sigma_Born_D[nqbins][nvbins], sigma_Rad_D[nqbins][nvbins];
Float_t sigma_BornIn_D[nqbins][nvbins], sigma_BornQE_D[nqbins][nvbins];
Float_t sigma_RadEl_D[nqbins][nvbins], sigma_RadQE_D[nqbins][nvbins], sigma_RadDIS_D[nqbins][nvbins];
Float_t C_corr_D[nqbins][nvbins];
Float_t C_corr_vD[nvbins];

Float_t sigma_Born_C[nqbins][nvbins], sigma_Rad_C[nqbins][nvbins];
Float_t sigma_BornIn_C[nqbins][nvbins], sigma_BornQE_C[nqbins][nvbins];
Float_t sigma_RadEl_C[nqbins][nvbins], sigma_RadQE_C[nqbins][nvbins], sigma_RadDIS_C[nqbins][nvbins];
Float_t C_corr_C[nqbins][nvbins];
Float_t C_corr_vC[nvbins];

Float_t sigma_Born_Fe[nqbins][nvbins], sigma_Rad_Fe[nqbins][nvbins];
Float_t sigma_BornIn_Fe[nqbins][nvbins], sigma_BornQE_Fe[nqbins][nvbins];
Float_t sigma_RadEl_Fe[nqbins][nvbins], sigma_RadQE_Fe[nqbins][nvbins], sigma_RadDIS_Fe[nqbins][nvbins];
Float_t C_corr_Fe[nqbins][nvbins];
Float_t C_corr_vFe[nvbins];

Float_t sigma_Born_Pb[nqbins][nvbins], sigma_Rad_Pb[nqbins][nvbins];
Float_t sigma_BornIn_Pb[nqbins][nvbins], sigma_BornQE_Pb[nqbins][nvbins];
Float_t sigma_RadEl_Pb[nqbins][nvbins], sigma_RadQE_Pb[nqbins][nvbins], sigma_RadDIS_Pb[nqbins][nvbins];
Float_t C_corr_Pb[nqbins][nvbins];
Float_t C_corr_vPb[nvbins];

Float_t RC_D[nqbins][nvbins],RC_C[nqbins][nvbins],RC_Fe[nqbins][nvbins],RC_Pb[nqbins][nvbins];
Float_t ratioRC_CD[nqbins][nvbins], ratioRC_FeD[nqbins][nvbins], ratioRC_PbD[nqbins][nvbins];

// here sigma_Born_v  are summed from 2D grid on (Q2,nu)
Float_t sigma_Born_vD[nvbins], sigma_Born_vC[nvbins],sigma_Born_vFe[nvbins],sigma_Born_vPb[nvbins];
Float_t sigma_CC_vD[nvbins], sigma_CC_vC[nvbins],sigma_CC_vFe[nvbins],sigma_CC_vPb[nvbins];
// here average 
Float_t sigmatot_Born_C, sigmatot_Born_Fe, sigmatot_Born_Pb;
Float_t sigmatot_CC_C, sigmatot_CC_Fe, sigmatot_CC_Pb;
Float_t C_Corrtot_C, C_Corrtot_Fe, C_Corrtot_Pb;

Float_t sigma_Rad_vsumD[nvbins], sigma_Rad_vsumC[nvbins],sigma_Rad_vsumFe[nvbins],sigma_Rad_vsumPb[nvbins];
Float_t RC_vsumD[nvbins], RC_vsumC[nvbins],RC_vsumFe[nvbins],RC_vsumPb[nvbins];
Float_t  ratioRC_vsumCD[nvbins],ratioRC_vsumFeD[nvbins],ratioRC_vsumPbD[nvbins];

TGraphErrors *g_RC_D[nqbins],*g_RC_C[nqbins],*g_RC_Fe[nqbins],*g_RC_Pb[nqbins];
TGraphErrors *g_Born_Pb[nqbins],*g_Born_C[nqbins],*g_Born_Fe[nqbins],*g_Born_D[nqbins];
TGraphErrors *g_Born_vD,*g_Born_vC,*g_Born_vFe,*g_Born_vPb;
TGraphErrors *g_ratioRC_CD[nqbins], *g_ratioRC_FeD[nqbins], *g_ratioRC_PbD[nqbins];

TGraphErrors *g_Coulomb_vD,*g_Coulomb_vC,*g_Coulomb_vFe,*g_Coulomb_vPb;

void rc_Coulomb(void){
 gROOT->SetStyle("Plain");
 gStyle->SetPalette(1);
 gStyle->SetOptFit(1);
 gStyle->SetOptStat(0);
 gStyle->SetOptTitle(1);
 gStyle->SetOptDate(0);

 ////////////////////// READ RC ////////////////////////
Float_t Ebeam, x, Q2; 
Float_t  x_v, Q2_v; 
 ifstream in;

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
    C_corr_D[nq][nv] = 1.0 *C_corr_D[nq][nv];

  // cout << " **************** D ********************** " << endl; 
  // cout << Ebeam << " " << Eprime[nq][nv]<< " " << Theta[nq][nv]<< " "<< x<< " " << Q2 << " " << sigma_Born_D[nq][nv]<< " " 
//	   << sigma_BornIn_D[nq][nv]<< " " << sigma_BornQE_D[nq][nv]<< " " << sigma_Rad_D[nq][nv]<< " " 
//	   << sigma_RadEl_D[nq][nv]<< " "<< sigma_RadQE_D[nq][nv]
//	   << " " << sigma_RadDIS_D[nq][nv]<< " " <<C_corr_D[nq][nv]<< endl;
   }
  }
 in.close();

 in.open("clasC12.out");
 for(int nq=0; nq<nqbins;nq++){
  for(int nv=0;  nv<nvbins; nv++){  
   in >> Ebeam >> Eprime[nq][nv] >> Theta[nq][nv] >> x >> Q2 >> sigma_Born_C[nq][nv] >> sigma_BornIn_C[nq][nv] >> sigma_BornQE_C[nq][nv]
      >> sigma_Rad_C[nq][nv] >>  sigma_RadEl_C[nq][nv] >> sigma_RadQE_C[nq][nv] >> sigma_RadDIS_C[nq][nv] >> C_corr_C[nq][nv];
   
  // cout << " **************** C ********************** " << endl; 
   //cout << Ebeam << " " << Eprime[nq][nv]<< " " << Theta[nq][nv]<< " "<< x<< " " << Q2 << " " << sigma_Born_C[nq][nv]<< " " 
//	   << sigma_BornIn_C[nq][nv]<< " " << sigma_BornQE_C[nq][nv]<< " " << sigma_Rad_C[nq][nv]<< " " 
//	   << sigma_RadEl_C[nq][nv]<< " "<< sigma_RadQE_C[nq][nv]
//	   << " " << sigma_RadDIS_C[nq][nv]<< " " <<C_corr_C[nq][nv]<< endl;
   }
  }
 in.close();

// cout << " **************** Fe ********************** " << endl; 
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


// cout << " **************** Pb ********************** " << endl; 
 in.open("clasPb208.out");
 for(int nq=0; nq<nqbins;nq++){
  for(int nv=0;  nv<nvbins; nv++){  
   in >> Ebeam >> Eprime[nq][nv] >> Theta[nq][nv] >> x >> Q2 >> sigma_Born_Pb[nq][nv] >> sigma_BornIn_Pb[nq][nv] >> sigma_BornQE_Pb[nq][nv]
      >> sigma_Rad_Pb[nq][nv] >>  sigma_RadEl_Pb[nq][nv] >> sigma_RadQE_Pb[nq][nv] >> sigma_RadDIS_Pb[nq][nv] >> C_corr_Pb[nq][nv];
   
   }
  }
 in.close();


////////////////////// RC observed/absolute ////////////////// 
 
 cout << " "<< endl;
 cout << " ******************************** RC *********************************"<< endl;
 float deltaq;
 for(int nq=0; nq<nqbins;nq++){
  for(int nv=0;  nv<nvbins; nv++){  
   /////////// 2D ///////////
   RC_D[nq][nv]= sigma_Rad_D[nq][nv]/sigma_Born_D[nq][nv];
   RC_C[nq][nv]= sigma_Rad_C[nq][nv]/sigma_Born_C[nq][nv];
   RC_Fe[nq][nv]= sigma_Rad_Fe[nq][nv]/sigma_Born_Fe[nq][nv];
   RC_Pb[nq][nv]= sigma_Rad_Pb[nq][nv]/sigma_Born_Pb[nq][nv];
   
   ratioRC_CD[nq][nv] = RC_C[nq][nv]/RC_D[nq][nv];
   ratioRC_FeD[nq][nv] = RC_Fe[nq][nv]/RC_D[nq][nv];
   ratioRC_PbD[nq][nv] = RC_Pb[nq][nv]/RC_D[nq][nv];
   
   ///// prep for 1D coulomb ////
   sigma_Born_vD[nv] += sigma_Born_D[nq][nv];
   sigma_Born_vC[nv] += sigma_Born_C[nq][nv];
   sigma_Born_vFe[nv] += sigma_Born_Fe[nq][nv];
   sigma_Born_vPb[nv] += sigma_Born_Pb[nq][nv];
   sigma_CC_vD[nv] += (sigma_Born_D[nq][nv]/C_corr_D[nq][nv]); 
   sigma_CC_vC[nv] += (sigma_Born_C[nq][nv]/C_corr_C[nq][nv]); 
   sigma_CC_vFe[nv] += (sigma_Born_Fe[nq][nv]/C_corr_Fe[nq][nv]); 
   sigma_CC_vPb[nv] += (sigma_Born_Pb[nq][nv]/C_corr_Pb[nq][nv]); 
   //prep for average overall
   sigmatot_Born_C += sigma_Born_C[nq][nv];
   sigmatot_Born_Fe += sigma_Born_Fe[nq][nv];
   sigmatot_Born_Pb += sigma_Born_Pb[nq][nv];
   sigmatot_CC_C += (sigma_Born_C[nq][nv]/C_corr_C[nq][nv]);
   sigmatot_CC_Fe += (sigma_Born_Fe[nq][nv]/C_corr_Fe[nq][nv]);
   sigmatot_CC_Pb += (sigma_Born_Pb[nq][nv]/C_corr_Pb[nq][nv]);
 //cout<< " check1: nv=" << nv << " nq=" <<nq << " sigma_Born_C[nq][nv]=" <<sigma_Born_C[nq][nv] 
//	<< " C_corr_C[nq][nv]=" << C_corr_C[nq][nv]<<endl; 
  }
 }
 
 cout<<endl<<endl;
 //////// Coulomb average//////
 C_Corrtot_C=sigmatot_Born_C / sigmatot_CC_C;
 C_Corrtot_Fe=sigmatot_Born_Fe / sigmatot_CC_Fe;
 C_Corrtot_Pb=sigmatot_Born_Pb / sigmatot_CC_Pb;

 cout<< " <ave_CC_C>=" << (C_Corrtot_C-1)*100 <<"%  "  
    << " <ave_CC_Fe>=" << (C_Corrtot_Fe-1)*100<<"% "	
    << " <ave_CC_Pb>=" << (C_Corrtot_Pb-1)*100<<"% "<<endl;   	
cout<<endl; 
///////// COulomb in 1D ///////
   for(int nv=0;  nv<nvbins; nv++){
   C_corr_vD[nv] = sigma_Born_vD[nv]/sigma_CC_vD[nv];
   C_corr_vC[nv] = sigma_Born_vC[nv]/sigma_CC_vC[nv];
   C_corr_vFe[nv] = sigma_Born_vFe[nv]/sigma_CC_vFe[nv];
   C_corr_vPb[nv] = sigma_Born_vPb[nv]/sigma_CC_vPb[nv];
  
   cout << " Coulomb 1D in v-bin = " << nv << " D=" <<C_corr_vD[nv]
	 <<" C=" <<C_corr_vC[nv] << " Fe="<<C_corr_vFe[nv] << " Pb="<<C_corr_vPb[nv]<<endl;  
   }
 
 float sigma_Born_avePb[nvbins];
 for(int nv=0; nv<nvbins;nv++){
     sigma_Born_avePb[nv] = sigma_Born_vPb[nv]/3.;
 }

//////////////////////////// Write ////////////////////////// 
 ofstream out;
 out.open("eCoulomb_qv_all.txt");
 for(int nq=0;  nq<nqbins; nq++){
  for(int nv=0; nv<nvbins;nv++){
   out << nq << " " << nv << " " <<  C_corr_D[nq][nv] <<" " << C_corr_C[nq][nv] << " "  
	                         << C_corr_Fe[nq][nv] <<" " <<  C_corr_Pb[nq][nv] << endl;
  }
 }
 out.close();
 out.open("eCoulomb_v_all.txt");
 for(int nv=0; nv<nvbins;nv++){
  out <<nv << " " <<  C_corr_vD[nv] <<" " << C_corr_vC[nv] << " "
	  << C_corr_vFe[nv] <<" " <<  C_corr_vPb[nv] << endl;
 }

 out.close();
 ////////////////////////// Plot ////////////////////////
 cout << "Plotting" <<endl;  
 gStyle->SetOptTitle(0); 
 
 ////////////////////// 2D /////////////////////////////
 TCanvas *c = new TCanvas("c","c", nqbins*500, 400);
 TH2F *h = new TH2F("h","h",default_binning, 2.2, 4.2, default_binning, 0.995, 1.04);
 h->GetXaxis()->SetTitle("#nu, GeV");
 h->GetYaxis()->SetTitle("C_corr");
 h->GetYaxis()->SetTitleOffset(1.4);
 h->GetXaxis()->SetTitleOffset(0.9);
 h->GetXaxis()->SetTitleSize(0.055);
 h->GetYaxis()->SetTitleSize(0.04);
 h->GetYaxis()->SetLabelSize(0.04);
 h->GetXaxis()->SetLabelSize(0.04);
 c->Divide(nqbins,1, 0.00001, 0.00001);
 for(int nq=0; nq<nqbins;nq++){
   c->cd(nq+1)->SetGridy();
  g_RC_D[nq] = new TGraphErrors(nvbins, center_v, C_corr_D[nq],0, 0);
  g_RC_C[nq] = new TGraphErrors(nvbins, center_v, C_corr_C[nq],0, 0);
  g_RC_Fe[nq] = new TGraphErrors(nvbins, center_v, C_corr_Fe[nq],0, 0);
  g_RC_Pb[nq] = new TGraphErrors(nvbins, center_v, C_corr_Pb[nq],0,0);
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
	 
  auto legend = new TLegend(0.12,0.89,0.24,0.69);
  legend->AddEntry("D","D","p");
  legend->AddEntry("C","C","p");
  legend->AddEntry("Fe","Fe","p");
  legend->AddEntry("Pb","Pb","p");
	        legend->Draw();
 }
 c->Print("~/secure/ExternalsRC/rc_Coulomb.png");
 c->Clear();

/////////////// 1D Nu Coulomb/////////////////
 TCanvas *c1 = new TCanvas("c1","c1", 500, 400);
 TH2F *h1 = new TH2F("h1","h1",default_binning, 2.2, 4.2, default_binning, 0.995, 1.04);
 h1->GetXaxis()->SetTitle("#nu, GeV");
 h1->GetYaxis()->SetTitle("C_Corr");
 h1->GetYaxis()->SetTitleOffset(1.35);
 h1->GetXaxis()->SetTitleOffset(0.9);
 h1->GetXaxis()->SetTitleSize(0.04);
 h1->GetYaxis()->SetTitleSize(0.04);
 h1->GetYaxis()->SetLabelSize(0.04);
 h1->GetXaxis()->SetLabelSize(0.04);
 c1->Divide(nvbins,1, 0.00001, 0.00001);
 c1->SetGridy();

 g_Coulomb_vD = new TGraphErrors(nvbins, center_v, C_corr_vD,0, 0); 
 g_Coulomb_vC = new TGraphErrors(nvbins, center_v, C_corr_vC,0, 0); 
 g_Coulomb_vFe = new TGraphErrors(nvbins, center_v, C_corr_vFe,0, 0); 
 g_Coulomb_vPb = new TGraphErrors(nvbins, center_v, C_corr_vPb,0,0);
 g_Coulomb_vD->SetName("vD");
 g_Coulomb_vC->SetName("vC");
 g_Coulomb_vFe->SetName("vFe");
 g_Coulomb_vPb->SetName("vPb");
 
 h1->Draw();
 g_Coulomb_vD->SetMarkerColor(1);
 g_Coulomb_vD->SetMarkerStyle(20);
 g_Coulomb_vD->SetMarkerSize(1.7);
 g_Coulomb_vD->Draw("Psame");
 g_Coulomb_vC->SetMarkerColor(2);
 g_Coulomb_vC->SetMarkerStyle(20);
 g_Coulomb_vC->SetMarkerSize(1.7);
 g_Coulomb_vC->Draw("Psame");
 g_Coulomb_vFe->SetMarkerColor(4);
 g_Coulomb_vFe->SetMarkerStyle(20);
 g_Coulomb_vFe->SetMarkerSize(1.7);
 g_Coulomb_vFe->Draw("Psame");
 g_Coulomb_vPb->SetMarkerColor(8);
 g_Coulomb_vPb->SetMarkerStyle(20);
 g_Coulomb_vPb->SetMarkerSize(1.7);
 g_Coulomb_vPb->Draw("Psame");
 TLatex a;
 a.SetNDC();
 a.SetTextSize(0.06);
 a.DrawLatex(0.33,0.92, "1.0 < Q^{2} < 4.1 GeV^{2}");

  auto legend = new TLegend(0.12,0.89,0.24,0.69); 
  legend->AddEntry("vD","D","p");
  legend->AddEntry("vC","C","p");
  legend->AddEntry("vFe","Fe","p");
  legend->AddEntry("vPb","Pb","p");
  legend->Draw();
  c1->Print("~/secure/ExternalsRC/rc_vCoulomb_integQ2.png");
  c1->Clear();  

//////////////// BORN cross section //////////////
 //TCanvas *c2 = new TCanvas("c2","c2", 2*500, 400);
 TCanvas *c2 = new TCanvas("c2","c2", 500, 400);
 TH2F *h2 = new TH2F("h2","h2",default_binning, 2.2, 4.2, default_binning, 0.,10000.);
 h2->GetXaxis()->SetTitle("#nu, GeV");
 h2->GetYaxis()->SetTitle("#sigma(Pb Born), [nb/sr/MeV]");
 h2->GetYaxis()->SetTitleOffset(1.35);
 h2->GetXaxis()->SetTitleOffset(0.9);
 h2->GetXaxis()->SetTitleSize(0.04);
 h2->GetYaxis()->SetTitleSize(0.04);
 h2->GetYaxis()->SetLabelSize(0.035);
 h2->GetXaxis()->SetLabelSize(0.04);
 //c2->Divide(2,1, 0.001, 0.001);
 c2->SetGridy();
 h2->Draw();
 auto legend = new TLegend(0.7,0.9,0.9,0.65);

 for(int nq=0; nq<nqbins;nq++){
 g_Born_Pb[nq] = new TGraphErrors(nvbins, center_v, sigma_Born_Pb[nq],0, 0);
 g_Born_Pb[nq]->SetMarkerColor(nq+1);
 g_Born_Pb[nq]->SetMarkerStyle(20);
 g_Born_Pb[nq]->SetMarkerSize(1.7);
 g_Born_Pb[nq]->Draw("Psame");
 if(nq==0){g_Born_Pb[nq]->SetName("BornPb0");
	   legend->AddEntry("BornPb0",Form("%3.2f<Q^{2}<%3.2f ", q_values[nq], q_values[nq+1]),"p");}
 if(nq==1){g_Born_Pb[nq]->SetName("BornPb1");
	   legend->AddEntry("BornPb1",Form("%3.2f<Q^{2}<%3.2f ", q_values[nq], q_values[nq+1]),"p");}
 if(nq==2){g_Born_Pb[nq]->SetName("BornPb2");
	   legend->AddEntry("BornPb2",Form("%3.2f<Q^{2}<%3.2f ", q_values[nq], q_values[nq+1]),"p");}
 }  
 legend->Draw();
 c2->Print("~/secure/ExternalsRC/rc_Born_Pb.png");
 c2->Clear();  

 c2->SetGridy();
 h2->Draw();
 g_Born_vPb = new TGraphErrors(nvbins, center_v, sigma_Born_avePb,0, 0); 
 g_Born_vPb->SetMarkerColor(kBlue);
 g_Born_vPb->SetMarkerStyle(20);
 g_Born_vPb->SetMarkerSize(1.7);
 g_Born_vPb->Draw("Psame");
 a.DrawLatex(0.27,0.92, "1.0 < Q^{2} < 4.1 GeV^{2}");

 c2->Print("~/secure/ExternalsRC/rc_Born_avePb.png");
 c2->Clear();  

}
