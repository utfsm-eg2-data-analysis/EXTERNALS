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


Float_t Eprime[nqbins][nvbins], Theta[nqbins][nvbins], Xb[nqbins][nvbins];

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

Float_t sigma_Born_vD[nvbins], sigma_Born_vC[nvbins],sigma_Born_vFe[nvbins],sigma_Born_vPb[nvbins];
Float_t sigma_Rad_vD[nvbins], sigma_Rad_vC[nvbins],sigma_Rad_vFe[nvbins],sigma_Rad_vPb[nvbins];
Float_t RC_vD[nvbins], RC_vC[nvbins],RC_vFe[nvbins],RC_vPb[nvbins];
Float_t  ratioRC_vCD[nvbins],ratioRC_vFeD[nvbins],ratioRC_vPbD[nvbins];

Float_t  RC_Pb0[nqbins], RC_Pb1[nqbins], RC_Pb2[nqbins];
Float_t Xb0[nqbins], Xb1[nqbins], Xb2[nqbins];
Float_t Q20[nqbins], Q21[nqbins], Q22[nqbins];

TGraphErrors *g_RC_D[nqbins],*g_RC_C[nqbins],*g_RC_Fe[nqbins],*g_RC_Pb[nqbins];
TGraphErrors *g_ratioRC_CD[nqbins], *g_ratioRC_FeD[nqbins], *g_ratioRC_PbD[nqbins];
TGraphErrors *g_rcQ2X_Pb0;
TGraphErrors *g_rcQ2X_Pb1;
TGraphErrors *g_rcQ2X_Pb2;

TGraphErrors *g_RC_vD,*g_RC_vC,*g_RC_vFe,*g_RC_vPb;
TGraphErrors *g_ratioRC_vCD, *g_ratioRC_vFeD, *g_ratioRC_vPbD;

void rc_xQ2(void){
 gROOT->SetStyle("Plain");
 gStyle->SetPalette(1);
 gStyle->SetOptFit(1);
 gStyle->SetOptStat(0);
 gStyle->SetOptTitle(1);
 gStyle->SetOptDate(0);

 ////////////////////// READ RC ////////////////////////
Float_t Ebeam, x, Q2[nqbins][nvbins]; 
 ifstream in;

 cout << " **************** D ********************** " << endl; 
 in.open("clasd2.out");
 for(int nq=0; nq<nqbins;nq++){
  for(int nv=0;  nv<nvbins; nv++){  
   in >> Ebeam
      >> Eprime[nq][nv]
      >> Theta[nq][nv]
      >> Xb[nq][nv] 
      >> Q2[nq][nv]
      >> sigma_Born_D[nq][nv]
      >> sigma_BornIn_D[nq][nv] 
      >> sigma_BornQE_D[nq][nv] 
      >> sigma_Rad_D[nq][nv]
      >> sigma_RadEl_D[nq][nv]
      >> sigma_RadQE_D[nq][nv] 
      >> sigma_RadDIS_D[nq][nv]
      >> C_corr_D[nq][nv];
   
   cout << " nq="<<nq<< " nv="<< nv << " Q2=" << Q2[nq][nv] 
        << " nu="<< Q2[nq][nv]/(2*0.938*Xb[nq][nv]) << endl; 
	 //  <<" "<< Ebeam << " " << Eprime[nq][nv]<< " " << Theta[nq][nv]<< " "<< x<< " " << Q2[nq][nv]
	 //  << " nu="<< Q2[nq][nv]/(2*0.938*Xb[nq][nv]) << endl;
        /* <<" Crosssection:"
	   << sigma_Born_D[nq][nv]<< " " 
	   << sigma_BornIn_D[nq][nv]<< " " << sigma_BornQE_D[nq][nv]<< " " << sigma_Rad_D[nq][nv]<< " " 
	   << sigma_RadEl_D[nq][nv]<< " "<< sigma_RadQE_D[nq][nv]
	   << " " << sigma_RadDIS_D[nq][nv]<< " " <<C_corr_D[nq][nv]<< endl;
  */
    }
  }
 in.close();

 cout << " **************** C ********************** " << endl; 
 in.open("clasC12.out");
for(int nq=0; nq<nqbins;nq++){
 for(int nv=0;  nv<nvbins; nv++){  
   in >> Ebeam >> Eprime[nq][nv] >> Theta[nq][nv] >> Xb[nq][nv] >> Q2[nq][nv] >> sigma_Born_C[nq][nv] >> sigma_BornIn_C[nq][nv] >> sigma_BornQE_C[nq][nv]
      >> sigma_Rad_C[nq][nv] >>  sigma_RadEl_C[nq][nv] >> sigma_RadQE_C[nq][nv] >> sigma_RadDIS_C[nq][nv] >> C_corr_C[nq][nv];
   
 /*  cout << Ebeam << " " << Eprime[nq][nv]<< " " << Theta[nq][nv]<< " "<< x<< " " << Q2[nq][nv] << " " << sigma_Born_C[nq][nv]<< " " 
	   << sigma_BornIn_C[nq][nv]<< " " << sigma_BornQE_C[nq][nv]<< " " << sigma_Rad_C[nq][nv]<< " " 
	   << sigma_RadEl_C[nq][nv]<< " "<< sigma_RadQE_C[nq][nv]
	   << " " << sigma_RadDIS_C[nq][nv]<< " " <<C_corr_C[nq][nv]<< endl;
   */
   }
  }
 in.close();

 cout << " **************** Fe ********************** " << endl; 
 in.open("clasFe56.out");
for(int nq=0; nq<nqbins;nq++){
 for(int nv=0;  nv<nvbins; nv++){  
   in >> Ebeam >> Eprime[nq][nv] >> Theta[nq][nv] >> Xb[nq][nv] >> Q2[nq][nv] >> sigma_Born_Fe[nq][nv] >> sigma_BornIn_Fe[nq][nv] >> sigma_BornQE_Fe[nq][nv]
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
   in >> Ebeam >> Eprime[nq][nv] >> Theta[nq][nv] >> Xb[nq][nv] >> Q2[nq][nv] >> sigma_Born_Pb[nq][nv] >> sigma_BornIn_Pb[nq][nv] >> sigma_BornQE_Pb[nq][nv]
      >> sigma_Rad_Pb[nq][nv] >>  sigma_RadEl_Pb[nq][nv] >> sigma_RadQE_Pb[nq][nv] >> sigma_RadDIS_Pb[nq][nv] >> C_corr_Pb[nq][nv];
   
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
   deltaq = (q_values[nq+1]-q_values[nq]); 
   sigma_Born_vD[nv] += sigma_Born_D[nq][nv]* deltaq;
   sigma_Rad_vD[nv] += sigma_Rad_D[nq][nv]* deltaq;
   sigma_Born_vC[nv] += sigma_Born_C[nq][nv]* deltaq;
   sigma_Rad_vC[nv] += sigma_Rad_C[nq][nv]* deltaq;
   sigma_Born_vFe[nv] += sigma_Born_Fe[nq][nv]* deltaq;
   sigma_Rad_vFe[nv] += sigma_Rad_Fe[nq][nv]* deltaq;
   sigma_Born_vPb[nv] += sigma_Born_Pb[nq][nv]* deltaq;
   sigma_Rad_vPb[nv] += sigma_Rad_Pb[nq][nv]* deltaq;
  
   RC_D[nq][nv]= 1/(sigma_Rad_D[nq][nv]/sigma_Born_D[nq][nv]);
   RC_C[nq][nv]= 1/(sigma_Rad_C[nq][nv]/sigma_Born_C[nq][nv]);
   RC_Fe[nq][nv]= 1/(sigma_Rad_Fe[nq][nv]/sigma_Born_Fe[nq][nv]);
   RC_Pb[nq][nv]= 1/(sigma_Rad_Pb[nq][nv]/sigma_Born_Pb[nq][nv]);

   ratioRC_CD[nq][nv] = RC_C[nq][nv]/RC_D[nq][nv];
   ratioRC_FeD[nq][nv] = RC_Fe[nq][nv]/RC_D[nq][nv];
   ratioRC_PbD[nq][nv] = RC_Pb[nq][nv]/RC_D[nq][nv];

   cout << nq <<" " <<nv<< " RC(D)" << RC_D[nq][nv] << " RC(Pb)" << RC_Pb[nq][nv]<< " "<<
	 "Ratio(Pb/D)" << ratioRC_PbD[nq][nv]<<endl; 
  }
 }

  for(int nv=0;  nv<nvbins; nv++){
   RC_vD[nv]= 1/(sigma_Rad_vD[nv]/sigma_Born_vD[nv]);
   RC_vC[nv]= 1/(sigma_Rad_vC[nv]/sigma_Born_vC[nv]);
   RC_vFe[nv]= 1/(sigma_Rad_vFe[nv]/sigma_Born_vFe[nv]);
   RC_vPb[nv]= 1/(sigma_Rad_vPb[nv]/sigma_Born_vPb[nv]);
   
   ratioRC_vCD[nv] = RC_vC[nv] /RC_vD[nv];
   ratioRC_vFeD[nv] = RC_vFe[nv] /RC_vD[nv];
   ratioRC_vPbD[nv] = RC_vPb[nv] /RC_vD[nv]; 

  // cout << nv << " " << " ratio C/D = " << ratioRC_vCD[nv] 
//	      << " ratio Fe/D = " << ratioRC_vFeD[nv]<< " ratio Pb/D = " << ratioRC_vPbD[nv]<<endl;
   }
 
 ofstream out;
 out.open("eRC_qv_DCFePb.txt");
 for(int nq=0;  nq<nqbins; nq++){
  for(int nv=0; nv<nvbins;nv++){
   out << nq << " " << nv << " " << RC_D[nq][nv] << " "
	                        <<  RC_C[nq][nv] << " "
				<<  RC_Fe[nq][nv] << " "
				<<  RC_Pb[nq][nv] << " " <<endl;
  }
 }
 out.close();
 
 out.open("eRC_qv_ratios2D.txt");
 for(int nq=0;  nq<nqbins; nq++){
  for(int nv=0; nv<nvbins;nv++){
   out << nq << " " << nv << " " << ratioRC_CD[nq][nv]
	                  << " " << ratioRC_FeD[nq][nv]
	                  << " " << ratioRC_PbD[nq][nv]<<endl;
  }
 }
 out.close();

 //////////////////////////////////////////////////////////////////////////////////////
/*
 out.open("eRC_v_D.txt");
 for(int nv=0; nv<nvbins;nv++){
  out << nv << " " << RC_vD[nv] <<endl;
  }
 out.close();

 out.open("eRC_v_Pb.txt");
 for(int nv=0; nv<nvbins;nv++){
  out << nv << " " << RC_vPb[nv] << " " << ratioRC_vPbD[nv] <<endl;
  }
 out.close();
*/


 ////////////////////////// Plot ////////////////////////
 cout << "Plotting" <<endl;  
 gStyle->SetOptTitle(0); 
 
 ////////////////////// 2D /////////////////////////////
 TCanvas *c = new TCanvas("c","c", nqbins*500, 400);
 TH2F *h = new TH2F("h","h",default_binning, 0.1, 0.4, default_binning, 0.7, 1.1);
 h->GetXaxis()->SetTitle("x");
 h->GetYaxis()->SetTitle("Rad Corr=#sigma^{Born}/#sigma^{Rad}");
 h->GetYaxis()->SetTitleOffset(0.99);
 h->GetXaxis()->SetTitleOffset(0.9);
 h->GetXaxis()->SetTitleSize(0.055);
 h->GetYaxis()->SetTitleSize(0.047);
 h->GetYaxis()->SetLabelSize(0.04);
 h->GetXaxis()->SetLabelSize(0.052);
 c->Divide(nqbins,1, 0.00001,0.00001);
 for(int nq=0; nq<nqbins;nq++){
  c->cd(nq+1)->SetGridy();
  g_RC_D[nq] = new TGraphErrors(nvbins, Xb[nq], RC_D[nq],0, 0);
  g_RC_C[nq] = new TGraphErrors(nvbins, Xb[nq], RC_C[nq],0, 0);
  g_RC_Fe[nq] = new TGraphErrors(nvbins, Xb[nq],RC_Fe[nq],0, 0);
  g_RC_Pb[nq] = new TGraphErrors(nvbins, Xb[nq],RC_Pb[nq],0,0);
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
	 
  auto legend = new TLegend(0.8,0.95,0.92,0.75);
  legend->AddEntry("D","D","p");
  legend->AddEntry("C","C","p");
  legend->AddEntry("Fe","Fe","p");
  legend->AddEntry("Pb","Pb","p");
	        legend->Draw();
 }
 c->Print("~/secure/ExternalsRC/rc_xQ2.png");
 c->Clear();

/////////////// see Q2 dependence /////////// 
 TCanvas *c = new TCanvas("c","c", 500, 400);
 TH2F *h1 = new TH2F("h1","h1",default_binning, 1., 2.4, default_binning, 0.7, 1.1);
 h1->GetXaxis()->SetTitle("Q^{2}");
 h1->GetYaxis()->SetTitle("Rad Corr=#sigma^{Born}/#sigma^{Rad}");
 h1->GetYaxis()->SetTitleOffset(1.2);
 h1->GetXaxis()->SetTitleOffset(0.9);
 h1->GetXaxis()->SetTitleSize(0.05);
 h1->GetYaxis()->SetTitleSize(0.04);
 h1->GetYaxis()->SetLabelSize(0.04);
 h1->GetXaxis()->SetLabelSize(0.04);
 c->SetGridy();

 for(int nq=0; nq<nqbins;nq++){
    RC_Pb0[nq]=RC_Pb[nq][0];
    RC_Pb1[nq]=RC_Pb[nq][1];
    RC_Pb2[nq]=RC_Pb[nq][2];
    Q20[nq]=Q2[nq][0];
    Q21[nq]=Q2[nq][1];
    Q22[nq]=Q2[nq][2];
 }
for(int nq=0; nq<nqbins;nq++){
 cout << " nq="<<nq <<  " RC_Pb0[nq]="<<RC_Pb0[nq]<<endl;
 cout <<" nq="<<nq << " RC_Pb1[nq]="<<RC_Pb1[nq]<<endl;
 cout << " nq="<<nq <<" RC_Pb2[nq]="<<RC_Pb2[nq]<<endl;
}
///// here: at higher Nu, correction goes UP!!! /////

 g_rcQ2X_Pb0 = new TGraphErrors(nvbins, Q20, RC_Pb0,0,0); // v=0
 g_rcQ2X_Pb1 = new TGraphErrors(nvbins, Q21,RC_Pb1,0,0); // v=1
 g_rcQ2X_Pb2 = new TGraphErrors(nvbins, Q22,RC_Pb2,0,0); //v=2
 
 auto legend = new TLegend(0.6,0.95,0.92,0.75);
  h1->Draw();
  g_rcQ2X_Pb0->SetMarkerColor(1);
  g_rcQ2X_Pb0->SetMarkerStyle(20);
  g_rcQ2X_Pb0->SetMarkerSize(1.7);
  g_rcQ2X_Pb0->Draw("Psame");
 
  g_rcQ2X_Pb0->SetName("Pbqx0");
  legend->AddEntry("Pbqx0", "<#nu>=2.88 GeV" ,"p");
  legend->Draw();

  g_rcQ2X_Pb1->SetMarkerColor(2);
  g_rcQ2X_Pb1->SetMarkerStyle(20);
  g_rcQ2X_Pb1->SetMarkerSize(1.7);
  g_rcQ2X_Pb1->Draw("Psame");
  g_rcQ2X_Pb1->SetName("Pbqx1");
  legend->AddEntry("Pbqx1","<#nu>=3.47 GeV","p");
  legend->Draw();

  g_rcQ2X_Pb2->SetMarkerColor(3);
  g_rcQ2X_Pb2->SetMarkerStyle(20);
  g_rcQ2X_Pb2->SetMarkerSize(1.7);
  g_rcQ2X_Pb2->Draw("Psame");
  g_rcQ2X_Pb2->SetName("Pbqx2");
  legend->AddEntry("Pbqx2", "<#nu>=3.99 GeV ","p");
  legend->Draw();
  
  cout<< " Q2[0]="<< Q20[0]<<endl;
  cout<< " Q2[1]="<< Q21[1]<<endl;
  cout<< " Q2[2]="<< Q22[2]<<endl;

 c->Print("~/secure/ExternalsRC/rc_Q2Nu_Pb.png");
 c->Clear();

}
