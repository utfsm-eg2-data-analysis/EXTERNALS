#define ec_gam_cxx
#include "ec_gam.h"
#include <TH2.h>
#include <TStyle.h>

void ec_gam::Begin(TTree *tree)
{
   TString option = GetOption();
   Init(tree);
   Nevta = 0;
   Nevt_found = 0;
   if(AvailableEvents==0)AvailableEvents = Int_t(tree->GetEntries());

   cout << "*******************************HELLO! ***********************************" << endl;
  
   cout << "Activating prefered styles... "; cout.flush();
   gROOT->SetStyle("Plain");
   gStyle->SetPalette(1);
  // gStyle->SetOptFit(1111);
  // gStyle->SetOptStat(11);
  // gStyle->SetOptTitle(1);
   cout << "Done !" << endl << "Allocations/Initializations"; cout.flush();
   cout << endl << "The time is now ";cout.flush(); gSystem->Exec("date");cout<<endl;
  
   NewFile = new TFile("NewFile.root","recreate");
   TreeGam = new TTree("ec_gam","Tree with electron variables");
   MakeBranch();
   cout<< "Branch made" <<endl;  

   BookHists();

   cout << "Done !" << endl;
   cout << "looping" << endl;
}

void ec_gam::SlaveBegin(TTree *tree)
{
   TString option = GetOption();
   Init(tree);
}

Bool_t ec_gam::Process(Long64_t entry)
{
   fChain->GetTree()->GetEntry(entry);
   TString option = GetOption();
   Nevta++;
  
   Int_t evtsPrint = Int_t(TMath::Floor(Nevta/(AvailableEvents/100)));
   evtsPrint = evtsPrint*Int_t(AvailableEvents/100);
   if(Nevta==evtsPrint){
   cout << "\r"<<TMath::FloorNint(Nevta/(AvailableEvents/100.))<<"\%";cout.flush();
   }
   
   Calculate();
   FillTables();
   TreeGam->Fill();

   return kTRUE;
}

void ec_gam::SlaveTerminate()
{

}

void ec_gam::Terminate()
{

  WriteHists();
  cout<<" hist written " << endl;
  WriteTxt();
  cout<< " Txt is written " << endl;
  NewFile->Close();
  cout<<"file is closed"<<endl;
  

  cout << "Over. Ciao, baby!" << endl;
}


void ec_gam::Calculate(void){
// here calculate Eprime and Theta
 Eprime = Ebeam - trig_nu;
 Theta = 2*(TMath::ASin(TMath::Sqrt(trig_Q2/(4*Ebeam*Eprime)))) ;  
 Theta = TMath::RadToDeg()*Theta;
 
// cout<< "trig_nu = "<< trig_nu<< " Eprime=" << Eprime << " Theta=" << Theta<<endl;

}

void ec_gam::BookHists(void){
 qbin = -1;
 vbin = -1;
 h_etheta = new TH2F("h_etheta","Theta_Eprime",200,0.5,3.,200,15,55.); 
 for(int  nq=0;  nq<nqbins; nq++){
 for(int  nv=0; nv<nvbins;nv++){
 h_etheta_b[nq][nv] = new TH2F(Form("h_etheta_q%dv%d",nq,nv),Form("h_etheta_q%dv%d",nq,nv),200,0.5,3.,200,15,55.);
   }
  }

}


void ec_gam::FillTables(void){
 qbin = Get_q_bin(trig_Q2);
 vbin =  Get_v_bin(trig_nu); 
 if(qbin>-1&&qbin<nqbins &&vbin>-1&&vbin<nvbins){
  h_etheta->Fill(Eprime, Theta);
  h_etheta_b[qbin][vbin]->Fill(Eprime, Theta);
 }

}

void ec_gam::WriteHists(void){
cout<<endl;
///////////////////////// Binned //////////////////////////
 gStyle->SetOptTitle(11111);
 TCanvas *c1 = new TCanvas("c1","c1",300*9,200*7);
 c1->Divide(9,7, 0.0001,0.0001);
 for(int  nq=0;  nq<nqbins; nq++){
   for(int  nv=0; nv<nvbins;nv++){
    c1->cd(nq+1+9*nv);
    h_etheta_b[nq][nv]->Draw("colz");
    h_etheta_b[nq][nv]->GetXaxis()->SetTitle("E'(GeV)");
    h_etheta_b[nq][nv]->GetYaxis()->SetTitle("#theta(deg)");
    
    mean_theta[nq][nv] = h_etheta_b[nq][nv]->GetMean(2);
    mean_Eprime[nq][nv] =h_etheta_b[nq][nv]-> GetMean(1);    
     
   // cout << nq << " " << nv << " <Q2>="<< center_q[nq] << " <v>="<<center_v[nv] << " <theta>="<< h_etheta_b[nq][nv]->GetMean(2)
      //                      << " <E'>="<< h_etheta_b[nq][nv]-> GetMean(1)<<endl;
    }
   } 
 c1->Print("~/secure/ExternalsRC/ThetaEprime_Q2Nu97_binned.png");
 c1->Clear();
}
////////////////// TREE /////////////////////////////////

void ec_gam::MakeBranch(void){
 if(TreeGam!=NULL){
  if(qbin>-1&&qbin<nqbins &&vbin>-1&&vbin<nvbins){
  TreeGam->Branch("trig_theta",&trig_theta,"trig_theta/F");
  TreeGam->Branch("trig_Q2",&trig_Q2,"trig_Q2/F");
  TreeGam->Branch("trig_nu",&trig_nu,"trig_nu/F");
  TreeGam->Branch("trig_mom",&trig_mom,"trig_mom/F");
  TreeGam->Branch("Theta", &Theta,"Theta/F");
  TreeGam->Branch("Eprime", &Eprime,"Eprime/F");
  }
 }
 
 return;
}

void ec_gam::WriteTxt(void){
ofstream outfile;
outfile.open("qv97_Eprime_Theta.txt");
for(int  nq=0;  nq<nqbins; nq++){
  for(int  nv=0; nv<nvbins;nv++){
     //outfile << nq << " " << nv<< " " << mean_theta[nq][nv]<< " "<<mean_Eprime[nq][nv]<<endl;
     outfile << "5.014 "<< mean_Eprime[nq][nv]<< " " << mean_theta[nq][nv]<<endl;
   }
  }
outfile.close();

} // writeTXT

/////////////////////////////////////////////////////////
int  ec_gam::Get_q_bin(float q){
  if(q<qmin||q>qmax)return -1;
   for(int nq=0; nq<nqbins; nq++){
    if(q>q_values[nq]&&q<q_values[nq+1]) return nq;
		      }
   return -1;
}

 int  ec_gam::Get_v_bin(float v){
  if(v<vmin||v>vmax)return -1;
   for(int nv=0; nv<nvbins; nv++){
    if(v>v_values[nv]&&v<v_values[nv+1]) return nv;
		       }
   return -1;
 }



