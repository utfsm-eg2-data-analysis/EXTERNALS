/*****************************************/
/*  GetCentroids.cxx                     */
/*                                       */
/*  Andrés Bórquez                       */
/*                                       */
/*****************************************/

// This program gets the centroids of (Nu, Q2) bins based on a (Nu, Q2) binning given by a CSV file
// June 2021

#include "Headers.hxx"
#include "UX.hxx"

int main(int argc, char **argv) {

  gProgram = "GetCentroids";

  /*** INPUT ***/

  parseCommandLine(argc, argv);
  printOptions();

  TChain *InputChain = new TChain();

  TCut CutDIS = "Q2 > 1. && W > 2. && Yb < 0.85";
  TCut CutVertex;
  if (gTargetOption == "D") {  // unified D
    InputChain->Add(gWorkDir + "/out/GetSimpleTuple/data/C/*.root/ntuple_e");
    InputChain->Add(gWorkDir + "/out/GetSimpleTuple/data/Fe/*.root/ntuple_e");
    InputChain->Add(gWorkDir + "/out/GetSimpleTuple/data/Pb/*.root/ntuple_e");
    CutVertex = "TargType == 1 && vyec > -1.4 && vyec < 1.4";  // GST format
  } else {
    InputChain->Add(gWorkDir + "/out/GetSimpleTuple/data/" + gTargetOption + "/*.root/ntuple_e");
    CutVertex = "TargType == 2 && vyec > -1.4 && vyec < 1.4";  // GST format
  }

  /*** READ CSV FILE ***/

  std::ifstream BinningFile;
  BinningFile.open("binning.csv", std::ios::in);
  std::vector<Double_t> Nu_low;
  std::vector<Double_t> Nu_up;
  std::vector<Double_t> Q2_low;
  std::vector<Double_t> Q2_up;

  // exit program if ifstream could not open file
  if (!BinningFile) {
    std::cerr << "ERROR: File could not be opened" << std::endl;
    exit(EXIT_FAILURE);
  }

  TString auxLine;
  TString auxString[4];
  while (BinningFile >> auxLine) {
    for (Int_t i = 0; i < 4; i++) {
      auxString[i] = ((TObjString *)(auxLine.Tokenize(",")->At(i)))->String();
    }
    Nu_low.push_back(auxString[0].Atof());
    Nu_up.push_back(auxString[1].Atof());
    Q2_low.push_back(auxString[2].Atof());
    Q2_up.push_back(auxString[3].Atof());
  }

#ifdef DEBUG
  for (Size_t i = 0; i < Nu_low.size(); i++) {
    std::cout << Nu_low[i] << " - " << Nu_up[i] << ", " << Q2_low[i] << " - " << Q2_up[i] << std::endl;
  }
  std::cout << std::endl;
#endif

  /*** HISTOGRAMS ***/

#ifdef DEBUG
  TFile *RootOutputFile = new TFile("debug.root", "RECREATE");
#endif

  Int_t NTotalBins = (Int_t)Nu_low.size();
  TCut CutBin;
  Double_t meanX[NTotalBins];
  Double_t meanY[NTotalBins];
  TH2D *h[NTotalBins];
  for (Int_t i = 0; i < NTotalBins; i++) {
    CutBin = Form("Q2 > %.2f && Q2 < %.2f && Nu > %.2f && Nu < %.2f", Q2_low[i], Q2_up[i], Nu_low[i], Nu_up[i]);
    InputChain->Draw(Form("Q2:Nu>>h_%d(200, 2.2, 4.2, 300, 1.0, 4.0)", i), CutDIS && CutVertex && CutBin, "goff");
    h[i] = (TH2D *)gDirectory->GetList()->FindObject(Form("h_%d", i));
    meanX[i] = h[i]->GetMean(1);
    meanY[i] = h[i]->GetMean(2);
#ifdef DEBUG
    h[i]->Write();
    std::cout << meanX[i] << ", " << meanY[i] << std::endl;
#endif
  }
#ifdef DEBUG
  RootOutputFile->Close();
#endif

  /*** OUTPUT FILE ***/

  std::ofstream OutFile("centroids_" + gTargetOption + ".txt", std::ios::out);
  OutFile << "('Nu','Q2') (";
  for (Int_t i = 0; i < NTotalBins; i++) {
    OutFile << "(" << meanX[i] << "," << meanY[i] << ")";
    if (i + 1 == NTotalBins)
      OutFile << ")";
    else
      OutFile << ",";
  }
  OutFile.close();
}