/*****************************************/
/*  GetBinning.cxx                       */
/*                                       */
/*  Andrés Bórquez                       */
/*                                       */
/*****************************************/

// This program prints the binning of (Nu, Q2) in a CSV file
// October 2021

#include "Binning.hxx"
#include "Headers.hxx"
#include "UX.hxx"

int main(int argc, char **argv) {

  gProgram = "GetBinning";

  /*** INPUT ***/

  parseCommandLine(argc, argv);
  printOptions();

  /*** WRITE OUTPUT FILE ***/

  TString OutputFilename = "binning_" + gParticleOption + ".csv";
  std::ofstream OutFile(OutputFilename, std::ios::out);

  // related to particle option - now, get edges
  Double_t EdgesQ2[NbinsQ2 + 1];
  Double_t EdgesNu[NbinsNu + 1];
  if (gParticleOption == "omega") {
    for (Int_t q = 0; q < NbinsQ2 + 1; q++) {
      EdgesQ2[q] = kEdgesQ2_Omega[q];
    }
    for (Int_t n = 0; n < NbinsNu + 1; n++) {
      EdgesNu[n] = kEdgesNu_Omega[n];
    }
  } else {  // "eta"
    for (Int_t q = 0; q < NbinsQ2 + 1; q++) {
      EdgesQ2[q] = kEdgesQ2_Eta[q];
    }
    for (Int_t n = 0; n < NbinsNu + 1; n++) {
      EdgesNu[n] = kEdgesNu_Eta[n];
    }
  }

  std::cout << "Binning : (Nu, Q2) = " << NbinsNu << " x " << NbinsQ2 << std::endl;

  for (Int_t n = 0; n < NbinsNu; n++) {
    for (Int_t q = 0; q < NbinsQ2; q++) {
      OutFile << EdgesNu[n] << "," << EdgesNu[n + 1] << "," << EdgesQ2[q] << "," << EdgesQ2[q + 1] << std::endl;
    }
  }

  // close output file
  OutFile.close();

  std::cout << "The following file has been created: " << OutputFilename << std::endl;

  return 0;
}
