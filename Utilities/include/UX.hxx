#ifndef UX_HXX
#define UX_HXX

#ifndef HEADERS_HXX
#include "Headers.hxx"
#endif

/*** Global variables ***/

TString gWorkDir = getenv("WORKDIR");

TString gTargetOption;
TString gRunNumber;

TString gProgram = "";

/*** Input-related functions ***/

void printUsage() {
  std::cout << gProgram << " program. Usage is:" << std::endl;
  std::cout << std::endl;
  std::cout << "./" << gProgram << " -h" << std::endl;
  std::cout << "    prints this message and exits program" << std::endl;
  std::cout << std::endl;
  if (gProgram == "GetCentroids") {
    std::cout << "./" << gProgram << " -t[D, C, Fe, Pb]" << std::endl;
    std::cout << "    filters respective target" << std::endl;
    std::cout << std::endl;
  }
}

void parseCommandLine(int argc, char* argv[]) {
  Int_t c;
  if (argc == 1) {
    std::cerr << "Empty command line. Execute ./" << gProgram << " -h to print help." << std::endl;
    exit(0);
  }
  while ((c = getopt(argc, argv, "ht:")) != -1) switch (c) {
      case 'h':
        printUsage();
        exit(0);
        break;
      case 't':
        gTargetOption = optarg;
        break;
      default:
        std::cerr << "Unrecognized argument. Execute ./" << gProgram << " -h to print help." << std::endl;
        exit(0);
        break;
    }
}

void printOptions() {
  std::cout << "Executing " << gProgram << " program. The chosen parameters are: " << std::endl;
  std::cout << "  gTargetOption   = " << gTargetOption << std::endl;
  std::cout << std::endl;
}

#endif
