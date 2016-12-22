// main69.cc is a part of the PYTHIA event generator.
// Copyright (C) 2016 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Main program to generate charged hadron spectra from photon-initiated 
// processes, by combining four sub-runs with direct or resolved photons.

#include "Pythia8/Pythia.h"

using namespace Pythia8;

int main() {

  // Generator.
  Pythia pythia("../share/Pythia8/xmldoc", false);

  // Decrease the output.
  pythia.readString("Init:showChangedSettings = off");
  pythia.readString("Init:showChangedParticleData = off");
  pythia.readString("Next:numberCount = 0");
  pythia.readString("Next:numberShowInfo = 0");
  pythia.readString("Next:numberShowProcess = 0");
  pythia.readString("Next:numberShowEvent = 0");

  // Shorthand for some public members of pythia (also static ones).
  Settings& settings = pythia.settings;
  Info& info = pythia.info;

  // Generate photon-photon events in leptonic or photon beams.
  bool photonsFromElectrons = true;
  bool biasSampling         = false;

  // Select gamma PDF: currently only one choice available.  
  pythia.readString("PDF:gammaSet = 1");
  pythia.readString("PDF:lepton2gammaSet = 1");

  // Optionally use different PDFs for hard process.
  // pythia.readString("PDF:useHard = on");
  // pythia.readString("PDF:GammaHardSet = LHAPDF5:SASG.LHgrid/5");
  // pythia.readString("PDF:GammaHardSet = LHAPDF5:GRVG0.LHgrid/1");

  // Beam parameters.
  pythia.readString("Beams:eCM = 200.");

  // Electron beam with photons inside it. 
  if ( photonsFromElectrons) {
    pythia.readString("Beams:idA = -11");
    pythia.readString("Beams:idB =  11");
    pythia.readString("PDF:lepton2gamma = on");

    // Cuts on photon virtuality and invariant mass of gamma-gamma pair.
    pythia.readString("Photon:Q2max = 1.0");
    pythia.readString("Photon:Wmin  = 10.0");
    // pythia.readString("Photon:Wmax  = 125.0");

  // Fixed-energy photon beams.
  } else {
    pythia.readString("Beams:idA = 22");
    pythia.readString("Beams:idB = 22");
  }

  // Number of events per run.
  int nEvent = 10000;

  // Kinematic limits and histogram limits.  
  double pTcut = 5.0;
  double pT0   = 0.0;
  double pTmin = 0.0;
  double pTmax = 40.0;
  int nBinsPT  = 40;

  Hist pTtot("Total charged hadron pT distribution", nBinsPT, pTmin, pTmax);
  Hist pTsoft("Soft QCD contribution from resolved", nBinsPT, pTmin, pTmax);
  Hist pThard("Hard QCD contribution from resolved", nBinsPT, pTmin, pTmax);
  Hist pTresdir("Resolved-direct contribution", nBinsPT, pTmin, pTmax);
  Hist pTdirres("Direct-resolved contribution", nBinsPT, pTmin, pTmax);
  Hist pTdirdir("Direct-direct contribution", nBinsPT, pTmin, pTmax);
  Hist pTiRun("Contribution from Run i", nBinsPT, pTmin, pTmax);

  // Loop over relevant processes.
  for ( int iRun = 0; iRun < 5; ++iRun) {

    // Set the type of gamma-gamma process:
    // 1 = resolved-resolved
    // 2 = resolved-direct
    // 3 = direct-resolved
    // 4 = direct-direct
    int photonMode = iRun == 0 ? 1 : iRun;
    settings.mode("Photon:ProcessType", photonMode);

    // First run: soft and QCD processes for resolved photons.
    if ( iRun == 0 ) {
      pythia.readString("SoftQCD:nonDiffractive = on");
      settings.parm("PhaseSpace:pTHatMin", pT0);
   
    // Second run: hard QCD processes for resolved photons.
    } else if ( iRun == 1 ) {
      pythia.readString("HardQCD:all = on");
      pythia.readString("SoftQCD:nonDiffractive = off");
      settings.parm("PhaseSpace:pTHatMin", pTcut);

      // Hard processes with pT-weight.
      if (biasSampling) {
        pythia.readString("PhaseSpace:bias2Selection = on");
        pythia.readString("PhaseSpace:bias2SelectionPow = 3.");
        settings.parm("PhaseSpace:bias2SelectionRef", pTcut);
      }

    // Third and fourth run: direct-resolved processes.
    } else if ( (iRun == 2) || (iRun == 3) ) {
      pythia.readString("HardQCD:all = off");
      pythia.readString("PhotonParton:all = on");
      settings.parm("PhaseSpace:pTHatMin", pT0);
      pythia.readString("PartonLevel:MPI = off");

    // Fifth run: direct-direct QCD processes.
    } else {
      pythia.readString("PhotonParton:all = off");
      pythia.readString("PhotonCollision:gmgm2qqbar = on");
      pythia.readString("PhotonCollision:gmgm2ccbar = on");
      pythia.readString("PhotonCollision:gmgm2bbbar = on");

      // No need for pT weight for direct-direct gm+gm collisions.
      if (!photonsFromElectrons)
        pythia.readString("PhaseSpace:bias2Selection = off");
    }

    // Initialize the generator.
    pythia.init();

    // Clear the histogram.
    pTiRun.null();

    // Begin event loop. Skip if fails.
    for (int iEvent = 0; iEvent < nEvent; ++iEvent) {

      // Generate next event.
      if (!pythia.next()) continue;

      // List the first process and event for each contribution.
      if (iEvent == 0) {
        pythia.process.list();
        pythia.event.list();
      }

      // Possible event weights.
      double weight = info.weight();

      // Store pThat and discard soft processes that overlaps with hard ones.
      double pThat = info.pTHat();
      if ( iRun == 0 && pThat > pTcut ) continue;

      // Loop over event record and find charged final state particles.
      for (int i = 0; i < pythia.event.size(); ++i){
        if ( pythia.event[i].isFinal() && pythia.event[i].isCharged() ){

          // Store pT and add pT to the histogram.
          double pTch = pythia.event[i].pT();
          pTiRun.fill(pTch, weight);
        }
      }

    }

    // Show statistics after each run (erorrs cumulate).
    pythia.stat();

    // Normalize to cross section [mb].
    double sigmaNorm = info.sigmaGen() / info.weightSum();
    double pTBin     = (pTmax - pTmin) / (1. * nBinsPT);
    pTiRun          *= sigmaNorm / pTBin;

    if (iRun == 0) pTsoft   = pTiRun;
    if (iRun == 1) pThard   = pTiRun;
    if (iRun == 2) pTresdir = pTiRun;
    if (iRun == 3) pTdirres = pTiRun;
    if (iRun == 4) pTdirdir = pTiRun;
    pTtot += pTiRun;

  // End of loop over runs.
  }

  // Print histograms.
  cout << pTsoft << pThard << pTresdir << pTdirres << pTdirdir << pTtot;
  
  // Done.
  return 0;
}
