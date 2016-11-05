// GammaKinematics.h is a part of the PYTHIA event generator.
// Copyright (C) 2016 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Header file for kinematics selection of photons from lepton beams.

#ifndef Pythia8_GammaKinematics_H
#define Pythia8_GammaKinematics_H

#include "Pythia8/Basics.h"
#include "Pythia8/BeamParticle.h"
#include "Pythia8/Info.h"
#include "Pythia8/Settings.h"

namespace Pythia8 {

//==========================================================================

// Class to sample the virtuality and transverse momentum of emitted photons.

class GammaKinematics {

public:

  // Constructor.
  GammaKinematics() {}

  // Sample the trial or final event kinematics.
  bool init(Info* infoPtrIn, Settings* settingsPtrIn, Rndm* rndmPtrIn,
    BeamParticle* beamAPtrIn, BeamParticle* beamBPtrIn);
  bool sampleKTgamma();
  bool finalize();

  // Methods to pass along the sampled values.
  double getQ2gamma1()   const {return Q2gamma1;}
  double getQ2gamma2()   const {return Q2gamma2;}
  double getQ2min1()     const {return Q2min1;}
  double getQ2min2()     const {return Q2min2;}
  double getPhi1()       const {return phi1;}
  double getPhi2()       const {return phi2;}
  double getKT1()        const {return kT1;}
  double getKT2()        const {return kT2;}
  double eCMsub()        const {return mGmGm;}

private:

  // Pointer to various information on the generation.
  Info*         infoPtr;

  // Pointer to the settings database.
  Settings*     settingsPtr;

  // Pointer to the random number generator.
  Rndm*         rndmPtr;

  // Pointers to incoming beams.
  BeamParticle* beamAPtr;
  BeamParticle* beamBPtr;

  // Kinematics variables.
  double Q2maxGamma, Wmin, Wmax, eCM, sCM, m2BeamA, m2BeamB, m2sA, m2sB,
         Q2min1, Q2min2, xGamma1, xGamma2, Q2gamma1, Q2gamma2, phi1, phi2,
         kT1, kT2, mGmGm, m2GmGm, theta1, theta2, theta1Max, theta2Max;
};

//==========================================================================

} // end namespace Pythia8

#endif // Pythia8_GammaKinematics_H
