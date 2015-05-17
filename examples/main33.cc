// main33.cc is a part of the PYTHIA event generator.
// Copyright (C) 2015 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.
// Author: Philip Ilten.

// An example where the HVQ POWHEGBOX matrix element binary is
// interfaced directly with PYTHIA. For this example to run correctly
// PYTHIA must be configured with the
// --with-powheg-bin=<path to directory containing only POWHEG binaries>
// option. This will build plugin libraries of the name
// libpythia8powheg<binary name>.so in the library directory.
// For these plugin libraries to build correctly, special compiler flags
// must have been used when building the POWHEGBOX binaries. These are
// "-rdynamic -fPIE -fPIC -pie". The following SED command will correctly
// insert them into the relevant POWHEGBOX Makefile:
//     sed -i "s/F77= gfortran/F77= gfortran -rdynamic -fPIE -fPIC -pie/g"
//     Makefile
// For this specific example the library libpythia8powheghvq.so must
// have been built.

#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/PowhegProcs.h"

using namespace Pythia8;

int main() {

  Pythia pythia;

  // The constructor PowhegProcs(process name, pythia, run directory,
  // PDF filename) where run directory by default is ./powhegrun and
  // the PDF filename is empty. If using native PDFs rather than
  // LHAPDF, the PDF filename is the full name of the PDF file to copy
  // to the run directory.
  PowhegProcs procs(&pythia, "hvq");

  // The commands readFile and readString are used to configure the
  // POWHEG matrix element. If a setting is repeated a warning is
  // issued and the most recent setting is used.
  procs.readFile("main33.cmnd");

  // This init call must be made before PYTHIA is initialized. It
  // copies the POWHEG input and PDF file to the POWHEG run directory.
  procs.init();

  // Further configuration of PYTHIA can still be performed. By
  // default the user hooks are set to a PowhegHooks instance and the
  // following configuration is passed to PYTHIA.
  //     POWHEG:nFinal = 2
  //     POWHEG:veto = 1
  //     POWHEG:vetoCount = 3
  //     POWHEG:pThard = 2
  //     POWHEG:pTemt = 0
  //     POWHEG:emitted = 0
  //     POWHEG:pTdef = 1
  //     POWHEG:MPIveto = 0
  //     POWHEG:QEDveto = 2
  //     Beams:frameType = 5
  // Note that POWHEG:nFinal should be changed to whatever is
  // appropriate for the matrix element being used.
  pythia.init();

  // Run PYTHIA. The random numbers are taken from the associated
  // PYTHIA random number generator.
  for (int iEvent = 0; iEvent < 100; ++iEvent) pythia.next();

  pythia.stat();
  return 0;
}
