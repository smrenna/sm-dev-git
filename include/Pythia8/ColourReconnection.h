// ColourReconnection.h is a part of the PYTHIA event generator.
// Copyright (C) 2014 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Header file for the Colour reconnection handling.
// Reconnect the colours between the partons before hadronization.
// It Contains the following classes:
// ColourDipole, ColourParticle, ColourJunction, ColourReconnection.

#ifndef Pythia8_ColourReconnection_H
#define Pythia8_ColourReconnection_H

#include "Pythia8/Basics.h"
#include "Pythia8/BeamParticle.h"
#include "Pythia8/Event.h"
#include "Pythia8/FragmentationFlavZpT.h"
#include "Pythia8/Info.h"
#include "Pythia8/ParticleData.h"
#include "Pythia8/StringFragmentation.h"
#include "Pythia8/PartonDistributions.h"
#include "Pythia8/PartonSystems.h"
#include "Pythia8/PythiaStdlib.h"
#include "Pythia8/Settings.h"
#include "Pythia8/StringLength.h"

namespace Pythia8 {

//==========================================================================

// Contain a single colour chain. It always start from a quark and goes to 
// an anti quark or from an anti-junction to at junction 
// (or possible combinations).

class ColourDipole {

public:

  // Constructor.
  ColourDipole( int colIn = 0, int iColIn = 0, int iAcolIn = 0, 
    int colReconnectionIn = 0, bool isJunIn = false, bool isAntiJunIn = false,
    bool isActiveIn = true, bool isRealIn = false) : col(colIn), iCol(iColIn), 
    iAcol(iAcolIn), colReconnection(colReconnectionIn), isJun(isJunIn),
    isAntiJun(isAntiJunIn),isActive(isActiveIn), isReal(isRealIn)
    {leftDip = 0; rightDip = 0; iColLeg = 0; iAcolLeg = 0; printed = false;}

  double mDip(Event & event) {
    if (isJun || isAntiJun) return 1E9;
    else return m(event[iCol].p(),event[iAcol].p());
  }

  // Members.
  int    col, iCol, iAcol, iColLeg, iAcolLeg, colReconnection;
  bool   isJun, isAntiJun, isActive, isReal, printed, inChain;
  ColourDipole *leftDip, *rightDip;
  vector<ColourDipole *> colDips, acolDips;
  double p1p2;

  // Printing function, mainly intended for debugging.
  void print();

};

//==========================================================================

// Junction class. In addition to the normal junction class, also contains a 
// list of dipoles connected to it.

class ColourJunction : public Junction {

public:

  ColourJunction(const Junction& ju) : Junction(ju) {
      for(int i = 0;i < 3;++i) {
	dips[i] = 0; dipsOrig[i] = 0;}
  }
  ColourJunction(const ColourJunction& ju) : Junction(Junction(ju)) {
    for(int i = 0;i < 3;++i) {
      dips[i] = ju.dips[i]; dipsOrig[i] = ju.dipsOrig[i];}
  }
  ColourJunction& operator=( const ColourJunction& ju) {
    this->Junction::operator=(ju);
    for(int i = 0;i < 3;++i) {
      dips[i] = ju.dips[i]; dipsOrig[i] = ju.dipsOrig[i];
    }
    return (*this);
  }
  
  ColourDipole * getColDip(int i) {return dips[i];}
  void setColDip(int i, ColourDipole * dip) {dips[i] = dip;}
  ColourDipole * dips[3]; 
  ColourDipole * dipsOrig[3];
  void print();

};

//==========================================================================

// ColourParticle class.

class ColourParticle : public Particle {

public:

 ColourParticle(const Particle& ju) : Particle(ju) {}

  vector<vector<ColourDipole *> > dips; 
  vector<bool> colEndIncluded, acolEndIncluded;
  vector<ColourDipole *> activeDips;
  bool isJun;
  int junKind;

  // Printing functions, intended for debugging.
  void  list();
  void listActiveDips();
  void print();

};

//==========================================================================

// The ColourReconnection class handles the colour reconnection.

//--------------------------------------------------------------------------

class ColourReconnection {

public:

  // Constructor 
  ColourReconnection() {}

  // Initialization.
  bool init( Info* infoPtrIn, Settings& settings, Rndm* rndmPtrIn,
    BeamParticle* beamAPtrIn, BeamParticle* beamBPtrIn,
    PartonSystems* partonSystemsPtrIn);

  // Do Colour reconnection for current event.
  bool next( Event & event, int oldSize);

private:
  
  // list of current dipoles.
  vector<ColourDipole*> dipoles;
  vector<ColourJunction> junctions;
  vector<ColourParticle> particles;

  vector<int> iColEnd, iAcolEnd, iColAndAcol;
  vector<vector<int> > iColJun;

  // Variables needed.
  int    nSys, nReconCols, swap1, swap2, reconnectMode, flipMode;
  bool   allowJunctions, sameNeighbourCol;
  double eCM, sCM, pT0, pT20Rec, pT0Ref, ecmRef, ecmPow, reconnectRange, 
         m0, m0sqr, m2Lambda, fracGluon, dLambdaCut;

  // Pointer to various information on the generation.
  Info*          infoPtr;

  // Pointer to the random number generator.
  Rndm*          rndmPtr;

  // Pointers to the two incoming beams.
  BeamParticle*  beamAPtr;
  BeamParticle*  beamBPtr;

  // Pointer to information on subcollision parton locations.
  PartonSystems* partonSystemsPtr;

  // This is only to access the function call junctionRestFrame.
  StringFragmentation stringFragmentation;

  // This class is used to calculate the string length.
  StringLength stringLength;

  // Do colour reconnection for the event using the new model.
  bool nextNew( Event & event, int oldSize);

  // Simple test swap between two dipoles.
  void swapDipoles(ColourDipole* dip1, ColourDipole* dip2, bool back = false);
  
  // Setup the dipoles.
  void setupDipoles( Event& event, int iFirst = 0);

  // Form pseuparticle of a given dipole (or junction system).
  void makePseudoParticle( ColourDipole* dip, int status, 
    bool setupDone = false);

  // Find Length of string neighbouring the given dipole assuming the dipole 
  // was collapsed.
  double neighbourLength(ColourDipole* dipNeigh, vector<ColourDipole*> &dips, 
    int iOld0, int iOld1, int iNew1, int iOld2 = -7, int iOld3 = -7, 
    int iNew2 = -7);
 
  // Find the indices in the particle list of the junction and also their
  // respectively leg numbers.
  bool getJunctionIndicies(ColourDipole* dip, int &iJun, int &i0, int &i1, 
    int &i2, int &junLeg0, int &junLeg1, int &junLeg2);
 
  // Make a test pseudo particle used to calculate string lengths.
  int makeTestParticle( ColourDipole* dip);
  
  // Form all possible pseudoparticles.
  void makeAllPseudoParticles(Event & event, int iFirst = 0);

  // Update all colours in the event.
  void updateEvent( Event& event, int iFirst = 0);

  double calculateStringLength( ColourDipole* dip, 
    vector<ColourDipole*> & dips);

  // Calculate the string length for two event indicies.
  double calculateStringLength( int i, int j);

  // Calculate the length of a single junction 
  // given the 3 entries in the particle list.
  double calculateJunctionLength(int i, int j, int k);

  // Calculate the length of a double junction,
  // given the 4 entries in the particle record.
  // First two are expected to be the quarks and second two the anti quarks.
  double calculateDoubleJunctionLength( int i, int j, int k, int l);

  // Find all the particles connected in the junction.
  // If a single junction, the size of iParticles should be 3.
  // For multiple junction structures, the size will increase.
  bool findJunctionParticles( int iJun, vector<int>& iParticles, 
    vector<bool> &usedJuns, int &nJuns, vector<ColourDipole*> &dips);

  // Do a single trial reconnection, return true if colour was changed.
  bool singleReconnection( int iDip1 = -1, int iDip2 = -1);

  // Do a single trial reconnection to form a junction, 
  // return true if junction is formed.
  bool singleJunction( Event& event, int iDip1, int iDip2);

  // Do a single trial reconnection to form a junction, 
  // return true if junction is formed.
  bool singleJunction( Event& event, int iDip1, int iDip2, int iDip3);

  // Calculate length of dipoles not original included.
  double calculateAdditionalLengths(vector<ColourDipole*> oldDips, 
    vector<ColourDipole*> newDips);

  // Print the chain containing the dipole.
  void listChain(ColourDipole* dip);

  // Print all the chains.
  void listAllChains();

  // Print dipoles, intended for debuggning purposes.
  void listDipoles( bool onlyActive = false, bool onlyReal = false);

  // Print particles, intended for debugging purposes.
  void listParticles();
  
  // Print junctions, intended for debugging purposes.
  void listJunctions();

  // Check that the current dipole setup is consistent.
  void checkDipoles();

  // Calculate the invariant mass of a dipole.
  double mDip(ColourDipole* dip);

  // The old MPI-based scheme.
  bool reconnectMPIs( Event& event, int oldSize);

  // Vectors and methods needed for the new gluon-move model.

  // Array of (indices of) all final coloured particles. 
  vector<int> iReduceCol, iExpandCol;

  // Array of all lambda distances between coloured partons.
  int nColMove;
  vector<double> lambdaijMove;

  // Function to return lambda value from array.
  double lambda12Move( int i, int j) {
    int iAC = iReduceCol[i]; int jAC = iReduceCol[j]; 
    return lambdaijMove[nColMove * min( iAC, jAC) + max( iAC, jAC)];
  }    

  // Function to return lambda(i,j) + lambda(i,k) - lambda(j,k).
  double lambda123Move( int i, int j, int k) {
    int iAC = iReduceCol[i]; int jAC = iReduceCol[j]; int kAC = iReduceCol[k];
    return lambdaijMove[nColMove * min( iAC, jAC) + max( iAC, jAC)]
         + lambdaijMove[nColMove * min( iAC, kAC) + max( iAC, kAC)]
         - lambdaijMove[nColMove * min( jAC, kAC) + max( jAC, kAC)];
  }     

  // The new gluon-move scheme.
  bool reconnectMove( Event& event, int oldSize);

};

//==========================================================================

} // end namespace Pythia8

#endif // Pythia8_ColourReconnection_H
