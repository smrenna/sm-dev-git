// ColosurReconnection.cc is a part of the PYTHIA event generator.
// Copyright (C) 2014 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Function definitions (not found in the header) for the
// ColourReconnection class.

#include "Pythia8/ColourReconnection.h"

namespace Pythia8 {

//==========================================================================

// The BeamDipole class is purely internal to reconnectMPIs.

class BeamDipole {

public:

  // Constructor.
  BeamDipole( int colIn = 0, int iColIn = 0, int iAcolIn = 0)
    : col(colIn), iCol(iColIn), iAcol(iAcolIn) {}

  // Members.
  int    col, iCol, iAcol;
  double p1p2;
 
};

//==========================================================================

// The ColourDipole class.

//--------------------------------------------------------------------------

// Printing function, inteded for debugging.

void ColourDipole::print() {
 
  cout << setw(10) << this << setw(6) << col << setw(3) << colReconnection 
       << setw(6) << iCol << setw(5) << iAcol << setw(6) << iColLeg << setw(5) 
       << iAcolLeg << setw(6) << isJun << setw(5) << isAntiJun  << setw(10) 
       << p1p2 << " colDips: ";
  for (int i = 0;i < int(colDips.size());++i)
    cout << setw(10) << colDips[i];
  cout  <<  " acolDips: ";
  for (int i = 0;i < int(acolDips.size());++i)
    cout << setw(10) << acolDips[i];
  cout << setw(3) << isActive << endl;

}

//==========================================================================

// The InfoGluonMove class is purely internal to reconnectMove.

class InfoGluonMove{

public:

  // Constructors.
  InfoGluonMove(int i1in, int col1in, int acol1in, int iCol1in, int iAcol1in,
    int col2in, int iCol2in, int iAcol2in, double lambdaRefIn, 
    double dLambdaIn) : i1(i1in), i2(0), col1(col1in), acol1(acol1in),
    iCol1(iCol1in), iAcol1(iAcol1in), col2(col2in), iCol2(iCol2in),
    iAcol2(iAcol2in), lambdaRef(lambdaRefIn), dLambda(dLambdaIn) {}
  InfoGluonMove(int i1in, int i2in, int iCol1in, int iAcol1in, int iCol2in,
    int iAcol2in, int dLambdaIn) : i1(i1in), i2(i2in), col1(0), acol1(0),
    iCol1(iCol1in), iAcol1(iAcol1in), col2(0), iCol2(iCol2in), 
    iAcol2(iAcol2in), lambdaRef(0.), dLambda(dLambdaIn) {}

  // Members.
  int i1, i2, col1, acol1, iCol1, iAcol1, col2, iCol2, iAcol2;
  double lambdaRef, dLambda;

}; 

//==========================================================================

// The ColourJunction class.

//--------------------------------------------------------------------------

// Printing function, inteded for debugging.

void ColourJunction::print() {

  cout << setw(6) << kind() << setw(6)
       << col(0) << setw(6) << col(1) << setw(6) << col(2) << setw(6)
       << endCol(0) << setw(6) << endCol(1) << setw(6) << endCol(2) << setw(6)
       << status(0) << setw(6) << status(1) << setw(6) << status(2) << setw(10)
       << dips[0] << setw(10) << dips[1] << setw(10) << dips[2] << setw(10) 
       << "\n";
  cout << "     " << setw(10) << dipsOrig[0] << setw(10) << dipsOrig[1] 
       << setw(10) << dipsOrig[2] << endl;

}

//==========================================================================

// The ColourParticle class.

//--------------------------------------------------------------------------

// Printing function, inteded for debugging.

void ColourParticle::list() {

  const Particle& pt = (*this);
  
  // Basic line for a particle, always printed.
  cout << setw(10) << pt.id() << "   " << left
       << setw(18) << pt.nameWithStatus(18) << right << setw(4)
       << pt.status() << setw(6) << pt.mother1() << setw(6)
       << pt.mother2() << setw(6) << pt.daughter1() << setw(6)
       << pt.daughter2() << setw(6) << pt.col() << setw(6) << pt.acol()
       << setprecision(3)
       << setw(11) << pt.px() << setw(11) << pt.py() << setw(11)
       << pt.pz() << setw(11) << pt.e() << setw(11) << pt.m() << "\n";

}

//--------------------------------------------------------------------------

// Printing function, inteded for debugging.

void ColourParticle::listActiveDips() {

  cout << "active dips: " << endl;
  for(int i = 0;i < int(activeDips.size()); ++i)
    activeDips[i]->print();

}

//--------------------------------------------------------------------------

// Printing function, inteded for debugging.

void ColourParticle::print() {

  cout << "---   Particle   ---" << endl;
  for(int i = 0;i < int(dips.size());++i) {
    cout << "(" <<colEndIncluded[i] << ") ";
    for(int j = 0;j < int(dips[i].size());++j) {
      cout << dips[i][j]->iCol << " (" << dips[i][j]->col << ") ";
      if( j == int(dips[i].size() - 1))
	cout << dips[i][j]->iAcol << " (" << acolEndIncluded[i] << ")" << endl;
    }
  }

}

//==========================================================================

// The ColourReconnection class.

//--------------------------------------------------------------------------

// Initialization.

bool ColourReconnection::init( Info* infoPtrIn, Settings& settings, 
  Rndm* rndmPtrIn,  BeamParticle* beamAPtrIn, BeamParticle* beamBPtrIn,
  PartonSystems* partonSystemsPtrIn) {

  // Save pointers.
  infoPtr             = infoPtrIn;
  rndmPtr             = rndmPtrIn;
  beamAPtr            = beamAPtrIn;
  beamBPtr            = beamBPtrIn;
  partonSystemsPtr    = partonSystemsPtrIn; 

  // Total and squared CM energy at nominal energy.
  eCM                 = infoPtr->eCM();
  sCM                 = eCM * eCM;

  // Choice of reconnection model.
  reconnectMode       = settings.mode("ColourReconnection:mode");

  // pT0 scale of MPI; used in the MPI-based reconnection model.
  pT0Ref              = settings.parm("MultipartonInteractions:pT0Ref");
  ecmRef              = settings.parm("MultipartonInteractions:ecmRef");
  ecmPow              = settings.parm("MultipartonInteractions:ecmPow");
  pT0                 = pT0Ref * pow(eCM / ecmRef, ecmPow);

  // Parameter of the MPI-based reconnection model.
  reconnectRange      = settings.parm("ColourReconnection:range");
  pT20Rec             = pow2(reconnectRange * pT0);
 
  // Parameters of the new reconnection model.
  m0                  = settings.parm("ColourReconnection:m0");
  m0sqr               = pow2(m0);
  allowJunctions      = settings.flag("ColourReconnection:allowJunctions");
  nReconCols          = settings.mode("ColourReconnection:nColours");
  sameNeighbourCol = settings.flag("ColourReconnection:sameNeighbourColours");

  // Parameters of gluon-move model.
  m2Lambda            = settings.parm("ColourReconnection:m2Lambda");
  fracGluon           = settings.parm("ColourReconnection:fracGluon");
  dLambdaCut          = settings.parm("ColourReconnection:dLambdaCut");
  flipMode            = settings.mode("ColourReconnection:flipMode");
  
  // Initialize StringLength class.
  stringLength.init(infoPtr, settings);
  
  // Done.
  return true;

}

//--------------------------------------------------------------------------

// Do colour reconnection for current event.

bool ColourReconnection::next( Event& event, int iFirst) {

  // MPI-based reconnection model.
  if (reconnectMode == 0) return reconnectMPIs(event, iFirst);

  // New reconnection model.
  else if (reconnectMode == 1) return nextNew(event, iFirst);

  // Gluon-move model.
  else if (reconnectMode == 2) return reconnectMove(event, iFirst);

  // Undefined.
  else {
    infoPtr->errorMsg("Warning in ColourReconnection::next: "
		      "Colour reconnecion mode not found");
    return true;
  }

}

//--------------------------------------------------------------------------

// Do new colour reconnection for current event.

bool ColourReconnection::nextNew( Event& event, int iFirst) {

  // Clear old records.
  while (!dipoles.empty()) {
    delete dipoles.back(); 
    dipoles.pop_back();
  }
  particles.clear();
  junctions.clear();

  // Setup dipoles and make pseudo particles.
  setupDipoles(event, iFirst);
  makeAllPseudoParticles(event, iFirst);
  
  // Check that everything is in order.
  checkDipoles();
  
  // Start big loop.
  int nTotalJun = 0, nTotalRecon = 0;
  for (int iOuterLoop = 0; iOuterLoop < 20; ++iOuterLoop) {
    bool finished = true;

    // Do inner loop for string reconnections.
    for (int iTry = 0;iTry < 20; ++iTry) {
      int  nRecons = 0;

      // Split dipoles into the 9 different "colours".
      vector<vector<int> > iDips;
      iDips.resize(nReconCols);
      for(int i = 0;i < int(iDips.size()); ++i)
	iDips[i] = vector<int>();
    
      for(int i = 0;i < int(dipoles.size()); ++i)
	if(dipoles[i]->isActive)
	  iDips[dipoles[i]->colReconnection].push_back(i);
      
      // Loop over each colour individually.
      for (int i = 0;i < int(iDips.size()); ++i)
	for (int j = 0; j < int(iDips[i].size()); ++j)
	  for (int k = j + 1; k < int(iDips[i].size()); ++k) 
	    if (singleReconnection(iDips[i][j],iDips[i][k])) 
	      nRecons++;  
      
      checkDipoles();
      // Store total number of reconnections.
      nTotalRecon += nRecons;
      
      // Keep going untill no reconnection happened.
      if (nRecons == 0) 
	break;
    }

    // Loop over list of dipoles to try and form junction structures.
    if (allowJunctions)
    for (int iTry = 0; iTry < 20; ++iTry) {
      int nJuncs = 0;
      // Split dipoles into three categories.
      vector<vector<int> > iDips;
      iDips.resize(3);
      for(int i = 0;i < int(iDips.size()); ++i)
	iDips[i] = vector<int>();
      
      for(int i = 0;i < int(dipoles.size()); ++i)
	if(dipoles[i]->isActive)
	  iDips[dipoles[i]->colReconnection % 3].push_back(i);
      
      // Loop over different "colours" (now only three different groups).
      for (int i = 0;i < int(iDips.size()); ++i)
	for (int j = 0; j < int(iDips[i].size()); ++j)
	  for (int k = j + 1; k < int(iDips[i].size()); ++k)
	    if(singleJunction(event,iDips[i][j],iDips[i][k]))
	      nJuncs++;
      
      // Loop over different "colours" (now only three different groups).
      for (int i = 0;i < int(iDips.size()); ++i)
	for (int j = 0; j < int(iDips[i].size()); ++j)
	  for (int k = j + 1; k < int(iDips[i].size()); ++k)
	    for (int l = k + 1; l < int(iDips[i].size()); ++l)
	      if(singleJunction(event,iDips[i][j],iDips[i][k], iDips[i][l]))
		nJuncs++;
      
      // Store total number of junctions.
      nTotalJun += nJuncs;
      
      // Break if no new junction structures were found.
      if (nJuncs == 0) 
	break;
      else
	finished = false;
    }
    // If no junctions were made, the overall loop is finished.
    if (finished)
      break;
  }
  
  updateEvent(event, iFirst);

  // Done.
  return true;  
}


//--------------------------------------------------------------------------

// Setup initial guess on dipoles, here all colours are assumed 
// to be found in the final state.

void ColourReconnection::setupDipoles( Event& event, int iFirst) {

  // Make vectors needed for storage of chains.
  vector< vector<int > > chains;
  vector<bool> isJun;
  vector<bool> isAntiJun;
  vector<bool> isGluonLoop;
  vector<bool> inChain(event.size(),false);
 
  // Find all quarks and follow untill no more colour.
  for (int i = iFirst; i < event.size(); ++i) {
    if (event[i].isFinal() && !inChain[i] && 
        event[i].isQuark() && event[i].id() > 0) {
      int curCol = event[i].col();
      inChain[i] = true;
      vector<int> chain;
      chain.push_back(i);
      isAntiJun.push_back(false);
      isJun.push_back(false);
      isGluonLoop.push_back(false);
      for (int iSteps = 0; curCol != 0 && iSteps < 1000; ++iSteps) {

	// Check for particles with correct anti colour.
	for (int j = iFirst; j < event.size(); j++) {
	  if (event[j].isFinal() && !inChain[j] && event[j].acol() == curCol) {
	    chain.push_back(j);
	    inChain[j] = true;
	    curCol = event[j].col();
	    break;
	  }
	}
	
	// Check for junction with correct colour.
	for (int iJun = 0; iJun < event.sizeJunction(); ++iJun) {
	  for (int j = 0; j < 3; ++j) {
	    if (event.colJunction(iJun,j) == curCol) {
	      isJun[isJun.size() -1] = true;
	      curCol = 0;
	      chain.push_back( -(10 + 10 * iJun + j) );
	    }
	  }
	}
      }   
      chains.push_back(chain);
    }
  }
  
  // Start from anti-junction and make chains.
  for (int i = 0; i < event.sizeJunction(); ++i) {
  
    // First check if junction belongs to the right diffractive system.
    int checkCol = event.colJunction(i,0);
    bool wrongSystem = false;
    for (int j = 0; j < iFirst; ++j)
      if (event[j].isFinal() && event[j].acol() == checkCol)
	wrongSystem = true;
    if (wrongSystem)
      continue;
  
    // Loop over legs of anti junction.
    if (event.kindJunction(i) == 2) 
    for (int jCol = 0; jCol < 3; ++jCol) {
      int curCol = event.colJunction(i,jCol);
      vector<int> chain;
      chain.push_back( -(10 + 10 * i + jCol));
      isAntiJun.push_back(true);
      isJun.push_back(false);
      isGluonLoop.push_back(false);
      for (int iSteps = 0; curCol != 0 && iSteps < 1000; ++iSteps) {

	// Check for particles with correct anti colour.
	for (int j = iFirst; j < event.size(); ++j) 
        if (event[j].isFinal() && !inChain[j] && 
	    event[j].acol() == curCol) {
	  chain.push_back(j);
	  inChain[j] = true;
	  curCol = event[j].col();
	  break;
	}
	  
	// Check for junction with correct colour.
	for (int iJun = 0; iJun < event.sizeJunction(); ++iJun) 
	if (event.kindJunction(iJun) == 1) 
	for (int j = 0; j < 3; ++j) 
	if (event.colJunction(iJun,j) == curCol) {
	  isJun[isJun.size() - 1] = true;
	  curCol = 0;
	  chain.push_back( -(10 + 10 * iJun + j));
	}
      } 
      chains.push_back(chain);
    }
  }
 
  // Find all gluon loops.
  for (int i = iFirst; i < event.size(); ++i) 
  if (event[i].isFinal() && !inChain[i] && event[i].col() != 0) {
    int curCol = event[i].col();
    inChain[i] = true;
    vector<int> chain;
    chain.push_back(i);
    isAntiJun.push_back(false);
    isJun.push_back(false);
    isGluonLoop.push_back(true);
    for (int iSteps = 0; curCol != 0 && iSteps < 1000; ++iSteps) {
      bool foundNext = false;
      for (int j = iFirst; j < event.size(); ++j) 
      if (event[j].isFinal() && !inChain[j] && event[j].acol() == curCol) {
	chain.push_back(j);
	inChain[j] = true;
	curCol = event[j].col();
	foundNext = true;
	break;
      }
      	
      if(!foundNext)
	break;
    }
    chains.push_back(chain);
  }
 
  // Form dipoles from chains.
  for (int i = 0; i < int(chains.size()); ++i) {
    if (chains[i].size() == 1) continue;
    int lastCol = -1;
    int firstCol = 0;

    // Start from the first and form the dipoles.
    // Two consececutive dipoles can not share the same colour.
    for (int j = 0; j < int(chains[i].size()); ++j) {  
      if (j != int(chains[i].size() - 1)) {

	// Start by picking new colour.
        int newCol = lastCol;
        while (newCol == lastCol && !sameNeighbourCol) {
          newCol = (int)(rndmPtr->flat() * nReconCols); 
        }

	// Need to check whether the quark comes from a g->qqbar split.
	// If that is the case, it can not have the same as qbar.
        if (j == 0 && !isAntiJun[i] && !isGluonLoop[i]) {

	  int iMother = event[event[ chains[i][j] ].iTopCopy()].mother1();
	  if ( event[iMother].idAbs() == 21) {
            vector<int> sisters = event[ chains[i][j] ].sisterList(true);
	    // Need to have only one sister and need to be the anti particle.
	    if (sisters.size() == 1 && event[ chains[i][j] ].id() 
		== - event[ sisters[0] ].id()) {

	      // Try to find dipole with sister.
	      int colSis = -1;
	      for (int k = 0; k < int(dipoles.size()); ++k) 
		if (dipoles[k]->iAcol == sisters[0]) {
		  colSis = dipoles[k]->colReconnection;
		  break;
		}
	    
	      // If the two colours are the same, pick a new.
	      while (colSis == newCol && !sameNeighbourCol)
		newCol = (int)(rndmPtr->flat() * nReconCols);
	    }
	  }
	}
      
        // Check if quark end comes from g->qqbar split.
	// If so check that the two quarks get different colours.
	if (j == int(chains[i].size() - 2) && !isJun[i] && !isGluonLoop[i]) {

	  int iMother = event[event[chains[i][j + 1]].iTopCopy()].mother1();
	  if (event[iMother].idAbs() == 21) {
            vector<int> sisters = event[ chains[i][j + 1] ].sisterList(true);
	    // Need to have only one sister and need to be the anti particle.
	    if (sisters.size() == 1 && event[ chains[i][j + 1] ].id() 
		== - event[ sisters[0] ].id()) {

	      // Try to find dipole with sister.
	      int colSis = -1;
	      for (int k = 0; k < int(dipoles.size()); ++k) 
		if (dipoles[k]->iCol == sisters[0]) {
		  colSis = dipoles[k]->colReconnection;
		  break;
		}
	    
	      // If the two colours are the same, pick a new.
	      while ((colSis == newCol || newCol == lastCol) 
		     && !sameNeighbourCol)
		newCol = (int)(rndmPtr->flat() * nReconCols);
	    }
	  }
	}

	// Special case for junction splitting if multiple gluons 
	// between the junctions.
	if ((chains[i][j] > 0 && event[chains[i][j]].status() == 65) || 
	    (chains[i][j + 1] > 0 && 
             event[ chains[i][j + 1] ].status() == 65) ) {

 	  // Find sisters.
	  vector<int> sisters;
	  if (chains[i][j] > 0 && event[ chains[i][j] ].status() == 65)
	    sisters = event[ chains[i][j] ].sisterList();
	  else
	    sisters = event[ chains[i][j + 1] ].sisterList();

          if (sisters.size() == 3 ) {
	    
	    // Find colour of sisters.
	    int acolSis1 = -1, acolSis2 = -1, acolSis3 = -1;
	    int colSis1 = -1, colSis2 = -1, colSis3 = -1;
	    for (int k = 0;k < int(dipoles.size()); ++k)  {
	      if( dipoles[k]->iAcol == sisters[0]) 
		acolSis1 = dipoles[k]->colReconnection;

	      if( dipoles[k]->iAcol == sisters[1]) 
		acolSis2 = dipoles[k]->colReconnection;
	    
	      if( dipoles[k]->iAcol == sisters[2]) 
		acolSis3 = dipoles[k]->colReconnection;

	      if( dipoles[k]->iCol == sisters[0]) 
		colSis1 = dipoles[k]->colReconnection;

	      if( dipoles[k]->iCol == sisters[1]) 
		colSis2 = dipoles[k]->colReconnection;
	    
	      if( dipoles[k]->iCol == sisters[2]) 
		colSis3 = dipoles[k]->colReconnection;
	    }

	    // If the two colours are the same, pick a new.
	    while ((colSis1 == newCol || colSis2 == newCol || 
                   colSis3 == newCol || acolSis1 == newCol || 
                   acolSis2 == newCol || acolSis3 == newCol) 
		   && !sameNeighbourCol)
	      newCol = (int)(rndmPtr->flat() * nReconCols);
	  }
	}

	// Update stored colours.
	if (j == 0) firstCol = newCol;
	lastCol = newCol;

	// Check if it is anti junction need special dipole.
	if (j == 0 && isAntiJun[i]) {
	  int col = event.colJunction( - int(chains[i][j]/10) - 1, 
				       -chains[i][j] % 10);
	  dipoles.push_back(new ColourDipole(col, chains[i][j], 
	    chains[i][j+1], newCol));
	  dipoles.back()->isAntiJun = true;
	}

	// Otherwise just make the dipole.
	else dipoles.push_back(new ColourDipole(event[ chains[i][j] ].col(), 
          chains[i][j], chains[i][j+1], newCol));

	// If the chain in end a junction mark it.
	if (j == int(chains[i].size() - 2) && isJun[i]) 
          dipoles.back()->isJun = true;

	// Update the links between dipoles.
	if (j > 0) {
	  dipoles[dipoles.size() - 1]->leftDip  = dipoles[dipoles.size() - 2];
	  dipoles[dipoles.size() - 2]->rightDip = dipoles[dipoles.size() - 1];
	}    
      } 

      // If last particle has anti colour it should be possible to connect it 
      // to the first particle in the chain. (e.g. gluon loop)
      else 
      if (isGluonLoop[i]) 
      if (event[ chains[i][j] ].col() == event[ chains[i][0] ].acol()) {
	int newCol = lastCol;
	while ( (newCol == lastCol || newCol == firstCol)
		&& !sameNeighbourCol) {
	  newCol = (int)(rndmPtr->flat() * nReconCols); 
	}
	dipoles.push_back(new ColourDipole(event[ chains[i][j] ].col(), 
	  chains[i][j], chains[i][0], newCol));
	
	// Update links between dipoles.
	dipoles[dipoles.size() - 1]->leftDip = dipoles[dipoles.size() - 2];
	dipoles[dipoles.size() - 2]->rightDip = dipoles[dipoles.size() - 1];
	dipoles[dipoles.size() - chains[i].size()]->leftDip = 
	  dipoles[dipoles.size() -1];
	dipoles[dipoles.size() - 1]->rightDip = 
	  dipoles[dipoles.size() - chains[i].size()];
	  
      }
    }
  }

  // Setup junction list.
  iColJun.clear();
  iColJun.resize(event.sizeJunction());
  for (int i = 0; i < int(iColJun.size()); ++i) iColJun[i] = vector<int>(3,0);
  
  // Loop over event and store indecies.
  for (int i = iFirst; i < event.size(); ++i)
  if (event[i].isFinal()) 
  for (int j = 0; j < event.sizeJunction(); ++j)
  for (int jLeg = 0; jLeg < 3; ++jLeg)
  if (event[i].col() == event.colJunction(j,jLeg) ||
      event[i].acol() == event.colJunction(j,jLeg))
    iColJun[j][jLeg] = i;
  
  // Loop over junction and store indecies.
  for (int i = 0;i < event.sizeJunction(); ++i)
  for (int iLeg = 0; iLeg < 3; ++iLeg)
  for (int j = i + 1;j < event.sizeJunction(); ++j)
  for (int jLeg = 0; jLeg < 3; ++jLeg)
  if (event.colJunction(i, iLeg) == event.colJunction(j, jLeg)) {
    iColJun[i][iLeg] = -(10*j + 10 + jLeg);
    iColJun[j][jLeg] = -(10*i + 10 + iLeg);
  }

  // Done.
}

//--------------------------------------------------------------------------

// Calculate the string length of a dipole.

double ColourReconnection::calculateStringLength(ColourDipole * dip, 
  vector<ColourDipole *> &dips) {

  // Check if the dipole has already been included.
  for (int i = 0;i < int(dips.size()); ++i)
    if (dips[i] == dip)
      return 0;

  // Handle simple strigs with no junction connections.
  if (!dip->isJun && !dip->isAntiJun) {
    // If dipole is larger than m0
    if (mDip(dip) >= m0) {
      dips.push_back(dip);
      double dipLength =  calculateStringLength(dip->iCol, dip->iAcol);
      return dipLength;
    }
    
    // Otherwise form a pseudo particle and calculate length of neighbours.
    else {
      // MIGHT NEED TO CONSIDER CASE, IF IT ACTUALLY PREFERS 
      // TO GO TOGETHER WITH NEIGHBOUR.
      dips.push_back(dip);
      int iNew = makeTestParticle(dip);
      double dipLength = 0;
      // First handle one side.
      for(int i = 0; i < int(particles[ dip->iCol ].activeDips.size()); ++i) 
        dipLength += neighbourLength(particles[ dip->iCol ].activeDips[i], dips,
         dip->iCol, dip->iAcol, iNew);
      // Then the other side.
      for(int i = 0; i < int(particles[dip->iAcol].activeDips.size()); ++i) 
         dipLength += neighbourLength(particles[dip->iAcol].activeDips[i], dips,
           dip->iCol, dip->iAcol, iNew);			      

      // Remove the test particle. 
      particles.pop_back();
      
      // Done.
      return dipLength;
    }
  } 

  // Handle Junctions.
  else {  
 
    // If the dipole has mass above m0.
    if (mDip(dip) >= m0) {
      // Start by finding all particles connected to the junction system.
      vector<int> iParticles;
      vector<bool> usedJuns(junctions.size(),false);
      int nJuns = 0;
      if(dip->isJun) {
	if (!findJunctionParticles( -int(dip->iAcol/10) - 1, iParticles, 
          usedJuns, nJuns, dips)) return 1E9;
      } else
	if (!findJunctionParticles(-int(dip->iCol/10) - 1, iParticles, 
          usedJuns, nJuns, dips)) return 1E9;
      
      // If it is a single junction.
      if(int(iParticles.size()) == 3) 
	return calculateJunctionLength(iParticles[0], iParticles[1],
          iParticles[2]);
      
      // If it is a junction pair.
      else if(int(iParticles.size()) == 4) 
        return calculateDoubleJunctionLength(iParticles[0], 
          iParticles[1], iParticles[2], iParticles[3]);

      // If any other number of junction legs return high number.
      else return 1E9;
    }
      
    // Otherwise form a pseudo particle and calculate length of neighbours.
    else {
   
      // Make test particle.
      int iNew = makeTestParticle(dip);
      double dipLength = 0;
      
      // Get indicies for dipole.
      int i0, i1, i2;
      int junLeg0, junLeg1, junLeg2;
      int iJun;
      
      // Should always return positive, mainly sanity check.
      if (!getJunctionIndicies(dip, iJun, i0, i1, i2, junLeg0, 
                               junLeg1, junLeg2)) {
	particles.pop_back();
	return 1E9;
      }
      
      // Add now passive dipoles to list of dipoles.
      dips.push_back(dip);
      dips.push_back(junctions[iJun].dips[junLeg1]);
      
      // First handle one leg.
      for (int i = 0;i < int(particles[i0].activeDips.size()); ++i) 
        dipLength += neighbourLength(particles[i0].activeDips[i], dips, 
          i0, i1, iNew);
      
      // Second leg.
      for (int i = 0;i < int(particles[i1].activeDips.size()); ++i) 
        dipLength += neighbourLength(particles[i1].activeDips[i], dips,
          i0, i1, iNew);
      
      // Third leg is special.
      // Need to first tell it to not see it is as a junction.
      if (dip->isJun) {
	junctions[iJun].dips[junLeg2]->isJun = false;
	junctions[iJun].dips[junLeg2]->iAcol = iNew;
      } else {
	junctions[iJun].dips[junLeg2]->isAntiJun = false;
	junctions[iJun].dips[junLeg2]->iCol = iNew;
      }

      // Calculate the length.
      dipLength += neighbourLength(junctions[iJun].dips[junLeg2], dips,
        i0, i1, iNew);
      
      // Restore correct setup.
      if (dip->isJun) {
	junctions[iJun].dips[junLeg2]->isJun = true;
	junctions[iJun].dips[junLeg2]->iAcol = -( (iJun+1)*10 + junLeg2);
      } else {
	junctions[iJun].dips[junLeg2]->isAntiJun = true;
	junctions[iJun].dips[junLeg2]->iCol = -( (iJun+1)*10 + junLeg2);
      }
      particles.pop_back();
      
      return dipLength;
    }
  }
  
}
//--------------------------------------------------------------------------
 
// Update all colours in the event.

void ColourReconnection::updateEvent( Event& event, int iFirst) {

  // Start by making a new copy of particles.
  int oldSize = event.size();
  for(int i = iFirst; i < oldSize;++i)
    if(event[i].isFinal()) event.copy(i,66);
  
  // Copy over junctions.
  event.clearJunctions();
  for(int i = 0;i < int(junctions.size()); ++i) {
    for(int j = 0;j < 3; ++j) {
      if ( junctions[i].dipsOrig[j] != 0) {
	junctions[i].col(j, junctions[i].dipsOrig[j]->col); 
      }
    }
    event.appendJunction(Junction(junctions[i]));
  }
  
  // Assign colour according to the real dipoles.
  for(int i = 0;i < int(dipoles.size()); ++i)
    if(dipoles[i]->isReal) {
      if(dipoles[i]->iCol >= 0)
	event[ event[ dipoles[i]->iCol ].daughter1() ].col(dipoles[i]->col);
      else
	event.colJunction(-(dipoles[i]->iCol/10 + 1), -dipoles[i]->iCol % 10,
          dipoles[i]->col);
      if(dipoles[i]->iAcol >= 0)
	event[ event[ dipoles[i]->iAcol ].daughter1() ].acol(dipoles[i]->col);
      else
	event.colJunction(-(dipoles[i]->iAcol/10 + 1), -dipoles[i]->iAcol % 10,
          dipoles[i]->col); 
    }
}

//--------------------------------------------------------------------------

// Find all the particles connected in the junction.
// If a single junction, the size of iParticles should be 3.
// For multiple junction structures, the size will increase.

bool ColourReconnection::findJunctionParticles(int iJun, 
  vector<int>& iParticles, vector<bool> &usedJuns, int &nJuns, 
  vector<ColourDipole*> &dips ) {
  
  // Mark current junction as used.
  usedJuns[iJun] = true;
  nJuns++;

  // It is not possible to handle junction structures larger than two.
  if (nJuns > 2) 
    return false;
  
  // Find particles connected to the 
  if (junctions[iJun].kind() % 2 == 1) 
    for (int i = 0; i < 3; ++i)
      iParticles.push_back(junctions[iJun].dips[i]->iCol);
  else
    for (int i = 0; i < 3; ++i)
      iParticles.push_back(junctions[iJun].dips[i]->iAcol);
  
  // Add dipoles if not already included.
  for (int i = 0; i < 3; ++i) {
    bool addDip = true;
    for (int j = 0; j < int(dips.size()); ++j) {
      if (dips[j] == junctions[iJun].dips[i]) {
	addDip = false;
	break;
      }
    }
    if (addDip) dips.push_back(junctions[iJun].dips[i]);
  }

  // Check whether it connects to any other junctions.
  for (int i = 0; i < int(iParticles.size()); ++i)
  if (iParticles[i] < 0) {
    int iNewJun = - int(iParticles[i] / 10) -1;
    iParticles.erase(iParticles.begin() + i);
    i--;
    if (!usedJuns[iNewJun] && !findJunctionParticles( iNewJun, iParticles, 
      usedJuns, nJuns, dips) )
      return false;
  }
  
  // Done.
  return true;
}

//--------------------------------------------------------------------------

// Calculate string length for two indices in the particles record.

double ColourReconnection::calculateStringLength(int i, int j) {
  
  // Find rest frame of particles.
  Vec4 p1 =  particles[i].p();
  Vec4 p2 =  particles[j].p();
  Vec4 pSum = p1 + p2;
  
  // Check that particles are not completely parallel.
  // This should not happen.
  if (REtaPhi(p1,p2)  < 1E-5) return 1E9;

  // Boost to rest frame.
  p1.bstback(pSum);
  p2.bstback(pSum);
  
  // Calculate string length.
  Vec4 p0(0,0,0,1.);
  return stringLength.getLength(p1, p0) + stringLength.getLength(p2, p0);

}

//--------------------------------------------------------------------------

// Calculate the length of a single junction given the 3 entries in the event.

double ColourReconnection::calculateJunctionLength(int i, 
  int j, int k) {
  
  // Need to be separate indices.
  if ( i == j || i == k || j == k) return 1E9;
  
  Vec4 p1 = particles[i].p();
  Vec4 p2 = particles[j].p();
  Vec4 p3 = particles[k].p();
  
  return stringLength.getJuncLength(p1, p2, p3);

}

//--------------------------------------------------------------------------

// Calculate the length of a double junction given the 4 particle entries.
// The first two are expected to be quarks, the second two to be antiquarks.

double ColourReconnection::calculateDoubleJunctionLength( int i, int j, 
  int k, int l) {
  
  // Need to be separate indices.
  if (i == j || i == k || i == l || j == k || j == l || k == l) return 1E9;
  
  // Calculate minimum length of new colour structure.
  double origLength = calculateStringLength( i, k) + 
    calculateStringLength( j, l);
  double minLength  = calculateStringLength( i, j) + 
    calculateStringLength( k, l); 

  // If minimum length is larger than original length, reject already here.
  if (origLength < minLength) return minLength;
  
  Vec4 p1 = particles[i].p();
  Vec4 p2 = particles[j].p();
  Vec4 p3 = particles[k].p();
  Vec4 p4 = particles[l].p();
  
  return stringLength.getJuncLength(p1, p2, p3, p4);

}

//--------------------------------------------------------------------------

// Do a single trial emission.

bool ColourReconnection::singleReconnection( int iDip1, int iDip2) {
  
  // Select trial dipoles.
  if(iDip1 < 0) iDip1 = (int)(rndmPtr->flat() * dipoles.size());
  if(iDip2 < 0) iDip2 = (int)(rndmPtr->flat() * dipoles.size());

  ColourDipole* dip1 = dipoles[iDip1];
  ColourDipole* dip2 = dipoles[iDip2];
  
  // Do nothing if it is the same dipole.
  if(iDip1 == iDip2)
    return false;

  // No colour reconnection if colours do not match.
  if (dip1->colReconnection != dip2->colReconnection)
    return false;

  // If it is not active return
  if(!dip1->isActive) return false;
  if(!dip2->isActive) return false;

  // Not possible to connect a gluon with itself.
  if(dip1->iCol == dip2->iAcol)
    return false;
  if(dip1->iAcol == dip2->iCol)
    return false;

  vector<ColourDipole*> oldDips, newDips;
  
  // Calculate old string length.
  double oldStringLength = calculateStringLength(dip1, oldDips)
    + calculateStringLength( dip2, oldDips);
 
  // Make test configuration.
  swapDipoles(dip1,dip2);
  
  // IF BOTH DIPS HAVE MASS BELOW M0 NOT TREATED CORRECTLY IF THEY ARE 
  // CONNECTED. FIRST DIPOLE IS MADE INTO A SINGLE PARTICLE, THEN REVERTED 
  // AND FOLLOWED BY A SECOND. PROBABLY NOT REALLY THAT LIKELY TO BE A 
  // SERIOUS PROBLEM, BUT SHOULD BE NOTED.
  
 // Calculate new string lengths
  double newStringLength = calculateStringLength(dip1, newDips)
    +  calculateStringLength(dip2, newDips);
 
  // Swap back.
  swapDipoles(dip1,dip2, true);
 
  // First check if new combination was not useable.
  if (newStringLength >= 0.5E9) return false; 
  
  // Add dipoles to old length that was included in the new length calculation.
  oldStringLength += calculateAdditionalLengths(oldDips, newDips);

  // Check if it is a favourable setup, if not return.
  if (oldStringLength <= newStringLength) return false;
 
  // Update to new configuration.
  
  // If both acols ends are normal particles.
  if (dip1->iAcol >= 0 && dip2->iAcol >= 0) {
    swap(particles[dip1->iAcol].dips[dip1->iAcolLeg].front()->iAcol,
         particles[dip2->iAcol].dips[dip2->iAcolLeg].front()->iAcol);
    swap(particles[dip1->iAcol].dips[dip1->iAcolLeg].front(),
         particles[dip2->iAcol].dips[dip2->iAcolLeg].front());
  // If only dip1 has normal acol end.
  } else if (dip1->iAcol >= 0) {
    swap(particles[dip1->iAcol].dips[dip1->iAcolLeg].front()->iAcol,
	 junctions[-(dip2->iAcol / 10 + 1)].dipsOrig[-dip2->iAcol % 10]->iAcol);
    swap(particles[dip1->iAcol].dips[dip1->iAcolLeg].front(),
	 junctions[-(dip2->iAcol / 10 + 1)].dipsOrig[-dip2->iAcol % 10]);
  // If only dip2 has normal acol end.
  } else if(dip2->iAcol >= 0) {
    swap(particles[dip2->iAcol].dips[dip2->iAcolLeg].front()->iAcol,
	 junctions[-(dip1->iAcol / 10 + 1)].dipsOrig[-dip1->iAcol % 10]->iAcol);
    swap(particles[dip2->iAcol].dips[dip2->iAcolLeg].front(),
	 junctions[-(dip1->iAcol / 10 + 1)].dipsOrig[-dip1->iAcol % 10]);
  // If both ends are junctions.
  } else {
    swap(junctions[ -(dip1->iAcol / 10 + 1) ].dipsOrig[
           -dip1->iAcol % 10 ]->iAcol,
	 junctions[ -(dip2->iAcol / 10 + 1) ].dipsOrig[ 
           -dip2->iAcol % 10 ]->iAcol); 
    swap(junctions[ -(dip1->iAcol / 10 + 1) ].dipsOrig[ 
           -dip1->iAcol % 10 ],
	 junctions[ -(dip2->iAcol / 10 + 1) ].dipsOrig[ 
           -dip2->iAcol % 10]);
  }

  // Swap the dipoles.
  swapDipoles(dip1, dip2);

  // If new particles are below treshhold, form pseudoParticles.
  if(mDip(dip1) < m0) makePseudoParticle(dip1, 68, true);
  if(mDip(dip2) < m0) makePseudoParticle(dip2, 68, true);
 
  return true;
}

//--------------------------------------------------------------------------

// Simple test swap between two dipoles.

void ColourReconnection::swapDipoles(ColourDipole* dip1, 
  ColourDipole* dip2, bool back) {
 
  // Swap the anti colour of the dipoles.
  swap(dip1->iAcol, dip2->iAcol);
  swap(dip1->isJun, dip2->isJun);
  swap(dip1->iAcolLeg, dip2->iAcolLeg);
  
  // Update the active dipoles.
  // Only change 1 active dipole, this is to avoid problems when switching back.
  if (dip1->iAcol != dip2->iAcol) {
    if (!back) {
      if (dip1->iAcol >= 0) 
      for (int i = 0; i < int(particles[dip1->iAcol].activeDips.size()); ++i)
      if (particles[dip1->iAcol].activeDips[i] == dip2) {
	particles[dip1->iAcol].activeDips[i] = dip1;
	swap1 = i;
	break;
      }
      if (dip2->iAcol >= 0) 
      for (int i = 0; i < int(particles[dip2->iAcol].activeDips.size()); ++i)
      if (particles[dip2->iAcol].activeDips[i] == dip1) {
	particles[dip2->iAcol].activeDips[i] = dip2;
	swap2 = i;
	break;	    
      }
    } else {
      if (dip1->iAcol >= 0) particles[dip1->iAcol].activeDips[swap2] = dip1;
      if (dip2->iAcol >= 0) particles[dip2->iAcol].activeDips[swap1] = dip2;
    }
  }

  // Update list of junctions (only junctions, anti junctions stay the same).
  for (int i = 0; i < int(junctions.size()); ++i) 
  if (junctions[i].kind() % 2 == 1) 
  for (int iLeg = 0; iLeg < 3; ++iLeg) {
    if (junctions[i].dips[iLeg] == dip1) {
      junctions[i].dips[iLeg] = dip2;
      continue;
    }
    if (junctions[i].dips[iLeg] == dip2) {
      junctions[i].dips[iLeg] = dip1;
      continue;
    }
  }

  // Done.
}

//--------------------------------------------------------------------------

// Calculate length of dipoles not original included.

double ColourReconnection::calculateAdditionalLengths(
  vector<ColourDipole*> oldDips, vector<ColourDipole*> newDips) {
  
  // Remove dipoles that was already included in both cases.
  for (int i = 0; i < int(oldDips.size()); ++i) 
  for (int j = 0; j < int(newDips.size()); ++j) 
  if ( oldDips[i] == newDips[j]) {
    oldDips.erase(oldDips.begin() + i);
    newDips.erase(newDips.begin() + j);
    i--;
    break;
  }

  // If any dipoles remain in the old case something went wrong.
  if (oldDips.size() != 0) {
    infoPtr->errorMsg("Warning in ColourReconnection::"
      "calculateAdditionalLengths: Old neighbour dipoles remained");
    return 1E9;
  }
  
  double additionalLength = 0;
  // Add missed lengths to the old length.
  while (newDips.size() > 0) {
    vector<ColourDipole *> dips;
    additionalLength += calculateStringLength(newDips[0], dips);
   
    // Remove the dipoles that are now included.
    for (int i = 0; i < int(dips.size()); ++i) {
      bool found = false;
      for (int j = 0; j < int(newDips.size()); ++j) {
	if ( dips[i] == newDips[j]) {
	newDips.erase(newDips.begin() + j);
	found = true;
	break;
	}
      }
      // If I add one not on the list something went wrong.
      if (!found) {
	infoPtr->errorMsg("Warning in ColourReconnection::"
	  "calculateAdditionalLengths: Could not find new neighbours");
	return 1E9;
      }
    }
  }
  
  return additionalLength;
}

//--------------------------------------------------------------------------

// Do a single trial emission.

bool ColourReconnection::singleJunction( Event& event,int iDip1, int iDip2) {

  // Select trial dipoles.
  if (iDip1 < 0) iDip1 = (int)(rndmPtr->flat() * dipoles.size());
  if (iDip2 < 0) iDip2 = (int)(rndmPtr->flat() * dipoles.size());

  // Do nothing if it is the same dipole.
  if(iDip1 == iDip2)
    return false;
  
  ColourDipole* dip1 = dipoles[iDip1];
  ColourDipole* dip2 = dipoles[iDip2];

  int iCol1  = dipoles[iDip1]->iCol;
  int iCol2  = dipoles[iDip2]->iCol;
  int iAcol1 = dipoles[iDip1]->iAcol;
  int iAcol2 = dipoles[iDip2]->iAcol;
  
  // Not possible to connet a gluon with itself.
  if (iCol1 == iCol2) return false;
  if (iAcol1 == iAcol2) return false;

  // Check that all dipoles are active.
  if (!dip1->isActive || !dip2->isActive) return false;

  // Do nothing if one of the dipoles is a junction or anti junction.
  if (dipoles[iDip1]->isJun || dipoles[iDip1]->isAntiJun) return false;
  if (dipoles[iDip2]->isJun || dipoles[iDip2]->isAntiJun) return false;

  // Only accept one third of the pairs.
  if ( (dipoles[iDip1]->colReconnection) % 3 != 
        dipoles[iDip2]->colReconnection % 3) return false;

  vector<ColourDipole*> oldDips, newDips;

  // Calculate old string length.
  double oldStringLength = calculateStringLength(dip1, oldDips)
    + calculateStringLength( dip2, oldDips);
  
  // Calculate new length and compare to old.
  double newStringLength = 0;

  // If both the new "dipoles" have mass below m0.
  if (m(particles[iCol1].p(), particles[iCol2].p()) < m0 &&
      m(particles[iAcol1].p(), particles[iAcol2].p()) < m0) {

    // Make test particles
    particles.push_back(particles[iCol1]);
    particles.back().p(particles[iCol1].p() + particles[iCol2].p());
    int iNew1 = particles.size() -1;

    particles.push_back(particles[iAcol1]);
    particles.back().p(particles[iAcol1].p() + particles[iAcol2].p());
    int iNew2 = particles.size() -1;
    
    // Calculate length of neighbours.
    newDips.push_back(dip1);
    newDips.push_back(dip2);
    newStringLength += calculateStringLength(iNew1, iNew2);
    for (int i = 0; i < int(particles[iCol1].activeDips.size()); ++i)
      newStringLength += neighbourLength(particles[iCol1].activeDips[i], 
        newDips, iCol1, iCol2, iNew1, iAcol1, iAcol2, iNew2);
    
    for (int i = 0; i < int(particles[iCol2].activeDips.size()); ++i)
      newStringLength += neighbourLength(particles[iCol2].activeDips[i], 
        newDips, iCol1, iCol2, iNew1, iAcol1, iAcol2, iNew2);
    
    for (int i = 0; i < int(particles[iAcol1].activeDips.size()); ++i)
      newStringLength += neighbourLength(particles[iAcol1].activeDips[i], 
        newDips, iCol1, iCol2, iNew1, iAcol1, iAcol2, iNew2);

    for (int i = 0; i < int(particles[iAcol2].activeDips.size()); ++i)
      newStringLength += neighbourLength(particles[iAcol2].activeDips[i], 
        newDips, iCol1, iCol2, iNew1, iAcol1, iAcol2, iNew2);
    
    // Remove test particles.
    particles.pop_back();
    particles.pop_back();
    
  } // If only the col end has mass below m0.
  else if (m(particles[dipoles[iDip1]->iCol].p(), 
	     particles[dipoles[iDip2]->iCol].p()) < m0) {

    // Make test particles
    particles.push_back(particles[iCol1]);
    particles.back().p(particles[iCol1].p() + particles[iCol2].p());
    int iNew1 = particles.size() -1;

    newDips.push_back(dip1);
    newDips.push_back(dip2);
    
    // Calculate string length of the junction left.
    newStringLength += calculateJunctionLength(iNew1, iAcol1, iAcol2);
    
    // Calculate stringlength of the neighbouring particles.
    for (int i = 0;i < int(particles[iCol1].activeDips.size()); ++i)
      newStringLength += neighbourLength(particles[iCol1].activeDips[i], 
        newDips, iCol1, iCol2, iNew1);
    
    for (int i = 0;i < int(particles[iCol2].activeDips.size()); ++i)
      newStringLength += neighbourLength(particles[iCol2].activeDips[i], 
        newDips, iCol1, iCol2, iNew1);

    // Remove test particle.
    particles.pop_back();
    
  // If only the anti colour end has mass below m0.
  } else if (m(particles[dipoles[iDip1]->iAcol].p(), 
	       particles[dipoles[iDip2]->iAcol].p()) < m0) {

    // Make test particles
    particles.push_back(particles[iAcol1]);
    particles.back().p(particles[iAcol1].p() + particles[iAcol2].p());
    int iNew1 = particles.size() -1;

    newDips.push_back(dip1);
    newDips.push_back(dip2);
    
    // Calulate string length of the junction left.
    newStringLength += calculateJunctionLength(iNew1, iCol1, iCol2);
    
    // Calculate string length of the neighbouring particles.
    for (int i = 0;i < int(particles[iAcol1].activeDips.size()); ++i)
      newStringLength += neighbourLength(particles[iAcol1].activeDips[i], 
        newDips, iAcol1, iAcol2, iNew1);
    
    for (int i = 0;i < int(particles[iAcol2].activeDips.size()); ++i)
      newStringLength += neighbourLength(particles[iAcol2].activeDips[i], 
        newDips, iAcol1, iAcol2, iNew1);

    // Remove test particle.
    particles.pop_back();

  // If no "dipoles" have mass below m0. 
  } else {
    newDips.push_back(dip1);
    newDips.push_back(dip2);
    newStringLength = calculateDoubleJunctionLength(iCol1, iCol2, iAcol1, 
      iAcol2);
  }
  
  // Add additional lengths to the old string lengths.
  oldStringLength += calculateAdditionalLengths(oldDips, newDips);

  // Do nothing if old configuration was better.
  if (newStringLength > oldStringLength) return false;

  // Append new junctions.
  // Start by getting colour tags.
  int oldCol1 = dipoles[iDip1]->col;
  int oldCol2 = dipoles[iDip2]->col;
  int newCol1  = event.nextColTag();
  int newCol2 = event.nextColTag();
  int newCol3 = event.nextColTag();
  
  // Need to make 3 new real dipoles and 3 active dipoles.

  // First make dipoles between junctions.

  // Choose the new reconnection colour to be different from the others.
  int colReconnection = int(nReconCols * rndmPtr->flat());
  while (colReconnection == dipoles[iDip1]->colReconnection 
    && colReconnection == dipoles[iDip2]->colReconnection)
    colReconnection = int(nReconCols * rndmPtr->flat());
  
  // Need one active and one real dipole.
  int iJun = junctions.size();
  int iAntiJun = junctions.size() + 1;
  dipoles.push_back(new ColourDipole(newCol1, -( iAntiJun * 10 + 10 + 2) , 
    -( iJun * 10 + 10 + 2), colReconnection, true, true, false, true));
  int iReal1 = dipoles.size() - 1;
  dipoles.push_back(new ColourDipole(newCol1, -( iAntiJun * 10 + 10 + 2) , 
    -( iJun * 10 + 10 + 2), colReconnection, true, true));
  int iActive1 = dipoles.size() - 1;
 
  // Now make dipole between anti junction and iAcol1.
  // Start by finding real iAcol.
  int iAcol1real 
    = particles[iAcol1].dips[dipoles[iDip1]->iAcolLeg].front()->iAcol;
 
  dipoles.push_back(new ColourDipole(newCol2, -( iAntiJun * 10 + 10), 
    iAcol1real, dipoles[iDip1]->colReconnection, false, true, false, true));
  int iReal2 = dipoles.size() - 1;
  particles[iAcol1].dips[dipoles[iDip1]->iAcolLeg].front() = dipoles.back();
  
  dipoles.push_back(new ColourDipole(newCol2, -( iAntiJun * 10 + 10), 
    iAcol1, dipoles[iDip1]->colReconnection, false, true));
  dipoles.back()->iAcolLeg = dipoles[iDip1]->iAcolLeg;
  int iActive2 = dipoles.size() - 1;
 
  // Now make dipole between anti junction and iAcol1.
  // Start by finding real iAcol.
  int iAcol2real 
    = particles[iAcol2].dips[dipoles[iDip2]->iAcolLeg].front()->iAcol;
  dipoles.push_back(new ColourDipole(newCol3, -( iAntiJun * 10 + 10 + 1), 
    iAcol2real, dipoles[iDip2]->colReconnection, false, true, false, true));
  int iReal3 = dipoles.size() - 1;
  particles[iAcol2].dips[dipoles[iDip2]->iAcolLeg].front() = dipoles.back();
  
  dipoles.push_back(new ColourDipole(newCol3, -( iAntiJun * 10 + 10 + 1), 
    iAcol2, dipoles[iDip2]->colReconnection, false, true));
  dipoles.back()->iAcolLeg = dipoles[iDip2]->iAcolLeg;
  int iActive3 = dipoles.size() - 1;

  // Update already existing dipoles.

  // Update real dipoles.
  particles[iCol1].dips[dipoles[iDip1]->iColLeg].back()->iAcol 
    = - (iJun * 10 + 10);
  particles[iCol2].dips[dipoles[iDip2]->iColLeg].back()->iAcol 
    = - (iJun * 10 + 10 + 1);
  particles[iCol1].dips[dipoles[iDip1]->iColLeg].back()->isJun = true;
  particles[iCol2].dips[dipoles[iDip2]->iColLeg].back()->isJun = true;

  // Update active dipoles.
  dipoles[iDip1]->isJun = true;
  dipoles[iDip2]->isJun = true;
  dipoles[iDip1]->iAcol = - (iJun * 10 + 10);
  dipoles[iDip2]->iAcol = - (iJun * 10 + 10 + 1);
  dipoles[iDip1]->iAcol = - (iJun * 10 + 10);
  dipoles[iDip2]->iAcol = - (iJun * 10 + 10 + 1);
  dipoles[iDip1]->iAcolLeg = 0;
  dipoles[iDip2]->iAcolLeg = 0;
  
  // Update active for anti particles.
  // Normally should only contain active dipoles once,
  // only problem is if the two dipole ends are the same particle.
  for (int i = 0; i < int(particles[iAcol1].activeDips.size()); ++i) 
    if (particles[iAcol1].activeDips[i] == dip1) {
      particles[iAcol1].activeDips[i] = dipoles[iActive2];
      break;
    }
  
  for (int i = 0; i < int(particles[iAcol2].activeDips.size()); ++i) 
    if (particles[iAcol2].activeDips[i] == dip2) {
      particles[iAcol2].activeDips[i] = dipoles[iActive3];
    break;
    }

  // Add the junctions to the event.
  junctions.push_back(Junction(1, oldCol1, oldCol2, newCol1));
  junctions.push_back(Junction(2, newCol2, newCol3, newCol1));
  
  // Set junction information.
  junctions[iJun].dipsOrig[0] = 
    particles[iCol1].dips[dipoles[iDip1]->iColLeg].back();
  junctions[iJun].dipsOrig[1] = 
    particles[iCol2].dips[dipoles[iDip2]->iColLeg].back();
  junctions[iJun].dipsOrig[2] = dipoles[iReal1];
  junctions[iJun].dips[0] = dip1;
  junctions[iJun].dips[1] = dip2;
  junctions[iJun].dips[2] = dipoles[iActive1];
  
  // Set anti junction information.
  junctions[iAntiJun].dips[0] = dipoles[iActive2];
  junctions[iAntiJun].dips[1] = dipoles[iActive3];
  junctions[iAntiJun].dips[2] = dipoles[iActive1];
  junctions[iAntiJun].dipsOrig[0] = dipoles[iReal2];
  junctions[iAntiJun].dipsOrig[1] = dipoles[iReal3];
  junctions[iAntiJun].dipsOrig[2] = dipoles[iReal1];

  // Make pseudo particles if needed.
  if (mDip(dip1) < m0) makePseudoParticle(dip1, 68, true);
  if (mDip(dipoles[iActive2]) < m0) 
    makePseudoParticle(dipoles[iActive2], 68, true);

  // Check dipoles.
  checkDipoles();
  
  // Done.
  return true;
}

//--------------------------------------------------------------------------

// Do a single trial emission.

bool ColourReconnection::singleJunction(Event& event, int iDip1, 
  int iDip2, int iDip3) {

  // Not possible to connect a gluon with itself.
  if(dipoles[iDip1]->iCol == dipoles[iDip2]->iAcol)
    return false;
  if(dipoles[iDip1]->iAcol == dipoles[iDip2]->iCol)
    return false;

  // NEED TO ALLOW FOR ONE JUNCTION AND ONE ANTI JUNCTION.
  // SHOULD BE MINOR PART, SO IGNORE FOR NOW.
  // Do nothing if one of the dipoles is a junction or anti junction.
  if (dipoles[iDip1]->isJun || dipoles[iDip1]->isAntiJun)
    return false;
  if (dipoles[iDip2]->isJun || dipoles[iDip2]->isAntiJun)
    return false;
  if (dipoles[iDip3]->isJun || dipoles[iDip3]->isAntiJun)
    return false;

  // Check that all dipoles are active.
  if (!dipoles[iDip1]->isActive || !dipoles[iDip2]->isActive || 
      !dipoles[iDip3]->isActive) return false;

  // Store particle indecies.
  int iCol1 = dipoles[iDip1]->iCol;
  int iCol2 = dipoles[iDip2]->iCol;
  int iCol3 = dipoles[iDip3]->iCol;
  int iAcol1 = dipoles[iDip1]->iAcol;
  int iAcol2 = dipoles[iDip2]->iAcol;
  int iAcol3 = dipoles[iDip3]->iAcol;
  
  // Only accept one third of the pairs.
  if ( dipoles[iDip1]->colReconnection % 3 
    != dipoles[iDip2]->colReconnection % 3 
    || dipoles[iDip1]->colReconnection % 3 
    != dipoles[iDip3]->colReconnection % 3) return false;
  
  // Need to be either 1-4-7 or 1-1-1.
  if ( !(dipoles[iDip1]->colReconnection == dipoles[iDip2]->colReconnection 
      && dipoles[iDip1]->colReconnection == dipoles[iDip3]->colReconnection)
    && !(dipoles[iDip1]->colReconnection != dipoles[iDip2]->colReconnection
      && dipoles[iDip1]->colReconnection != dipoles[iDip3]->colReconnection 
      && dipoles[iDip2]->colReconnection != dipoles[iDip3]->colReconnection) )
    return false;
  
  // Store dipoles.
  ColourDipole* dip1 = dipoles[iDip1];
  ColourDipole* dip2 = dipoles[iDip2];
  ColourDipole* dip3 = dipoles[iDip3];

  vector<ColourDipole*> oldDips, newDips;

  // Calculate the old string length.
  double oldStringLength = 0;
  oldStringLength += calculateStringLength(dip1,oldDips);
  oldStringLength += calculateStringLength(dip2,oldDips);
  oldStringLength += calculateStringLength(dip3,oldDips);
  
  // Calculate the new string length.
  double newStringLength = 0;
  newDips.push_back(dip1);
  newDips.push_back(dip2);
  newDips.push_back(dip3);
  
  // Calcualte the new string length, start by doing the colour junction.
  // Different cases depending on whether the new particle forms 
  // a pseudo particle.
  
  // If all particles should be combined to a single pseudoparticle.
  if ( (particles[iCol1].p() + particles[iCol2].p() 
       + particles[iCol3].p() ).mCalc() < m0) {

    // Make test particles
    particles.push_back(particles[iCol1]);
    particles.back().p( particles[iCol1].p() + particles[iCol2].p() 
		      + particles[iCol3].p() );
    int iNew1 = particles.size() -1;
    
    // Calculate neighbour lengths.
    for (int i = 0;i < int(particles[iCol1].activeDips.size()); ++i)
      newStringLength += neighbourLength(particles[iCol1].activeDips[i], 
        newDips, iCol1, iCol2, iNew1, iCol3, -7, iNew1);

    for (int i = 0;i < int(particles[iCol2].activeDips.size()); ++i)
      newStringLength += neighbourLength(particles[iCol2].activeDips[i], 
        newDips, iCol1, iCol2, iNew1, iCol3, -7, iNew1);

    for (int i = 0;i < int(particles[iCol3].activeDips.size()); ++i)
      newStringLength += neighbourLength(particles[iCol3].activeDips[i], 
        newDips, iCol1, iCol2, iNew1, iCol3, -7, iNew1);

    // Remove test particle
    particles.pop_back();
    
  // If a single pseudo particle needs to be made.
  } else if ((particles[iCol1].p() + particles[iCol2].p()).mCalc() < m0 ||
	   (particles[iCol1].p() + particles[iCol3].p()).mCalc() < m0 ||
	   (particles[iCol2].p() + particles[iCol3].p()).mCalc() < m0) {

    // Find the pair to that should be combined.
    int i1 = iCol1,i2 = iCol2, i3 = iCol3;
    if ((particles[iCol1].p() + particles[iCol3].p()).mCalc() < 
	min( (particles[iCol1].p() + particles[iCol2].p()).mCalc(), 
	     (particles[iCol2].p() + particles[iCol3].p()).mCalc())) {
      i2 = iCol3;
      i3 = iCol2;
    } else if ((particles[iCol2].p() + particles[iCol3].p()).mCalc() < 
	min( (particles[iCol1].p() + particles[iCol2].p()).mCalc(), 
	     (particles[iCol1].p() + particles[iCol3].p()).mCalc())) {
      i1 = iCol3;
      i2 = iCol1;
    }

    // Make test particles.
    particles.push_back(particles[i1]);
    particles.back().p(particles[i1].p() + particles[i2].p());
    int iNew1 = particles.size() -1;

    // Calculate new neighbour lengths.
    for (int i = 0; i < int(particles[i1].activeDips.size()); ++i)
      newStringLength += neighbourLength(particles[i1].activeDips[i], 
	 newDips, i1, i2, iNew1);

    for (int i = 0; i < int(particles[i2].activeDips.size()); ++i)
      newStringLength += neighbourLength(particles[i2].activeDips[i], 
        newDips, i1, i2, iNew1);
    
    // Add the internal junction string.
    newStringLength += calculateStringLength(iNew1,i3);
    
    // Remove test particle.
    particles.pop_back();

  // No pseudo particles needed.
  } else {
    newStringLength += calculateJunctionLength(iCol1,iCol2,iCol3);
  }
  
  // Do it for the anti junction.

  // If all particles should be combined to a single pseudoparticle.
  if ( (particles[iAcol1].p() + particles[iAcol2].p() 
	+ particles[iAcol3].p() ).mCalc() < m0) {

    // Make test particles
    particles.push_back(particles[iAcol1]);
    particles.back().p(particles[iAcol1].p() + particles[iAcol2].p() 
		       + particles[iAcol3].p());
    int iNew1 = particles.size() -1;
    
    // Calculate new neightbour lengths.
    for (int i = 0;i < int(particles[iAcol1].activeDips.size()); ++i)
      newStringLength += neighbourLength(particles[iAcol1].activeDips[i], 
        newDips, iAcol1, iAcol2, iNew1, iAcol3, -7, iNew1);

    for (int i = 0;i < int(particles[iAcol2].activeDips.size()); ++i)
      newStringLength += neighbourLength(particles[iAcol2].activeDips[i], 
        newDips, iAcol1, iAcol2, iNew1, iAcol3, -7, iNew1);

    for (int i = 0;i < int(particles[iAcol3].activeDips.size()); ++i)
      newStringLength += neighbourLength(particles[iAcol3].activeDips[i], 
        newDips, iAcol1, iAcol2, iNew1, iAcol3, -7, iNew1);

    // Remove test particle.
    particles.pop_back();
    
  // If a single pseudo particle needs to be made.
  } else if((particles[iAcol1].p() + particles[iAcol2].p()).mCalc() < m0 ||
	    (particles[iAcol1].p() + particles[iAcol3].p()).mCalc() < m0 ||
	    (particles[iAcol2].p() + particles[iAcol3].p()).mCalc() < m0) {

    // Find the pair to that should be combined.
    int i1 = iAcol1,i2 = iAcol2, i3 = iAcol3;
    if ((particles[iAcol1].p() + particles[iAcol3].p()).mCalc() < 
	min( (particles[iAcol1].p() + particles[iAcol2].p()).mCalc(), 
	     (particles[iAcol2].p() + particles[iAcol3].p()).mCalc())) {
      i2 = iAcol3;
      i3 = iAcol2;
    } else if((particles[iAcol2].p() + particles[iAcol3].p()).mCalc() < 
	min( (particles[iAcol1].p() + particles[iAcol2].p()).mCalc(), 
	     (particles[iAcol1].p() + particles[iAcol3].p()).mCalc())) {
      i1 = iAcol3;
      i2 = iAcol1;
    }

    // Make test particles.
    particles.push_back(particles[i1]);
    particles.back().p(particles[i1].p() + particles[i2].p());
    int iNew1 = particles.size() -1;

    // Calculate new neighbour lengths.
    for (int i = 0;i < int(particles[i1].activeDips.size()); ++i)
      newStringLength += neighbourLength(particles[i1].activeDips[i], 
	 newDips, i1, i2, iNew1);

    for (int i = 0;i < int(particles[i2].activeDips.size()); ++i)
      newStringLength += neighbourLength(particles[i2].activeDips[i], 
        newDips, i1, i2, iNew1);
    
    // Add the internal junction string.
    newStringLength += calculateStringLength(iNew1,i3);
    
    // Remove test particle.
    particles.pop_back();

  // No pseudo particles needed.
  } else {
    newStringLength += calculateJunctionLength(iAcol1,iAcol2,iAcol3);
  }

  // Add dipoles to old length that was included in the new length calculation.
  oldStringLength += calculateAdditionalLengths(oldDips, newDips);
    
  if (newStringLength > oldStringLength)
    return false;
 
  // Append new junctions.
  int oldCol1 = dipoles[iDip1]->col;
  int oldCol2 = dipoles[iDip2]->col;
  int oldCol3 = dipoles[iDip3]->col;
  int newCol1 = event.nextColTag();
  int newCol2 = event.nextColTag();
  int newCol3 = event.nextColTag();
  
  // Store new junction indices.
  int iJun = junctions.size();
  int iAntiJun = junctions.size() + 1;
  
  // Now make dipole between anti junction and iAcol1.
  // Start by finding real iAcol.
  int iAcol1real 
    = particles[iAcol1].dips[dipoles[iDip1]->iAcolLeg].front()->iAcol;
  dipoles.push_back(new ColourDipole(newCol1, -( iAntiJun * 10 + 10), 
    iAcol1real, dipoles[iDip1]->colReconnection, false, true, false, true));
  int iReal1 = dipoles.size() - 1;
  particles[iAcol1].dips[dipoles[iDip1]->iAcolLeg].front() = dipoles.back();

  dipoles.push_back(new ColourDipole(newCol1, -( iAntiJun * 10 + 10), 
    iAcol1, dipoles[iDip1]->colReconnection, false, true));
  dipoles.back()->iAcolLeg = dipoles[iDip1]->iAcolLeg;
  int iActive1 = dipoles.size() - 1;
 
  // Now make dipole between anti junction and iAcol2.
  // Start by finding real iAcol2.
  int iAcol2real 
    = particles[iAcol2].dips[dipoles[iDip2]->iAcolLeg].front()->iAcol;
  dipoles.push_back(new ColourDipole(newCol2, -( iAntiJun * 10 + 10 + 1), 
    iAcol2real, dipoles[iDip2]->colReconnection, false, true, false, true));
  int iReal2 = dipoles.size() - 1;
  particles[iAcol2].dips[dipoles[iDip2]->iAcolLeg].front() = dipoles.back();
  
  dipoles.push_back(new ColourDipole(newCol2, -( iAntiJun * 10 + 10 + 1), 
    iAcol2, dipoles[iDip2]->colReconnection, false, true));
  dipoles.back()->iAcolLeg = dipoles[iDip2]->iAcolLeg;
  int iActive2 = dipoles.size() - 1;

  // Now make dipole between anti junction and iAcol3.
  // Start by finding real iAcol3.
  int iAcol3real 
    = particles[iAcol3].dips[dipoles[iDip3]->iAcolLeg].front()->iAcol;
  dipoles.push_back(new ColourDipole(newCol3, -( iAntiJun * 10 + 10 + 2), 
    iAcol3real, dipoles[iDip3]->colReconnection, false, true, false, true));
  int iReal3 = dipoles.size() - 1;
  particles[iAcol3].dips[dipoles[iDip3]->iAcolLeg].front() = dipoles.back();
  
  dipoles.push_back(new ColourDipole(newCol3, -( iAntiJun * 10 + 10 + 2), 
    iAcol3, dipoles[iDip3]->colReconnection, false, true));
  dipoles.back()->iAcolLeg = dipoles[iDip3]->iAcolLeg;
  int iActive3 = dipoles.size() - 1;

  // Update already existing dipoles.

  // Update real dipoles.
  particles[iCol1].dips[dipoles[iDip1]->iColLeg].back()->iAcol 
    = - (iJun * 10 + 10);
  particles[iCol2].dips[dipoles[iDip2]->iColLeg].back()->iAcol 
    = - (iJun * 10 + 10 + 1);
  particles[iCol3].dips[dipoles[iDip3]->iColLeg].back()->iAcol 
    = - (iJun * 10 + 10 + 2);
  particles[iCol1].dips[dipoles[iDip1]->iColLeg].back()->isJun = true;
  particles[iCol2].dips[dipoles[iDip2]->iColLeg].back()->isJun = true;
  particles[iCol3].dips[dipoles[iDip3]->iColLeg].back()->isJun = true;

  // Update active dipoles.
  dipoles[iDip1]->isJun = true;
  dipoles[iDip2]->isJun = true;
  dipoles[iDip3]->isJun = true;
  dipoles[iDip1]->iAcol = - (iJun * 10 + 10);
  dipoles[iDip2]->iAcol = - (iJun * 10 + 10 + 1);
  dipoles[iDip3]->iAcol = - (iJun * 10 + 10 + 2);
  dipoles[iDip1]->iAcolLeg = 0;
  dipoles[iDip2]->iAcolLeg = 0;
  dipoles[iDip3]->iAcolLeg = 0;
  
  // Update active dipoles for anti particles.
  for (int i = 0; i < int(particles[iAcol1].activeDips.size()); ++i) 
    if (particles[iAcol1].activeDips[i] == dip1)
      particles[iAcol1].activeDips[i] = dipoles[iActive1];
  for (int i = 0; i < int(particles[iAcol2].activeDips.size()); ++i) 
    if (particles[iAcol2].activeDips[i] == dip2)
      particles[iAcol2].activeDips[i] = dipoles[iActive2];
  for (int i = 0; i < int(particles[iAcol3].activeDips.size()); ++i) 
    if (particles[iAcol3].activeDips[i] == dip3)
      particles[iAcol3].activeDips[i] = dipoles[iActive3];
  
  // Add the junctions to the event.
  junctions.push_back(Junction(1, oldCol1, oldCol2, oldCol3));
  junctions.push_back(Junction(2, newCol1, newCol3, newCol3));

  // Update junction ends.
  junctions[iJun].dipsOrig[0] = 
    particles[iCol1].dips[dipoles[iDip1]->iColLeg].back();
  junctions[iJun].dipsOrig[1] = 
    particles[iCol2].dips[dipoles[iDip2]->iColLeg].back();
  junctions[iJun].dipsOrig[2] =
    particles[iCol3].dips[dipoles[iDip3]->iColLeg].back();
  junctions[iJun].dips[0] = dip1;
  junctions[iJun].dips[1] = dip2;
  junctions[iJun].dips[2] = dip3;

  // Update the anti junction.
  junctions[iAntiJun].dips[0] = dipoles[iActive1];
  junctions[iAntiJun].dips[1] = dipoles[iActive2];
  junctions[iAntiJun].dips[2] = dipoles[iActive3];
  junctions[iAntiJun].dipsOrig[0] = dipoles[iReal1];
  junctions[iAntiJun].dipsOrig[1] = dipoles[iReal2];
  junctions[iAntiJun].dipsOrig[2] = dipoles[iReal3];

  // Make pseudo particles if needed.
  if (dip1->isActive && mDip(dip1) < m0)
    makePseudoParticle(dip1, 68, true);
  if (dip2->isActive && mDip(dip2) < m0)
    makePseudoParticle(dip2, 68, true);
  if (dip3->isActive && mDip(dip3) < m0)
    makePseudoParticle(dip3, 68, true);
  
  if (dipoles[iActive1]->isActive && mDip(dipoles[iActive1]) < m0)
    makePseudoParticle(dipoles[iActive1], 68, true);
  if (dipoles[iActive2]->isActive && mDip(dipoles[iActive2]) < m0)
    makePseudoParticle(dipoles[iActive2], 68, true);
  if (dipoles[iActive3]->isActive && mDip(dipoles[iActive3]) < m0)
    makePseudoParticle(dipoles[iActive3], 68, true);
  
  // Check dipoles.
  checkDipoles();
  
  // Done.
  return true;
}

// ------------------------------------------------------------------

void ColourReconnection::makePseudoParticle(ColourDipole* dip , int status, 
  bool setupDone) {
  
  // If it is a normal dipole that needs to be combined.
  if (!dip->isJun && !dip->isAntiJun) {
    // Start by storing variables for easier use.
    int iCol = dip->iCol;
    int iAcol = dip->iAcol;
    int iColLeg = dip->iColLeg;
    int iAcolLeg = dip->iAcolLeg;
   
    // Make new pseudo particle.
    int iNew = particles.size();
    particles.push_back(particles[iCol]);
    particles[iNew].acol(particles[iCol].acol());
    particles[iNew].col(particles[iAcol].col());
    particles[iNew].mother1(iCol);
    particles[iNew].mother2(iAcol);
    particles[iNew].id(99);
    particles[iNew].status(status);
    particles[iNew].isJun = false;
    particles[iAcol].statusNeg();
    particles[iAcol].daughter1(iNew);
    particles[iCol].statusNeg();
    particles[iCol].daughter1(iNew);
    if (iCol != iAcol)
      particles[iNew].p(particles[iCol].p() + particles[iAcol].p());
    else
      particles[iNew].p(particles[iCol].p());

    // Add all the dipoles from the old pseudo particle.
    // First from particle 1.
    particles[iNew].dips = particles[dip->iCol].dips;
    particles[iNew].colEndIncluded = particles[dip->iCol].colEndIncluded;
    particles[iNew].acolEndIncluded = particles[dip->iCol].acolEndIncluded;
    
    // Then particle 2.
    if(iCol != iAcol) {
      for(int i = 0;i < int(particles[dip->iAcol].dips.size()); ++i)  {
        if (i != dip->iAcolLeg) {
          // If it is not the same leg, add as separate vector.
	  particles[iNew].dips.push_back(particles[dip->iAcol].dips[i]);
	  particles[iNew].colEndIncluded.push_back(
            particles[dip->iAcol].colEndIncluded[i]);
	  particles[iNew].acolEndIncluded.push_back(
            particles[dip->iAcol].acolEndIncluded[i]);
        } // If it is the same leg, at at the end of the vector.
        else {
	  particles[iNew].acolEndIncluded[iColLeg] = 
            particles[iAcol].acolEndIncluded[i];
	  particles[iNew].dips[iColLeg].pop_back();
	  particles[iNew].dips[iColLeg].insert(
	    particles[iNew].dips[iColLeg].end(), 
            particles[iAcol].dips[i].begin(), particles[iAcol].dips[i].end() );
        }
      }
    }
    if (iCol != iAcol) {
      // Update the dipole legs to the new particle.
      for (int i = 0; i < int(particles[iAcol].activeDips.size()); ++i) {
        if ( particles[iAcol].activeDips[i]->iAcol == iAcol) {
          if (particles[iAcol].activeDips[i]->iAcolLeg < iAcolLeg) 
	    particles[iAcol].activeDips[i]->iAcolLeg += 
              particles[iCol].dips.size();
	  else if (particles[iAcol].activeDips[i]->iAcolLeg == iAcolLeg)
	    particles[iAcol].activeDips[i]->iAcolLeg = iColLeg;
	  else if (particles[iAcol].activeDips[i]->iAcolLeg > iAcolLeg)
	    particles[iAcol].activeDips[i]->iAcolLeg += 
              particles[iCol].dips.size() - 1;
        }
        if (particles[iAcol].activeDips[i]->iCol == iAcol) {
	  if (particles[iAcol].activeDips[i]->iColLeg < iAcolLeg)
	    particles[iAcol].activeDips[i]->iColLeg += 
              particles[iCol].dips.size();
	  else if (particles[iAcol].activeDips[i]->iColLeg == iAcolLeg)
	    particles[iAcol].activeDips[i]->iColLeg = iColLeg;
	  else if (particles[iAcol].activeDips[i]->iColLeg > iAcolLeg)
	    particles[iAcol].activeDips[i]->iColLeg += 
              particles[iCol].dips.size() - 1;
        }
      }
    }

    // Update list of active dipoles.  
    particles[iNew].activeDips.clear();
    particles[iNew].activeDips.insert(particles[iNew].activeDips.end(),
      particles[iCol].activeDips.begin(), particles[iCol].activeDips.end());
    if (iCol != iAcol)
      particles[iNew].activeDips.insert(particles[iNew].activeDips.end(),
        particles[iAcol].activeDips.begin(), particles[iAcol].activeDips.end());
    
    // Remove the now inactive dipole.
    for (int i = 0; i < int(particles[iNew].activeDips.size()); ++i)
      if (particles[iNew].activeDips[i] == dip) {
	particles[iNew].activeDips.erase(
          particles[iNew].activeDips.begin() + i);
	i--;
      }
    
    // Update the indecies in the active dipoles.
    for (int i = 0; i < int(particles[iNew].activeDips.size()); ++i) {
      if (particles[iNew].activeDips[i]->iCol == iAcol)
	particles[iNew].activeDips[i]->iCol = iNew;
      if (particles[iNew].activeDips[i]->iCol == iCol)
	particles[iNew].activeDips[i]->iCol = iNew;
      if (particles[iNew].activeDips[i]->iAcol == iAcol)
	particles[iNew].activeDips[i]->iAcol = iNew;
      if (particles[iNew].activeDips[i]->iAcol == iCol)
	particles[iNew].activeDips[i]->iAcol = iNew;
      particles[iNew].activeDips[i]->p1p2 
        = mDip(particles[iNew].activeDips[i]);  
    }
    
    // If it is a combination of the same particle,
    // check if any double active dipoles
    if (iCol == iAcol)
    for (int i = 0; i < int(particles[iNew].activeDips.size()); ++i) 
    for (int j = i + 1; j < int(particles[iNew].activeDips.size()); ++j) 
    if (particles[iNew].activeDips[i] == particles[iNew].activeDips[j]) {
      particles[iNew].activeDips.erase(particles[iNew].activeDips.begin() + j);
      j--;
    }
    
    // mark the internal dipole as not active.
    dip->isActive = false;
    
    // Done.
    return;
  }
  
  // If both ends are connected to a junction something went wrong!
  else if (dip->isJun && dip->isAntiJun) {
    return;
  }
  else {

    // Find junction index and first leg to combine.
    int iJun, i0, i1, i2, junLeg0, junLeg1, junLeg2;
    getJunctionIndicies(dip, iJun, i0, i1, i2, junLeg0, junLeg1, junLeg2);    
    ColourDipole* dip2 = junctions[iJun].dips[junLeg1];
    ColourDipole* dip3 = junctions[iJun].dips[junLeg2];

    // Add new particle.
    int iNew = particles.size();
    particles.push_back(ColourParticle( Particle( 99, status, i0, i1, 0, 0, 0, 
      0, particles[i0].p() + particles[i1].p() ) ) ); 
    particles[iNew].isJun = true;
    particles[iNew].junKind = junctions[iJun].kind();
    if (i0 == i1) particles[iNew].p(particles[i0].p());
    
    // Update old particles.
    particles[i0].statusNeg();
    particles[i0].daughter1(iNew);
    particles[i1].statusNeg();
    particles[i1].daughter1(iNew);
    
    // Update list of internal dipoles.
    particles[iNew].dips.clear();
    particles[iNew].dips.insert(particles[iNew].dips.end(),
      particles[i0].dips.begin(),particles[i0].dips.end());
    if (i0 != i1)
      particles[iNew].dips.insert(particles[iNew].dips.end(),
        particles[i1].dips.begin(),particles[i1].dips.end());

    // Update list of whether colour ending is included.
    particles[iNew].colEndIncluded.clear();
    particles[iNew].colEndIncluded.insert(
      particles[iNew].colEndIncluded.end(),
      particles[i0].colEndIncluded.begin(), 
      particles[i0].colEndIncluded.end() );
    if (i0 != i1)
      particles[iNew].colEndIncluded.insert(
        particles[iNew].colEndIncluded.end(),
        particles[i1].colEndIncluded.begin(),
        particles[i1].colEndIncluded.end() );
   
    // Update list of whether anti colour ending is included.
    particles[iNew].acolEndIncluded.clear();
    particles[iNew].acolEndIncluded.insert(
      particles[iNew].acolEndIncluded.end(),
      particles[i0].acolEndIncluded.begin(), 
      particles[i0].acolEndIncluded.end() );
    if (i0 != i1)
      particles[iNew].acolEndIncluded.insert(
        particles[iNew].acolEndIncluded.end(),
        particles[i1].acolEndIncluded.begin(), 
        particles[i1].acolEndIncluded.end() );
   
    // Third particle just need to add one to list of dipoles.
    if (dip->isJun && i2 >= 0 && i2 != i0 && i2 != i1) {
      particles[iNew].dips.push_back(particles[i2].dips[dip3->iColLeg]);
      particles[iNew].dips.back().erase(particles[iNew].dips.back().begin(), 
        particles[iNew].dips.back().end() - 1);
      
    }
    if (dip->isAntiJun && i2 >= 0 && i2 != i0 && i2 != i1) {
      particles[iNew].dips.push_back(particles[i2].dips[dip3->iAcolLeg]);
      particles[iNew].dips.back().erase(
        particles[iNew].dips.back().begin() + 1, 
        particles[iNew].dips.back().end() );
    }

    // Add endings for the third particle.
    if (i2 != i0 && i2 != i1) {
      particles[iNew].acolEndIncluded.push_back(false);
      particles[iNew].colEndIncluded.push_back(false);
    }

    // Special case if it is J-J connection.
    if (i2 < 0) {
      particles[iNew].dips.push_back(vector<ColourDipole *>());
      
      // Find the real dipole to add to dipole list.
      for (int i = 0; i < int(dipoles.size()); ++i)
	if (dipoles[i]->isReal && dipoles[i]->iCol == dip3->iCol &&
	    dipoles[i]->iAcol == dip3->iAcol)
	  particles[iNew].dips.back().push_back(dipoles[i]);
      
      // Change ending.
      particles[iNew].acolEndIncluded.back() = true;
      particles[iNew].colEndIncluded.back()  = true;
    }
    
    // The endings need to reflect the new junction structure.
    if (dip->isJun)
    for (int i = 0; i < int(particles[iNew].acolEndIncluded.size()); ++i)
      particles[iNew].acolEndIncluded[i] = true;
    else
    for (int i = 0; i < int(particles[iNew].colEndIncluded.size()); ++i)
      particles[iNew].colEndIncluded[i] = true;
    
    // Update active dipoles, first junction case.
    // Set the now internal dipoles as inactive.    
    dip->isActive = false;
    dip2->isActive = false;
    dip3->isActive = true;

    // Update the dipole legs to the new particle.
    // Only need to do it for the iAcol particle, 
    // since nothing changes for the iCol particle.
    if (i0 != i1) 
    for (int i = 0; i < int(particles[i1].activeDips.size()); ++i) {
      if (particles[i1].activeDips[i]->iAcol == i1) 
	particles[i1].activeDips[i]->iAcolLeg += particles[i0].dips.size();

      if (particles[i1].activeDips[i]->iCol == i1) 
	particles[i1].activeDips[i]->iColLeg += particles[i0].dips.size();	
    }
    
    // Update list of active dipoles.  
    particles[iNew].activeDips.clear();
    particles[iNew].activeDips.insert(particles[iNew].activeDips.end(),
      particles[i0].activeDips.begin(), particles[i0].activeDips.end());
    if (i0 != i1)
      particles[iNew].activeDips.insert(particles[iNew].activeDips.end(),
        particles[i1].activeDips.begin(), particles[i1].activeDips.end());
    if (i2 != i0 && i2 != i1)
      particles[iNew].activeDips.push_back(dip3);
    
    // Remove the now inactive dipoles.
    for (int i = 0; i < int(particles[iNew].activeDips.size()); ++i) {
      if (particles[iNew].activeDips[i] == dip) {
	particles[iNew].activeDips.erase(
          particles[iNew].activeDips.begin() + i);
	i--;
	continue;
      }
      if (particles[iNew].activeDips[i] == dip2) {
	particles[iNew].activeDips.erase(
          particles[iNew].activeDips.begin() + i);
	i--;
	continue;
      }
    }

    // Update the indecies in the active dipoles.
    for (int i = 0; i < int(particles[iNew].activeDips.size()); ++i) {
      if (particles[iNew].activeDips[i]->iCol == i1)
	particles[iNew].activeDips[i]->iCol = iNew;
      if (particles[iNew].activeDips[i]->iCol == i0)
	particles[iNew].activeDips[i]->iCol = iNew;
      if (particles[iNew].activeDips[i]->iAcol == i1)
	particles[iNew].activeDips[i]->iAcol = iNew;
      if (particles[iNew].activeDips[i]->iAcol == i0)
	particles[iNew].activeDips[i]->iAcol = iNew;
      particles[iNew].activeDips[i]->p1p2 
        = mDip(particles[iNew].activeDips[i]);
    }

    // The third dip is no longer connected to a junction.
    if (dip->isJun) {
      dip3->isJun = false;
      dip3->iAcol = iNew;
      if (i2 != i0 && i2 != i1)
        dip3->iAcolLeg = particles[iNew].dips.size() - 1;
    }
    else  {
      dip3->isAntiJun = false;
      dip3->iCol = iNew;
      if (i2 != i0 && i2 != i1)
        dip3->iColLeg = particles[iNew].dips.size() - 1;
    }
    
    // Possible for the new dip to have a low m0.
    if (setupDone && mDip(dip3) < m0)
      makePseudoParticle(dip3,status, true);
  }
  
  // Done.

}

// ------------------------------------------------------------------

// Help function to sort dipoles in right order.

bool sortFunc(ColourDipole* a, ColourDipole* b) {

    return (a->p1p2 < b->p1p2);
  
}

// ------------------------------------------------------------------

// Form all pseudoparticles below m0.

void ColourReconnection::makeAllPseudoParticles( Event & event, int iFirst) {
  
  // Make junctions.
  for (int i = 0; i < event.sizeJunction(); ++i) 
    junctions.push_back(event.getJunction(i));

  // Make new copy of all the dipoles.
  int oldSize = int(dipoles.size());
  for (int i = 0; i < oldSize; ++i) {
    dipoles.push_back(new ColourDipole(*dipoles[i]));
    dipoles[i + oldSize]->iColLeg = 0;
    dipoles[i + oldSize]->iAcolLeg = 0;
    dipoles[i]->iColLeg = 0;
    dipoles[i]->iAcolLeg = 0;
    dipoles[i]->isActive = false;
    dipoles[i]->isReal = true;
    dipoles[i + oldSize]->isReal = false;

    // Store original dipoles connected to junctions.
    if (dipoles[i]->iCol < 0) {
      junctions[-(dipoles[i]->iCol / 10 + 1)].dipsOrig[(-dipoles[i]->iCol) 
        % 10] = dipoles[i];
    }
    if (dipoles[i]->iAcol < 0) {
      junctions[-(dipoles[i]->iAcol / 10 + 1)].dipsOrig[-(dipoles[i]->iAcol 
        % 10)] = dipoles[i];
    }
  }
   
  // Set up the coldDips and acolDips.
  for (int i = 0; i < oldSize; ++i) {
    if (dipoles[i]->leftDip != 0) 
    for (int j = 0; j < oldSize; ++j)
    if (dipoles[i]->leftDip == dipoles[j]) {
      dipoles[i + oldSize]->colDips.push_back(dipoles[j + oldSize]);
      break;
    }
   
    if (dipoles[i]->rightDip != 0) 
    for (int j = 0; j < oldSize; ++j)
    if (dipoles[i]->rightDip == dipoles[j]) {
      dipoles[i + oldSize]->acolDips.push_back(dipoles[j + oldSize]);
      break;
    }
  }

  // Start by copying event record to make pseudoparticles.
  // The pseudoparticles also need to gain
  for (int i = iFirst; i < event.size(); ++i) 
  if (event[i].isFinal()) {
    particles.push_back(ColourParticle(event[i]));
    particles.back().dips.resize(1,vector<ColourDipole *>());
    
    // Set up dipoles.
    for (int j = 0; j < int(dipoles.size()); ++j) {
      if (dipoles[j]->iCol == i) { 
	if (dipoles[j]->isActive) {
	  dipoles[j]->iCol = particles.size() - 1;
	  particles.back().activeDips.push_back(dipoles[j]);
	}
	else particles.back().dips[0].push_back(dipoles[j]);
      }

      if (dipoles[j]->iAcol == i) {
	if (dipoles[j]->isActive) {
	  dipoles[j]->iAcol = particles.size() - 1;
	  particles.back().activeDips.push_back(dipoles[j]);
	}
	else particles.back().dips[0].insert(particles.back().dips[0].begin(),
          dipoles[j]);
      }
    }
    
    // Tell whether dipoles are connected to other dipoles.
    if (event[i].isQuark() && event[i].id() > 0)
      particles.back().colEndIncluded.push_back(true);
    else particles.back().colEndIncluded.push_back(false);
    
    if (event[i].isQuark() && event[i].id() < 0)
      particles.back().acolEndIncluded.push_back(true);
    else particles.back().acolEndIncluded.push_back(false);
  }
  
  // Inserting a copy of the event record, but now with full 
  // pseudo particle setup.
  // This is mainly to avoid having to distinguish between combining
  // original particles and pseudoparticles.
  
  // Set right dipole connections in junctions.
  for (int i = 0; i < int(dipoles.size()); ++i) {
    if (dipoles[i]->iCol < 0) {
      int j = (- dipoles[i]->iCol / 10) - 1;
      int jLeg = - dipoles[i]->iCol % 10;
      junctions[j].setColDip(jLeg, dipoles[i]);
    }
    if (dipoles[i]->iAcol < 0) {
      int j = (- dipoles[i]->iAcol / 10) - 1;
      int jLeg = - dipoles[i]->iAcol % 10;
      junctions[j].setColDip(jLeg, dipoles[i]);
    }
  }

  // Make sure all dipoles masses are set correctly.
  for (int i = 0; i < int(dipoles.size()); ++i) {
    if (dipoles[i]->isActive)
      dipoles[i]->p1p2 = mDip(dipoles[i]);
    else
      dipoles[i]->p1p2 = 1E9;
  }

  // Keep making pseudo particles until they are above the threshold.
  while (true) {
    sort(dipoles.begin(),dipoles.end(),sortFunc);
    bool finished = true;
    for (int i = 0; i < int(dipoles.size()); ++i) {
      if (!dipoles[i]->isActive) continue;
      if (dipoles[i]->p1p2 < m0) {
	makePseudoParticle( dipoles[i], 68);
	finished = false;
	break;
      }
      else break;
    }
    if (finished) break;
  }
  // Sort the dipoles.
  sort(dipoles.begin(),dipoles.end(),sortFunc);
  
  // Done.
  return;

}

// ------------------------------------------------------------------

// Print statements if something is wrong in dipole setup.
// Does not have a return statement.

void ColourReconnection::checkDipoles() {

  for (int i = 0;i < int(dipoles.size()); ++i) {
    if (dipoles[i] == 0) { cout << "dipole empty" << endl;}
    if (dipoles[i]->isActive) {
      if (dipoles[i]->iCol >= 0) {
	bool foundMyself = false;
	for (int j = 0; j < int(particles[ dipoles[i]->iCol ].
          activeDips.size()); ++j) {
	  if (!particles[dipoles[i]->iCol].activeDips[j]->isActive) {
	    infoPtr->errorMsg("Error in ColourReconnection::checkDipoles:"
	      "Found inactive dipole, where only active was expected");
	  }
	  if (particles[dipoles[i]->iCol].activeDips[j] == dipoles[i])
	    foundMyself = true;
	}

	if (!foundMyself) {
	  infoPtr->errorMsg("Error in ColourReconnection::checkDipoles:"
	    "Linking between active dipoles and particles is wrong");
	}
	if (dipoles[i]->iColLeg 
          >= int(particles[dipoles[i]->iCol].dips.size())) {
	  infoPtr->errorMsg("Error in ColourReconnection::checkDipoles:"
	    "Original dipoles not stored correct");
	}

	// Check that linking to old dipoles work.
	if (dipoles[i]->col != 
	   particles[dipoles[i]->iCol].dips[dipoles[i]->iColLeg].back()->col) {
	   infoPtr->errorMsg("Error in ColourReconnection::checkDipoles:"
	    "Original dipoles do not match in");
	}
      }

      if (dipoles[i]->iAcol >= 0) {
	bool foundMyself = false;
	for (int j = 0;j < int(particles[ dipoles[i]->iAcol ].
	  activeDips.size()); ++j) {
	 
	  if (!particles[dipoles[i]->iAcol].activeDips[j]->isActive) {
	    infoPtr->errorMsg("Error in ColourReconnection::checkDipoles:"
	      "Found inactive dipole, where only active was expected");
	  }
	   if (particles[dipoles[i]->iAcol].activeDips[j] == dipoles[i])
	    foundMyself = true;
	}

	if (!foundMyself) {
	   infoPtr->errorMsg("Error in ColourReconnection::checkDipoles:"
	    "Linking between active dipoles and particles is wrong");
	}
	if (dipoles[i]->iAcolLeg >= int(particles[dipoles[i]->iAcol].
          dips.size() )) {
	  infoPtr->errorMsg("Error in ColourReconnection::checkDipoles:"
	    "Original dipoles not stored correct");
	}

	// Check that linking to old dipoles work
      	if (dipoles[i]->col != particles[dipoles[i]->iAcol].
	    dips[dipoles[i]->iAcolLeg].front()->col) {
	   infoPtr->errorMsg("Error in ColourReconnection::checkDipoles:"
	    "Original dipoles do not match in");
	}
      }
    }
  }
}

// ------------------------------------------------------------------

// Print all the chains.

void ColourReconnection::listAllChains() {

  cout << "  ----- PRINTING CHIANS -----  " << endl;
  for(int i = 0;i < int(dipoles.size());i++)
    dipoles[i]->printed = false;

  for(int i = 0;i < int(dipoles.size()); ++i)
    if(!dipoles[i]->printed)
      listChain(dipoles[i]);
  cout << "  ----- PRINTED CHIANS -----  " << endl;

}

// ------------------------------------------------------------------
  
// Print the chain containing the dipole.

void ColourReconnection::listChain(ColourDipole *dip) {

  // THIS LINE INVALIDATES THE REST.
  return ;

  // Make sure not an empty pointer.
  if(dip == 0) return;

  // If chain is not active, just print it.
  if (!dip->isActive) {
    cout << "inactive:  " << dip->iCol << " (" << dip->p1p2 << ", " 
	 << dip->col << ") " << dip->iAcol << endl;
    return;
  }

  ColourDipole * colDip = dip;
  bool isGluonLoop = false;
  for(int i =0;i < int(dipoles.size());++i)
    dipoles[i]->inChain = false;
  
  // Try to reach one end of the chain.
  while (int(colDip->colDips.size()) == 1) {
    isGluonLoop = true;
    colDip->inChain = true;
    
    bool reachedEnd = false;
    for (int i = 0; i < int(colDip->colDips[0]->colDips.size()); ++i) {

      // Check if end is reached.
      if (colDip == colDip->colDips[0]->colDips[i]) {
	reachedEnd = true;
	if (colDip->iCol == colDip->colDips[0]->colDips[i]->iAcol) 
	  isGluonLoop = true;	
      }
    }
    
    if (reachedEnd) break;
    
    colDip = colDip->colDips[0];   
    if (colDip->inChain) break;
  }
  
  // Set all dipoles to not used.
  for(int i =0;i < int(dipoles.size());++i)
    dipoles[i]->inChain = false;
  
  // Start the printing.
  while (colDip->acolDips.size() == 1) {
    cout << colDip->iCol << " (" << colDip->p1p2 <<  ", " << colDip->col 
	 << ") (" << colDip->isActive << ")"; 
    colDip->printed = true;
    colDip->inChain = true;
    
    bool reachedEnd = false;
    for (int i = 0; i < int(colDip->acolDips[0]->acolDips.size()); ++i) 
    if (colDip == colDip->acolDips[0]->acolDips[i]) {
      reachedEnd = true;
      if (colDip->iCol == colDip->acolDips[0]->acolDips[i]->iAcol) 
	isGluonLoop = true;
    }
    
    // If the end is reached.
    if (colDip->acolDips[0]->inChain || (reachedEnd && !isGluonLoop)) {
      cout << " " << colDip->iAcol << endl;
      break;
    }
    
    colDip = colDip->acolDips[0];
  }

  // Done.
}

// ------------------------------------------------------------------

// Find length of string neighbouring the given dipole assuming the dipole 
// was collapsed.

double ColourReconnection::neighbourLength(ColourDipole* dipNeigh, 
  vector<ColourDipole*>& dips, int iOld0, int iOld1, int iNew1, 
  int iOld2, int iOld3, int iNew2) {
 
  // First check that I did not already include dipole.
  for (int j = 0;j < int(dips.size()); ++j)
    if (dips[j] == dipNeigh)
      return 0;
  
  // Find correct indicies of neighbouring dipoles.
  int iCol = dipNeigh->iCol;
  int iAcol = dipNeigh->iAcol;
  if (iCol == iOld0 || iCol == iOld1)
    iCol = iNew1;
  if (iAcol == iOld0 || iAcol == iOld1)
    iAcol = iNew1;

  if (iCol == iOld2 || iCol == iOld2)
    iCol = iNew2;
  if (iAcol == iOld3 || iAcol == iOld3)
    iAcol = iNew2;

  // If neighbour is a simple string.
  if (!dipNeigh->isJun && !dipNeigh->isAntiJun) {
    dips.push_back(dipNeigh);
    return calculateStringLength(iCol, iAcol);
  }

  // If neighbour is a junction / anti junction.
  vector<int> iParticles;
  vector<bool> usedJuns(junctions.size(),false);
  int nJuns = 0;

  // Find particles in junction/anti-junction structures.
  if (dipNeigh->isJun) {
    if (!findJunctionParticles( -int(dipNeigh->iAcol/10) - 1, iParticles, 
      usedJuns, nJuns, dips)) return 1E9;
  } else {
    if (!findJunctionParticles( -int(dipNeigh->iCol/10) - 1, iParticles, 
      usedJuns, nJuns, dips)) return 1E9;
  }

  // Update list of particles to reflect the combination.
  for (int j = 0; j < int(iParticles.size()); ++j) {
    if (iParticles[j] == iOld0 || iParticles[j] == iOld1)
      iParticles[j] = iNew1;
    
    if (iParticles[j] == iOld2 || iParticles[j] == iOld3)
      iParticles[j] = iNew2;
  }
  
  // If it is a single junction.
  if (int(iParticles.size()) == 3) 
    return calculateJunctionLength(iParticles[0], iParticles[1],
      iParticles[2]);

  // If it is a junction pair.
  else if (int(iParticles.size()) == 4) 
    return calculateDoubleJunctionLength(iParticles[0], 
      iParticles[1], iParticles[2], iParticles[3]);
  
  // Something went wrong.
  else return 1E9;
}

// ------------------------------------------------------------------

bool ColourReconnection::getJunctionIndicies(ColourDipole * dip, int &iJun, 
  int &i0, int &i1, int &i2, int &junLeg0, int &junLeg1, int &junLeg2) {
  
  // Find junction index and first leg to combine.
  int indxJun = dip->iCol;
  if (dip->iAcol < 0)
      indxJun = dip->iAcol;
  iJun = (- indxJun / 10) - 1;
  junLeg0 = -(indxJun % 10);
  junLeg1 = 1; 
  junLeg2 = 2;
  if (junLeg0 == 1) junLeg1 = 0;
  else if (junLeg0 == 2) junLeg2 = 0;  
 
  if (dip->iCol < 0) {
    i0 = dip->iAcol;
    i1 = junctions[iJun].dips[junLeg1]->iAcol;
    i2 = junctions[iJun].dips[junLeg2]->iAcol;
  }
  else {
    i0 = dip->iCol;
    i1 = junctions[iJun].dips[junLeg1]->iCol;
    i2 = junctions[iJun].dips[junLeg2]->iCol;
  }
  
  // It is not possible to form a pseudoparticle if only a single particle is 
  // connected to the junction.
  if (i1 < 0 && i2 < 0) return false;
  
  // Check which two particle should form the pseudoparticle.
  double m1 = 1E9, m2 = 1E9;
  if (i1 >= 0)
    m1 = m(particles[i0].p(),particles[i1].p());
  if (i2 >= 0)
    m2 = m(particles[i0].p(),particles[i2].p());
  
  if (m1 > m2) {
    swap(i1,i2);
    swap(junLeg1,junLeg2);
  }
  // Force switch if i0 == i2
  if (i0 == i2) {
    swap(i1,i2);
    swap(junLeg1,junLeg2);
  }
  
  return true;  
}
  
// ------------------------------------------------------------------

// Make a test pseudo particle used to calculate string lengths.

int ColourReconnection::makeTestParticle( ColourDipole * dip) {
  
  // If simple dipole.
  if (dip->iCol >= 0 && dip->iAcol >= 0) {
    particles.push_back(particles[dip->iCol]);
    particles.back().p(particles[dip->iCol].p() + particles[dip->iAcol].p());
    return particles.size() -1;
  }

  // If junction structure.
  int i1 = -1,i2 = -1,i3 = -1;
  if(dip->iCol >= 0) {
    particles.push_back(particles[dip->iCol]);
    i1 = dip->iCol;
    
    // Find possible pseudocandidates.
    for (int i = 0;i < 3 ;i++) {
      if (junctions[- (dip->iAcol / 10 + 1)].dips[i] != dip) {
	if (i2 == -1)
	  i2 = junctions[- (dip->iAcol / 10 + 1)].dips[i]->iCol;
	else
	  i3 = junctions[- (dip->iAcol / 10 + 1)].dips[i]->iCol;
      }
    }
  } else {
    particles.push_back(particles[dip->iAcol]);
    i1 = dip->iAcol;
    
    // Find possible pseudocandidates.
    for (int i = 0;i < 3 ;i++) {
      if (junctions[- (dip->iCol / 10 + 1)].dips[i] != dip) {
	if (i2 == -1)
	  i2 = junctions[- (dip->iCol / 10 + 1)].dips[i]->iAcol;
	else
	  i3 = junctions[- (dip->iCol / 10 + 1)].dips[i]->iAcol;
      }
    }
  }
  // Calculate the invariant masses to find the minimum.
  double m1 = 1E9,m2 = 1E9;
  if (i2 >= 0)
    m1 = m(particles[i1].p(),particles[i2].p());
  if (i3 >= 0)
    m2 = m(particles[i1].p(),particles[i3].p());
  if(m1 < m2)
    particles.back().p(particles[i1].p() + particles[i2].p());
  else
    particles.back().p(particles[i1].p() + particles[i3].p());

  // Return the index of the new pseudo test particle.
  return particles.size() - 1;

}

// ------------------------------------------------------------------

// Calculate the invariant mass of a dipole.

double ColourReconnection::mDip(ColourDipole* dip) {

  // If double junction no invariant mass is given.
  if (dip->isJun && dip->isAntiJun) return 1E9;
  // If it has a single junction end.
  else if (dip->isJun || dip->isAntiJun) {
    int iJun, i0, i1, i2, junLeg0, junLeg1, junLeg2;
    getJunctionIndicies(dip, iJun, i0, i1, i2, junLeg0, junLeg1, junLeg2);
    if (i0 == i1)
      return particles[i0].m();
    if (i1 < 0)
      return 1E9;
    return m(particles[i0].p(),particles[i1].p());
  } // No junction ends.
  else {
    if(dip->iCol == dip->iAcol)
      return particles[dip->iCol].m();
    else
      return m(particles[dip->iCol].p(),particles[dip->iAcol].p());
  }

}

// ------------------------------------------------------------------

// Print dipoles, intended for debuggning purposes.

void ColourReconnection::listDipoles(bool onlyActive, bool onlyReal) {

  cout << " --- listing dipoles ---" << endl;
  for(int i = 0;i < int(dipoles.size()); ++i) {
    if(onlyActive && !dipoles[i]->isActive)
      continue;
    if(onlyReal && !dipoles[i]->isReal)
      continue;
    dipoles[i]->print();
  }
  cout << " --- finished listing ---" << endl;

}

// ------------------------------------------------------------------

// Print particles, intended for debugging purposes.

void ColourReconnection::listParticles() {

  for (int i = 0; i < int(particles.size()); ++i) {
    const ColourParticle& pt = particles[i];

    // Basic line for a particle, always printed.
    cout << setw(6) << i << setw(10) << pt.id() << "   " << left
       << setw(18) << pt.nameWithStatus(18) << right << setw(4)
       << pt.status() << setw(6) << pt.mother1() << setw(6)
       << pt.mother2() << setw(6) << pt.daughter1() << setw(6)
       << pt.daughter2() << setw(6) << pt.col() << setw(6) << pt.acol()
       << setprecision(3)
       << setw(11) << pt.px() << setw(11) << pt.py() << setw(11)
       << pt.pz() << setw(11) << pt.e() << setw(11) << pt.m();
    for(int j = 0;j < int(pt.activeDips.size());++j)
      cout << setw(10) << pt.activeDips[j];
    cout << "\n";
  }

}

// ------------------------------------------------------------------

// Print junctions, intended for debugging purposes.

void ColourReconnection::listJunctions() {

  cout << " --- listing junctions ---" << endl;
  for(int i = 0;i < int(junctions.size()); ++i)
    junctions[i].print();
  cout << " --- finished listing ---" << endl;

}

// ------------------------------------------------------------------

// Allow colour reconnections by mergings of MPI collision subsystems.
// iRec is system that may be reconnected, by moving its gluons to iSys,
// where minimal pT (or equivalently Lambda) is used to pick location.
// Therefore all dipoles in iSys have to be found, and all gluons in iRec.
// Matching q-qbar pairs are treated by analogy with gluons.
// Note: owing to rescatterings some outgoing partons must be skipped.

bool ColourReconnection::reconnectMPIs( Event&  event, int oldSize) {

  // References to beams to simplify indexing.
  BeamParticle& beamA = *beamAPtr;
  BeamParticle& beamB = *beamBPtr;

  // Prepare record of which systems should be merged onto another.
  // The iSys system must have colour in final state to attach to it.
  nSys = partonSystemsPtr->sizeSys();
  vector<int>  iMerge(nSys);
  vector<bool> hasColour(nSys);
  for (int iSys = 0; iSys < nSys; ++iSys) {
    iMerge[iSys] = iSys;
    bool hasCol = false;
    for (int iMem = 0; iMem < partonSystemsPtr->sizeOut(iSys); ++iMem) {
      int iNow = partonSystemsPtr->getOut( iSys, iMem);
      if (event[iNow].isFinal() && (event[iNow].col() > 0
        || event[iNow].acol() > 0) ) {
        hasCol = true;
        break;
      }
    }
    hasColour[iSys] = hasCol;
  }

  // Loop over systems to decide which should be reconnected.
  for (int iRec = nSys - 1; iRec > 0; --iRec) {

    // Determine reconnection strength from pT scale of system.
    double pT2Rec  = pow2( partonSystemsPtr->getPTHat(iRec) );
    double probRec = pT20Rec / (pT20Rec + pT2Rec);

    // Loop over other systems iSys at higher pT scale and
    // decide whether to reconnect the iRec gluons onto one of them.
    for (int iSys = iRec - 1; iSys >= 0; --iSys)
    if (hasColour[iSys] && probRec > rndmPtr->flat()) {

      // The iRec system and all merged with it to be merged with iSys.
      iMerge[iRec] = iSys;
      for (int iRec2 = iRec + 1; iRec2 < nSys; ++iRec2)
      if (iMerge[iRec2] == iRec) iMerge[iRec2] = iSys;

      // Once a system has been merged do not test it anymore.
      break;
    }
  }

  // Loop over systems. Check whether other systems to be merged with it.
  for (int iSys = 0; iSys < nSys; ++iSys) {
    int nMerge = 0;
    for (int iRec = iSys + 1; iRec < nSys; ++iRec)
    if (iMerge[iRec] == iSys) ++nMerge;
    if (nMerge == 0) continue;

    // Incoming partons not counted if rescattered.
    int  iInASys = partonSystemsPtr->getInA(iSys);
    bool hasInA  = (beamA[iSys].isFromBeam());
    int  iInBSys = partonSystemsPtr->getInB(iSys);
    bool hasInB  = (beamB[iSys].isFromBeam());

    // Begin find dipoles in iSys system.
    vector<BeamDipole> bmdipoles;
    int sizeOut = partonSystemsPtr->sizeOut(iSys);
    for (int iMem = 0; iMem < sizeOut; ++iMem) {

      // Find colour dipoles to beam remnant.
      int iNow = partonSystemsPtr->getOut( iSys, iMem);
      if (!event[iNow].isFinal()) continue;
      int col = event[iNow].col();
      if (col > 0) {
        if      (hasInA && event[iInASys].col() == col)
          bmdipoles.push_back( BeamDipole( col, iNow, iInASys ) );
        else if (hasInB && event[iInBSys].col() == col)
          bmdipoles.push_back( BeamDipole( col, iNow, iInBSys ) );
 
        // Find colour dipole between final-state partons.
        else for (int iMem2 = 0; iMem2 < sizeOut; ++iMem2)
        if (iMem2 != iMem) {
          int iNow2 = partonSystemsPtr->getOut( iSys, iMem2);
          if (!event[iNow2].isFinal()) continue;
          if (event[iNow2].acol() == col) {
            bmdipoles.push_back( BeamDipole( col, iNow, iNow2) );
            break;
          }
        }
      }

      // Find anticolour dipoles to beam remnant.
      int acol = event[iNow].acol();
      if (acol > 0) {
        if      (hasInA && event[iInASys].acol() == acol)
          bmdipoles.push_back( BeamDipole( acol, iInASys, iNow ) );
        else if (hasInB && event[iInBSys].acol() == acol)
          bmdipoles.push_back( BeamDipole( acol, iInBSys, iNow ) );
      }
    }
   
    // Skip mergings if no dipoles found.
    if (bmdipoles.size() == 0) continue;

    // Find dipole sizes.
    for (int iDip = 0; iDip < int(bmdipoles.size()); ++iDip)
      bmdipoles[iDip].p1p2 = event[bmdipoles[iDip].iCol].p()
                           * event[bmdipoles[iDip].iAcol].p();
    
    // Loop over systems iRec to be merged with iSys.
    for (int iRec = iSys + 1; iRec < nSys; ++iRec) {
      if (iMerge[iRec] != iSys) continue;

      // Information on iRec. Vectors for gluons and anything else.
      int sizeRec = partonSystemsPtr->sizeOut(iRec);
      int iInARec = partonSystemsPtr->getInA(iRec);
      int iInBRec = partonSystemsPtr->getInB(iRec);
      int nGluRec = 0;
      vector<int>    iGluRec;
      vector<double> pT2GluRec;
      int nAnyRec = 0;
      vector<int>    iAnyRec;
      vector<bool>   freeAnyRec;

      // Copy of gluon positions in descending order.
      for (int iMem = 0; iMem < sizeRec; ++iMem) {
        int iNow = partonSystemsPtr->getOut( iRec, iMem);
        if (!event[iNow].isFinal()) continue;
        if (event[iNow].isGluon()) {
          ++nGluRec;
          iGluRec.push_back( iNow );
          pT2GluRec.push_back( event[iNow].pT2() );
          for (int i = nGluRec - 1; i > 1; --i) {
            if (pT2GluRec[i - 1] > pT2GluRec[i]) break;
            swap(   iGluRec[i - 1],   iGluRec[i] );
            swap( pT2GluRec[i - 1], pT2GluRec[i] );
          }
        // Copy of anything else, mainly quarks, in no particular order.
        } else {
          ++nAnyRec;
          iAnyRec.push_back( iNow );
          freeAnyRec.push_back( true );
        }
      }

      // For each gluon in iRec now find the dipole that gives the smallest
      // (pGlu * pI) (pGlu * pJ) / (pI * pJ), i.e. minimal pT (and Lambda).
      for (int iGRec = 0; iGRec < nGluRec; ++iGRec) {
        int    iGlu      = iGluRec[iGRec];
        Vec4   pGlu      = event[iGlu].p();
        int    iDipMin   = 0;
        double pT2DipMin = sCM;
        for (int iDip = 0; iDip < int(bmdipoles.size()); ++iDip) {
          double pT2Dip = (pGlu * event[bmdipoles[iDip].iCol].p())
            * (pGlu * event[bmdipoles[iDip].iAcol].p()) / bmdipoles[iDip].p1p2;
          if (pT2Dip < pT2DipMin) {
            iDipMin   = iDip;
            pT2DipMin = pT2Dip;
          }
        }

        // Attach the gluon to the dipole, i.e. split the dipole in two.
        int colGlu   = event[iGlu].col();
        int acolGlu  = event[iGlu].acol();
        int colDip   = bmdipoles[iDipMin].col;
        int iColDip  = bmdipoles[iDipMin].iCol;
        int iAcolDip = bmdipoles[iDipMin].iAcol;
        event[iGlu].acol( colDip );
        if (event[iAcolDip].acol() == colDip)
             event[iAcolDip].acol( colGlu );
        else event[iAcolDip].col(  colGlu );
        bmdipoles[iDipMin].iAcol = iGlu;
        bmdipoles[iDipMin].p1p2 = event[iColDip].p() * pGlu;
        bmdipoles.push_back( BeamDipole( colGlu, iGlu, iAcolDip ) );
        bmdipoles.back().p1p2 = pGlu * event[iAcolDip].p();
     
        // Remove gluon from old system: reconnect colours.
        for (int i = oldSize; i < event.size(); ++i)
        if (i != iGlu && i != iAcolDip) {
          if (event[i].isFinal()) {
            if (event[i].acol() == colGlu) event[i].acol( acolGlu );
          } else {
              if (event[i].col()  == colGlu) event[i].col( acolGlu );
          }
        }

        // Update any junction legs that match reconnected dipole.
        for (int iJun = 0; iJun < event.sizeJunction(); ++iJun) {

          // Only junctions need to be updated, not antijunctions.
          if (event.kindJunction(iJun) % 2 == 0) continue;
          for (int leg = 0; leg < 3; ++leg) {
            int col = event.colJunction(iJun, leg);
            if (col == colDip)
              event.colJunction(iJun, leg, colGlu);
          }
        }
        
      }

      // See if any matching quark-antiquark pairs among the rest.
      for (int iQRec = 0; iQRec < nAnyRec; ++iQRec) {
        int iQ  = iAnyRec[iQRec];
        int idQ = event[iQ].id();
        if (freeAnyRec[iQRec] && idQ > 0 && idQ < 6)
        for (int iQbarRec = 0; iQbarRec < nAnyRec; ++iQbarRec) {
          int iQbar  = iAnyRec[iQbarRec];
          if (freeAnyRec[iQbarRec] && event[iQbar].id() == -idQ) {

            // Check that these can be traced back to same gluon splitting.
            // For now also avoid qqbar pairs produced in rescatterings.??
            int iTopQ    = event[iQ].iTopCopyId();
            int iTopQbar = event[iQbar].iTopCopyId();
            int iMother  = event[iTopQ].mother1();
            if (event[iTopQbar].mother1() == iMother
              && event[iMother].isGluon() && event[iMother].status() != -34
              && event[iMother + 1].status() != -34 ) {

              // Now find the dipole that gives the smallest
              // ((pQ + pQbar) * pI) ((pQ + pQbar) * pJ) / (pI * pJ).
              Vec4   pGlu      = event[iQ].p() + event[iQbar].p();
              int    iDipMin   = 0;
              double pT2DipMin = sCM;
              for (int iDip = 0; iDip < int(bmdipoles.size()); ++iDip) {
                double pT2Dip = (pGlu * event[bmdipoles[iDip].iCol].p())
                  * (pGlu * event[bmdipoles[iDip].iAcol].p())
                  / bmdipoles[iDip].p1p2;
                if (pT2Dip < pT2DipMin) {
                  iDipMin   = iDip;
                  pT2DipMin = pT2Dip;
                }
              }

              // Attach the q-qbar pair to the dipole, i.e. split the dipole.
              int colGlu   = event[iQ].col();
              int acolGlu  = event[iQbar].acol();
              int colDip   = bmdipoles[iDipMin].col;
              int iColDip  = bmdipoles[iDipMin].iCol;
              int iAcolDip = bmdipoles[iDipMin].iAcol;
              event[iQbar].acol( colDip );
              if (event[iAcolDip].acol() == colDip)
                   event[iAcolDip].acol( colGlu );
              else event[iAcolDip].col(  colGlu );
              bmdipoles[iDipMin].iAcol = iQbar;
              bmdipoles[iDipMin].p1p2 = event[iColDip].p() * event[iQbar].p();
              bmdipoles.push_back( BeamDipole( colGlu, iQ, iAcolDip ) );
              bmdipoles.back().p1p2 = event[iQ].p() * event[iAcolDip].p();
     
              // Remove q-qbar pair from old system: reconnect colours.
              freeAnyRec[iQRec]    = false;
              freeAnyRec[iQbarRec] = false;
              for (int i = oldSize; i < event.size(); ++i)
              if (i != iQRec && i != iQbarRec && i != iColDip
                && i != iAcolDip) {
                if (event[i].isFinal()) {
                  if (event[i].acol() == colGlu) event[i].acol( acolGlu );
                } else {
                    if (event[i].col()  == colGlu) event[i].col( acolGlu );
                }
              }
               
              // Update any junction legs that match reconnected dipole.
              for (int iJun = 0; iJun < event.sizeJunction(); ++iJun) {
                
                // Only junctions need to be updated, not antijunctions.
                if (event.kindJunction(iJun) % 2 == 0) continue;
                for (int leg = 0; leg < 3; ++leg) {
                  int col = event.colJunction(iJun, leg);
                  if (col == colDip)
                    event.colJunction(iJun, leg, colGlu);
                }
              }
              
              // Done with processing of q-qbar pairs.
            }
          }
        }
      }

      // If only two beam gluons left of system, set their colour = anticolour.
      // Used by BeamParticle::remnantColours to skip irrelevant gluons.
      if ( event[iInARec].isGluon() && !event[iInARec].isRescatteredIncoming()
        && event[iInBRec].isGluon() && !event[iInBRec].isRescatteredIncoming()
        && event[iInARec].col() == event[iInBRec].acol()
        && event[iInARec].acol() == event[iInBRec].col() ) {
          event[iInARec].acol( event[iInARec].col() );
          event[iInBRec].acol( event[iInBRec].col() );
      }

    // End of loops over iRec and iSys systems.
    }
  }

  // Done.
  return true;
}

//--------------------------------------------------------------------------

// Allow colour reconnections by moving gluons from their current location
// to another colour line. Also optionally flip two colour chains.

bool ColourReconnection::reconnectMove( Event&  event, int oldSize) {

  // Create or reset arrays to prepare for the new event analysis.
  vector<int> iGlu; 
  iReduceCol.resize( event.size() );
  iExpandCol.clear();
  map<int, int> colMap, acolMap;
  map<int, int>::iterator colM, acolM;
  vector<InfoGluonMove> infoGM; 

  // Temporary variables.
  int iNow, colNow, acolNow, iColNow, iAcolNow, col2Now, iCol2Now, iAcol2Now;
  double lambdaRefNow, dLambdaNow;

  // Loop over all final particles. Store (fraction of) gluons to move.
  for (int i = oldSize; i < event.size(); ++i) if (event[i].isFinal()) {
    if (event[i].id() == 21 && rndmPtr->flat() < fracGluon) 
      iGlu.push_back(i);

    // Store location of all colour and anticolour particles and indices.
    if (event[i].col() > 0 || event[i].acol() > 0) {
      iReduceCol[i] = iExpandCol.size();
      iExpandCol.push_back(i); 
      if (event[i].col() > 0) colMap[event[i].col()] = i;
      if (event[i].acol() > 0) acolMap[event[i].acol()] = i;
    }
  }

  // Erase (anti)colours for (anti)junctions and skip adjacent gluons. 
  for (int iJun = 0; iJun < event.sizeJunction(); ++iJun) {
    if (event.kindJunction(iJun) == 1) {
      for (int j = 0; j < 3; ++j) {
	int jCol = event.colJunction( iJun, j);
	for (colM = colMap.begin(); colM != colMap.end(); ++colM) 
	if (colM->first == jCol) {
	  colMap.erase( colM);
	  break;
	}
	for (int iG = 0; iG < int(iGlu.size()); ++iG)
	if (event[iGlu[iG]].col() == jCol) {
	  iGlu.erase(iGlu.begin() + iG);
	  break;
	} 
      }   
    } else if (event.kindJunction(iJun) == 2) {
      for (int j = 0; j < 3; ++j) {
	int jCol = event.colJunction( iJun, j);
	for (acolM = acolMap.begin(); acolM != acolMap.end(); ++acolM) 
	if (acolM->first == jCol) {
	  acolMap.erase( acolM);
	  break;
	}
	for (int iG = 0; iG < int(iGlu.size()); ++iG)
	if (event[iGlu[iG]].acol() == jCol) {
	  iGlu.erase(iGlu.begin() + iG);
	  break;
	} 
      }   
    }
  }

  // Error checks.
  int nGlu = iGlu.size(); 
  int nCol = colMap.size();   
  if (int(acolMap.size()) != nCol) {
    infoPtr->errorMsg("Error in MBReconUserHooks: map sizes do not match"); 
    return false;
  }
  colM  = colMap.begin();
  acolM = acolMap.begin();
  for (int iCol = 0; iCol < nCol; ++iCol) {
    if (colM->first != acolM->first) {
      infoPtr->errorMsg("Error in MBReconUserHooks: map elements"
	" do not match"); 
      return false;
    }
    ++colM;
    ++acolM; 
  }

  // Calculate and tabulate lambda between any pair of coloured partons. 
  nColMove = iExpandCol.size();   
  lambdaijMove.resize( pow2(nColMove) );
  for (int iAC = 0; iAC < nColMove - 1; ++iAC) {
    int i = iExpandCol[iAC];
    for (int jAC = iAC + 1; jAC < nColMove; ++jAC) {
      int j = iExpandCol[jAC];
      lambdaijMove[nColMove * iAC + jAC] 
	= log(1. + m2( event[i], event[j]) / m2Lambda);
    }
  }

  // Set up initial possible gluon moves with lambda gains/losses.
  for (int iG = 0; iG < nGlu; ++iG) {

    // Gluon and its current neighbours.
    iNow     = iGlu[iG];
    colNow   = event[iNow].col();
    acolNow  = event[iNow].acol();
    iColNow  = acolMap[colNow];
    iAcolNow = colMap[acolNow];

    // Addition to Lambda of gluon in current position.
    lambdaRefNow = lambda123Move( iNow, iColNow, iAcolNow);  

    // Loop over all colour lines where gluon could be inserted.
    for (colM = colMap.begin(); colM != colMap.end(); ++colM) {
      col2Now   = colM->first;  
      iCol2Now  = colMap[col2Now];
      iAcol2Now = acolMap[col2Now];

      // Addition to total Lambda if gluon moved to be inserted on line.
      dLambdaNow = (iCol2Now == iNow || iAcol2Now == iNow
	|| iColNow == iAcolNow) ? 2e4
	: lambda123Move( iNow, iCol2Now, iAcol2Now) - lambdaRefNow;

      // Add new container for gluon and colour line information.  
      infoGM.push_back( InfoGluonMove( iNow, colNow, acolNow, iColNow, 
	iAcolNow, col2Now, iCol2Now, iAcol2Now, lambdaRefNow, dLambdaNow ));
    }
  }
  int nPair = infoGM.size();

  // Keep on looping over moves until no further negative dLambda. 
  for ( int iMove = 0; iMove < nGlu; ++iMove) {
    int    iPairMin   = -1;
    double dLambdaMin = 1e4; 

    // Find lowest dLambda.
    for (int iPair = 0; iPair < nPair; ++iPair)
    if (infoGM[iPair].dLambda < dLambdaMin) {
      iPairMin   = iPair;
      dLambdaMin = infoGM[iPair].dLambda;
    } 

    // Break if no shift below upper limit found.
    if (dLambdaMin > -dLambdaCut) break;

    // Partons and colours involved in move.
    InfoGluonMove& selSM = infoGM[iPairMin]; 
    int i1Sel     = selSM.i1;
    int iCol1Sel  = selSM.iCol1; 
    int iAcol1Sel = selSM.iAcol1; 
    int iCol2Sel  = selSM.iCol2; 
    int iAcol2Sel = selSM.iAcol2; 
    int iCol2Mod[3]  = { iAcol1Sel , i1Sel     , iCol2Sel    };
    int col2Mod[3]   = { selSM.col1, selSM.col2, selSM.acol1};      

    // Remove gluon from old colour line and insert on new colour line.
    for (int i = 0; i < 3; ++i) {
      event[ iCol2Mod[i] ].col( col2Mod[i] );
      colMap[ col2Mod[i] ] = iCol2Mod[i];
    }      

    // Update info for partons with new colors.
    int  i1Now    = 0;
    bool doUpdate = false;
    for (int iPair = 0; iPair < nPair; ++iPair) {
      InfoGluonMove& tmpSM = infoGM[iPair];
      if (tmpSM.i1 != i1Now) {
	i1Now = tmpSM.i1;
	doUpdate = false;
	if (i1Now == i1Sel || i1Now == iCol1Sel || i1Now == iAcol1Sel 
	  || i1Now == iCol2Sel || i1Now == iAcol2Sel) {
	  colNow       = event[i1Now].col();
	  acolNow      = event[i1Now].acol();
	  iColNow      = acolMap[colNow];
	  iAcolNow     = colMap[acolNow];
	  lambdaRefNow = lambda123Move( i1Now, iColNow, iAcolNow);
	  doUpdate     = true;
	}
      }
      if (doUpdate) {
	tmpSM.col1      = colNow;
	tmpSM.acol1     = acolNow;
	tmpSM.iCol1     = iColNow;
	tmpSM.iAcol1    = iAcolNow;
	tmpSM.lambdaRef = lambdaRefNow;
      }
    }

    // Update info on dLambda for affected particles and colour lines.
    for (int iPair = 0; iPair < nPair; ++iPair) {
      InfoGluonMove& tmpSM = infoGM[iPair];
      int iMod = -1;
      for (int i = 0; i < 3; ++i) if (tmpSM.col2 == col2Mod[i]) iMod = i;
      if (iMod > -1) tmpSM.iCol2 = iCol2Mod[iMod];
      if (tmpSM.i1 == i1Sel || tmpSM.i1 == iCol1Sel || tmpSM.i1 == iAcol1Sel 
	|| tmpSM.i1 == iCol2Sel || tmpSM.i1 == iAcol2Sel || iMod > -1) 
	tmpSM.dLambda = (tmpSM.iCol2 == tmpSM.i1 || tmpSM.iAcol2 == tmpSM.i1
	  || tmpSM.iCol1 == tmpSM.iAcol1) ? 2e4
	  : lambda123Move( tmpSM.i1, tmpSM.iCol2, tmpSM.iAcol2)
	  - tmpSM.lambdaRef;
    }

  // End of loop over gluon shifting.
  }

  // Done if no flip.
  if (flipMode == 0) return true;

  // Array with colour lines, and where each line begins and ends.
  vector<int> iTmpFlip, iVecFlip, iBegFlip, iEndFlip; 

  // Variables for minimum search.
  int i1c, i1a, i2c, i2a, i1cMin, i1aMin, i2cMin, i2aMin, iSMin;
  double dLambdaFlip, dLambdaFlipMin;
  vector<InfoGluonMove> flipMin;

  // Grab all colour ends. 
  for (int i = oldSize; i < event.size(); ++i) 
  if (event[i].isFinal() && event[i].col() > 0 && event[i].acol() == 0) {
    iTmpFlip.clear();
    iTmpFlip.push_back( i);   

    // Step through colour neighbours to catch system. 
    iNow = i;
    acolM = acolMap.find( event[iNow].col() );
    bool foundEnd = false;
    while (acolM != acolMap.end()) {
      iNow = acolM->second;
      iTmpFlip.push_back( iNow);   
      if (event[iNow].col() == 0) {
	foundEnd = true;
	break;
      } 
      acolM = acolMap.find( event[iNow].col() );
    }  

    // Store acceptable system, optionally including junction legs.
    if (foundEnd || flipMode == 2) {
      iBegFlip.push_back( iVecFlip.size());
      for (int j = 0; j < int(iTmpFlip.size()); ++j) 
	iVecFlip.push_back( iTmpFlip[j]);
      iEndFlip.push_back( iVecFlip.size());
    }
  }

  // Optionally search for antijunction legs: grab all anticolour ends. 
  if (flipMode == 2) for (int i = oldSize; i < event.size(); ++i) 
  if (event[i].isFinal() && event[i].acol() > 0 && event[i].col() == 0) {
    iTmpFlip.clear();
    iTmpFlip.push_back( i);   

    // Step through anticolour neighbours to catch system. 
    iNow = i;
    colM = colMap.find( event[iNow].acol() );
    bool foundEnd = false;
    while (colM != colMap.end()) {
      iNow = colM->second;
      iTmpFlip.push_back( iNow);   
      if (event[iNow].acol() == 0) {
	foundEnd = true;
	break;
      } 
      colM = colMap.find( event[iNow].acol() );
    }  

    // Store acceptable system, but do not doublecount q - (n g) - qbar.
    if (!foundEnd) {
      iBegFlip.push_back( iVecFlip.size());
      for (int j = 0; j < int(iTmpFlip.size()); ++j) 
	iVecFlip.push_back( iTmpFlip[j]);
      iEndFlip.push_back( iVecFlip.size());
    }
  }

  // Loop through all system pairs.
  int nSysFlip = iBegFlip.size();
  for (int iSys1 = 0; iSys1 < nSysFlip - 1; ++iSys1) 
  if (iBegFlip[iSys1] >= 0) 
  for (int iSys2 = iSys1 + 1; iSys2 < nSysFlip; ++iSys2) 
  if (iBegFlip[iSys2] >= 0) {
    i1cMin     = 0;
    i1aMin     = 0;
    i2cMin     = 0;
    i2aMin     = 0; 
    dLambdaFlipMin = 1e4;

    // Loop through all possible flip locations for a pair.
    for (int j1 = iBegFlip[iSys1]; j1 < iEndFlip[iSys1] - 1; ++j1) 
    for (int j2 = iBegFlip[iSys2]; j2 < iEndFlip[iSys2] - 1; ++j2) {
      i1c = iVecFlip[j1];
      i1a = iVecFlip[j1 + 1];  
      i2c = iVecFlip[j2];
      i2a = iVecFlip[j2 + 1];  
      dLambdaFlip = lambda12Move( i1c, i2a) + lambda12Move( i2c, i1a)  
		  - lambda12Move( i1c, i1a) - lambda12Move( i2c, i2a);
      if (dLambdaFlip < dLambdaFlipMin) {
	i1cMin = i1c;
	i1aMin = i1a;
	i2cMin = i2c;
	i2aMin = i2a;
	dLambdaFlipMin = dLambdaFlip;
      }
    }

    // Store possible flips if low enough dLambdaMin.
    if (dLambdaFlipMin < -dLambdaCut) flipMin.push_back( InfoGluonMove( 
      iSys1, iSys2, i1cMin, i1aMin, i2cMin, i2aMin, dLambdaFlipMin) );
  } 
  int flipSize = flipMin.size();

  // Search for lowest possible flip among unused systems.
  for (int iFlip = 0; iFlip < min( nSysFlip / 2, flipSize); ++iFlip) { 
    iSMin = -1; 
    dLambdaFlipMin  = 1e4;
    for (int iSys12 = 0; iSys12 < flipSize; ++iSys12) 
    if (flipMin[iSys12].i1 >= 0 && flipMin[iSys12].dLambda < dLambdaFlipMin) {
      iSMin   = iSys12;
      dLambdaFlipMin = flipMin[iSys12].dLambda; 
    }

    // Do flip. Mark flipped systems.
    if (iSMin >= 0) {
      InfoGluonMove& flipNow = flipMin[iSMin];
      int iS1 = flipNow.i1;
      int iS2 = flipNow.i2;
      event[ flipNow.iAcol1 ].acol( event[flipNow.iCol2].col() );    
      event[ flipNow.iAcol2 ].acol( event[flipNow.iCol1].col() );
       for (int iSys12 = 0; iSys12 < flipSize; ++iSys12) 
      if ( flipMin[iSys12].i1 == iS1 || flipMin[iSys12].i1 == iS2
	|| flipMin[iSys12].i2 == iS1 || flipMin[iSys12].i2 == iS2)
	flipMin[iSys12].i1 = -1;
    }
    else break;
  }

  // Done.
  return true;

}

//==========================================================================

} // end namespace Pythia8
