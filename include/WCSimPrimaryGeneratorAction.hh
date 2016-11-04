#ifndef WCSimPrimaryGeneratorAction_h
#define WCSimPrimaryGeneratorAction_h

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"
#include "TTree.h"
#include "TChain.h"
#include "TRandom3.h"
#include <fstream>

class WCSimDetectorConstruction;
class G4ParticleGun;
class G4GeneralParticleSource;
class G4Event;
class WCSimPrimaryGeneratorMessenger;

class WCSimPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
  WCSimPrimaryGeneratorAction(WCSimDetectorConstruction*, G4String infile);
  ~WCSimPrimaryGeneratorAction();
  
public:
  void GeneratePrimaries(G4Event* anEvent);

  // Normal gun setting calls these functions to fill jhfNtuple and Root tree
  void SetVtx(G4ThreeVector i)     { vtx = i; };
  void SetBeamEnergy(G4double i)   { beamenergy = i; };
  void SetBeamDir(G4ThreeVector i) { beamdir = i; };
  void SetBeamPDG(G4int i)         { beampdg = i; };

  // These go with jhfNtuple
  G4int GetVecRecNumber(){return vecRecNumber;}
  G4int GetMode() {return mode;};
  G4int GetVtxVol() {return vtxvol;};
  G4ThreeVector GetVtx() {return vtx;}
  G4int GetNpar() {return npar;};
  G4int GetBeamPDG() {return beampdg;};
  G4double GetBeamEnergy() {return beamenergy;};
  G4ThreeVector GetBeamDir() {return beamdir;};
  G4int GetTargetPDG() {return targetpdg;};
  G4double GetTargetEnergy() {return targetenergy;};
  G4ThreeVector GetTargetDir() {return targetdir;};

  // older ...
  G4double GetNuEnergy() {return nuEnergy;};
  G4double GetEnergy() {return energy;};
  G4double GetXPos() {return xPos;};
  G4double GetYPos() {return yPos;};
  G4double GetZPos() {return zPos;};
  G4double GetXDir() {return xDir;};
  G4double GetYDir() {return yDir;};
  G4double GetZDir() {return zDir;};

private:
  WCSimDetectorConstruction*      myDetector;
  G4ParticleGun*                  particleGun;
  G4GeneralParticleSource*        MyGPS;  //T. Akiri: GPS to run Laser
  WCSimPrimaryGeneratorMessenger* messenger;
  G4String                        fInputFileName;
  
  // Variables set by the messenger
  G4bool   useMulineEvt;
  G4bool   useNormalEvt;
  G4bool   useLaserEvt;  //T. Akiri: Laser flag
  G4bool   useBeamEvt;
  std::fstream inputFile;
  G4String vectorFileName;
  G4bool   GenerateVertexInRock;

  // These go with jhfNtuple
  G4int mode;
  G4int vtxvol;
  G4ThreeVector vtx;
  G4int npar;
  G4int beampdg, targetpdg;
  G4ThreeVector beamdir, targetdir;
  G4double beamenergy, targetenergy;
  G4int vecRecNumber;

  G4double nuEnergy;
  G4double energy;
  G4double xPos, yPos, zPos;
  G4double xDir, yDir, zDir;

  G4int    _counterRock; 
  G4int    _counterCublic; 
  
  TChain* inputdata;
	
	Int_t inputEntry;
	Int_t entriesInThisTree;
	Int_t treeNumber;
	TBranch* runBranch=0, *vtxxBranch=0, *vtxyBranch=0, *vtxzBranch=0, *vtxtBranch=0, *pxBranch=0, *pyBranch=0, *pzBranch=0, *EBranch=0, *KEBranch=0, *pdgBranch=0, *nTankBranch=0, *nupdgBranch=0, *nuvtxxBranch=0, *nuvtxyBranch=0, *nuvtxzBranch=0, *nuvtxtBranch=0, *nuPVBranch=0, *nuvtxmatBranch=0, *geniePrimaryBranch=0;
	Int_t runbranchval, entrybranchval, ntankbranchval, nupdgval;
	Int_t* pdgbranchval=0, *genieprimarybranchval=0;
	Int_t pdgval, genieprimaryval;
	Double_t *vtxxbranchval=0, *vtxybranchval=0, *vtxzbranchval=0, *vtxtbranchval=0, *pxbranchval=0, *pybranchval=0, *pzbranchval=0, *ebranchval=0, *kebranchval=0;
	Double_t vtxxval, vtxyval, vtxzval, vtxtval, pxval, pyval, pzval, eval, keval, nuvtxxval, nuvtxyval, nuvtxzval, nuvtxtval;
	Char_t nupvval[100];
	Char_t numatval[100];
	
public:

  inline void SetMulineEvtGenerator(G4bool choice) { useMulineEvt = choice; }
  inline G4bool IsUsingMulineEvtGenerator() { return useMulineEvt; }

  inline void SetNormalEvtGenerator(G4bool choice) { useNormalEvt = choice; }
  inline G4bool IsUsingNormalEvtGenerator()  { return useNormalEvt; }

  //T. Akiri: Addition of function for the laser flag
  inline void SetLaserEvtGenerator(G4bool choice) { useLaserEvt = choice; }
  inline G4bool IsUsingLaserEvtGenerator()  { return useLaserEvt; }
  
  inline void SetBeamEvtGenerator(G4bool choice) { useBeamEvt = choice; }
  inline G4bool IsUsingBeamEvtGenerator()  { return useBeamEvt; }

  inline void OpenVectorFile(G4String fileName) 
  {
    if ( inputFile.is_open() ) 
      inputFile.close();

    vectorFileName = fileName;
    inputFile.open(vectorFileName, std::fstream::in);

    if ( !inputFile.is_open() ) {
      G4cout << "Vector file " << vectorFileName << " not found" << G4endl;
      exit(-1);
    }
  }
  inline G4bool IsGeneratingVertexInRock() { return GenerateVertexInRock; }
  inline void SetGenerateVertexInRock(G4bool choice) { GenerateVertexInRock = choice; }

};

#endif


