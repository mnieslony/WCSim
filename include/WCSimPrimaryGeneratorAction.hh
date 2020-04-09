#ifndef WCSimPrimaryGeneratorAction_h
#define WCSimPrimaryGeneratorAction_h

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ThreeVector.hh"
#include "G4SPSAngDistribution.hh"
#include "G4SPSPosDistribution.hh"
#include "G4RandomDirection.hh"
#include "globals.hh"
#include "TTree.h"
#include "TChain.h"
#include "TRandom3.h"
#include "jhfNtuple.h"
#include <vector>
#include <fstream>

// GENIE headers
#ifndef NO_GENIE
#include "GHEP/GHepParticle.h"
#include "GHEP/GHepRecord.h"
#include "Ntuple/NtpMCTreeHeader.h"
#include "Ntuple/NtpMCEventRecord.h"
#include "Interaction/Interaction.h"
#endif
#include "WCSimRootOptions.hh"

class WCSimDetectorConstruction;
class G4ParticleGun;
class G4GeneralParticleSource;
class G4Event;
class WCSimPrimaryGeneratorMessenger;
class G4Generator;

class WCSimPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
  WCSimPrimaryGeneratorAction(WCSimDetectorConstruction*);
  ~WCSimPrimaryGeneratorAction();
  
public:
  void GeneratePrimaries(G4Event* anEvent);

  // Normal gun setting calls these functions to fill jhfNtuple and Root tree
  void SetVtx(G4ThreeVector i)     { vtxs[0] = i; nvtxs = 1; };
  void SetBeamEnergy(G4double i, G4int n = 0)   { beamenergies[n] = i;};
  void SetBeamDir(G4ThreeVector i, G4int n = 0) { beamdirs[n] = i;};
  void SetBeamPDG(G4int i, G4int n = 0)         { beampdgs[n] = i;};
  void SetNvtxs(G4int i)     { nvtxs = i; };
  void SetVtxs(G4int i, G4ThreeVector v)     { vtxs[i] = v; };

  G4double ShootEnergyPositronCustom();
  G4double ShootEnergyNeutron();
  std::string CalculateResNucleus(int resNuclA, int resNuclZ);
  void LoadDeexcitationProb();
  void GenerateDeexcitation(std::vector<int> *Talys_pdg, std::vector<G4ThreeVector> *Talys_momdir, std::vector<double> *Talys_energy, int resNuclA, int resNuclZ);

  // These go with jhfNtuple
  G4int GetVecRecNumber(){return vecRecNumber;}
  G4int GetMode() {return mode;};
  G4int GetNvtxs() {return nvtxs;};
  G4int GetVtxVol(G4int n = 0) {return vtxsvol[n];};
  G4ThreeVector GetVtx(G4int n = 0) {return vtxs[n];}
  G4int GetNpar() {return npar;};
  G4int GetBeamPDG(G4int n = 0) {return beampdgs[n];};
  G4double GetBeamEnergy(G4int n = 0) {return beamenergies[n];};
  G4ThreeVector GetBeamDir(G4int n = 0) {return beamdirs[n];};
  G4int GetTargetPDG(G4int n = 0) {return targetpdgs[n];};
  G4double GetTargetEnergy(G4int n = 0) {return targetenergies[n];};
  G4ThreeVector GetTargetDir(G4int n = 0) {return targetdirs[n];};
  // ANNIE: trace upstream sources
  G4String  GetDirtFileName(){return dirtFileName;}
  G4String GetGenieFileName(){return genieFileName;}
  G4int    GetDirtEntryNum(){return dirtEntryNum;}
  G4int   GetGenieEntryNum(){return genieEntryNum;}

  // older ...
  G4double GetNuEnergy() {return nuEnergy;};
  G4double GetEnergy() {return energy;};
  G4double GetXPos() {return xPos;};
  G4double GetYPos() {return yPos;};
  G4double GetZPos() {return zPos;};
  G4double GetXDir() {return xDir;};
  G4double GetYDir() {return yDir;};
  G4double GetZDir() {return zDir;};

  G4String GetGeneratorTypeString();
  
  void SaveOptionsToOutput(WCSimRootOptions * wcopt);

private:
  WCSimDetectorConstruction*      myDetector;
  G4ParticleGun*                  particleGun;
  G4GeneralParticleSource*        MyGPS;  //T. Akiri: GPS to run Laser
  WCSimPrimaryGeneratorMessenger* messenger;

  // Angular and positional generators (AntiNu generator)
  G4SPSAngDistribution *theSPSAng = nullptr;
  G4SPSPosDistribution *theSPSPos = nullptr;
  G4ThreeVector thePosition;
  G4ThreeVector theDirection;
  G4bool isFirstEvent;

  // Variables set by the messenger
  G4bool   useMulineEvt;
  G4bool   useGunEvt;
  G4bool   useLaserEvt;  //T. Akiri: Laser flag
  G4bool   useBeamEvt;
  G4bool   useGPSEvt;
  G4bool   useAntiNuEvt;
  G4bool   useGenieEvt;
  std::fstream inputFile;
  std::fstream inputSpecFile;
  G4String vectorFileName;
  G4String spectrumFileName;
  G4bool   GenerateVertexInRock;
  
  G4String dirtFileName;
  G4String genieFileName;
  G4int dirtEntryNum;
  G4int genieEntryNum;

  // These go with jhfNtuple
  G4int mode;
  G4int nvtxs;
  G4int vtxsvol[MAX_N_PRIMARIES];
  G4ThreeVector vtxs[MAX_N_PRIMARIES];
  G4int npar;
  G4int beampdgs[MAX_N_PRIMARIES], targetpdgs[MAX_N_PRIMARIES];
  G4ThreeVector beamdirs[MAX_N_PRIMARIES], targetdirs[MAX_N_PRIMARIES];
  G4double beamenergies[MAX_N_PRIMARIES], targetenergies[MAX_N_PRIMARIES];
  G4int vecRecNumber;

  G4double nuEnergy;
  G4double energy;
  G4double xPos, yPos, zPos;
  G4double xDir, yDir, zDir;

  G4int    _counterRock; 
  G4int    _counterCublic; 
  
  TChain* inputdata;
  TChain* metadata;
  TChain* geniedata;

  Long64_t localEntry;
  Int_t inputEntry;
  Int_t entriesInThisTree;
  Int_t treeNumber;
  TBranch* runBranch=0, *vtxxBranch=0, *vtxyBranch=0, *vtxzBranch=0, *vtxtBranch=0, *pxBranch=0, *pyBranch=0, *pzBranch=0, *EBranch=0, *KEBranch=0, *pdgBranch=0, *nTankBranch=0, *nupdgBranch=0, *nuvtxxBranch=0, *nuvtxyBranch=0, *nuvtxzBranch=0, *nuvtxtBranch=0, *nuPVBranch=0, *nuvtxmatBranch=0, *nuprimaryBranch=0, *nufluxfilenameBranch=0, *genieentryBranch=0, *genierecordBranch=0;
  Int_t runbranchval, entrybranchval, ntankbranchval, nupdgval, genieentrybranchval, pdgval, nuprimaryval;
  Double_t vtxxval, vtxyval, vtxzval, vtxtval, pxval, pyval, pzval, eval, keval, nuvtxxval, nuvtxyval, nuvtxzval, nuvtxtval;
  Int_t* pdgbranchval=0, *nuprimarybranchval=0;
  Double_t* vtxxbranchval=0, *vtxybranchval=0, *vtxzbranchval=0, *vtxtbranchval=0, *pxbranchval=0, *pybranchval=0, *pzbranchval=0, *ebranchval=0, *kebranchval=0;
  Char_t nupvval[100];
  Char_t numatval[100];
  Char_t nufluxfilenameval[100];
#ifndef NO_GENIE
  genie::NtpMCEventRecord* genierecordval;
#endif

  //GENIE gst tree variables
  TBranch* genie_entry_branch=0, *genie_neutrinoflav_branch=0, *genie_cc_branch=0, *genie_nc_branch=0, *genie_Z_branch=0, *genie_A_branch=0, *genie_hitnuc_branch=0, *genie_qel_branch=0, *genie_res_branch=0, *genie_dis_branch=0, *genie_coh_branch=0, *genie_imd_branch=0, *genie_El_branch=0, *genie_pxl_branch=0, *genie_pyl_branch=0, *genie_pzl_branch=0, *genie_fspl_branch=0, *genie_final_n_branch=0, *genie_final_p_branch=0, *genie_final_pip_branch=0, *genie_final_pim_branch=0, *genie_final_pi0_branch=0, *genie_final_kp_branch=0, *genie_final_km_branch=0, *genie_final_k0_branch=0, *genie_neutrinoE_branch=0, *genie_neutrinopx_branch=0, *genie_neutrinopy_branch=0, *genie_neutrinopz_branch=0, *genie_vtxx_branch=0, *genie_vtxy_branch=0, *genie_vtxz_branch=0, *genie_pdg_final_branch=0, *genie_E_final_branch=0, *genie_px_final_branch=0, *genie_py_final_branch=0, *genie_pz_final_branch=0, *genie_vtxt_branch=0, *genie_num_final_branch=0;
  Int_t genie_entry, genie_neutrinoflav, genie_Z, genie_A, genie_hitnuc, genie_final_n, genie_final_p, genie_final_pip, genie_final_pim, genie_final_pi0, genie_final_kp, genie_final_km, genie_final_k0, genie_num_final, genie_fspl;
  Bool_t genie_cc, genie_nc, genie_qel, genie_res, genie_dis, genie_coh, genie_imd;
  Double_t genie_neutrinoE, genie_neutrinopx, genie_neutrinopy, genie_neutrinopz, genie_vtxx, genie_vtxy, genie_vtxz, genie_vtxt, genie_El, genie_pxl, genie_pyl, genie_pzl;
  Int_t* genie_pdg_final;
  Double_t* genie_E_final, genie_px_final, genie_py_final, genie_pz_final;

  //TALYS tree variables
  TFile *f_O15 = 0, *f_N15 = 0,*f_N14 = 0,*f_Li9 = 0,*f_Li7 = 0,*f_C14 = 0,*f_C13 = 0,*f_C11 = 0,*f_C10 = 0,*f_Be10 = 0,*f_Be9 = 0, *f_B11 = 0,*f_B10 = 0,*f_B9 = 0;
  TTree *talys_O15 = 0, *talys_N15 = 0, *talys_N14 = 0, *talys_Li9 = 0, *talys_Li7 = 0, *talys_C14 = 0, *talys_C13 = 0, *talys_C11 = 0, *talys_C10 = 0, *talys_Be10 = 0, *talys_Be9 = 0, *talys_B11 = 0, *talys_B10 = 0, *talys_B9 = 0;
  Int_t talys_channel;
  std::vector<double> talys_gammaEnergy, talys_neutronEnergy, talys_protonEnergy, talys_deuteronEnergy, talys_tritiumEnergy, talys_heliumEnergy, talys_alphaEnergy;
  TBranch *branch_talys_channel = 0, *branch_talys_gammaE = 0, *branch_talys_neutronE = 0, *branch_talys_protonE = 0, *branch_talys_deuteronE = 0, *branch_talys_tritiumE = 0, *branch_talys_heliumE = 0, *branch_talys_alphaE = 0;  

  TFile *f_O15gamma = 0, *f_N15gamma = 0,*f_N14gamma = 0,*f_Li9gamma = 0,*f_Li7gamma = 0,*f_C14gamma = 0,*f_C13gamma = 0,*f_C11gamma = 0,*f_C10gamma = 0,*f_Be10gamma = 0,*f_Be9gamma = 0, *f_B11gamma = 0,*f_B10gamma = 0,*f_B9gamma = 0;
  TTree *talys_O15gamma = 0, *talys_N15gamma = 0, *talys_N14gamma = 0, *talys_Li9gamma = 0, *talys_Li7gamma = 0, *talys_C14gamma = 0, *talys_C13gamma = 0, *talys_C11gamma = 0, *talys_C10gamma = 0, *talys_Be10gamma = 0, *talys_Be9gamma = 0, *talys_B11gamma = 0, *talys_B10gamma = 0, *talys_B9gamma = 0;
  Int_t resnuclZ, resnuclA;
  std::vector<int> resnuclLevel;
  std::vector<double> resnuclEnergy, resnuclPop;
  TBranch *branch_resnuclZ = 0, *branch_resnuclA = 0, *branch_resnuclLevel = 0, *branch_resnuclEnergy = 0, *branch_resnuclPop = 0; 

  //TALYS additional variables
  double ExcitationProb = 0.25;	//probability to be in an excited state for O16 after nucleon was knocked out
  std::map<std::string,TTree*> talys_treemap;
  std::map<std::string,TTree*> talys_gammatreemap;
  TTree *talys_current = 0;
  TTree *talys_currentgamma = 0;

  std::vector<int> talys_pdg;       //additional de-excitation final state pdgs
  std::vector<G4ThreeVector> talys_momdir; //additional de-excitation final state momenta
  std::vector<double> talys_energy; //de-excitation final state energies

  std::vector<std::vector<int>> deex_part_pdg;
  std::vector<std::vector<double>> deex_part_energy;
  std::vector<int> deex_channel;
  std::vector<int> deex_resnuclZ;
  std::vector<int> deex_resnuclA;
  std::vector<std::vector<int>> deex_resnuclLevel;
  std::vector<std::vector<double>> deex_resnuclEnergy;
  std::vector<std::vector<double>> deex_resnuclPopulation;

  //Strings containing different file directories
  G4String primariesDirectory;
  G4String neutrinosDirectory;
  G4String genieDirectory;
  G4String talysDirectory;
  G4bool loadNewPrimaries;
  G4bool loadNewGenie;
  G4int primariesoffset;	

  // antinu read-in of the energy spectrum
  std::vector<G4double> Espectrum;
  std::vector<G4double> Espectrum_positron;
  std::vector<G4double> ProbabilitySpec;
	
public:

  inline void SetMulineEvtGenerator(G4bool choice) { useMulineEvt = choice; }
  inline G4bool IsUsingMulineEvtGenerator() { return useMulineEvt; }

  inline void SetGunEvtGenerator(G4bool choice) { useGunEvt = choice; }
  inline G4bool IsUsingGunEvtGenerator()  { return useGunEvt; }

  //T. Akiri: Addition of function for the laser flag
  inline void SetLaserEvtGenerator(G4bool choice) { useLaserEvt = choice; }
  inline G4bool IsUsingLaserEvtGenerator()  { return useLaserEvt; }
  
  inline void SetBeamEvtGenerator(G4bool choice) { useBeamEvt = choice; }
  inline G4bool IsUsingBeamEvtGenerator()  { return useBeamEvt; }

  inline void SetAntiNuEvtGenerator(G4bool choice) { useAntiNuEvt = choice; }
  inline G4bool IsUsingAntiNuEvtGenerator()  { return useAntiNuEvt; }
  
  inline void SetGenieEvtGenerator(G4bool choice) { useGenieEvt = choice; }
  inline G4bool IsUsingGenieEvtGenerator()  { return useGenieEvt; }
  
  inline void SetGPSEvtGenerator(G4bool choice) { useGPSEvt = choice; }
  inline G4bool IsUsingGPSEvtGenerator()  { return useGPSEvt; }

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

  inline void OpenSpectrumFile(G4String spectrumFile)
  {
    if ( inputSpecFile.is_open() )
      inputSpecFile.close();

    spectrumFileName = spectrumFile;
    inputSpecFile.open(spectrumFileName, std::fstream::in);

    if (!inputSpecFile.is_open() ) {
      G4cout << "Energy spectrum file "<< spectrumFileName << "not found" << G4endl;
      exit(-1);
    }
  }
 
  void SetPosition(G4ThreeVector position) {theSPSPos->SetCentreCoords(position);}
  void SetRadius(G4double radius) {theSPSPos->SetRadius(radius);}
  void SetHalfZ(G4double height) {theSPSPos->SetHalfZ(height);}
  void SetRot1(G4ThreeVector rot) {theSPSPos->SetPosRot1(rot);}
  void SetRot2(G4ThreeVector rot) {theSPSPos->SetPosRot2(rot);}

  inline void SetPrimaryFilesDirectory(G4String directoryName) { primariesDirectory = directoryName; }
  inline void SetNeutrinoFilesDirectory(G4String directoryName) { neutrinosDirectory = directoryName; }
  inline void SetGenieFilesDirectory(G4String directoryName) { genieDirectory = directoryName; }
  inline void SetTalysFilesDirectory(G4String directoryName) { talysDirectory = directoryName; }
  inline void SetNewPrimariesFlag(G4bool flagin){ loadNewPrimaries=flagin; }
  void LoadNewPrimaries();
  void LoadNewGENIEFile();
  void LoadTalysFiles();
  void SetPrimariesOffset(G4int offset){ primariesoffset=offset; }
  inline G4bool IsGeneratingVertexInRock() { return GenerateVertexInRock; }
  inline void SetGenerateVertexInRock(G4bool choice) { GenerateVertexInRock = choice; }
#ifndef NO_GENIE
  genie::NtpMCEventRecord* GetGenieRecord() { return genierecordval; }
#endif

};

#endif


