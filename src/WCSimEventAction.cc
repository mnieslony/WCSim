#include "WCSimEventAction.hh"
#include "WCSimTrajectory.hh"
#include "WCSimRunAction.hh"
#include "WCSimPrimaryGeneratorAction.hh"
#include "WCSimWCHit.hh"
#include "WCSimWCDigi.hh"
#include "WCSimWCDigitizer.hh"
#include "WCSimWCTrigger.hh"
#include "WCSimWCAddDarkNoise.hh"
#include "WCSimWCPMT.hh"
#include "WCSimWCLAPPD.hh"
#include "WCSimDetectorConstruction.hh"

#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4EventManager.hh"
#include "G4UImanager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4VVisManager.hh"
#include "G4ios.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4TransportationManager.hh" 
#include "G4Navigator.hh" 
#include "G4SDManager.hh"
#include "G4DigiManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"

#include <set>
#include <iomanip>
#include <string>
#include <vector>
#include <climits>

#include "jhfNtuple.h"
#include "TTree.h"
#include "TFile.h"
#include "WCSimRootEvent.hh"
#include "TStopwatch.h"

// GENIE headers
#ifndef NO_GENIE
#include "GHEP/GHepParticle.h"
#include "Ntuple/NtpMCTreeHeader.h"
#include "Interaction/Interaction.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGUtils.h"
#include "PDG/PDGLibrary.h"
#endif
#include "NCVSD.hh"

#ifndef _SAVE_RAW_HITS
#define _SAVE_RAW_HITS
#ifndef _SAVE_RAW_HITS_VERBOSE
//#define _SAVE_RAW_HITS_VERBOSE
#endif
#endif
#ifndef SAVE_DIGITS_VERBOSE
//#define SAVE_DIGITS_VERBOSE
#endif
#ifndef TIME_DAQ_STEPS
//#define TIME_DAQ_STEPS
#endif

#ifndef NPMTS_VERBOSE
//#define NPMTS_VERBOSE 10
#endif

#ifndef HYPER_VERBOSITY
//#define HYPER_VERBOSITY
#endif

WCSimEventAction::WCSimEventAction(WCSimRunAction* myRun, 
				   WCSimDetectorConstruction* myDetector, 
				   WCSimPrimaryGeneratorAction* myGenerator)
  :runAction(myRun), generatorAction(myGenerator), 
   detectorConstructor(myDetector),
   ConstructedDAQClasses(false), LAPPDfile(0)
{
  DAQMessenger = new WCSimWCDAQMessenger(this);

  G4DigiManager* DMman = G4DigiManager::GetDMpointer();

  //create PMT response module
  WCSimWCPMT* WCDMPMT = new WCSimWCPMT( "WCReadoutPMT", myDetector, "tank");
  DMman->AddNewModule(WCDMPMT);

  //create dark noise module
  WCSimWCAddDarkNoise* WCDNM = new WCSimWCAddDarkNoise("WCDarkNoise", detectorConstructor, "tank");
  DMman->AddNewModule(WCDNM);

  // Repeat for MRD
#ifdef HYPER_VERBOSITY
  G4cout<<"WCSimEventAction::WCSimEventAction ☆ making new WCSimWCPMT for mrd with name WCReadoutPMT_MRD"<<G4endl;
#endif
  WCSimWCPMT* WCDMPMT_MRD = new WCSimWCPMT( "WCReadoutPMT_MRD", myDetector, "mrd");
  DMman->AddNewModule(WCDMPMT_MRD);
#ifdef HYPER_VERBOSITY
  G4cout<<"WCSimEventAction::WCSimEventAction ☆ making new WCSimWCAddDarkNoise for mrd with name WCDarkNoise_MRD"<<G4endl;
#endif
  WCSimWCAddDarkNoise* WCDNM_MRD = new WCSimWCAddDarkNoise("WCDarkNoise_MRD", detectorConstructor, "mrd");
  DMman->AddNewModule(WCDNM_MRD);
  // Repeat for FACC
  WCSimWCAddDarkNoise* WCDNM_FACC = new WCSimWCAddDarkNoise("WCDarkNoise_FACC", detectorConstructor, "facc");
  DMman->AddNewModule(WCDNM_FACC);
  WCSimWCPMT* WCDMPMT_FACC = new WCSimWCPMT( "WCReadoutPMT_FACC", myDetector, "facc");
  DMman->AddNewModule(WCDMPMT_FACC);
  
  WCSimWCLAPPD* WCDMLAPPD = new WCSimWCLAPPD( "WCReadoutLAPPD", myDetector);
  DMman->AddNewModule(WCDMLAPPD);
  
 //----------------------
 // ****************
 //  LAPPD TRUTH HITS
 // ****************
  lappdhit_x = new G4double[klappdhitnmax]; // where/when was the hit?
  lappdhit_y = new G4double[klappdhitnmax];
  lappdhit_z = new G4double[klappdhitnmax];
  lappdhit_t = new G4double[klappdhitnmax];
  lappdhit_process = new G4int[klappdhitnmax];    // what was the interaction process?
  lappdhit_particleID = new G4int[klappdhitnmax]; // what was the particle type interacting?
  lappdhit_trackID = new G4int[klappdhitnmax];    // what was the track ID
  lappdhit_edep = new G4double[klappdhitnmax];    // how much energy was deposited?
  lappdhit_objnum = new G4int[klappdhitnmax];     // which geometry object was hit?
  //lappdhit_copynum = new G4int[klappdhitnmax];    // which copy was hit?

}

WCSimEventAction::~WCSimEventAction()
{
  delete DAQMessenger;
  
  // Actions to cleanup LAPPD truth hits:
  // ====================================
   if(LAPPDfile){LAPPDfile->Close(); delete LAPPDfile; LAPPDfile=0;}
   
   delete[] lappdhit_x;
   delete[] lappdhit_y;
   delete[] lappdhit_z;
   delete[] lappdhit_t;
   delete[] lappdhit_process;
   delete[] lappdhit_particleID;
   delete[] lappdhit_trackID;
   delete[] lappdhit_edep;
   delete[] lappdhit_objnum;
}

void WCSimEventAction::CreateDAQInstances()
{
  if(ConstructedDAQClasses) {
    G4cerr << "WCSimEventAction::CreateDAQInstances() has already been called. Exiting..." << G4endl;
    return;
    exit(-1);
  }

  G4cout << "Creating digitizer and trigger class instances in WCSimEventAction::CreateDAQInstances()" << G4endl;

  G4DigiManager* DMman = G4DigiManager::GetDMpointer();

  //create your choice of digitizer module
  if(DigitizerChoice == "SKI") {
    WCSimWCDigitizerSKI* WCDM = new WCSimWCDigitizerSKI("WCReadoutDigits", detectorConstructor, DAQMessenger, "tank");
    DMman->AddNewModule(WCDM);
  }
  else {
    G4cerr << "Unknown DigitizerChoice " << DigitizerChoice << G4endl;
    exit(-1);
  }

  //create your choice of trigger module
  if(TriggerChoice == "NDigits") {
    WCSimWCTriggerNDigits* WCTM = new WCSimWCTriggerNDigits("WCReadout", detectorConstructor, DAQMessenger, "tank");
    DMman->AddNewModule(WCTM);
  }
  else if(TriggerChoice == "NDigits2") {
    WCSimWCTriggerNDigits2* WCTM = new WCSimWCTriggerNDigits2("WCReadout", detectorConstructor, DAQMessenger, "tank");
    DMman->AddNewModule(WCTM);
  }
  else if(TriggerChoice == "NoTrigger") {
    WCSimWCTriggerNoTrigger* WCTM = new WCSimWCTriggerNoTrigger("WCReadout", detectorConstructor, DAQMessenger, "tank");
    DMman->AddNewModule(WCTM);
  }
  else {
    G4cerr << "Unknown TriggerChoice " << TriggerChoice << G4endl;
    exit(-1);
  }

  // Repeat for MRD & FACC //
  if(DigitizerChoice=="SKI"){
    // Repeat for MRD
#ifdef HYPER_VERBOSITY
    G4cout<<"WCSimEventAction::CreateDAQInstances ☆ making new WCSimWCDigitizerSKI for mrd with name WCReadoutDigits_MRD"<<G4endl;
#endif
    WCSimWCDigitizerSKI* WCDM_MRD = new WCSimWCDigitizerSKI("WCReadoutDigits_MRD", detectorConstructor, DAQMessenger, "mrd");
    DMman->AddNewModule(WCDM_MRD);
    // repeat for FACC
    WCSimWCDigitizerSKI* WCDM_FACC = new WCSimWCDigitizerSKI("WCReadoutDigits_FACC", detectorConstructor, DAQMessenger, "facc");
    DMman->AddNewModule(WCDM_FACC);
  }
  // whatever the choice of tank trigger style is, MRD and FACC read out whatever digits are within tank trigger windows
#ifdef HYPER_VERBOSITY
  G4cout<<"WCSimEventAction::CreateDAQInstances ☆ making new WCSimWCTriggerOnTankDigits for mrd with name WCReadout_MRD"<<G4endl;
#endif
  WCSimWCTriggerOnTankDigits* WCTM_MRD = new WCSimWCTriggerOnTankDigits("WCReadout_MRD", detectorConstructor, DAQMessenger, "mrd");
  DMman->AddNewModule(WCTM_MRD);
  // repeat for facc
  WCSimWCTriggerOnTankDigits* WCTM_FACC = new WCSimWCTriggerOnTankDigits("WCReadout_FACC", detectorConstructor, DAQMessenger, "facc");
  DMman->AddNewModule(WCTM_FACC);

  ConstructedDAQClasses = true;
}


void WCSimEventAction::BeginOfEventAction(const G4Event* evt)
{
  if(!ConstructedDAQClasses)
    CreateDAQInstances();
  
  // equivalent of beginrunaction using event id
  if(evt->GetEventID()==0){ CreateNewLAPPDFile(); }
}

void WCSimEventAction::EndOfEventAction(const G4Event* evt)
{
  G4cout<<G4endl;
  G4cout<<"############# WCSIM  BEGIN END OF EVENT ACTION  ################"<<G4endl;
  // ----------------------------------------------------------------------
  //  Get Particle Table
  // ----------------------------------------------------------------------
  
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  
  // ----------------------------------------------------------------------
  //  Get Trajectory Container
  // ----------------------------------------------------------------------
  
  G4TrajectoryContainer* trajectoryContainer = evt->GetTrajectoryContainer();
  
  G4int n_trajectories = 0;
  if (trajectoryContainer) n_trajectories = trajectoryContainer->entries();
  
  // ----------------------------------------------------------------------
  //  Get Event Information
  // ----------------------------------------------------------------------
  
  G4int         event_id = evt->GetEventID();
  G4int         mode     = generatorAction->GetMode();

  G4int         nvtxs   = generatorAction->GetNvtxs();
  G4ThreeVector vtxs[MAX_N_PRIMARIES];
  G4int         vtxsvol[MAX_N_PRIMARIES];
  for( Int_t u=0; u<nvtxs; u++ ){
    vtxs[u]      = generatorAction->GetVtx(u);
    vtxsvol[u]   = WCSimEventFindStartingVolume(vtxs[u]);
  }
  G4int         vecRecNumber = generatorAction->GetVecRecNumber();
#ifndef NO_GENIE
  genie::NtpMCEventRecord* genierecordntpl = generatorAction->GetGenieRecord();
  
  if(genierecordntpl){
    genie::EventRecord* gevtRec = genierecordntpl->event;
    genie::Interaction* genieint = gevtRec->Summary();
    G4cout<<"This event was a "<<genieint->ProcInfo().AsString()
    //genieint->ScatteringTypeAsString()<<", "<<genieint->InteractionTypeAsString()
          <<" interaction of a "<<genieint->InitState().ProbeE(genie::kRfLab)<<"GeV "
          << genieint->InitState().Probe()->GetName() << " on a ";
    int nuc_pdgc = genieint->InitState().Tgt().HitNucPdg();
    if ( genie::pdg::IsNeutronOrProton(nuc_pdgc) ) {
    TParticlePDG * p = genie::PDGLibrary::Instance()->Find(nuc_pdgc);
    G4cout<< p->GetName();} else { G4cout<<"PDC-Code = " << nuc_pdgc;}
    G4cout<<" in ";
    TParticlePDG * tgt = genie::PDGLibrary::Instance()->Find( genieint->InitState().Tgt().Pdg() );
    if(tgt){G4cout<<tgt->GetName();} 
    else {G4cout<<"[Z="<<genieint->InitState().Tgt().Z()<<", A="<<genieint->InitState().Tgt().A()<<"]";}
    //<< ", PDG-Code = " << fTgt->Pdg();
    G4cout<<" producing a ";
    int ipos = gevtRec->RemnantNucleusPosition(); 
    if(ipos>-1){ G4cout<<gevtRec->Particle(ipos)->Energy()<<"GeV "<<gevtRec->Particle(ipos)->Name(); }
    G4cout<<" and a ";
    ipos = gevtRec->FinalStatePrimaryLeptonPosition();
    if(ipos>-1){ G4cout<<gevtRec->Particle(ipos)->Energy()<<"GeV "<<gevtRec->Particle(ipos)->Name(); }
    G4cout<<G4endl<<"Q^2 was "<<genieint->Kine().Q2()<<"XX, ";
    TLorentzVector& k1 = *(gevtRec->Probe()->P4());
    TLorentzVector& k2 = *(gevtRec->FinalStatePrimaryLepton()->P4());
    double costhl = TMath::Cos( k2.Vect().Angle(k1.Vect()) ); 
    G4cout<<"with final state lepton ejected at Cos(θ)="<<costhl<<G4endl;
    G4cout<<"Additional final state particles included "<<G4endl;
    G4cout<< " N(p) = "       << genieint->ExclTag().NProtons()
          << " N(n) = "       << genieint->ExclTag().NNeutrons()
          << G4endl
          << " N(pi^0) = "    << genieint->ExclTag().NPi0()
          << " N(pi^+) = "    << genieint->ExclTag().NPiPlus()
          << " N(pi^-) = "    << genieint->ExclTag().NPiMinus()
          <<G4endl;
    }
    
// use like:
// if (interaction.ProcInfo().IsQuasiElastic()) { ... }
// double Q2 = interaction.Kine().Q2();
// double Ev = interaction.InitState().ProbeE(kRfLab);
// int Z = interaction.InitState().Tgt().Z();
// double Ethr = interaction.PhaseSpace().Threshold();
    
    
/////////////////////////////
//    nupdg
//    nuvtxval (TLorentzVector)
//    nuvtx volumeName
//    nuvtx material
//		target energy
/////////////////////////////
#endif
    
  // ----------------------------------------------------------------------
  //  Get WC Hit Collection
  // ----------------------------------------------------------------------
    
  G4SDManager* SDman = G4SDManager::GetSDMpointer();

  // Get Hit collections of this event
  G4HCofThisEvent* HCE         = evt->GetHCofThisEvent();
  WCSimWCHitsCollection* WCHC = 0;

  G4String WCIDCollectionName = detectorConstructor->GetIDCollectionName();
  G4int collectionID;
  WCSimWCHitsCollection* WCHClappd = 0;
  G4String WCIDCollectionName2 = detectorConstructor->GetIDCollectionName2();
  G4int collectionID2;
  std::vector<int> objnumv;
  G4cout<<"-------------------------> I'm in event: event_id: "<<event_id<<G4endl;
  if (HCE)
  { 
    collectionID = SDman->GetCollectionID(WCIDCollectionName);
    WCHC = (WCSimWCHitsCollection*)HCE->GetHC(collectionID);
    //G4cout<<WCIDCollectionName<<" has "<<WCHC->entries()<<" entries"<<G4endl;
    collectionID2 = SDman->GetCollectionID(WCIDCollectionName2);
    WCHClappd = (WCSimWCHitsCollection*)HCE->GetHC(collectionID2);
    //G4cout<<"-----------name= "<<name<<" name2= "<<name2<<"--------"<<G4endl;
  } else {G4cout<<"Could not find hit colletion of event!"<<G4endl;}
  /** used if save raw hits is enabled in preprocessor - WCHC is a member variable used in FillRootEvent() */
  
  // To use Do like This:
  // --------------------
  //   if (WCHC)
  //     for (G4int i=0; i< WCHC->entries() ;i++)
  //       G4cout << (*WCHC)[i]->GetTotalPe() << G4endl;
  
  // ----------------------------------------------------------------------
  //  Get Digitized Hit Collection
  // ----------------------------------------------------------------------
  //G4cout<<"WCIDCollectionName= "<<WCIDCollectionName<<" WCIDCollectionName2= "<<WCIDCollectionName2<<G4endl;
  //G4cout<<"____________ WCHC->entries()= "<<int(WCHC->entries())<<G4endl;
  //G4cout<<"____________ WCHClappd->entries()= "<<int(WCHClappd->entries())<<G4endl;
  //if(WCHClappd->entries()>0.){ G4cout<<"GOTIT!!!!!!"<<G4endl; }
  for (G4int ii=0; ii< WCHClappd->entries() ;ii++){
    //G4cout<<"total pe @ LAPPDs: "<< (*WCHClappd)[ii]->GetTotalPe() << G4endl;
    G4int   lappdID         = (*WCHClappd)[ii]->GetTubeID();
    objnumv.push_back((*WCHClappd)[ii]->GetTubeID());
    G4float timelappd           = (*WCHClappd)[ii]->GetTime(ii);
    //G4cout<<"ii= "<<ii<<" @ "<<timelappd<<" lappdID "<<lappdID<<" @ "<<timelappd<<G4endl; 
  }
  
  //--------- STORE lappd HITS -------------
 //  //=========================================
  lappdevt=evt->GetEventID();
  G4double totalEnergy2 = 0.;
  G4int numberOfHits2 = WCHClappd->GetSize();
  G4cout << "&&&&  A total of " << numberOfHits2 << " hits on lappd were recorded!" << G4endl;
  if (WCHClappd!=0) {
    G4int totalpes_perevt = 0; G4int lappd_numhits0=0;
    for (G4int hitnum=0; hitnum<numberOfHits2; hitnum++) {
       G4cout<<"retrieving details for lappd hit "<<hitnum<<G4endl;
       lappd_numhits0++;
       WCSimWCHit* aHit = (*WCHClappd)[hitnum];
       G4ThreeVector hitPos = aHit->GetPos();
       hitPosx=hitPos.x();
       hitPosy=hitPos.y();
       hitPosz=hitPos.z();
       hitTime = aHit->GetTime(hitnum);
       //G4cout<<"hitTime= "<<hitTime<<G4endl;
       //hitParticleName = aHit->GetParticleName();
       hitTrackID = aHit->GetTrackID();
       for (G4int ip =0; ip < (*WCHClappd)[hitnum]->GetTotalPe(); ip++){
         G4cout<<"retrieving pe "<<ip<<" for lappd hit "<<"hitnum"<<G4endl;
       //hitPartCode = aHit->GetParticleID();
         lappdhit_process[hitnum] = hitProcessCode;
         double strip_coorx = ((*WCHClappd)[hitnum]->GetStripPosition(ip).x());
         double strip_coory = ((*WCHClappd)[hitnum]->GetStripPosition(ip).y());
         double strip_coorz = ((*WCHClappd)[hitnum]->GetStripPosition(ip).z());
         //G4cout<<"totalpes_perevt= "<<totalpes_perevt<<"--->GetStripPosition= "<<(*WCHClappd)[hitnum]->GetStripPosition(ip)<<" strip_coorx= "<<strip_coorx<<" strip_coory= "<<strip_coory<<G4endl;
         lappdhit_stripcoorx.push_back(strip_coorx);
         lappdhit_stripcoory.push_back(strip_coory);
         lappdhit_stripcoorz.push_back(strip_coorz);
         totalpes_perevt++;
       }
       //hitPartCode = ConvertParticleNameToCode(hitParticleName); // or use hitParticleID - from PDGEncoding?
       //hitProcessName = aHit->GetProcessName();
       //hitProcessCode = ConvertProcessNameToCode(hitProcessName);
       //if(hitProcessCode<0){G4cout<<"lappd hit process unaccounted"<<G4endl;}
       hitEdep = (*WCHClappd)[hitnum]->GetTotalPe();  //aHit->GetEdeposit();
       //hitCopyNum = aHit->GetCopyNum();
       //hitPhysical = (std::string)aHit->GetPhysical();
       // aHit->Print();
       lappdhit_x[hitnum] = hitPosx;
       lappdhit_y[hitnum] = hitPosy;
       lappdhit_z[hitnum] = hitPosz;
       lappdhit_t[hitnum] = hitTime;
       lappdhit_particleID[hitnum] = hitPartCode;
       lappdhit_trackID[hitnum] = hitTrackID;
       lappdhit_edep[hitnum] = hitEdep;
       //lappdhit_copynum[hitnum] = hitCopyNum;
       lappdhit_objnum[hitnum] = (objnumv[hitnum]);
       //G4cout<<"lappdhit_objnum[hitnum] = "<<lappdhit_objnum[hitnum]<<G4endl;
    }
    lappd_numhits = lappd_numhits0;
    lappdhit_totalpes_perevt = totalpes_perevt;
    //G4cout<<"lappd_numhits= "<<lappd_numhits<<" lappdhit_totalpes_perevt= "<<lappdhit_totalpes_perevt<<G4endl;
    //LAPPDtree ->Fill();
    /*for(int m=0; m<totalpes_perevt; m++){
      G4cout<<"tellme: "<<lappdhit_stripcoorx[m]<<","<<lappdhit_stripcoory[m]<<","<<lappdhit_stripcoorz[m]<<G4endl;
     }*/
  }
  G4cout<<"done processing lappd hits"<<G4endl;

  // Get a pointer to the Digitizing Module Manager
  G4DigiManager* DMman = G4DigiManager::GetDMpointer();

  // Get a pointer to the WC PMT module
  WCSimWCPMT* WCDMPMT =
    (WCSimWCPMT*)DMman->FindDigitizerModule("WCReadoutPMT");
 
  // new MFechner, aug 2006
  // need to clear up the old info inside PMT
  WCDMPMT->ReInitialize();
 
  WCSimWCLAPPD* WCDMLAPPD =
    (WCSimWCLAPPD*)DMman->FindDigitizerModule("WCReadoutLAPPD");
  G4cout<<"Reinitializing WCDLAPPD"<<G4endl;
  WCDMLAPPD->ReInitialize();
  
#ifdef TIME_DAQ_STEPS
  TStopwatch* ms = new TStopwatch();
  ms->Start();
#endif

  //Convert the hits to PMT pulse
  G4cout<<"Dititizing WCDMPMT"<<G4endl;
  WCDMPMT->Digitize();
  G4cout<<"Digizing WCDMLAPPD"<<G4endl;
  WCDMLAPPD->Digitize();
  
  // Do the Dark Noise, then Digitization, then Trigger
  // First, add Dark noise hits before digitizing
    
  //Get a pointer to the WC Dark Noise Module
  WCSimWCAddDarkNoise* WCDNM =
    (WCSimWCAddDarkNoise*)DMman->FindDigitizerModule("WCDarkNoise");
  
  //Add the dark noise
  G4cout<<"Adding dark noise with WCDNM"<<G4endl;
  WCDNM->AddDarkNoise();

  // Next, do the digitization
  
  //Get a pointer to the WC Digitizer Module
  WCSimWCDigitizerBase* WCDM =
    (WCSimWCDigitizerBase*)DMman->FindDigitizerModule("WCReadoutDigits");

  //Digitize the hits
  G4cout<<"Digitizing WCDM"<<G4endl;
  WCDM->Digitize();

  // Finally, apply the trigger
  
  //Get a pointer to the WC Trigger Module
  WCSimWCTriggerBase* WCTM =
    (WCSimWCTriggerBase*)DMman->FindDigitizerModule("WCReadout");
  
  //tell it the dark noise rate (for calculating the average dark occupancy -> can adjust the NDigits threshold)
  G4cout<<"Setting dark rate for WCDNM"<<G4endl;
  WCTM->SetDarkRate(WCDNM->GetDarkRate());
  
  //Apply the trigger
  // This takes the digits, and places them into trigger gates
  // Also throws away digits not contained in an trigger gate
  G4cout<<"Digitizing WCTM"<<G4endl;
  WCTM->Digitize();

#ifdef TIME_DAQ_STEPS
  ms->Stop();
  G4cout << " Digtization :  Real = " << ms->RealTime() 
    	    << " ; CPU = " << ms->CpuTime() << "\n";  
#endif

   // Get the post-noise hit collection for the WC
   G4int WCDChitsID = DMman->GetDigiCollectionID("WCRawPMTSignalCollection");
   WCSimWCDigitsCollection * WCDC_hits = (WCSimWCDigitsCollection*) DMman->GetDigiCollection(WCDChitsID);
  
   // Get the digitized collection for the WC
   G4int WCDCID = DMman->GetDigiCollectionID("WCDigitizedCollection");
   WCSimWCTriggeredDigitsCollection * WCDC = (WCSimWCTriggeredDigitsCollection*) DMman->GetDigiCollection(WCDCID);

   //----- for lappds -----
   G4int WCDChitsIDlappd = DMman->GetDigiCollectionID("WCRawLAPPDSignalCollection");
   WCSimWCDigitsCollection * WCDC_hitslappd = (WCSimWCDigitsCollection*) DMman->GetDigiCollection(WCDChitsIDlappd);
  
   // Get the digitized collection for the WC
   G4int WCDCIDlappd = DMman->GetDigiCollectionID("WCDigitizedCollectionLAPPD");
   WCSimWCTriggeredDigitsCollection * WCDClappd = (WCSimWCTriggeredDigitsCollection*) DMman->GetDigiCollection(WCDCIDlappd);
  //-------------

  //--------- LAPPDs -----------
  if (WCDC_hitslappd)
  {
    //add the truth raw hits
    // Both the pre- and post-PMT smearing hit times are accessible
    // Choose to save just the pre-smeared times for now
   #ifdef _SAVE_RAW_HITS_VERBOSE
    G4cout<<"RAW HITS"<<G4endl;
   #endif
   // wcsimrootevent->SetNumTubesHit(WCDC_hitslappd->entries());
  // std::vector<float> truetime2, smeartime2;
  // std::vector<int>   primaryParentID2;
   double hit_time_smear, hit_time_true;
   int hit_parentid; G4int id0=0;
   //loop over the DigitsCollection
   for(int idigi = 0; idigi < WCDC_hitslappd->entries(); idigi++) {
      G4cout<<"Adding LAPPD raw hit "<<idigi<<G4endl;
      int digi_tubeid = (*WCDC_hitslappd)[idigi]->GetTubeID();
      //G4cout<<"digi_lappdid= "<<digi_tubeid<<"/"<<(*WCDC_hitslappd)[idigi]->GetTotalPe()<<G4endl;
      lappdhit_totalpes_perlappd2.push_back((*WCDC_hitslappd)[idigi]->GetTotalPe());

      for(G4int id = 0; id < (*WCDC_hitslappd)[idigi]->GetTotalPe(); id++){
        G4cout<<"Adding pe "<<id<<"for LAPPD digi "<<idigi<<G4endl;
        id0++;
        hit_time_true  = (*WCDC_hitslappd)[idigi]->GetPreSmearTime(id);
        hit_parentid = (*WCDC_hitslappd)[idigi]->GetParentID(id);
        //G4cout<<"0___LAPPD idigi= "<<idigi<<" id= "<<id<<"/"<<(*WCDC_hitslappd)[idigi]->GetTotalPe()<<G4endl;
        //G4cout<<"id0= "<<id0<<" hit_time_true= "<<hit_time_true<<" hit_parentid= "<<hit_parentid<<G4endl;
        lappdhit_truetime2.push_back(hit_time_true);
        lappdhit_primaryParentID2.push_back(hit_parentid);
        ////---strip number and digitised hits-----
        int stripno = (*WCDC_hitslappd)[idigi]->GetStripNo(id);
        lappdhit_stripnum.push_back(stripno);
        //G4cout<<"LAPPD idigi= "<<idigi<<" id= "<<id<<" stripno= "<<stripno<<G4endl;
        std::map<int,double> stripno_peak=(*WCDC_hitslappd)[idigi]->GetNeighStripNo(id);
        int stiphit=0;
        G4cout<<"processing stripno_peaks"<<G4endl;
        for(std::map<int,double>::iterator m1=stripno_peak.begin(); m1!=stripno_peak.end(); ++m1){
           //G4cout<<"DIGImap____"<<(m1)->first<<","<<(m1)->second<<G4endl;
           stiphit++;
           lappdhit_neighstripnum.push_back( (m1)->first );
           lappdhit_neighstrippeak.push_back( (m1)->second );
        }
        stripno_peak.clear();
        //G4cout<<"We had "<<stiphit<<" neighbouring strips hit!"<<G4endl;
        lappdhit_NoOfneighstripsHit.push_back(stiphit);
       
        std::map<int,double> stripno_time=(*WCDC_hitslappd)[idigi]->GetNeighStripTime(id);
        G4cout<<"processing stripno_times"<<G4endl;
        for(std::map<int,double>::iterator m2=stripno_time.begin(); m2!=stripno_time.end(); ++m2){
          //G4cout<<"DIGItime____"<<(m2)->first<<","<<(m2)->second<<G4endl;
          lappdhit_neighstrip_time.push_back( (m2)->second );
        }
        stripno_time.clear();
        std::map<int,double> stripno_lefttime=(*WCDC_hitslappd)[idigi]->GetNeighStripLeftTime(id);
        G4cout<<"processing stripno_lefttimes"<<G4endl;
        for(std::map<int,double>::iterator m3=stripno_lefttime.begin(); m3!=stripno_lefttime.end(); ++m3){
          //G4cout<<"DIGIlefttime____"<<(m3)->first<<","<<(m3)->second<<G4endl;
          lappdhit_neighstrip_lefttime.push_back( (m3)->second );
        }
        stripno_lefttime.clear();
        std::map<int,double> stripno_righttime=(*WCDC_hitslappd)[idigi]->GetNeighStripRightTime(id);
        G4cout<<"processing stripno_righttimes"<<G4endl;
        for(std::map<int,double>::iterator m4=stripno_righttime.begin(); m4!=stripno_righttime.end(); ++m4){
          //G4cout<<"DIGIrighttime____"<<(m4)->first<<","<<(m4)->second<<G4endl;
          lappdhit_neighstrip_righttime.push_back( (m4)->second );
        }
        stripno_righttime.clear();
        //-------
        ////#ifdef _SAVE_RAW_HITS_VERBOSE
        hit_time_smear = (*WCDC_hitslappd)[idigi]->GetTime(id);
        lappdhit_smeartime2.push_back(hit_time_smear);
        //G4cout<<"hit_time_smear= "<<hit_time_smear<<G4endl;
       //#endif
      }//id
   #ifdef _SAVE_RAW_HITS_VERBOSE
      if(digi_tubeid < NPMTS_VERBOSE) {
        G4cout << "Adding " << truetime2.size()
               << " Cherenkov hits in tube " << digi_tubeid
               << " with truetime:smeartime:primaryparentID";
        for(G4int id = 0; id < truetime2.size(); id++) {
           G4cout << " " << truetime[id]
                  << ":" << smeartime[id]
                  << ":" << primaryParentID[id];
        }//id
       G4cout << G4endl;
      }
#endif
     
      /*smeartime2.clear();
      truetime2.clear();
      primaryParentID2.clear();*/
    }//idigi
   }
  G4cout<<"Filling LAPPDtree"<<G4endl;
  LAPPDtree->Fill();
  LAPPDfile->cd();
  G4cout<<"Writing LAPPD file"<<G4endl;
  LAPPDfile->Write("",TObject::kOverwrite);
  G4cout<<"clearing up"<<G4endl;
  lappdhit_NoOfneighstripsHit.clear();
  lappdhit_stripcoorx.clear();
  lappdhit_stripcoory.clear();
  lappdhit_stripcoorz.clear();
  lappdhit_totalpes_perlappd2.clear();
  objnumv.clear();
  lappdhit_smeartime2.clear();
  lappdhit_truetime2.clear();
  lappdhit_primaryParentID2.clear();
  lappdhit_stripnum.clear();
  lappdhit_neighstripnum.clear();
  lappdhit_neighstrippeak.clear();
  lappdhit_neighstrip_time.clear();
  lappdhit_neighstrip_lefttime.clear();
  lappdhit_neighstrip_righttime.clear();
  
   //___________________________
   /*   
   // To use Do like This:
   // --------------------
   if(WCDC) 
     for (G4int i=0; i < WCDC->entries(); i++) 
       {
	 G4int   tubeID         = (*WCDC)[i]->GetTubeID();
	 G4float photoElectrons = (*WCDC)[i]->GetPe(i);
	 G4float time           = (*WCDC)[i]->GetTime(i);
	 //	 G4cout << "time " << i << " " <<time << G4endl; 
	 //	 G4cout << "tubeID " << i << " " <<tubeID << G4endl; 
	 //	 G4cout << "Pe " << i << " " <<photoElectrons << G4endl; 
	 //   (*WCDC)[i]->Print();
       }
   */
   G4cout<<"moving to MRD & FACC"<<G4endl;
  // Repeat the steps for the MRD and FACC
  G4cout<<G4endl<<G4endl;
  G4String WCMRDCollectionName = detectorConstructor->GetMRDCollectionName();
  if(HCE){
  collectionID = SDman->GetCollectionID(WCMRDCollectionName);
  WCSimWCHitsCollection* WCHC_MRD = (WCSimWCHitsCollection*)HCE->GetHC(collectionID);
#ifdef HYPER_VERBOSITY
    G4cout<<"WCSimEventAction::EndOfEventAction ☆ (WCSimWCHitsCollection*)"<<WCMRDCollectionName
          <<" has "<<WCHC_MRD->entries()<<" entries"<<G4endl;
#endif
  }
  WCSimWCPMT* WCDMPMT_MRD = (WCSimWCPMT*)DMman->FindDigitizerModule("WCReadoutPMT_MRD");
  if(WCDMPMT_MRD==0){G4cout<<"WCReadoutPMT_MRD digitzer module not found!"<<G4endl;}
  WCDMPMT_MRD->ReInitialize();
#ifdef HYPER_VERBOSITY
  G4cout<<"WCSimEventAction::EndOfEventAction ☆ Calling Digitize on (WCSimWCPMT*)WCReadoutPMT_MRD"<<G4endl;
#endif
  WCDMPMT_MRD->Digitize();
  WCSimWCAddDarkNoise* WCDNM_MRD = (WCSimWCAddDarkNoise*)DMman->FindDigitizerModule("WCDarkNoise_MRD");
  if(WCDNM_MRD==0){G4cout<<"WCDarkNoise_MRD dark noise module not found!"<<G4endl;}
#ifdef HYPER_VERBOSITY
  G4cout<<"WCSimEventAction::EndOfEventAction ☆ Calling AddDarkNoise on (WCSimWCAddDarkNoise*)WCDarkNoise_MRD"<<G4endl;
#endif
  WCDNM_MRD->AddDarkNoise();
  WCSimWCDigitizerBase* WCDM_MRD = (WCSimWCDigitizerBase*)DMman->FindDigitizerModule("WCReadoutDigits_MRD");
  if(WCDM_MRD==0){G4cout<<"WCReadoutDigits_MRD digitizer module not found!"<<G4endl;}
#ifdef HYPER_VERBOSITY
  G4cout<<"WCSimEventAction::EndOfEventAction ☆ Calling Digitize on (WCSimWCDigitizerBase*)WCReadoutDigits_MRD"<<G4endl;
#endif
  WCDM_MRD->Digitize();
  WCSimWCTriggerBase* WCTM_MRD = (WCSimWCTriggerBase*)DMman->FindDigitizerModule("WCReadout_MRD");
  if(WCTM_MRD==0){G4cout<<"WCReadout_MRD trigger module not found!"<<G4endl;}
  WCTM_MRD->SetDarkRate(WCDNM_MRD->GetDarkRate());
#ifdef HYPER_VERBOSITY
  G4cout<<"WCSimEventAction::EndOfEventAction ☆ Calling Digitize on (WCSimWCTriggerBase*)WCReadout_MRD"<<G4endl;
#endif
  WCTM_MRD->Digitize();
  /** these are retrieved to pass to FillRootEvent() */
#ifdef HYPER_VERBOSITY
  G4cout<<"WCSimEventAction::EndOfEventAction ☆ retrieving raw hits (WCSimWCDigitsCollection*)WCRawPMTSignalCollection_MRD for FillRootEvent, which has ";
#endif
  G4int WCDChitsID_MRD = DMman->GetDigiCollectionID("WCRawPMTSignalCollection_MRD");
  WCSimWCDigitsCollection * WCDC_hits_MRD = (WCSimWCDigitsCollection*) DMman->GetDigiCollection(WCDChitsID_MRD);
#ifdef HYPER_VERBOSITY
  if(WCDC_hits_MRD){G4cout<<WCDC_hits_MRD->entries();} else {G4cout<<"no";} G4cout<<" entries"<<G4endl;
  G4cout<<"WCSimEventAction::EndOfEventAction ☆ retrieving readout hits (WCSimWCTriggeredDigitsCollection*)WCDigitizedCollection_MRD for FillRootEvent, which has ";
#endif
  G4int WCDCID_MRD = DMman->GetDigiCollectionID("WCDigitizedCollection_MRD");
  WCSimWCTriggeredDigitsCollection * WCDC_MRD = (WCSimWCTriggeredDigitsCollection*) DMman->GetDigiCollection(WCDCID_MRD);
#ifdef HYPER_VERBOSITY
  if(WCDC_hits_MRD){G4cout<<WCDC_MRD->entries();} else {G4cout<<"no";} G4cout<<" entries"<<G4endl;
#endif
  ///////////////////////////////
  // Repeat for FACC
  G4cout<<G4endl<<G4endl;
  G4String WCFACCCollectionName = detectorConstructor->GetFACCCollectionName();
  if(HCE){
  collectionID = SDman->GetCollectionID(WCFACCCollectionName);
  WCSimWCHitsCollection* WCHC_FACC = (WCSimWCHitsCollection*)HCE->GetHC(collectionID);
  }
  WCSimWCPMT* WCDMPMT_FACC = (WCSimWCPMT*)DMman->FindDigitizerModule("WCReadoutPMT_FACC");
  WCDMPMT_FACC->ReInitialize();
  WCDMPMT_FACC->Digitize();
  WCSimWCAddDarkNoise* WCDNM_FACC = (WCSimWCAddDarkNoise*)DMman->FindDigitizerModule("WCDarkNoise_FACC");
  WCDNM_FACC->AddDarkNoise();
  WCSimWCDigitizerBase* WCDM_FACC = (WCSimWCDigitizerBase*)DMman->FindDigitizerModule("WCReadoutDigits_FACC");
  WCDM_FACC->Digitize();
  WCSimWCTriggerBase* WCTM_FACC = (WCSimWCTriggerBase*)DMman->FindDigitizerModule("WCReadout_FACC");
  WCTM_FACC->SetDarkRate(WCDNM_FACC->GetDarkRate());
  WCTM_FACC->Digitize();
  G4int WCDChitsID_FACC = DMman->GetDigiCollectionID("WCRawPMTSignalCollection_FACC");
  WCSimWCDigitsCollection * WCDC_hits_FACC = (WCSimWCDigitsCollection*) DMman->GetDigiCollection(WCDChitsID_FACC);
  G4int WCDCID_FACC = DMman->GetDigiCollectionID("WCDigitizedCollection_FACC");
  WCSimWCTriggeredDigitsCollection * WCDC_FACC = (WCSimWCTriggeredDigitsCollection*) DMman->GetDigiCollection(WCDCID_FACC);
  /////////////////////////////////////////////////////

  // ----------------------------------------------------------------------
  //  Fill Ntuple
  // ----------------------------------------------------------------------

   jhfNtuple.mode   = mode;         // interaction mode
   jhfNtuple.nvtxs = nvtxs;       // number of vertices
   for( Int_t u=0; u<nvtxs; u++ ){
     jhfNtuple.vtxsvol[u] = vtxsvol[u];       // volume of vertex
     // unit mismatch between geant4 and reconstruction, M Fechner
     jhfNtuple.vtxs[u][0] = vtxs[u][0]/cm; // interaction vertex
     jhfNtuple.vtxs[u][1] = vtxs[u][1]/cm; // interaction vertex
     jhfNtuple.vtxs[u][2] = vtxs[u][2]/cm; // interaction vertex
   }
   jhfNtuple.vecRecNumber = vecRecNumber; //vectorfile record number
   
   // mustop, pstop, npar will be filled later
   
   // Next in the ntuple is an array of tracks.
   // We will keep count with npar
   
   G4int npar = 0;
   
   // First two tracks for each vertex are special: beam and target

   for( Int_t u=0; u<nvtxs; u++ ){
     /////////////////////////////////
     // npar = 0        NEUTRINO /////
     /////////////////////////////////
   
     G4int         beampdg;
     G4double      beamenergy;
     G4ThreeVector beamdir;
     
     beampdg    = generatorAction->GetBeamPDG(u);
     beamenergy = generatorAction->GetBeamEnergy(u);
     beamdir    = generatorAction->GetBeamDir(u);
  
     jhfNtuple.ipnu[npar]    = beampdg;               // id
     jhfNtuple.flag[npar]    = -1;                    // incoming neutrino
     jhfNtuple.m[npar]       = 0.0;                   // mass (always a neutrino)
     jhfNtuple.p[npar]       = beamenergy;            // momentum magnitude
     jhfNtuple.E[npar]       = beamenergy;            // energy 
     jhfNtuple.startvol[npar]= -1;                    // starting volume, vtxvol should be referred
     jhfNtuple.stopvol[npar] = -1;                    // stopping volume 
     jhfNtuple.dir[npar][0]  = beamdir[0];            // direction 
     jhfNtuple.dir[npar][1]  = beamdir[1];            // direction 
     jhfNtuple.dir[npar][2]  = beamdir[2];            // direction 
     jhfNtuple.pdir[npar][0] = beamenergy*beamdir[0]; // momentum-vector 
     jhfNtuple.pdir[npar][1] = beamenergy*beamdir[1]; // momentum-vector 
     jhfNtuple.pdir[npar][2] = beamenergy*beamdir[2]; // momentum-vector 
     // M Fechner, same as above
     jhfNtuple.stop[npar][0] = vtxs[u][0]/cm;  // stopping point (not meaningful)
     jhfNtuple.stop[npar][1] = vtxs[u][1]/cm;  // stopping point (not meaningful)
     jhfNtuple.stop[npar][2] = vtxs[u][2]/cm;  // stopping point (not meaningful)
     jhfNtuple.parent[npar] = 0;

     npar++;
     /////////////////////////////////
     // npar = 1        TARGET ///////
     /////////////////////////////////
     
     G4double      targetpmag = 0.0, targetmass = 0.0;
     G4int         targetpdg    = generatorAction->GetTargetPDG(u);
     G4double      targetenergy = generatorAction->GetTargetEnergy(u);
     G4ThreeVector targetdir    = generatorAction->GetTargetDir(u);

     if (targetpdg!=0) {            // protects against seg-fault
       if (targetpdg > 999)         // 16O nucleus not in pdg table
	 targetmass = targetenergy; // 16O is at rest, so E = m
       else
	 targetmass = particleTable->FindParticle(targetpdg)->GetPDGMass();
       if (targetenergy > targetmass) 
	 //      targetpmag = sqrt(targetenergy*targetenergy - targetmass*targetenergy);
	 // MF : bug fix
	 targetpmag = sqrt(targetenergy*targetenergy - targetmass*targetmass);
       else // protect against NaN
	 targetpmag = 0.0;
     }
     
     jhfNtuple.ipnu[npar]     = targetpdg;    // id
     jhfNtuple.flag[npar]    = -2;            // target
     jhfNtuple.m[npar]       = targetmass;    // mass (always a neutrino)
     jhfNtuple.p[npar]       = targetpmag;    // momentum magnitude
     jhfNtuple.E[npar]       = targetenergy;  // energy (total!) 
     jhfNtuple.startvol[npar] = -1;           // starting volume 
     jhfNtuple.stopvol[npar] = -1;            // stopping volume 
     jhfNtuple.dir[npar][0]  = targetdir[0];  // direction 
     jhfNtuple.dir[npar][1]  = targetdir[1];  // direction 
     jhfNtuple.dir[npar][2]  = targetdir[2];  // direction 
     // MF feb9,2006 : we want the momentum, not the energy...
     //  jhfNtuple.pdir[npar][0] = targetenergy*targetdir[0];  // momentum-vector 
     //  jhfNtuple.pdir[npar][1] = targetenergy*targetdir[1];  // momentum-vector 
     //  jhfNtuple.pdir[npar][2] = targetenergy*targetdir[2];  // momentum-vector 
     jhfNtuple.pdir[npar][0] = targetpmag*targetdir[0];  // momentum-vector 
     jhfNtuple.pdir[npar][1] = targetpmag*targetdir[1];  // momentum-vector 
     jhfNtuple.pdir[npar][2] = targetpmag*targetdir[2];  // momentum-vector 
     // M Fechner, same as above
     jhfNtuple.stop[npar][0] = vtxs[u][0]/cm;  // stopping point (not meaningful)
     jhfNtuple.stop[npar][1] = vtxs[u][1]/cm;  // stopping point (not meaningful)
     jhfNtuple.stop[npar][2] = vtxs[u][2]/cm;  // stopping point (not meaningful)
     jhfNtuple.parent[npar] = 0;

     npar++;
   }

  ////////////////////////
  // npar > nvertices  ///
  ////////////////////////

  // Draw Charged Tracks
  G4int PDG_e=11,PDG_v_e=12,PDG_gam=22;
  for( Int_t u=0; u<nvtxs; u++ ){
    G4int trkid_e=INT_MAX,trkid_v_e=INT_MAX,trkid_gam=INT_MAX;
    G4int idx_e=INT_MAX,idx_v_e=INT_MAX,idx_gam=INT_MAX;
    for (G4int i=0; i < n_trajectories; i++) 
      {
	WCSimTrajectory* trj = 
	  (WCSimTrajectory*)((*(evt->GetTrajectoryContainer()))[i]);
	if(u==0){       // only need to draw the tracks the first time round
	        if (trj->GetCharge() != 0.) trj->DrawTrajectory(50);
	}
	// now also search for any decay products of the primaries. 
	if(abs(trj->GetPDGEncoding()) == PDG_e && trj->GetParentID() == u+1 && trj->GetTrackID() < trkid_e) {
	  trkid_e = trj->GetTrackID();
	  idx_e = i;
	}
	if(abs(trj->GetPDGEncoding()) == PDG_v_e && trj->GetParentID() == u+1 && trj->GetTrackID() < trkid_v_e) {
	  trkid_v_e = trj->GetTrackID();
	  idx_v_e = i;
	}
	if(abs(trj->GetPDGEncoding()) == PDG_gam && trj->GetParentID() == u+1 && trj->GetTrackID() < trkid_gam) {
	  trkid_gam = trj->GetTrackID();
	  idx_gam = i;
	}
      }

    if(idx_e != INT_MAX) {
      WCSimTrajectory* trj =
	(WCSimTrajectory*)((*(evt->GetTrajectoryContainer()))[idx_e]);
      jhfNtuple.ipnu[npar]     = trj->GetPDGEncoding();    // id
      jhfNtuple.flag[npar]    = 0;            // target
      jhfNtuple.m[npar]       = particleTable->FindParticle(trj->GetPDGEncoding())->GetPDGMass();    // mass (always a neutrino)
      jhfNtuple.p[npar]       = trj->GetInitialMomentum().mag();    // momentum magnitude
      jhfNtuple.E[npar]       = sqrt(jhfNtuple.m[npar]*jhfNtuple.m[npar]+jhfNtuple.p[npar]*jhfNtuple.p[npar]);  // energy (total!) 
      jhfNtuple.startvol[npar] = -1;           // starting volume 
      jhfNtuple.stopvol[npar] = -1;            // stopping volume 
      jhfNtuple.dir[npar][0]  = trj->GetInitialMomentum().unit().getX();  // direction 
      jhfNtuple.dir[npar][1]  = trj->GetInitialMomentum().unit().getY();  // direction 
      jhfNtuple.dir[npar][2]  = trj->GetInitialMomentum().unit().getZ();  // direction 
      // MF feb9,2006 : we want the momentum, not the energy...
      //  jhfNtuple.pdir[npar][0] = targetenergy*targetdir[0];  // momentum-vector 
      //  jhfNtuple.pdir[npar][1] = targetenergy*targetdir[1];  // momentum-vector 
      //  jhfNtuple.pdir[npar][2] = targetenergy*targetdir[2];  // momentum-vector 
      jhfNtuple.pdir[npar][0] = trj->GetInitialMomentum().getX();  // momentum-vector 
      jhfNtuple.pdir[npar][1] = trj->GetInitialMomentum().getY();  // momentum-vector 
      jhfNtuple.pdir[npar][2] = trj->GetInitialMomentum().getZ();  // momentum-vector 
      // M Fechner, same as above
      jhfNtuple.stop[npar][0] = vtxs[u][0]/cm;  // stopping point (not meaningful)
      jhfNtuple.stop[npar][1] = vtxs[u][1]/cm;  // stopping point (not meaningful)
      jhfNtuple.stop[npar][2] = vtxs[u][2]/cm;  // stopping point (not meaningful)
      jhfNtuple.parent[npar] = 0;
      
      npar++; 
    }
    
    if(idx_v_e != INT_MAX) {
      WCSimTrajectory* trj =
	(WCSimTrajectory*)((*(evt->GetTrajectoryContainer()))[idx_v_e]);
      jhfNtuple.ipnu[npar]     = trj->GetPDGEncoding();    // id
      jhfNtuple.flag[npar]    = 0;            // target
      jhfNtuple.m[npar]       = particleTable->FindParticle(trj->GetPDGEncoding())->GetPDGMass();    // mass (always a neutrino)
      jhfNtuple.p[npar]       = trj->GetInitialMomentum().mag();    // momentum magnitude
      jhfNtuple.E[npar]       = sqrt(jhfNtuple.m[npar]*jhfNtuple.m[npar]+jhfNtuple.p[npar]*jhfNtuple.p[npar]);  // energy (total!) 
      jhfNtuple.startvol[npar] = -1;           // starting volume 
      jhfNtuple.stopvol[npar] = -1;            // stopping volume 
      jhfNtuple.dir[npar][0]  = trj->GetInitialMomentum().unit().getX();  // direction 
      jhfNtuple.dir[npar][1]  = trj->GetInitialMomentum().unit().getY();  // direction 
      jhfNtuple.dir[npar][2]  = trj->GetInitialMomentum().unit().getZ();  // direction 
      // MF feb9,2006 : we want the momentum, not the energy...
      //  jhfNtuple.pdir[npar][0] = targetenergy*targetdir[0];  // momentum-vector 
      //  jhfNtuple.pdir[npar][1] = targetenergy*targetdir[1];  // momentum-vector 
      //  jhfNtuple.pdir[npar][2] = targetenergy*targetdir[2];  // momentum-vector 
      jhfNtuple.pdir[npar][0] = trj->GetInitialMomentum().getX();  // momentum-vector 
      jhfNtuple.pdir[npar][1] = trj->GetInitialMomentum().getY();  // momentum-vector 
      jhfNtuple.pdir[npar][2] = trj->GetInitialMomentum().getZ();  // momentum-vector 
      // M Fechner, same as above
      jhfNtuple.stop[npar][0] = vtxs[u][0]/cm;  // stopping point (not meaningful)
      jhfNtuple.stop[npar][1] = vtxs[u][1]/cm;  // stopping point (not meaningful)
      jhfNtuple.stop[npar][2] = vtxs[u][2]/cm;  // stopping point (not meaningful)
      jhfNtuple.parent[npar] = 0;
      
      npar++; 
    }
    
    
    if(idx_gam != INT_MAX) {
      WCSimTrajectory* trj =
	(WCSimTrajectory*)((*(evt->GetTrajectoryContainer()))[idx_gam]);
      jhfNtuple.ipnu[npar]     = trj->GetPDGEncoding();    // id
      jhfNtuple.flag[npar]    = 0;            // target
      jhfNtuple.m[npar]       = particleTable->FindParticle(trj->GetPDGEncoding())->GetPDGMass();    // mass (always a neutrino)
      jhfNtuple.p[npar]       = trj->GetInitialMomentum().mag();    // momentum magnitude
      jhfNtuple.E[npar]       = sqrt(jhfNtuple.m[npar]*jhfNtuple.m[npar]+jhfNtuple.p[npar]*jhfNtuple.p[npar]);  // energy (total!) 
      jhfNtuple.startvol[npar] = -1;           // starting volume 
      jhfNtuple.stopvol[npar] = -1;            // stopping volume 
      jhfNtuple.dir[npar][0]  = trj->GetInitialMomentum().unit().getX();  // direction 
      jhfNtuple.dir[npar][1]  = trj->GetInitialMomentum().unit().getY();  // direction 
      jhfNtuple.dir[npar][2]  = trj->GetInitialMomentum().unit().getZ();  // direction 
      // MF feb9,2006 : we want the momentum, not the energy...
      //  jhfNtuple.pdir[npar][0] = targetenergy*targetdir[0];  // momentum-vector 
      //  jhfNtuple.pdir[npar][1] = targetenergy*targetdir[1];  // momentum-vector 
      //  jhfNtuple.pdir[npar][2] = targetenergy*targetdir[2];  // momentum-vector 
      jhfNtuple.pdir[npar][0] = trj->GetInitialMomentum().getX();  // momentum-vector 
      jhfNtuple.pdir[npar][1] = trj->GetInitialMomentum().getY();  // momentum-vector 
      jhfNtuple.pdir[npar][2] = trj->GetInitialMomentum().getZ();  // momentum-vector 
      // M Fechner, same as above
      jhfNtuple.stop[npar][0] = vtxs[u][0]/cm;  // stopping point (not meaningful)
      jhfNtuple.stop[npar][1] = vtxs[u][1]/cm;  // stopping point (not meaningful)
      jhfNtuple.stop[npar][2] = vtxs[u][2]/cm;  // stopping point (not meaningful)
      jhfNtuple.parent[npar] = 0;
      
      npar++; 
    }
  }

   //G4cout << "end________________ Filling Root Event: "<<event_id<<" lappdevt: "<<lappdevt<< G4endl;
   G4cout<<"Filling Root Event: "<<event_id<<G4endl;

   //   G4cout << "event_id: " << &event_id << G4endl;
   // G4cout << "jhfNtuple: " << &jhfNtuple << G4endl;
   //  G4cout << "WCHC: " << &WCHC << G4endl;
   //  G4cout << "WCDC: " << &WCDC << G4endl;
   //  G4cout << "WCFVHC: " << &WCFVHC << G4endl;
   //  G4cout << "WCFVDC: " << &WCFVDC << G4endl;
   // G4cout << "lArDHC: " << &lArDHC << G4endl;
   // G4cout << "FGDxHC: " << &FGDxHC << G4endl;
   // G4cout << "FGDyHC: " << &FGDyHC << G4endl;
   // G4cout << "MRDxHC: " << &MRDxHC << G4endl;
   // G4cout << "MRDyHC: " << &MRDyHC << G4endl;
   
  jhfNtuple.npar = npar;
  
  FillRootEvent(event_id,
		jhfNtuple,
		trajectoryContainer,
		WCDC_hits,
		WCDC,
		"tank");
G4cout<<"Filling MRD Root Event"<<G4endl;
  FillRootEvent(event_id,
		jhfNtuple,
		trajectoryContainer,
		WCDC_hits_MRD,
		WCDC_MRD,
		"mrd");
G4cout<<"Filling FACC Root Event"<<G4endl;
  FillRootEvent(event_id,
		jhfNtuple,
		trajectoryContainer,
		WCDC_hits_FACC,
		WCDC_FACC,
		"facc");
  
  TTree* tree = GetRunAction()->GetTree();
  TBranch* tankeventbranch = tree->GetBranch("wcsimrootevent");
  tree->SetEntries(tankeventbranch->GetEntries());
  //tree->SetEntries(GetRunAction()->GetNumberOfEventsGenerated());
  TFile* hfile = tree->GetCurrentFile();
  hfile->cd();
  // MF : overwrite the trees -- otherwise we have as many copies of the tree
  // as we have events. All the intermediate copies are incomplete, only the
  // last one is useful --> huge waste of disk space.
  hfile->Write("",TObject::kOverwrite);
  
  G4cout<<"events generated so far: "<<(GetRunAction()->GetNumberOfEventsGenerated())<<G4endl;
  if(event_id%1000==0&&event_id!=0){
    GetRunAction()->CreateNewOutputFile();
    CreateNewLAPPDFile(); // this must always come *after* the runaction version
  }
  
  G4cout<<"############# WCSIM FINISH END OF EVENT ACTION  ################"<<G4endl;

}
//TODO: Starting and Stopping Volume finders also need to be modified to add MRD and FACC volumes
G4int WCSimEventAction::WCSimEventFindStartingVolume(G4ThreeVector vtx)
{
  // Get volume of starting point (see GEANT4 FAQ)

  G4int vtxvol = -1;

  G4Navigator* tmpNavigator = 
    G4TransportationManager::GetTransportationManager()->
    GetNavigatorForTracking();

  G4VPhysicalVolume* tmpVolume = tmpNavigator->LocateGlobalPointAndSetup(vtx);
  G4String       vtxVolumeName = tmpVolume!=0 ? tmpVolume->GetName() : "";



  if ( vtxVolumeName == "outerTube" ||
	    vtxVolumeName == "innerTube" ||
	    vtxVolumeName == "rearEndCap"|| 
	    vtxVolumeName == "frontEndCap" )
    vtxvol = 10;

  else if ( vtxVolumeName.contains("WC") && !vtxVolumeName.contains("FV") )
  {
// aah original line  ->if (vtxVolumeName.contains("WCBarrel"))
	  if ((vtxVolumeName.contains("WCBarrel"))|| (vtxVolumeName.contains("Tank")))	//aah I added "Tank" as MB equivalent of Barrel
      vtxvol = 10;
    else if (vtxVolumeName == "WCBox")
      vtxvol = -2;
    else if (vtxVolumeName.contains("PMT") ||
	     vtxVolumeName.contains("LAPPD") ||
	     vtxVolumeName.contains("Cap") ||
	     vtxVolumeName.contains("Cell"))
      vtxvol = 11;
    else if (vtxVolumeName.contains("OD"))
      vtxvol = 12;
    else
    {
      G4cout << vtxVolumeName << " unkown vtxVolumeName " << G4endl;
      vtxvol = -3;
    }
  }
  else if ( vtxVolumeName == "expHall" )
    vtxvol = 0;
  else if ( vtxVolumeName == "catcher" )
    vtxvol = 40;
  
  if(vtxvol<0){
    //G4cout<<"############# unkown vertex volume: "<<vtxVolumeName<<" ################"<<G4endl;
  }
  return vtxvol;
}

G4int WCSimEventAction::WCSimEventFindStoppingVolume(G4String stopVolumeName)
{
  G4int stopvol = -1;

  if ( stopVolumeName.contains("WC") && !stopVolumeName.contains("FV") )
  {
//aah Original line->    if (stopVolumeName.contains("WCBarrel"))
	  if ((stopVolumeName.contains("WCBarrel"))|| (stopVolumeName.contains("Tank")))	// aah I added "Tank" as MB equivalent of Barrel
      stopvol = 10;
    else if (stopVolumeName == "WCBox")
      stopvol = -2;
    else if (stopVolumeName.contains("PMT") ||
	     stopVolumeName.contains("LAPPD") ||
	     stopVolumeName.contains("Cap") ||
	     stopVolumeName.contains("Cell"))
      stopvol = 11;
    else if (stopVolumeName.contains("OD"))
      stopvol = 12;
    else
    {
      G4cout << stopVolumeName << " unkown stopVolumeName " << G4endl;
      stopvol = -3;
    }
  }

  else if ( stopVolumeName.contains("FV") )
  {
    if (stopVolumeName == "WCFVBarrel" ||
	stopVolumeName == "WCFVAnnulus" ||
	stopVolumeName == "WCFVRing" )
      stopvol = 10;
    else if (stopVolumeName.contains("FVPMT"))
      stopvol = 13;
    else
    {
      G4cout << stopVolumeName << " unkown stopVolumeName " << G4endl;
      stopvol = -3;
    }
  }
  else if ( stopVolumeName == "expHall" )
    stopvol = 0;
  else if ( stopVolumeName == "catcher" )
    stopvol = 40;

  if(stopvol<0){
    //G4cout<<"############# unkown vertex volume: "<<stopVolumeName<<" ################"<<G4endl;
  }
  return stopvol;
}

void WCSimEventAction::FillRootEvent(G4int event_id, 
				     const struct ntupleStruct& jhfNtuple,
				     G4TrajectoryContainer* TC,
				     WCSimWCDigitsCollection* WCDC_hits,
				     WCSimWCTriggeredDigitsCollection* WCDC,
				     G4String detectorElement
				     )
{
  // Fill up a Root event with stuff from the ntuple

  WCSimRootEvent* wcsimrootsuperevent = GetRunAction()->GetRootEvent(detectorElement);

  // start with the first "sub-event"
  // if the WC digitization requires it, we will add another subevent
  // for the WC.
  // all the rest goes into the first "sub-event".
  WCSimRootTrigger* wcsimrootevent = wcsimrootsuperevent->GetTrigger(0);
  G4DigiManager* DMman = G4DigiManager::GetDMpointer();
  WCSimWCTriggerBase* WCTM;
  
  if(detectorElement=="tank"){
    WCTM = (WCSimWCTriggerBase*)DMman->FindDigitizerModule("WCReadout");
  } else if(detectorElement=="mrd"){
    WCTM = (WCSimWCTriggerBase*)DMman->FindDigitizerModule("WCReadout_MRD");
  } else if(detectorElement=="facc"){
    WCTM = (WCSimWCTriggerBase*)DMman->FindDigitizerModule("WCReadout_FACC");
  }
  
  // get number of gates
  int ngates = WCTM->NumberOfGatesInThisEvent(); 

  G4cout << "ngates "<<detectorElement<<" =  " << ngates << "\n";
  for (int index = 0 ; index < ngates ; index++) 
    {
      if (index >=1 ) {
	wcsimrootsuperevent->AddSubEvent();
	wcsimrootevent = wcsimrootsuperevent->GetTrigger(index);
	wcsimrootevent->SetHeader(event_id,0,
				   0,index+1); // date & # of subevent 
	wcsimrootevent->SetMode(jhfNtuple.mode);
      }
      wcsimrootevent->SetTriggerInfo(WCTM->GetTriggerType(index),
				     WCTM->GetTriggerInfo(index));
    }
  

  // Fill the header
  // Need to add run and date
  wcsimrootevent = wcsimrootsuperevent->GetTrigger(0);
  wcsimrootevent->SetHeader(event_id,0,0); // will be set later.

  // Fill other info for this event

  wcsimrootevent->SetMode(jhfNtuple.mode);

  wcsimrootevent->SetNvtxs(jhfNtuple.nvtxs);
  for( Int_t u=0; u<jhfNtuple.nvtxs; u++ ){
    wcsimrootevent->SetVtxsvol(u,jhfNtuple.vtxsvol[u]);
    for (int j=0;j<3;j++)
      {
	wcsimrootevent->SetVtxs(u,j,jhfNtuple.vtxs[u][j]);
      }
  }
  wcsimrootevent->SetJmu(jhfNtuple.jmu);
  wcsimrootevent->SetJp(jhfNtuple.jp);
  wcsimrootevent->SetNpar(jhfNtuple.npar);
  wcsimrootevent->SetVecRecNumber(jhfNtuple.vecRecNumber);

  // Add the tracks with the particle information
  // First two tracks come from jhfNtuple, as they are special

  int k;
  //Modify to add decay products
  for (k=0;k<jhfNtuple.npar;k++) // should be just 2
  {
    float dir[3];
    float pdir[3];
    float stop[3];
    float start[3];
    for (int l=0;l<3;l++)
    {
      dir[l]=jhfNtuple.dir[k][l];
      pdir[l]=jhfNtuple.pdir[k][l];
      stop[l]=jhfNtuple.stop[k][l];
      start[l]=jhfNtuple.start[k][l];
	//G4cout<< "start[" << k << "][" << l <<"]: "<< jhfNtuple.start[k][l] <<G4endl;
    }

    // Add the track to the TClonesArray
    wcsimrootevent->AddTrack(jhfNtuple.ipnu[k], 
			      jhfNtuple.flag[k], 
			      jhfNtuple.m[k], 
			      jhfNtuple.p[k], 
			      jhfNtuple.E[k], 
			      jhfNtuple.startvol[k], 
			      jhfNtuple.stopvol[k], 
			      dir, 
			      pdir, 
			      stop,
			      start,
			      jhfNtuple.parent[k],
			     jhfNtuple.time[k],0); 
  }

  // the rest of the tracks come from WCSimTrajectory

  std::set<int> pizeroList;
  // added by M Fechner, dec 16th, 2004
  std::set<int> muonList;
  std::set<int> antimuonList;
  // same, april 7th 2005
  std::set<int> pionList;
  std::set<int> antipionList;
  std::set<int> primaryList;

  // Pi0 specific variables
  Float_t pi0Vtx[3];
  Int_t   gammaID[2];
  Float_t gammaE[2];
  Float_t gammaVtx[2][3];
  Int_t   r = 0;

  G4int n_trajectories = 0;
  if (TC)
    n_trajectories = TC->entries();

  // M Fechner : removed this limit to get to the primaries...
  //if (n_trajectories>50)  // there is no need for this limit, but it has
  //n_trajectories=50;    // existed in previous versions of the code.  It also
                          // makes the ROOT file smaller.  

  for (int i=0; i <n_trajectories; i++) 
  {
    WCSimTrajectory* trj = (WCSimTrajectory*)(*TC)[i];

    // If this track is a pizero remember it for later
    if ( trj->GetPDGEncoding() == 111)
      pizeroList.insert(trj->GetTrackID());
    // If it is a mu+/mu- also remember it
    if ( trj->GetPDGEncoding() == 13 ) muonList.insert(trj->GetTrackID());
    if ( trj->GetPDGEncoding() == -13 ) antimuonList.insert(trj->GetTrackID());
    if ( trj->GetPDGEncoding() == 211 ) pionList.insert(trj->GetTrackID());
    if ( trj->GetPDGEncoding() == -211 ) antipionList.insert(trj->GetTrackID());

    if( trj->GetParentID() == 0 ) primaryList.insert(trj->GetTrackID());

    // Process primary tracks or the secondaries from pizero or muons...

    if ( trj->GetSaveFlag() )
    {
      // initial point of the trajectory
      G4TrajectoryPoint* aa =   (G4TrajectoryPoint*)trj->GetPoint(0) ;   
      if(detectorElement=="tank") runAction->incrementEventsGenerated();
	
      G4int         ipnu   = trj->GetPDGEncoding();
      G4int         id     = trj->GetTrackID();
      G4int         flag   = 0;    // will be set later
      G4double      mass   = trj->GetParticleDefinition()->GetPDGMass();
      G4ThreeVector mom    = trj->GetInitialMomentum();
      G4double      mommag = mom.mag();
      G4double      energy = sqrt(mom.mag2() + mass*mass);
      G4ThreeVector Stop   = trj->GetStoppingPoint();
      G4ThreeVector Start  = aa->GetPosition();

      G4String stopVolumeName = trj->GetStoppingVolume()->GetName();
      G4int    stopvol     = WCSimEventFindStoppingVolume(stopVolumeName);
      G4int    startvol    = WCSimEventFindStartingVolume(Start);

      G4double ttime = trj->GetGlobalTime(); 

      G4int parentType;

     
      // Right now only secondaries whose parents are pi0's are stored
      // This may change later
      // M Fechner : dec 16, 2004 --> added decay e- from muons
      if (trj->GetParentID() == 0){
	parentType = 0;
      } else if (pizeroList.count(trj->GetParentID())   ) {
	parentType = 111;
      } else if (muonList.count(trj->GetParentID())     ) {
	parentType = 13;
      } else if (antimuonList.count(trj->GetParentID()) ) {
	parentType = -13;
      } else if (antipionList.count(trj->GetParentID()) ) {
	parentType = -211;
      } else if (pionList.count(trj->GetParentID()) ) {
	parentType = 211;
      } else if (primaryList.count(trj->GetParentID()) ) {
	parentType = 1;
      } else {  // no identified parent, but not a primary
	parentType = 999;
      }


      // G4cout << parentType << " " << ipnu << " " 
      //	     << id << " " << energy << "\n";

      // fill ntuple
      float dir[3];
      float pdir[3];
      float stop[3];
      float start[3];
      for (int l=0;l<3;l++)
      {
	dir[l]= mom[l]/mommag; // direction 
	pdir[l]=mom[l];        // momentum-vector 
	stop[l]=Stop[l]/CLHEP::cm; // stopping point 
	start[l]=Start[l]/CLHEP::cm; // starting point 
	//G4cout<<"part 2 start["<<l<<"]: "<< start[l] <<G4endl;
      }

      // Add the track to the TClonesArray, watching out for times
      if ( ! ( (ipnu==22)&&(parentType==999))  ) {
	int choose_event=0;

	if (ngates)
	{
	  if ( ttime > WCTM->GetTriggerTime(0)+950. && WCTM->GetTriggerTime(1)+950. > ttime ) choose_event=1; 
	  if ( ttime > WCTM->GetTriggerTime(1)+950. && WCTM->GetTriggerTime(2)+950. > ttime ) choose_event=2; 
	  if (choose_event >= ngates) choose_event = ngates-1; // do not overflow the number of events
	
	}
	wcsimrootevent= wcsimrootsuperevent->GetTrigger(choose_event);
	wcsimrootevent->AddTrack(ipnu, 
				  flag, 
				  mass, 
				  mommag, 
				  energy,
				  startvol, 
				  stopvol, 
				  dir, 
				  pdir, 
				  stop,
				  start,
				  parentType,
				 ttime,id); 
      }
      
      if (detectorConstructor->SavePi0Info() == true)
      {
	G4cout<<"Pi0 parentType: " << parentType <<G4endl;
	if (parentType == 111)
	{
	  if (r>1)
	    G4cout<<"WARNING: more than 2 primary gammas found"<<G4endl;
	  else
	  {

	    for (int y=0;y<3;y++)
	    {
	      pi0Vtx[y] = start[y];
	      gammaVtx[r][y] = stop[y];
	    }

	    gammaID[r] = id;
	    gammaE[r] = energy;
	    r++;
	
	    //amb79
		G4cout<<"Pi0 data: " << id <<G4endl;
		wcsimrootevent->SetPi0Info(pi0Vtx, gammaID, gammaE, gammaVtx);
	  }
	}
      }
    }
  }
  // Add the Cherenkov hits
  wcsimrootevent = wcsimrootsuperevent->GetTrigger(0);
  wcsimrootevent->SetNumTubesHit(jhfNtuple.numTubesHit);
#ifdef _SAVE_RAW_HITS

  if (WCDC_hits) 
  {
    //add the truth raw hits
    // Both the pre- and post-PMT smearing hit times are accessible
    // Choose to save just the pre-smeared times for now
#ifdef _SAVE_RAW_HITS_VERBOSE
    G4cout<<"RAW HITS"<<G4endl;
#endif
    wcsimrootevent->SetNumTubesHit(WCDC_hits->entries());
    std::vector<float> truetime, smeartime;
    std::vector<int>   primaryParentID;
    double hit_time_smear, hit_time_true;
    int hit_parentid;
    //loop over the DigitsCollection
    for(int idigi = 0; idigi < WCDC_hits->entries(); idigi++) {
      int digi_tubeid = (*WCDC_hits)[idigi]->GetTubeID();
      for(G4int id = 0; id < (*WCDC_hits)[idigi]->GetTotalPe(); id++){
	hit_time_true  = (*WCDC_hits)[idigi]->GetPreSmearTime(id);
	hit_parentid = (*WCDC_hits)[idigi]->GetParentID(id);
	truetime.push_back(hit_time_true);
	primaryParentID.push_back(hit_parentid);
#ifdef _SAVE_RAW_HITS_VERBOSE
	hit_time_smear = (*WCDC_hits)[idigi]->GetTime(id);
	smeartime.push_back(hit_time_smear);
#endif
      }//id
#ifdef _SAVE_RAW_HITS_VERBOSE
      if(digi_tubeid < NPMTS_VERBOSE) {
	G4cout << "Adding " << truetime.size()
	       << " Cherenkov hits in tube " << digi_tubeid
	       << " with truetime:smeartime:primaryparentID";
	for(G4int id = 0; id < truetime.size(); id++) {
	  G4cout << " " << truetime[id]
		 << ":" << smeartime[id]
		 << ":" << primaryParentID[id];
	}//id
	G4cout << G4endl;
      }
#endif
      wcsimrootevent->AddCherenkovHit(digi_tubeid,
				      truetime,
				      primaryParentID);
      smeartime.clear();
      truetime.clear();
      primaryParentID.clear();
    }//idigi
  }//if(WCDC_hits)
#endif //_SAVE_RAW_HITS

  // Add the digitized hits
  if (WCDC) 
  {
#ifdef SAVE_DIGITS_VERBOSE
    G4cout << "DIGI HITS" << G4endl;
#endif

    G4float sumq_tmp = 0.;
    
    for ( int index = 0 ; index < ngates ; index++)
      {	
	sumq_tmp = 0.0;	
	G4float gatestart;
	int countdigihits = 0;
	wcsimrootevent = wcsimrootsuperevent->GetTrigger(index);
	for (k=0;k<WCDC->entries();k++)
	  {
	    if ( (*WCDC)[k]->HasHitsInGate(index)) {
	      std::vector<float> vec_pe                  = (*WCDC)[k]->GetPe(index);
	      std::vector<float> vec_time                = (*WCDC)[k]->GetTime(index);
	      std::vector<std::vector<int> > vec_digicomp = (*WCDC)[k]->GetDigiCompositionInfo(index);
	      const int tubeID                           = (*WCDC)[k]->GetTubeID();
	      assert(vec_pe.size() == vec_time.size());
	      assert(vec_pe.size() == vec_digicomp.size());
	      for(unsigned int iv = 0; iv < vec_pe.size(); iv++) {
#ifdef SAVE_DIGITS_VERBOSE
		if(tubeID < NPMTS_VERBOSE) {
		  G4cout << "Adding digit " << iv 
			 << " for PMT " << tubeID
			 << " pe "   << vec_pe[iv]
			 << " time " << vec_time[iv]
			 << " digicomp";
		  for(unsigned int ivv = 0; ivv < vec_digicomp[iv].size(); ivv++)
		    G4cout << " " << vec_digicomp[iv][ivv];
		  G4cout << G4endl;
		}
#endif
		assert(vec_digicomp[iv].size() > 0);
		wcsimrootevent->AddCherenkovDigiHit(vec_pe[iv], vec_time[iv],
						    tubeID, vec_digicomp[iv]);
		sumq_tmp += vec_pe[iv];
		countdigihits++;
	      }//iv
	    }//Digit exists in Gate
	  }//k
	wcsimrootevent->SetNumDigitizedTubes(countdigihits);
	wcsimrootevent->SetSumQ(sumq_tmp);

#ifdef SAVE_DIGITS_VERBOSE
	G4cout << "checking digi hits ...\n";
	G4cout << "hits collection size (number of PMTs hit) =  " << 
	  wcsimrootevent->GetCherenkovHits()->GetEntries() << "\n";
	G4cout << "hits collection size (number of true photon + dark noise hits) =  " << 
	  wcsimrootevent->GetCherenkovHitTimes()->GetEntries() << "\n";
	G4cout << "digihits collection size =  " << 
	  wcsimrootevent->GetCherenkovDigiHits()->GetEntries() << "\n";
	G4cout << "tracks collection size =  " << 
	  wcsimrootevent->GetTracks()->GetEntries() 
	       <<" get ntracks = " <<  wcsimrootevent->GetNtrack() << "\n";
#endif
	gatestart = WCTM->GetTriggerTime(index);
	WCSimRootEventHeader*HH = wcsimrootevent->GetHeader();
	HH->SetDate(int(gatestart));
      }//index (loop over ngates)
    
    // end of loop over WC trigger gates --> back to the main sub-event
    wcsimrootevent = wcsimrootsuperevent->GetTrigger(0);
    
  }

  for (int i = 0 ; i < wcsimrootsuperevent->GetNumberOfEvents(); i++) {
    wcsimrootevent = wcsimrootsuperevent->GetTrigger(i);
    //G4cout << ">>>Root event "
	  // <<std::setw(5)<<wcsimrootevent->GetHeader()->GetEvtNum()<<"\n";
    //   if (WCDC){
    // G4cout <<"WC digi:"<<std::setw(4)<<wcsimrootevent->GetNcherenkovdigihits()<<"  ";
    // G4cout <<"WC digi sumQ:"<<std::setw(4)<<wcsimrootevent->GetSumQ()<<"  ";
  }
  
#ifdef _SAVE_RAW_HITS
  //if (WCHC)
  //     G4cout <<"WC:"<<std::setw(4)<<wcsimrootevent->GetNcherenkovhits()<<"  ";
  //    if (WCFVHC)
  //G4cout <<"WCFV:"<<std::setw(4)<<wcsimrootevent->GetINcherenkovhits()<<"  ";
#endif
  
  //  if (WCFVDC){
  //G4cout <<"WCFV digi:"<<std::setw(4)<<wcsimrootevent->GetNcherenkovdigihits()<<"  ";
  //G4cout <<"WCFV digi sumQ:"<<std::setw(4)<<wcsimrootevent->GetSumQ()<<"  ";
  //  }
  TTree* tree = GetRunAction()->GetTree();
  TBranch* branch = GetRunAction()->GetBranch(detectorElement);
  branch->Fill();
  //tree->Fill();
//  TFile* hfile = tree->GetCurrentFile();
//  // MF : overwrite the trees -- otherwise we have as many copies of the tree
//  // as we have events. All the intermediate copies are incomplete, only the
//  // last one is useful --> huge waste of disk space.
//  hfile->Write("",TObject::kOverwrite);
  
  // M Fechner : reinitialize the super event after the writing is over
  wcsimrootsuperevent->ReInitialize();
  
}

void WCSimEventAction::CreateNewLAPPDFile(){
  if(LAPPDfile){LAPPDfile->Close(); delete LAPPDfile; LAPPDfile=0;}
  // LAPPD file
  LAPPDRootFileName = GetRunAction()->GetRootFileNameBase() + "_lappd_" + std::to_string(GetRunAction()->GetOutputFileNum()) + ".root";
  LAPPDfile = new TFile(LAPPDRootFileName.c_str(),"RECREATE","WCSim LAPPD file");
  LAPPDfile->SetCompressionLevel(2);
  LAPPDfile->cd();
  LAPPDtree = new TTree("LAPPDTree","LAPPDTree");
  // Create the branches
  LAPPDtree->Branch("lappdevt",&lappdevt);
  LAPPDtree->Branch("lappd_numhits",&lappd_numhits, "lappd_numhits/I");
  LAPPDtree->Branch("lappdhit_totalpes_perevt", &lappdhit_totalpes_perevt, "lappdhit_totalpes_perevt/I");
  LAPPDtree->Branch("lappdhit_totalpes_perlappd2", &lappdhit_totalpes_perlappd2);

  LAPPDtree->Branch("lappdhit_x",lappdhit_x,"lappdhit_x[lappd_numhits]/D");
  LAPPDtree->Branch("lappdhit_y",lappdhit_y,"lappdhit_y[lappd_numhits]/D");
  LAPPDtree->Branch("lappdhit_z",lappdhit_z,"lappdhit_z[lappd_numhits]/D");
  LAPPDtree->Branch("lappdhit_t",lappdhit_t,"lappdhit_t[lappd_numhits]/D");
  LAPPDtree->Branch("lappdhit_stripcoorx", &lappdhit_stripcoorx);
  LAPPDtree->Branch("lappdhit_stripcoory", &lappdhit_stripcoory);
  LAPPDtree->Branch("lappdhit_stripcoorz", &lappdhit_stripcoorz);
  LAPPDtree->Branch("lappdhit_process",lappdhit_process,"lappdhit_process[lappd_numhits]/I");
  LAPPDtree->Branch("lappdhit_particleID",lappdhit_particleID,"lappdhit_particleID[lappd_numhits]/I");
  LAPPDtree->Branch("lappdhit_trackID",lappdhit_trackID,"lappdhit_trackID[lappd_numhits]/I");
  LAPPDtree->Branch("lappdhit_edep",lappdhit_edep,"lappdhit_edep[lappd_numhits]/D");
  LAPPDtree->Branch("lappdhit_objnum",lappdhit_objnum,"lappdhit_objnum[lappd_numhits]/I");
  //LAPPDtree->Branch("lappdhit_copynum",lappdhit_copynum,"lappdhit_copynum[lappd_numhits]/I");
  LAPPDtree->Branch("lappdhit_stripnum", &lappdhit_stripnum);
  LAPPDtree->Branch("lappdhit_truetime2",&lappdhit_truetime2);
  LAPPDtree->Branch("lappdhit_smeartime2", &lappdhit_smeartime2);
  LAPPDtree->Branch("lappdhit_primaryParentID2",&lappdhit_primaryParentID2);
  LAPPDtree->Branch("lappdhit_NoOfneighstripsHit", &lappdhit_NoOfneighstripsHit);
  LAPPDtree->Branch("lappdhit_neighstripnum", &lappdhit_neighstripnum);
  LAPPDtree->Branch("lappdhit_neighstrippeak", &lappdhit_neighstrippeak);
  LAPPDtree->Branch("lappdhit_neighstrip_time", &lappdhit_neighstrip_time);
  LAPPDtree->Branch("lappdhit_neighstrip_lefttime", &lappdhit_neighstrip_lefttime);
  LAPPDtree->Branch("lappdhit_neighstrip_righttime", &lappdhit_neighstrip_righttime);
}
