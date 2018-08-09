#include "WCSimTrackingAction.hh"
#include "WCSimTrajectory.hh"
#include "G4ParticleTypes.hh"
#include "G4TrackingManager.hh"
#include "G4Track.hh"
#include "G4ios.hh"
#include "G4VProcess.hh"
#include "WCSimTrackInformation.hh"
#include "G4TransportationManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include <algorithm>

WCSimTrackingAction::WCSimTrackingAction(){
  ProcessList.insert("Decay") ;
  ProcessList.insert("nCapture");
  ProcessList.insert("MuonMinusCaptureAtRest");
  //ProcessList.insert("conv");
  ParticleList.insert(111);  // pi0
  ParticleList.insert(211);  // pion+
  ParticleList.insert(-211); // pion-
  ParticleList.insert(321);  // kaon+
  ParticleList.insert(-321); // kaon-
  ParticleList.insert(311);  // kaon0
  ParticleList.insert(-311); // kaon0 bar
  ParticleList.insert(12);   // nu_e
  ParticleList.insert(-12);  // nubar_e
  ParticleList.insert(13);   // mu-
  ParticleList.insert(-13);  // mu+
  ParticleList.insert(14);   // nu_mu
  ParticleList.insert(-14);  // nubar_mu
  ParticleList.insert(2112); // neutron
  ParticleList.insert(2212); // proton
//  ParticleList.insert(11);   // e-    // do not save electrons unless they are from Decay process (mu decay)
//  ParticleList.insert(-11);  // e+
//  Don't put gammas there or there'll be too many
  
}

WCSimTrackingAction::~WCSimTrackingAction(){;}

void WCSimTrackingAction::PreUserTrackingAction(const G4Track* aTrack){
  G4float percentageOfCherenkovPhotonsToDraw = 0.0;
  
  if (aTrack->GetDefinition() != G4OpticalPhoton::OpticalPhotonDefinition()
       || G4UniformRand() < percentageOfCherenkovPhotonsToDraw){
      WCSimTrajectory* thisTrajectory = new WCSimTrajectory(aTrack);
      fpTrackingManager->SetTrajectory(thisTrajectory);
      fpTrackingManager->SetStoreTrajectory(true);
  } else {
      fpTrackingManager->SetStoreTrajectory(false);
  }
  
  /*
  // implemented to allow photon tracks to be drawn during photon debugging, 
  // but interferes with saving of primary information.
  WCSimTrackInformation* anInfo = new WCSimTrackInformation();
  G4Track* theTrack = (G4Track*)aTrack;
  anInfo->WillBeSaved(false);
  theTrack->SetUserInformation(anInfo);
  */
}

void WCSimTrackingAction::PostUserTrackingAction(const G4Track* aTrack){
  
  // added by M Fechner
  const G4VProcess* creatorProcess = aTrack->GetCreatorProcess();
  //  if ( creatorProcess )
  //    G4cout << "process name " << creatorProcess->GetProcessName() << G4endl;
  G4int thispdg = aTrack->GetDefinition()->GetPDGEncoding();
  if(aTrack->GetDefinition()==G4OpticalPhoton::OpticalPhotonDefinition()) thispdg=100;

  WCSimTrackInformation* anInfo;
  if (aTrack->GetUserInformation()){
    anInfo = (WCSimTrackInformation*)(aTrack->GetUserInformation());
  } else {
    anInfo = new WCSimTrackInformation();
  }

  // is it a primary ?
  // is the process in the set ? 
  // is the particle in the set ?
  // is it a gamma 
  // due to lazy evaluation of the 'or' in C++ the order is important
//  if( aTrack->GetParentID()==0 || 
//      ((creatorProcess!=0) && ProcessList.count(creatorProcess->GetProcessName()) ) || 
//      (ParticleList.count(thispdg) )
//      || (thispdg==22 && aTrack->GetTotalEnergy() > 50.0*CLHEP::MeV)
//      )
  // save more information. code from wcsim github issue 197.
  if( ( aTrack->GetParentID()==0 ) ||
      ( (creatorProcess!=0) && ProcessList.count(creatorProcess->GetProcessName()) ) ||
      ( ParticleList.count(thispdg) ) ||
      ( thispdg==22 && aTrack->GetTotalEnergy()>1.0*MeV ) ||
      ( creatorProcess->GetProcessName()=="muMinusCaptureAtRest" && aTrack->GetTotalEnergy()>1.0*MeV )
    ){
    // G4cout << "track # " << aTrack->GetTrackID() << " is worth saving\n";
    // G4cout << "It is a " <<aTrack->GetDefinition()->GetParticleName() << G4endl;
    anInfo->WillBeSaved(true);
  } else {
    anInfo->WillBeSaved(false);
  }

  if( thispdg==111 ){
    pi0List.insert(aTrack->GetTrackID()); // list of all pi0-s 
  }

  if((aTrack->GetParentID()==0 && aTrack->GetDefinition()!=G4OpticalPhoton::OpticalPhotonDefinition()) // primary
      || ( thispdg==22 && pi0List.count(aTrack->GetParentID()) ) // primary gamma from a pi0
    ){
    anInfo->SetPrimaryParentID(aTrack->GetTrackID());
  }

  G4Track* theTrack = (G4Track*)aTrack;
  theTrack->SetUserInformation(anInfo);
  
  // pass primary parent ID to children
  G4TrackVector* secondaries = fpTrackingManager->GimmeSecondaries();
  if(secondaries){
      for(size_t i=0;i<secondaries->size();i++){
        WCSimTrackInformation* infoSec = new WCSimTrackInformation(anInfo);
                 infoSec->WillBeSaved(false); // ADDED BY MFECHNER, temporary, 30/8/06
        infoSec->SetParentPdg(thispdg);
        (*secondaries)[i]->SetUserInformation(infoSec);
      }
  }

  if ( aTrack->GetDefinition() != G4OpticalPhoton::OpticalPhotonDefinition()){
    WCSimTrajectory *currentTrajectory = (WCSimTrajectory*)fpTrackingManager->GimmeTrajectory();
  
    G4ThreeVector currentPosition      = aTrack->GetPosition();
    G4VPhysicalVolume* currentVolume   = aTrack->GetVolume();
    G4double currentTime = aTrack->GetGlobalTime();
    G4ThreeVector currentMomentum = aTrack->GetMomentum();

    currentTrajectory->SetStoppingPoint(currentPosition);
    currentTrajectory->SetStoppingVolume(currentVolume);
    currentTrajectory->SetStoppingTime(currentTime);
    currentTrajectory->SetStoppingMomentum(currentMomentum);
    currentTrajectory->SetParentPdg(anInfo->GetParentPdg());

    // if UserTrackInformation is tagged to be saved, mark trajectory for WCSimEventAction;
    if (anInfo->isSaved())
      currentTrajectory->SetSaveFlag(true);
    else
      currentTrajectory->SetSaveFlag(false);
  }
  
  static int line=0;
  if(line%100000==0){ //100000
    G4cout<<"  PostUserTrackingAction call number: "<<line
          <<", "<<aTrack->GetDefinition()->GetParticleName();
    if(creatorProcess) G4cout<<" from "<<creatorProcess->GetProcessName();
    else G4cout<<"primary";
    G4cout<<" in "<<aTrack->GetVolume()->GetName()<<G4endl; 
  }
  line++;
  
  /*
  if( (aTrack->GetParentID()==0) // primary particle
      && (abs(thispdg)==13) ){ // is a muon
      G4ThreeVector endpos = aTrack->GetPosition();
      WCSimTrajectory* trj = (WCSimTrajectory*)fpTrackingManager->GimmeTrajectory();
      G4TrajectoryPoint* startpnt = (G4TrajectoryPoint*)trj->GetPoint(0);
      G4ThreeVector startpos  = startpnt->GetPosition();
      G4cout<<"Primary muon started at ("<<startpos.x()<<", "<<startpos.y()<<", "<<startpos.z()<<")"
            <<"and ended at ("<<endpos.x()<<", "<<endpos.y()<<", "<<endpos.z()<<")"<<G4endl;
  }
  */
}





