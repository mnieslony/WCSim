// ====================================================================
//   LAPPDSD.cc
//
//   26/11/15 M. O'Flaherty (based on SBsimLAPPDSD by 2006/03/03 K. Hiraide)
// ====================================================================
#include "G4VPhysicalVolume.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4VTouchable.hh"
#include "G4TouchableHistory.hh"
#include "G4RunManager.hh"

#include "LAPPDSD.hh"
#include "LAPPDHit.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4OpBoundaryProcess.hh"

//#include <CLHEP/Random/Randomize.h>

//////////////////////////////////////////////////
LAPPDSD::LAPPDSD(const G4String& name) : G4VSensitiveDetector(name) {

  collectionName.insert("lappdHitsCollection");
  G4cout<<"Constructing LAPPD Sensitive Detector."<<G4endl;
  
  hcid = -1;

  //G4RunManager *runMgr = G4RunManager::GetRunManager();
  //detector = (WCLiteDetectorConstruction*)runMgr->GetUserDetectorConstruction();
  //mrdresp = new SBsimLAPPDResponse(detector);
  //mrddb   = detector->GetLAPPDDB();
}

///////////////////////////////
LAPPDSD::~LAPPDSD(){}

//////////////////////////////////////////////////
void LAPPDSD::Initialize(G4HCofThisEvent* HCTE){

  // create hit collection(s)
  G4cout<<"Initializing LAPPD hit collection"<<G4endl;
  hitsCollection = new LAPPDHitsCollection(SensitiveDetectorName, collectionName[0]); 
  
  // push H.C. to "Hit Collection of This Event"
  
  hcid= GetCollectionID(0); 
  //if(hcid < 0) { hcid = G4SDManager::GetSDMpointer()->GetCollectionID(0); }
  HCTE-> AddHitsCollection(hcid, hitsCollection);

  // initialization
  nhit=0;

}

//////////////////////////////////////////////////////////////
G4bool LAPPDSD::ProcessHits(G4Step* aStep, G4TouchableHistory*) { // weird syntax; type declared but no variable??

    G4Track* aTrack = aStep->GetTrack();
    G4ParticleDefinition* particleType = aTrack->GetDefinition();
    // don't make hits from optical absorption...
    if(particleType==G4OpticalPhoton::OpticalPhotonDefinition()){return true;}
    // seem to have StepStatus = fGeomBoundary, ProcessDefinedStep = Transportation. not absorption..??
    	
    // don't make hits if no energy deposited...
    G4double eDeposited = aStep->GetTotalEnergyDeposit();
    if (eDeposited == 0) return true;

    ++nhit;    
    newHit= new LAPPDHit(aStep);
    newHit->SetHitID(nhit);      
    hitsCollection->insert(newHit);
    //G4cout<<G4endl<<"An LAPPD hit!"<<G4endl;
    //newHit->Print();

  // Scintillator response (quenching, attenuation)
  // This calls in code from SBSimLAPPDResponse which applies attenuation and Birk's constant to modify the deposited enegy. 
  // Probably good to implement, but depends on the geometry.
  //mrdresp->ApplyScintiResponse(&edepbuf,aTrack,preStepPoint);	//pass in ref to edepbuf (eDeposited from track) which is modified.

  /*
  if(particle->GetPDGCharge()!=0 && kineticE>0) {
    G4double dedx = emcal.GetDEDX(kineticE, particle, material)/(MeV/mm);
    steplength[id] += (aStep->GetStepLength())/mm;
    dedxtimesstep[id] += dedx*((aStep->GetStepLength())/mm);
  }
  */
  return true;
}

//////////////////////////////////////////////////////
void LAPPDSD::EndOfEvent(G4HCofThisEvent* /*HCTE*/)
//////////////////////////////////////////////////////
{
	//SciBooNE code adds noise at this point. Good idea?
	G4cout << "**LAPPDSD end of event action**" << G4endl;
	G4cout << "Number of LAPPD hits: " << nhit << G4endl;
	/*for (G4int i=0;i<5;i++){
		LAPPDHit* ahit = (LAPPDHit*)hitsCollection->GetHit(i);
		ahit->Print();
	}*/
	//hitsCollection->PrintAllHits();
}

/////////////////////////////
void LAPPDSD::DrawAll() {
  hitsCollection-> DrawAllHits();
}

//////////////////////////////
void LAPPDSD::PrintAll() {
  if(nhit>0){
    G4cout << "##### LAPPD Hit #####" << G4endl;
    hitsCollection-> PrintAllHits();
    G4cout << G4endl;
  }
}
