#include "WCSimWCPMT.hh"
#include "WCSimWCDigi.hh"
#include "WCSimWCHit.hh"

#include "G4EventManager.hh"
#include "G4Event.hh"
#include "G4SDManager.hh"
#include "G4DigiManager.hh"
#include "G4ios.hh"
#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"

#include "WCSimDetectorConstruction.hh"
#include "WCSimPmtInfo.hh"
#include "WCSimPMTObject.hh"


#include <vector>
// for memset
#include <cstring>


extern "C" void skrn1pe_(float* );
//extern "C" void rn1pe_(float* ); // 1Kton

WCSimWCPMT::WCSimWCPMT(G4String name,
				   WCSimDetectorConstruction* myDetector, G4String detectorElement)
  :G4VDigitizerModule(name), detectorElement(detectorElement)
{
  //G4String colName = "WCRawPMTSignalCollection";
  this->myDetector = myDetector;
  if(detectorElement=="tank"){
  	collectionName.push_back("WCRawPMTSignalCollection");	// ☆
  } else if(detectorElement=="mrd"){
  	collectionName.push_back("WCRawMRDSignalCollection");
  } else if(detectorElement=="facc"){
  	collectionName.push_back("WCRawFACCSignalCollection");
  }
  DigiHitMapPMT.clear();
  

}

WCSimWCPMT::~WCSimWCPMT(){
 
}

G4double WCSimWCPMT::rn1pe(){
  //G4String WCIDCollectionName = myDetector->GetIDCollectionName();	// ☆
  WCSimPMTObject * PMT;
  if(detectorElement=="tank"){
  	PMT = myDetector->GetPMTPointer(myDetector->GetIDCollectionName());	// ☆ arg WCIDCollectionName
  } else if(detectorElement=="mrd"){ 
  	PMT = myDetector->GetPMTPointer(myDetector->GetMRDCollectionName());
  } else if(detectorElement=="facc"){
  	PMT = myDetector->GetPMTPointer(myDetector->GetFACCCollectionName());
  }
  G4int i;
  G4double random = G4UniformRand();
  G4double random2 = G4UniformRand(); 
  G4float *qpe0;
  qpe0 = PMT->Getqpe();
  for(i = 0; i < 501; i++){
    
    if (random <= *(qpe0+i)) break;
  }
  if(i==500)
    random = G4UniformRand();
  
  return (G4double(i-50) + random2)/22.83;
  
}


void WCSimWCPMT::Digitize()
{
  // Create a DigitCollection and retrieve the appropriate hitCollection name based on detectorElement
  G4String WCCollectionName;
  if(detectorElement=="tank"){
  	DigitsCollection = new WCSimWCDigitsCollection ("WCDigitizedCollection",collectionName[0]); //☆
  	// Get the Associated Hit collection ID
  	WCCollectionName = myDetector->GetIDCollectionName(); //☆
  } else if(detectorElement=="mrd"){
  	DigitsCollection = new WCSimWCDigitsCollection ("WCDigitizedCollection_MRD",collectionName[0]);
  	// Get the Associated Hit collection ID
  	WCCollectionName = myDetector->GetMRDCollectionName();
  } else if(detectorElement=="facc"){
  	DigitsCollection = new WCSimWCDigitsCollection ("WCDigitizedCollection_FACC",collectionName[0]);
  	// Get the Associated Hit collection ID
  	WCCollectionName = myDetector->GetFACCCollectionName();
  }

  G4DigiManager* DigiMan = G4DigiManager::GetDMpointer();
  // Get the Associated Hit collection ID
  G4int WCHCID = DigiMan->GetHitsCollectionID(WCCollectionName);	// ☆
  // The Hits collection
  WCSimWCHitsCollection* WCHC = (WCSimWCHitsCollection*)(DigiMan->GetHitsCollection(WCHCID));
  
  if (WCHC) {
    MakePeCorrection(WCHC);
  }

  StoreDigiCollection(DigitsCollection);

}


void WCSimWCPMT::MakePeCorrection(WCSimWCHitsCollection* WCHC)
{ 

  //Get the PMT info for hit time smearing
  G4String WCCollectionName;	// ☆
  if(detectorElement=="tank"){
  	WCCollectionName = myDetector->GetIDCollectionName();
  } else if(detectorElement=="mrd"){
  	WCCollectionName = myDetector->GetMRDCollectionName();
  } else if(detectorElement=="facc"){
  	WCCollectionName = myDetector->GetFACCCollectionName();
  }
  
  WCSimPMTObject * PMT = myDetector->GetPMTPointer(WCCollectionName);	//☆

  for (G4int i=0; i < WCHC->entries(); i++)
    {

      //G4double peCutOff = .3;
      // MF, based on S.Mine's suggestion : global scaling factor applied to
      // all the smeared charges.
      // means that we need to increase the collected light by
      // (efficiency-1)*100% to
      // match K2K 1KT data  : maybe due to PMT curvature ?

      //G4double efficiency = 0.985; // with skrn1pe (AP tuning) & 30% QE increase in stacking action

      // Get the information from the hit
      G4int   tube         = (*WCHC)[i]->GetTubeID();
      G4double peSmeared = 0.0;
      double time_PMT, time_true;

	  for (G4int ip =0; ip < (*WCHC)[i]->GetTotalPe(); ip++){
	    time_true = (*WCHC)[i]->GetTime(ip);
	    peSmeared = rn1pe();
	    int parent_id = (*WCHC)[i]->GetParentID(ip);

	    //apply time smearing
	    float Q = (peSmeared > 0.5) ? peSmeared : 0.5;
	    time_PMT = time_true + PMT->HitTimeSmearing(Q);

	    if ( DigiHitMapPMT[tube] == 0) {
	      WCSimWCDigi* Digi = new WCSimWCDigi();
	      Digi->SetLogicalVolume((*WCHC)[0]->GetLogicalVolume());
	      Digi->AddPe(time_PMT);	
	      Digi->SetTubeID(tube);
	      Digi->SetPe(ip,peSmeared);
	      Digi->SetTime(ip,time_PMT);
	      Digi->SetPreSmearTime(ip,time_true);
	      Digi->SetParentID(ip,parent_id);
	      DigiHitMapPMT[tube] = DigitsCollection->insert(Digi);
	    }	
	    else {
	      (*DigitsCollection)[DigiHitMapPMT[tube]-1]->AddPe(time_PMT);
	      (*DigitsCollection)[DigiHitMapPMT[tube]-1]->SetLogicalVolume((*WCHC)[0]->GetLogicalVolume());
	      (*DigitsCollection)[DigiHitMapPMT[tube]-1]->SetTubeID(tube);
	      (*DigitsCollection)[DigiHitMapPMT[tube]-1]->SetPe(ip,peSmeared);
	      (*DigitsCollection)[DigiHitMapPMT[tube]-1]->SetTime(ip,time_PMT);
	      (*DigitsCollection)[DigiHitMapPMT[tube]-1]->SetPreSmearTime(ip,time_true);
	      (*DigitsCollection)[DigiHitMapPMT[tube]-1]->SetParentID(ip,parent_id);
	    }
      
	  } // Loop over hits in each PMT
    }// Loop over PMTs
}


