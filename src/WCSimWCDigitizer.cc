#include "WCSimWCDigitizer.hh"
#include "WCSimWCPMT.hh"
#include "WCSimWCLAPPD.hh"
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
#include "WCSimDarkRateMessenger.hh"

#include <vector>
// for memset
#include <cstring>
#include <exception>

#ifndef NPMTS_VERBOSE
//#define NPMTS_VERBOSE 10
#endif

#ifndef HYPER_VERBOSITY
//#define HYPER_VERBOSITY
#endif


// *******************************************
// BASE CLASS
// *******************************************

#ifndef WCSIMWCDIGITIZER_VERBOSE
//#define WCSIMWCDIGITIZER_VERBOSE
#endif

WCSimWCDigitizerBase::WCSimWCDigitizerBase(G4String name,
					   WCSimDetectorConstruction* inDetector,
					   WCSimWCDAQMessenger* myMessenger,
					   DigitizerType_t digitype, G4String detectorElementin)
  :G4VDigitizerModule(name), myDetector(inDetector), DAQMessenger(myMessenger), DigitizerType(digitype), detectorElement(detectorElementin), DigitizerClassName("")
{
  G4String colName;
  if(detectorElement=="tank"){
  	colName = "WCDigitizedStoreCollection";
  } else if(detectorElement=="mrd"){
  	colName = "WCDigitizedStoreCollection_MRD";
  } else if(detectorElement=="facc"){
  	colName = "WCDigitizedStoreCollection_FACC";
  }
  collectionName.push_back(colName);
  ReInitialize();
#ifdef HYPER_VERBOSITY
  if(detectorElement=="mrd"){
  G4cout<<"WCSimWCDigitizerBase::WCSimWCDigitizerBase ☆ recording collection name "<<colName<<" for "<<detectorElement<<G4endl;}
#endif
  if(DAQMessenger == NULL) {
    G4cerr << "WCSimWCDAQMessenger pointer is NULL when passed to WCSimWCDigitizerBase constructor. Exiting..." 
	   << G4endl;
    exit(-1);
  }
}

WCSimWCDigitizerBase::~WCSimWCDigitizerBase(){
}

void WCSimWCDigitizerBase::GetVariables()
{
  //set the options to digitizer-specific defaults
  DigitizerDeadTime          = GetDefaultDeadTime();
  DigitizerIntegrationWindow = GetDefaultIntegrationWindow();
  ExtendDigitizerIntegrationWindow = GetDefaultExtendIntegrationWindow();
  DoPhotonIntegration        = GetDefaultDoPhotonIntegration();

  //read the .mac file to override them
  if(DAQMessenger != NULL) {
    DAQMessenger->TellMeAboutTheDigitizer(this);
    DAQMessenger->SetDigitizerOptions(detectorElement);
  }
  else {
    G4cerr << "WCSimWCDAQMessenger pointer is NULL when used in WCSimWCDigitizerBase::GetVariables(). Exiting..." 
	   << G4endl;
    exit(-1);
  }

  G4cout << "Using digitizer deadtime " << DigitizerDeadTime << " ns" << G4endl;
  G4cout << "Using digitizer integration window " << DigitizerIntegrationWindow << " ns" << G4endl;
  if(ExtendDigitizerIntegrationWindow) G4cout<<"Will extend the digitizer integration window "
     << "when new hits arrive within an existing window" <<G4endl;
  else G4cout<<"Will not extend the digitizer integration window "
     << "when new hits arrive within an existing window" <<G4endl;
  if(DoPhotonIntegration) G4cout<<"Will integrate photons in the digitizer integration window into one digit" << G4endl;
  else G4cout<<"Will create a digit from each individual photon" <<G4endl;
}

void WCSimWCDigitizerBase::Digitize()
{
  //Input is WCSimWCDigitsCollection with raw PMT hits (photon + dark noise)
  //Output is WCSimWCDigitsCollection with digitied PMT hits

  //Clear the DigiStoreHitMap
  ReInitialize();
  //Temporary Storage of Digitized hits which is passed to the trigger
  DigiStore = new WCSimWCDigitsCollection(collectionName[0],collectionName[0]);

  G4DigiManager* DigiMan = G4DigiManager::GetDMpointer();
  
  // Get the PMT collection ID
  G4String rawcollectionName;
  if(detectorElement=="tank"){
  	rawcollectionName = "WCRawPMTSignalCollection";
  } else if(detectorElement=="mrd"){
  	rawcollectionName = "WCRawPMTSignalCollection_MRD";
  } else if(detectorElement=="facc"){
  	rawcollectionName = "WCRawPMTSignalCollection_FACC";
  }
  G4int WCHCID = DigiMan->GetDigiCollectionID(rawcollectionName);

  // Get the PMT Digits collection
  WCSimWCDigitsCollection* WCHCPMT = (WCSimWCDigitsCollection*)(DigiMan->GetDigiCollection(WCHCID));
  
#ifdef HYPER_VERBOSITY
  if(detectorElement=="mrd"){
  G4cout<<"WCSimWCDigitizerBase::Digitize ☆ making digits collection (WCSimWCDigitsCollection*)"<<collectionName[0]
  <<" for "<<detectorElement<<" and calling DigitizeHits on "<<rawcollectionName<<" to fill it"<<G4endl;}
#endif
  
  if (WCHCPMT) {
    DigitizeHits(WCHCPMT);
  } else {G4cout<<"WCSimWCDigitizerBase::Digitize didn't find hit collection for "<<detectorElement<<G4endl;}
  
  StoreDigiCollection(DigiStore);

}

bool WCSimWCDigitizerBase::AddNewDigit(int tube, int gate, float digihittime, float peSmearedin, std::vector<int> digi_comp)
{

  //gate is not a trigger, but just the position of the digit in the array
  //inside the WCSimWCDigi object
#ifdef WCSIMWCDIGITIZER_VERBOSE
  if(tube < NPMTS_VERBOSE) {
    G4cout<<"Adding hit "<<gate<<" in tube "<<tube
	  << " with time " << digihittime << " charge " << peSmearedin
	  << " (made of " << digi_comp.size() << " raw hits with IDs ";
    for(unsigned int iv = 0; iv < digi_comp.size(); iv++)
      G4cout << " " << digi_comp[iv] << ",";
    G4cout << ")";
  }
#endif

  if (peSmearedin > 0.0) {
      if ( DigiStoreHitMap[tube] == 0) {
	WCSimWCDigi* Digi = new WCSimWCDigi();
	Digi->SetTubeID(tube);
	Digi->SetPe(gate,peSmearedin);
	Digi->AddPe(digihittime);
	Digi->SetTime(gate,digihittime);
	Digi->AddDigiCompositionInfo(digi_comp);
	DigiStoreHitMap[tube] = DigiStore->insert(Digi);
#ifdef WCSIMWCDIGITIZER_VERBOSE
	if(tube < NPMTS_VERBOSE)
	  G4cout << " NEW HIT" << G4endl;
#endif
      }
      else {
	(*DigiStore)[DigiStoreHitMap[tube]-1]->SetPe(gate,peSmearedin);
	(*DigiStore)[DigiStoreHitMap[tube]-1]->SetTime(gate,digihittime);
	(*DigiStore)[DigiStoreHitMap[tube]-1]->AddPe(digihittime);
	(*DigiStore)[DigiStoreHitMap[tube]-1]->AddDigiCompositionInfo(digi_comp);
#ifdef WCSIMWCDIGITIZER_VERBOSE
	if(tube < NPMTS_VERBOSE)
	  G4cout << " DEJA VU" << G4endl;
#endif
      }
      return true;
  }//peSmearedin > 0
  else {
#ifdef WCSIMWCDIGITIZER_VERBOSE
    if(tube < NPMTS_VERBOSE)
      G4cout << "DIGIT REJECTED with charge " << peSmearedin
	     << " time " << digihittime << G4endl;
#endif
    return false;
  }
}

void WCSimWCDigitizerBase::SaveOptionsToOutput(WCSimRootOptions * wcopt)
{
  wcopt->SetDigitizerClassName(DigitizerClassName);
  wcopt->SetDigitizerDeadTime(DigitizerDeadTime);
  wcopt->SetDigitizerIntegrationWindow(DigitizerIntegrationWindow);
  wcopt->SetExtendIntegrationWindow(ExtendDigitizerIntegrationWindow);
  wcopt->SetDoPhotonIntegration(DoPhotonIntegration);
}


// *******************************************
// DERIVED CLASS
// *******************************************

WCSimWCDigitizerSKI::WCSimWCDigitizerSKI(G4String name,
					 WCSimDetectorConstruction* myDetector,
					 WCSimWCDAQMessenger* myMessenger, G4String detectorElement)
  : WCSimWCDigitizerBase(name, myDetector, myMessenger, kDigitizerSKI, detectorElement)
{
  DigitizerClassName = "SKI";
  GetVariables();
}

WCSimWCDigitizerSKI::~WCSimWCDigitizerSKI(){
}

void WCSimWCDigitizerSKI::DigitizeHits(WCSimWCDigitsCollection* WCHCPMT) {

#ifdef HYPER_VERBOSITY
  if(detectorElement=="mrd"){G4cout<<"WCSimWCDigitizerBase::DigitizeHits ☆ digitizing "<<WCHCPMT->entries()<<" entries"<<G4endl;}
#endif
  G4cout << "WCSimWCDigitizerSKI::DigitizeHits START ";
  if(detectorElement=="tank"){ G4cout<<"WCHCPMT->entries() = ";}
  if(detectorElement=="mrd"){  G4cout<<"HCMRD  ->entries() = ";}
  if(detectorElement=="facc"){ G4cout<<"HCFACC ->entries() = ";}
  G4cout<< WCHCPMT->entries() << G4endl;
  
  //loop over entires in WCHCPMT, each entry corresponds to
  //the photons on one PMT
  int absoluteindex=0;
  for (G4int i = 0 ; i < WCHCPMT->entries() ; i++)
    {

      //We must first sort hits by PMT in time.  This is very important as the code
      //assumes that each hit is in time order from lowest to highest.
      (*WCHCPMT)[i]->SortDigiMapsByHitTime();
      int tube = (*WCHCPMT)[i]->GetTubeID();
#ifdef WCSIMWCDIGITIZER_VERBOSE
      if(tube < NPMTS_VERBOSE) {
	G4cout << "tube " << tube
	       << " totalpe = " << (*WCHCPMT)[i]->GetTotalPe()
	       << " times";
	for(int ip = 0; ip < (*WCHCPMT)[i]->GetTotalPe(); ip++)
	  G4cout << " " << (*WCHCPMT)[i]->GetTime(ip);
	/*
	  G4cout<<" parents =\t";
	  for( G4int ip = 0 ; ip < (*WCHCPMT)[i]->GetTotalPe() ; ip++)
	  G4cout << " " << (*WCHCPMT)[i]->GetParentID(ip);
	*/
	G4cout <<G4endl;
      }
#endif

      //Sorting done.  Now we integrate the charge on each PMT.
      // Integration occurs for DigitizerIntegrationWindow ns (user set)
      // Digitizer is then dead for DigitizerDeadTime ns (user set)

      //look over all hits on the PMT
      //integrate charge and start digitizing
      float intgr_start=0;
      float upperlimit=0;
      G4double efficiency = 0.985; // with skrn1pe (AP tuning) & 30% QE increase in stacking action

      // Variables to store photon uniqueid that make up a digit
      int digi_unique_id   = 0;
      int photon_unique_id = 0;
      std::vector<int> digi_comp; 
      std::vector<float> digi_times;

      //loop over the hits on this PMT
#ifdef WCSIMWCDIGITIZER_VERBOSE
      if(detectorElement=="tank"){
        G4cout<<"tank pmt "<<i<<" had "<<(*WCHCPMT)[i]->GetTotalPe()<<" photon hits"<<G4endl;
      }
      int numdigitsrequested=0;
      int numdigitsrejectedthreshold=0;
      int numdigitsrejectedpe=0;
      int numdigitscreated=0;
#endif
      for( G4int ip = 0 ; ip < (*WCHCPMT)[i]->GetTotalPe() ; ip++)
	{
	  float time=0.;
	  try{
	    time = (*WCHCPMT)[i]->GetTime(ip);
	  }
	  catch (...){
	    G4cout<<"Exception in WCSimWCDigitizerSKI::DigitizeHits call to WCSimWCDigi::GetTime "
	          <<G4endl<<"Attempt to retreive time from pe "<<ip<<" in WCHCPMT entry "<<i<<G4endl;
	    G4cout<<"This digi had "<<(*WCHCPMT)[i]->GetTotalPe()<<" total pes"<<G4endl;
	    assert(false);
	  }
          float pe = (*WCHCPMT)[i]->GetPe(ip);

	  //start the integration time as the time of the first hit
	  //Hits must be sorted in time
	  if( (DoPhotonIntegration==false) || (ip==0) ){
	    intgr_start=time;
	    peSmeared = 0;
	    //Set the limits of the integration window [intgr_start,upperlimit]
	    upperlimit = intgr_start + DigitizerIntegrationWindow;
	  }
	  
#ifdef WCSIMWCDIGITIZER_VERBOSE
	  if(tube < NPMTS_VERBOSE)
	    G4cout << "ip "    << ip
		   << " pe "   << pe
		   << " time " << time
		   << " intgr_start " << intgr_start
		   << " upperlimit "  << upperlimit
		   << G4endl;
#endif

	  bool MakeDigit = false;
	  if(time >= intgr_start && time <= upperlimit) {
	    peSmeared += pe;
	    photon_unique_id = ip+absoluteindex;
	    digi_comp.push_back(photon_unique_id);
	    digi_times.push_back(time);
	    // extend the integration window, if enabled
	    if(ExtendDigitizerIntegrationWindow) upperlimit = time + DigitizerIntegrationWindow;
      
#ifdef WCSIMWCDIGITIZER_VERBOSE
	    if(tube < NPMTS_VERBOSE)
	      G4cout<<"INFO: time "<<time<<" digi_id "<<digi_unique_id<<" p_id "<<photon_unique_id<<G4endl;
#endif
	    //if this is the last digit, make sure to make the digit
	    if( (DoPhotonIntegration==false) || (ip + 1 == (*WCHCPMT)[i]->GetTotalPe()) ){
	      MakeDigit = true;
	    }
	    
	  }
	  //if ensures we don't append the same digit multiple times while in the integration window
	  else if(digi_comp.size()) {
	    //this hit is outside the integration time window.
	    //Charge integration is over.  The is now a DigitizerDeadTime ns dead
	    //time period where no hits can be recorded
	    MakeDigit = true;
	  }
	  
	  //Make digit here
	  if(MakeDigit) {
	    int iflag;
	    WCSimWCDigitizerSKI::Threshold(peSmeared,iflag);

	    //Check if previous hit passed the threshold.  If so we will digitize the hit
#ifdef WCSIMWCDIGITIZER_VERBOSE
	    numdigitsrequested++;
#endif
	    if(iflag == 0) {
	      //digitize hit
	      peSmeared *= efficiency;
	      // use the median time as the time of the digit
	      int median_index=std::min(1,int(double(digi_times.size())/2.));
	      float median_time = digi_times.at(median_index);
	      bool accepted = WCSimWCDigitizerBase::AddNewDigit(tube, digi_unique_id, median_time, peSmeared, digi_comp);
	      if(accepted) {
		digi_unique_id++;
#ifdef WCSIMWCDIGITIZER_VERBOSE
		numdigitscreated++;
	      } else {
	        numdigitsrejectedpe++;
#endif
	      }
	      assert(digi_comp.size());
	      digi_comp.clear();
	      digi_times.clear();
	    }
	    else {
	      //reject hit
#ifdef WCSIMWCDIGITIZER_VERBOSE
	      numdigitsrejectedthreshold++;
	      if(tube < NPMTS_VERBOSE)
		G4cout << "DIGIT REJECTED with time " << intgr_start << G4endl;
#endif
	      digi_comp.clear();
	      digi_times.clear();
	    }
	  }
	  
	  //Now try and deal with the next hit
	  if(time > upperlimit && time <= upperlimit + DigitizerDeadTime) {
	    //Now we need to reject hits that are after the integration
	    //period to the end of the veto signal
	    continue;
	  }
	  else if(time > upperlimit + DigitizerDeadTime){
#ifdef WCSIMWCDIGITIZER_VERBOSE
	    if(tube < NPMTS_VERBOSE)
	      G4cout<<"*** PREPARING FOR >1 DIGI ***"<<G4endl;
#endif
	    //we now need to start integrating from the hit
	    intgr_start=time;
	    peSmeared = pe;
	    //Set the limits of the integration window [intgr_start,upperlimit]
	    upperlimit = intgr_start + DigitizerIntegrationWindow;

	    //store the digi composition information
	    photon_unique_id = ip+absoluteindex;
            digi_comp.push_back(photon_unique_id);
            digi_times.push_back(time);

	    //if this is the last hit we must handle the creation of the digit 
	    //as the loop will not evaluate again
	    if(ip+1 == (*WCHCPMT)[i]->GetTotalPe()) {
	      int iflag;
	      WCSimWCDigitizerSKI::Threshold(peSmeared,iflag);
	      if(iflag == 0) {
		//digitize hit
		peSmeared *= efficiency;
		// use the median time as the time of the digit
		int median_index=std::min(1,int(double(digi_times.size())/2.));
		float median_time = digi_times.at(median_index);
		bool accepted = WCSimWCDigitizerBase::AddNewDigit(tube, digi_unique_id, median_time, peSmeared, digi_comp);
		if(accepted) {
		  digi_unique_id++;
		}
		assert(digi_comp.size());
		digi_comp.clear();
		digi_times.clear();
	      }
	      else {
		//reject hit
#ifdef WCSIMWCDIGITIZER_VERBOSE
		if(tube < NPMTS_VERBOSE)
		  G4cout << "DIGIT REJECTED with time " << intgr_start << G4endl;
#endif
		digi_comp.clear();
		digi_times.clear();
	      }
	    }
	  }
	}//ip (totalpe)
	absoluteindex+=(*WCHCPMT)[i]->GetTotalPe();
#ifdef WCSIMWCDIGITIZER_VERBOSE
	G4cout<<"Requested the creation of "<<numdigitsrequested<<" digits; "
	      <<numdigitscreated<<" were created, "
	      <<numdigitsrejectedthreshold<<" were rejected by digitizer threshold, "
	      <<numdigitsrejectedpe<<" were rejected by pe limit"<<G4endl;
#endif
    }//i (WCHCPMT->entries())
  G4cout<<"WCSimWCDigitizerSKI::DigitizeHits END DigiStore->entries() " << DigiStore->entries() << "\n";
  
#ifdef WCSIMWCDIGITIZER_VERBOSE
  G4cout<<"\n\n\nCHECK DIGI COMP:"<<G4endl;
  for (G4int idigi = 0 ; idigi < DigiStore->entries() ; idigi++){
    int tubeid = (*DigiStore)[idigi]->GetTubeID();
    if(tubeid < NPMTS_VERBOSE) {
      std::map< int, std::vector<int> > comp = (*DigiStore)[idigi]->GetDigiCompositionInfo();
      for(size_t i = 0; i < comp.size(); i++){
	G4cout << "tube "  << tubeid
	       << " gate " << i << " p_id";
	for(size_t iv = 0; iv < comp[i].size(); iv++) {
	  G4cout << " " << comp[i][iv];
	}//iv
	G4cout << G4endl;
      }//i
    }
  }//idigi
#endif
}
