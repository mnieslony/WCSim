#include "WCSimWCDAQMessenger.hh"
#include "WCSimEventAction.hh"
#include "WCSimWCDigitizer.hh"
#include "WCSimWCTrigger.hh"

#include "G4UIdirectory.hh"
#include "G4UIcommand.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithAString.hh"

#include <string>

WCSimWCDAQMessenger::WCSimWCDAQMessenger(WCSimEventAction* eventaction) :
  WCSimEvent(eventaction), SetDetectorElement(0)
{
  //G4cout<<"constructing an instance"<<G4endl;
  initialised = false;
}

WCSimWCDAQMessenger* WCSimWCDAQMessenger::iInstance = NULL;

WCSimWCDAQMessenger* WCSimWCDAQMessenger::GetInstance(){
  
  if(iInstance==NULL){
    //G4cout<<"making the static instance"<<G4endl;
    iInstance = new WCSimWCDAQMessenger(0);
  }
  //G4cout<<"returning static instance at "<<iInstance<<G4endl;
  return iInstance;
}

void WCSimWCDAQMessenger::Initialize(G4String detectorElementin)
{
  detectorElement = detectorElementin;
  G4cout<<"Initializing DAQ Messenger for detector element "<<detectorElement<<G4endl;
  initialiseString = " (this is a default set; it may be overwritten by user commands)";
  
  G4String defaultDigitizer = "SKI";
  G4String defaultTrigger = "NDigits";
  bool defaultMultiDigitsPerTrigger = false;
  int defaultDigitizerDeadTime = -99;
  int defaultDigitizerIntegrationWindow = -99;
  bool defaultExtendDigitizerIntegrationWindow = false;
  bool defaultDoPhotonIntegration = true;
  int defaultSaveFailuresTriggerMode = 0;
  double defaultSaveFailuresTriggerTime = 100;
  int defaultSaveFailuresPreTriggerWindow = -400;
  int defaultSaveFailuresPostTriggerWindow = 950;
  int defaultNDigitsTriggerThreshold = -99;
  int defaultNDigitsTriggerWindow = -99;
  bool defaultNDigitsTriggerAdjustForNoise = true;
  int defaultNDigitsPreTriggerWindow = -99;
  int defaultNDigitsPostTriggerWindow = -99;
  bool defaultPromptTrigger = false;
  int defaultPromptPreTriggerWindow = -400;
  int defaultPromptPostTriggerWindow = 950;

  if(not initialised){
    WCSimDAQDir = new G4UIdirectory("/DAQ/");
    WCSimDAQDir->SetGuidance("Commands to select DAQ options");

    DigitizerChoice = new G4UIcmdWithAString("/DAQ/Digitizer", this);
    DigitizerChoice->SetGuidance("Set the Digitizer type");
    DigitizerChoice->SetGuidance("Available choices are:\n"
			       "SKI\n"
			       );
    DigitizerChoice->SetParameterName("Digitizer", false);
    DigitizerChoice->SetCandidates(
				 "SKI "
				 );
    DigitizerChoice->AvailableForStates(G4State_PreInit, G4State_Idle);
    DigitizerChoice->SetDefaultValue(defaultDigitizer);

    TriggerChoice = new G4UIcmdWithAString("/DAQ/Trigger", this);
    TriggerChoice->SetGuidance("Set the Trigger type");
    TriggerChoice->SetGuidance("Available choices are:\n"
			     "NDigits\n"
			     "NDigits2\n"
			     "NoTrigger\n"
			     "TankDigits\n"
			     );
    TriggerChoice->SetParameterName("Trigger", false);
    TriggerChoice->SetCandidates(
			       "NDigits "
			       "NDigits2 "
			       "NoTrigger "
			       "TankDigits "
			       );
    TriggerChoice->AvailableForStates(G4State_PreInit, G4State_Idle);
    TriggerChoice->SetDefaultValue(defaultTrigger);

    MultiDigitsPerTrigger = new G4UIcmdWithABool("/DAQ/MultiDigitsPerTrigger", this);
    MultiDigitsPerTrigger->SetGuidance("Allow the number of digits per PMT per trigger to be > 1?");
    MultiDigitsPerTrigger->SetParameterName("MultiDigitsPerTrigger",true);
    MultiDigitsPerTrigger->SetDefaultValue(defaultMultiDigitsPerTrigger);

  //don't SetNewValue -> defaults class-specific and taken from GetDefault*()


    //Generic digitizer specific options
    DigitizerDir = new G4UIdirectory("/DAQ/DigitizerOpt/");
    DigitizerDir->SetGuidance("Generic commands for digitizers");

    DigitizerDeadTime = new G4UIcmdWithAnInteger("/DAQ/DigitizerOpt/DeadTime", this);
    DigitizerDeadTime->SetGuidance("The deadtime for the digitizer (in ns)");
    DigitizerDeadTime->SetParameterName("DigitizerDeadTime",true);
    DigitizerDeadTime->SetDefaultValue(defaultDigitizerDeadTime);
    //don't SetNewValue -> defaults class-specific and taken from GetDefault*()

    DigitizerIntegrationWindow = new G4UIcmdWithAnInteger("/DAQ/DigitizerOpt/IntegrationWindow", this);
    DigitizerIntegrationWindow->SetGuidance("The integration window for the digitizer (in ns)");
    DigitizerIntegrationWindow->SetParameterName("DigitizerIntegrationWindow",true);
    DigitizerIntegrationWindow->SetDefaultValue(defaultDigitizerIntegrationWindow);
    //don't SetNewValue -> defaults class-specific and taken from GetDefault*()

    ExtendDigitizerIntegrationWindow = new G4UIcmdWithABool("/DAQ/DigitizerOpt/ExtendIntegrationWindow", this);
    ExtendDigitizerIntegrationWindow->SetGuidance("Extend the digitizer integration window when new hits arrive in an existing one");
    ExtendDigitizerIntegrationWindow->SetParameterName("ExtendDigitizerIntegrationWindow",true);
    ExtendDigitizerIntegrationWindow->SetDefaultValue(defaultExtendDigitizerIntegrationWindow);

    DoPhotonIntegration = new G4UIcmdWithABool("/DAQ/DigitizerOpt/DoPhotonIntegration", this);
    DoPhotonIntegration->SetGuidance("Integrate photons within a window into one pulse (true) or convert each photon into a digit (false)");
    DoPhotonIntegration->SetParameterName("DoPhotonIntegration",true);
    DoPhotonIntegration->SetDefaultValue(defaultDoPhotonIntegration);

    //Save failure trigger specific options
    SaveFailuresTriggerDir = new G4UIdirectory("/DAQ/TriggerSaveFailures/");
    SaveFailuresTriggerDir->SetGuidance("Commands specific to the Save Failures trigger");

    SaveFailuresTriggerMode = new G4UIcmdWithAnInteger("/DAQ/TriggerSaveFailures/Mode", this);
    SaveFailuresTriggerMode->SetGuidance("0: save only triggered events; 1: save both triggered and failed events; 2: save only failed events");
    SaveFailuresTriggerMode->SetParameterName("SaveFailuresMode",true);
    SaveFailuresTriggerMode->SetDefaultValue(defaultSaveFailuresTriggerMode);

    SaveFailuresTriggerTime = new G4UIcmdWithADouble("/DAQ/TriggerSaveFailures/TriggerTime", this);
    SaveFailuresTriggerTime->SetGuidance("The trigger time for the events which failed other triggers");
    SaveFailuresTriggerTime->SetParameterName("SaveFailuresTime",true);
    SaveFailuresTriggerTime->SetDefaultValue(defaultSaveFailuresTriggerTime);

    SaveFailuresPreTriggerWindow = new G4UIcmdWithAnInteger("/DAQ/TriggerSaveFailures/PreTriggerWindow", this);
    SaveFailuresPreTriggerWindow->SetGuidance("Set the SaveFailures pretrigger window (in ns)");
    SaveFailuresPreTriggerWindow->SetParameterName("SaveFailuresPreTriggerWindow",false);
    SaveFailuresPreTriggerWindow->SetDefaultValue(defaultSaveFailuresPreTriggerWindow);

    SaveFailuresPostTriggerWindow = new G4UIcmdWithAnInteger("/DAQ/TriggerSaveFailures/PostTriggerWindow", this);
    SaveFailuresPostTriggerWindow->SetGuidance("Set the SaveFailures posttrigger window (in ns)");
    SaveFailuresPostTriggerWindow->SetParameterName("SaveFailuresPostTriggerWindow",false);
    SaveFailuresPostTriggerWindow->SetDefaultValue(defaultSaveFailuresPostTriggerWindow);


    //NDigits trigger specifc options
    NDigitsTriggerDir = new G4UIdirectory("/DAQ/TriggerNDigits/");
    NDigitsTriggerDir->SetGuidance("Commands specific to the NDigits trigger");
    
    NDigitsTriggerThreshold = new G4UIcmdWithAnInteger("/DAQ/TriggerNDigits/Threshold", this);
    NDigitsTriggerThreshold->SetGuidance("Set the NDigits trigger threshold");
    NDigitsTriggerThreshold->SetParameterName("NDigitsThreshold",false);
    NDigitsTriggerThreshold->SetDefaultValue(defaultNDigitsTriggerThreshold);
    //don't SetNewValue -> defaults class-specific and taken from GetDefault*()

    NDigitsTriggerWindow = new G4UIcmdWithAnInteger("/DAQ/TriggerNDigits/Window", this);
    NDigitsTriggerWindow->SetGuidance("Set the NDigits trigger window (in ns)");
    NDigitsTriggerWindow->SetParameterName("NDigitsWindow",false);
    NDigitsTriggerWindow->SetDefaultValue(defaultNDigitsTriggerWindow);
    //don't SetNewValue -> defaults class-specific and taken from GetDefault*()

    NDigitsTriggerAdjustForNoise = new G4UIcmdWithABool("/DAQ/TriggerNDigits/AdjustForNoise", this);
    NDigitsTriggerAdjustForNoise->SetGuidance("Adjust the NDigits trigger threshold automatically dependent on the average noise rate");
    NDigitsTriggerAdjustForNoise->SetParameterName("NDigitsAdjustForNoise",true);
    NDigitsTriggerAdjustForNoise->SetDefaultValue(defaultNDigitsTriggerAdjustForNoise);

    NDigitsPreTriggerWindow = new G4UIcmdWithAnInteger("/DAQ/TriggerNDigits/PreTriggerWindow", this);
    NDigitsPreTriggerWindow->SetGuidance("Set the NDigits pretrigger window (in ns)");
    NDigitsPreTriggerWindow->SetParameterName("NDigitsPreTriggerWindow",false);
    NDigitsPreTriggerWindow->SetDefaultValue(defaultNDigitsPreTriggerWindow);
    //don't SetNewValue -> defaults class-specific and taken from GetDefault*()

    NDigitsPostTriggerWindow = new G4UIcmdWithAnInteger("/DAQ/TriggerNDigits/PostTriggerWindow", this);
    NDigitsPostTriggerWindow->SetGuidance("Set the NDigits posttrigger window (in ns)");
    NDigitsPostTriggerWindow->SetParameterName("NDigitsPostTriggerWindow",false);
    NDigitsPostTriggerWindow->SetDefaultValue(defaultNDigitsPostTriggerWindow);
    //don't SetNewValue -> defaults class-specific and taken from GetDefault*()


    // Prompt (beam) Trigger specfic options
    PromptTriggerEnable = new G4UIcmdWithABool("/DAQ/PromptTrigger/Enable", this);
    PromptTriggerEnable->SetGuidance("Record all digits in a prompt window at the beginning of the event");
    PromptTriggerEnable->SetParameterName("PromptTriggerEnable",true);
    PromptTriggerEnable->SetDefaultValue(defaultPromptTrigger);

/// removed: this doesn't really make sense if the trigger is prompt; i.e. occurs at time 0
//    PromptPreTriggerWindow = new G4UIcmdWithAnInteger("/DAQ/PromptTrigger/PreTriggerWindow", this);
//    PromptPreTriggerWindow->SetGuidance("Set the prompt pretrigger window (in ns)");
//    PromptPreTriggerWindow->SetParameterName("PromptPreTriggerWindow",false);
//    PromptPreTriggerWindow->SetDefaultValue(defaultPromptPreTriggerWindow);
//    //don't SetNewValue -> defaults class-specific and taken from GetDefault*()

    PromptPostTriggerWindow = new G4UIcmdWithAnInteger("/DAQ/PromptTrigger/PostTriggerWindow", this);
    PromptPostTriggerWindow->SetGuidance("Set the prompt posttrigger window (in ns)");
    PromptPostTriggerWindow->SetParameterName("PromptPostTriggerWindow",false);
    PromptPostTriggerWindow->SetDefaultValue(defaultPromptPostTriggerWindow);
    //don't SetNewValue -> defaults class-specific and taken from GetDefault*()


    SetDetectorElement = new G4UIcmdWithAString("/DAQ/SetDetectorElement", this);
    SetDetectorElement->SetParameterName("detectorElement", false);
    SetDetectorElement->AvailableForStates(G4State_PreInit, G4State_Idle);
    SetDetectorElement->SetDefaultValue(detectorElement);
    SetDetectorElement->SetGuidance("Set the detector element for which DAQ settings will be applied");
    SetDetectorElement->SetGuidance("Available choices are:\n");
    
    initialised = true;
  }
  
  // 'SetGuidance' appends to the existing guidance
  SetDetectorElement->SetGuidance(G4String(detectorElement+"\n"));
  // 'SetCandidates' does not
  G4String thecandidates= SetDetectorElement->GetParameter(0)->GetParameterCandidates();
  thecandidates+=" "; thecandidates+=detectorElement; thecandidates+=" ";
  SetDetectorElement->SetCandidates(thecandidates);
  
  // set the silent digitizer options
  OptionsStore theseOptions;
  theseOptions.StoreDigitizerDeadTime = defaultDigitizerDeadTime;
  theseOptions.StoreDigitizerIntegrationWindow = defaultDigitizerIntegrationWindow;
  theseOptions.StoreExtendDigitizerIntegrationWindow = defaultExtendDigitizerIntegrationWindow;
  theseOptions.StoreDoPhotonIntegration = defaultDoPhotonIntegration;
  theseOptions.StoreDigitizerChoice = defaultDigitizer;
  theseOptions.StorePromptTrigger = defaultPromptTrigger;
  theseOptions.StorePromptPreWindow = defaultPromptPreTriggerWindow;
  theseOptions.StorePromptPostWindow = defaultPromptPostTriggerWindow;
  // set the silent trigger options
  theseOptions.StoreTriggerChoice = defaultTrigger;
  theseOptions.StoreMultiDigitsPerTrigger = defaultMultiDigitsPerTrigger;
  theseOptions.MultiDigitsPerTriggerSet = false; //this variable is bool & defaults are class specfic; use this to know if the default is overidden
  theseOptions.StoreNDigitsThreshold = defaultNDigitsTriggerThreshold;
  theseOptions.StoreNDigitsWindow = defaultNDigitsTriggerWindow;
  theseOptions.StoreNDigitsPreWindow = defaultNDigitsPreTriggerWindow;
  theseOptions.StoreNDigitsPostWindow = defaultNDigitsPostTriggerWindow;
  theseOptions.StoreSaveFailuresMode = defaultSaveFailuresTriggerMode;
  theseOptions.StoreSaveFailuresTime = defaultSaveFailuresTriggerTime;
  theseOptions.StoreSaveFailuresPreWindow = defaultSaveFailuresPreTriggerWindow;
  theseOptions.StoreSaveFailuresPostWindow = defaultSaveFailuresPostTriggerWindow;
  theseOptions.StoreNDigitsAdjustForNoise = defaultNDigitsTriggerAdjustForNoise;
  
  // Create the trigger + digitizer options in the map
  StoredOptions.emplace(detectorElement,theseOptions);
  
  // Set the verbose options: the 'SetNewValue' calls set the Stored values,
  // (which could more directly be set above) but also result in printout to the std::out
  // for brevity, we can comment these out if we don't want the output.
  SetNewValue(DigitizerChoice, defaultDigitizer);
  SetNewValue(TriggerChoice, defaultTrigger);
  SetNewValue(SaveFailuresTriggerMode, G4UIcommand::ConvertToString(defaultSaveFailuresTriggerMode));
  SetNewValue(SaveFailuresTriggerTime, G4UIcommand::ConvertToString(defaultSaveFailuresTriggerTime));
  SetNewValue(SaveFailuresPreTriggerWindow, G4UIcommand::ConvertToString(defaultSaveFailuresPreTriggerWindow));
  SetNewValue(SaveFailuresPostTriggerWindow, G4UIcommand::ConvertToString(defaultSaveFailuresPostTriggerWindow));
  SetNewValue(NDigitsTriggerAdjustForNoise, G4UIcommand::ConvertToString(defaultNDigitsTriggerAdjustForNoise));
  
  
  initialiseString = "";
}

WCSimWCDAQMessenger::~WCSimWCDAQMessenger()
{
  if(StoredOptions.size()){
    //G4cout<<"Cleaning up DAQ Messenger"<<G4endl;
    while(StoredOptions.size()){
      G4String nextelement = StoredOptions.begin()->first;
      //G4cout<<"Removing instance for detectorElement "<<nextelement<<G4endl;
      this->RemoveDAQMessengerInstance(nextelement);
    }
  }
  
  // clean up UI stuff
  if(SaveFailuresTriggerDir){ delete SaveFailuresTriggerDir; SaveFailuresTriggerDir=nullptr; }
  if(SaveFailuresTriggerMode){ delete SaveFailuresTriggerMode; SaveFailuresTriggerMode=nullptr; }
  if(SaveFailuresTriggerTime){ delete SaveFailuresTriggerTime; SaveFailuresTriggerTime=nullptr; }
  if(SaveFailuresPreTriggerWindow){ delete SaveFailuresPreTriggerWindow; SaveFailuresPreTriggerWindow=nullptr; }
  if(SaveFailuresPostTriggerWindow){ delete SaveFailuresPostTriggerWindow; SaveFailuresPostTriggerWindow=nullptr; }

  if(NDigitsTriggerDir){ delete NDigitsTriggerDir; NDigitsTriggerDir=nullptr; }
  if(NDigitsTriggerThreshold){ delete NDigitsTriggerThreshold; NDigitsTriggerThreshold=nullptr; }
  if(NDigitsTriggerWindow){ delete NDigitsTriggerWindow; NDigitsTriggerWindow=nullptr; }
  if(NDigitsTriggerAdjustForNoise){ delete NDigitsTriggerAdjustForNoise; NDigitsTriggerAdjustForNoise=nullptr; }
  if(NDigitsPreTriggerWindow){ delete NDigitsPreTriggerWindow; NDigitsPreTriggerWindow=nullptr; }
  if(NDigitsPostTriggerWindow){ delete NDigitsPostTriggerWindow; NDigitsPostTriggerWindow=nullptr; }

  if(PromptTriggerEnable){ delete PromptTriggerEnable; PromptTriggerEnable=nullptr; }
//  if(PromptPreTriggerWindow){ delete PromptPreTriggerWindow; PromptPreTriggerWindow=nullptr; }
  if(PromptPostTriggerWindow){ delete PromptPostTriggerWindow; PromptPostTriggerWindow=nullptr; }

  if(DigitizerDir){ delete DigitizerDir; DigitizerDir=nullptr; }
  if(DigitizerDeadTime){ delete DigitizerDeadTime; DigitizerDeadTime=nullptr; }
  if(DigitizerIntegrationWindow){ delete DigitizerIntegrationWindow; DigitizerIntegrationWindow=nullptr; }
  if(ExtendDigitizerIntegrationWindow){ delete ExtendDigitizerIntegrationWindow; ExtendDigitizerIntegrationWindow=nullptr; }
  if(DoPhotonIntegration){ delete DoPhotonIntegration; DoPhotonIntegration=nullptr; }

  if(DigitizerChoice){ delete DigitizerChoice; DigitizerChoice=nullptr; }
  if(TriggerChoice){ delete TriggerChoice; TriggerChoice=nullptr; }
  if(MultiDigitsPerTrigger){ delete MultiDigitsPerTrigger; MultiDigitsPerTrigger=nullptr; }
  if(WCSimDAQDir){ delete WCSimDAQDir; WCSimDAQDir=nullptr; }
  
  if(SetDetectorElement){ delete SetDetectorElement; SetDetectorElement=nullptr; }
  
}

void WCSimWCDAQMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{
  //Because this Messenger class contains options for classes that don't exist when options are
  // read in (Trigger and Digitizer class options) we need to store each options' value
  // for use in the Tell*() methods later

  if (command == DigitizerChoice) {
    G4cout << "Digitizer choice set to " << newValue << initialiseString.c_str() << G4endl;
    WCSimEvent->SetDigitizerChoice(newValue, detectorElement);
    StoredOptions.at(detectorElement).StoreDigitizerChoice = newValue;
  }
  else if (command == TriggerChoice) {
    G4cout << "Trigger choice set to " << newValue << initialiseString.c_str() << G4endl;
    WCSimEvent->SetTriggerChoice(newValue, detectorElement);
    StoredOptions.at(detectorElement).StoreTriggerChoice = newValue;
  }
  else if (command == MultiDigitsPerTrigger) {
    StoredOptions.at(detectorElement).StoreMultiDigitsPerTrigger = MultiDigitsPerTrigger->GetNewBoolValue(newValue);
    if(!StoredOptions.at(detectorElement).StoreMultiDigitsPerTrigger)
      G4cout << "Will restrict number of digits per PMT per trigger to <= 1" << initialiseString.c_str() << G4endl;
    else
      G4cout << "Will allow number of digits per PMT per trigger to go > 1" << initialiseString.c_str() << G4endl;
    if(initialised)
      StoredOptions.at(detectorElement).MultiDigitsPerTriggerSet = true;
  }

  //Generic digitizer options
  else if (command == DigitizerDeadTime) {
    G4cout << "Digitizer deadtime set to " << newValue << " ns" << initialiseString.c_str() << G4endl;
    StoredOptions.at(detectorElement).StoreDigitizerDeadTime = DigitizerDeadTime->GetNewIntValue(newValue);
  }
  else if (command == DigitizerIntegrationWindow) {
    G4cout << "Digitizer integration window set to " << newValue << " ns" << initialiseString.c_str() << G4endl;
    StoredOptions.at(detectorElement).StoreDigitizerIntegrationWindow = 
         DigitizerIntegrationWindow->GetNewIntValue(newValue);
  }
  else if (command == ExtendDigitizerIntegrationWindow) {
    G4cout<< "Digitizer integration window set to " << ((newValue=="true") ? "" : "not " ) 
          << "extend when new hits arrive in an existing window" << initialiseString.c_str() << G4endl;
    StoredOptions.at(detectorElement).StoreExtendDigitizerIntegrationWindow = 
         ExtendDigitizerIntegrationWindow->GetNewBoolValue(newValue);
  }
  else if (command == DoPhotonIntegration) {
    G4cout<< "Digitizer set to ";
    if(newValue=="true"){
      G4cout << "integrate photons within the digitizer integration window into Digits ";
    } else {
      G4cout << "convert individual photons into Digits ";
    } 
    G4cout << initialiseString.c_str() << G4endl;
    StoredOptions.at(detectorElement).StoreDoPhotonIntegration = 
         DoPhotonIntegration->GetNewBoolValue(newValue);
  }

  //Save failures "trigger"
  else if (command == SaveFailuresTriggerMode) {
    StoredOptions.at(detectorElement).StoreSaveFailuresMode = SaveFailuresTriggerMode->GetNewIntValue(newValue);
    std::string failuremode;
    if(StoredOptions.at(detectorElement).StoreSaveFailuresMode == 0)
      failuremode = "Saving only triggered events";
    else if(StoredOptions.at(detectorElement).StoreSaveFailuresMode == 1)
      failuremode = "Saving both triggered and failed events";
    else if(StoredOptions.at(detectorElement).StoreSaveFailuresMode == 2)
      failuremode = "Saving only failed events";
    else {
      G4cerr << "Unknown value of /DAQ/TriggerSaveFailures/Mode " << 
        StoredOptions.at(detectorElement).StoreSaveFailuresMode << " Exiting..." << G4endl;
      exit(-1);
    }
    G4cout << failuremode << initialiseString.c_str() << G4endl;
  }
  else if (command == SaveFailuresTriggerTime) {
    G4cout << "Trigger time for events which fail all triggers will be set to " << newValue << " ns" << initialiseString.c_str() << G4endl;
    StoredOptions.at(detectorElement).StoreSaveFailuresTime = SaveFailuresTriggerTime->GetNewDoubleValue(newValue);
  }
  else if (command == SaveFailuresPreTriggerWindow) {
    G4cout << "SaveFailures pretrigger window set to " << newValue << " ns" << initialiseString.c_str() << G4endl;
    StoredOptions.at(detectorElement).StoreSaveFailuresPreWindow = SaveFailuresPreTriggerWindow->GetNewIntValue(newValue);
  }
  else if (command == SaveFailuresPostTriggerWindow) {
    G4cout << "SaveFailures posttrigger window set to " << newValue << " ns" << initialiseString.c_str() << G4endl;
    StoredOptions.at(detectorElement).StoreSaveFailuresPostWindow = SaveFailuresPostTriggerWindow->GetNewIntValue(newValue);
  }

  //NDigits trigger
  else if (command == NDigitsTriggerThreshold) {
    G4cout << "NDigits trigger threshold set to " << newValue << initialiseString.c_str() << G4endl;
    StoredOptions.at(detectorElement).StoreNDigitsThreshold = NDigitsTriggerThreshold->GetNewIntValue(newValue);
  }
  else if (command == NDigitsTriggerAdjustForNoise) {
    StoredOptions.at(detectorElement).StoreNDigitsAdjustForNoise = NDigitsTriggerAdjustForNoise->GetNewBoolValue(newValue);
    if(StoredOptions.at(detectorElement).StoreNDigitsAdjustForNoise)
      G4cout << "Will adjust NDigits trigger threshold using average dark noise rate" << initialiseString.c_str() << G4endl;
  }
  else if (command == NDigitsTriggerWindow) {
    G4cout << "NDigits trigger window set to " << newValue << " ns" << initialiseString.c_str() << G4endl;
    StoredOptions.at(detectorElement).StoreNDigitsWindow = NDigitsTriggerWindow->GetNewIntValue(newValue);
  }
  else if (command == NDigitsPreTriggerWindow) {
    G4cout << "NDigits pretrigger window set to " << newValue << " ns" << initialiseString.c_str() << G4endl;
    StoredOptions.at(detectorElement).StoreNDigitsPreWindow = NDigitsPreTriggerWindow->GetNewIntValue(newValue);
  }
  else if (command == NDigitsPostTriggerWindow) {
    G4cout << "NDigits posttrigger window set to " << newValue << " ns" << initialiseString.c_str() << G4endl;
    StoredOptions.at(detectorElement).StoreNDigitsPostWindow = NDigitsPostTriggerWindow->GetNewIntValue(newValue);
  }
  
  //Prompt trigger
  else if (command == PromptTriggerEnable) {
    StoredOptions.at(detectorElement).StorePromptTrigger = PromptTriggerEnable->GetNewBoolValue(newValue);
    G4cout << "Will " << ((StoredOptions.at(detectorElement).StorePromptTrigger) ? "" : "not ") 
           << "add a prompt trigger window at the start of the event" << initialiseString.c_str() << G4endl;
  }
//  else if (command == PromptPreTriggerWindow) {
//    G4cout << "Prompt pretrigger window set to " << newValue << " ns" << initialiseString.c_str() << G4endl;
//    StoredOptions.at(detectorElement).StorePromptPreWindow = PromptPreTriggerWindow->GetNewIntValue(newValue);
//  }
  else if (command == PromptPostTriggerWindow) {
    G4cout << "Prompt posttrigger window set to " << newValue << " ns" << initialiseString.c_str() << G4endl;
    StoredOptions.at(detectorElement).StorePromptPostWindow = PromptPostTriggerWindow->GetNewIntValue(newValue);
  }
  
  // Detector Element to handle multiple digitizer and trigger classes
  else if(command == SetDetectorElement){
    if(StoredOptions.count(newValue)>0){
      detectorElement=newValue;
      G4cout << "Setting detectorElement value " << newValue << initialiseString.c_str() << G4endl;
    }
  }
}

void WCSimWCDAQMessenger::SetTriggerOptions(G4String detectorElementin)
{
  detectorElement = detectorElementin;
  G4cout << "Passing Trigger options to the trigger class instance " << detectorElement << G4endl;

  if(StoredOptions.at(detectorElement).MultiDigitsPerTriggerSet) {
    WCSimTrigger->SetMultiDigitsPerTrigger(StoredOptions.at(detectorElement).StoreMultiDigitsPerTrigger);
    if(!StoredOptions.at(detectorElement).StoreMultiDigitsPerTrigger)
      G4cout << "\tWill restrict number of digits per PMT per trigger to <= 1" << G4endl;
    else
      G4cout << "\tWill allow number of digits per PMT per trigger to go > 1" << initialiseString.c_str() << G4endl;
  }

  WCSimTrigger->SetSaveFailuresMode(StoredOptions.at(detectorElement).StoreSaveFailuresMode);
  std::string failuremode;
  if(StoredOptions.at(detectorElement).StoreSaveFailuresMode == 0)
    failuremode = "Saving only triggered events";
  else if(StoredOptions.at(detectorElement).StoreSaveFailuresMode == 1)
    failuremode = "Saving both triggered and failed events";
  else if(StoredOptions.at(detectorElement).StoreSaveFailuresMode == 2)
    failuremode = "Saving only failed events";
  G4cout << "\t" << failuremode << G4endl;
  WCSimTrigger->SetSaveFailuresTime(StoredOptions.at(detectorElement).StoreSaveFailuresTime);
  G4cout << "\tTrigger time for events which fail all triggers will be set to " 
         << StoredOptions.at(detectorElement).StoreSaveFailuresTime << " ns" << G4endl;
  if(StoredOptions.at(detectorElement).StoreSaveFailuresPreWindow >= -1E6) {
    WCSimTrigger->SetSaveFailuresPreTriggerWindow(StoredOptions.at(detectorElement).StoreSaveFailuresPreWindow);
    G4cout << "\tSaveFailures pretrigger window set to " 
           << StoredOptions.at(detectorElement).StoreSaveFailuresPreWindow << " ns" << G4endl;
  }
  if(StoredOptions.at(detectorElement).StoreSaveFailuresPostWindow >= 0) {
    WCSimTrigger->SetSaveFailuresPostTriggerWindow(StoredOptions.at(detectorElement).StoreSaveFailuresPostWindow);
    G4cout << "\tSaveFailures posttrigger window set to " 
           << StoredOptions.at(detectorElement).StoreSaveFailuresPostWindow << " ns" << G4endl;
  }

  if(StoredOptions.at(detectorElement).StoreNDigitsThreshold >= 0) {
   WCSimTrigger->SetNDigitsThreshold(StoredOptions.at(detectorElement).StoreNDigitsThreshold);
    G4cout << "\tNDigits trigger threshold set to " 
           << StoredOptions.at(detectorElement).StoreNDigitsThreshold << G4endl;
  }
  WCSimTrigger->SetNDigitsAdjustForNoise(StoredOptions.at(detectorElement).StoreNDigitsAdjustForNoise);
  if(StoredOptions.at(detectorElement).StoreNDigitsAdjustForNoise)
    G4cout << "\tWill adjust NDigits trigger threshold using average dark noise rate" << G4endl;
  if(StoredOptions.at(detectorElement).StoreNDigitsWindow >= 0) {
    WCSimTrigger->SetNDigitsWindow(StoredOptions.at(detectorElement).StoreNDigitsWindow);
    G4cout << "\tNDigits trigger window set to " 
           << StoredOptions.at(detectorElement).StoreNDigitsWindow << " ns" << G4endl;
  }
  if(StoredOptions.at(detectorElement).StoreNDigitsPreWindow >= 0) {
    WCSimTrigger->SetNDigitsPreTriggerWindow(StoredOptions.at(detectorElement).StoreNDigitsPreWindow);
    G4cout << "\tNDigits pretrigger window set to " 
           << StoredOptions.at(detectorElement).StoreNDigitsPreWindow << " ns" << G4endl;
  }
  if(StoredOptions.at(detectorElement).StoreNDigitsPostWindow >= 0) {
    WCSimTrigger->SetNDigitsPostTriggerWindow(StoredOptions.at(detectorElement).StoreNDigitsPostWindow);
    G4cout << "\tNDigits posttrigger window set to " 
           << StoredOptions.at(detectorElement).StoreNDigitsPostWindow << " ns" << G4endl;
  }

  WCSimTrigger->SetRecordPromptWindow(StoredOptions.at(detectorElement).StorePromptTrigger);
  G4cout << "\tWill "<<((StoredOptions.at(detectorElement).StorePromptTrigger) ? "" : "not ") 
         << "add a prompt trigger window at the start of the event" << G4endl;
  if(StoredOptions.at(detectorElement).StoreNDigitsPreWindow >= 0) {
    WCSimTrigger->SetPromptPreTriggerWindow(StoredOptions.at(detectorElement).StorePromptPreWindow);
    G4cout << "\tPrompt pretrigger window set to " 
           << StoredOptions.at(detectorElement).StorePromptPreWindow << " ns" << G4endl;
  }
  if(StoredOptions.at(detectorElement).StoreNDigitsPostWindow >= 0) {
    WCSimTrigger->SetPromptPostTriggerWindow(StoredOptions.at(detectorElement).StorePromptPostWindow);
    G4cout << "\tPrompt posttrigger window set to " 
           << StoredOptions.at(detectorElement).StorePromptPostWindow << " ns" << G4endl;
  }
}

void WCSimWCDAQMessenger::SetDigitizerOptions(G4String detectorElementin)
{
  detectorElement = detectorElementin;
  G4cout << "Passing Digitizer options to the digitizer class instance " << detectorElement << G4endl;
  if(StoredOptions.at(detectorElement).StoreDigitizerDeadTime >= 0) {
    WCSimDigitize->SetDigitizerDeadTime(StoredOptions.at(detectorElement).StoreDigitizerDeadTime);
    G4cout << "\tDigitizer deadtime set to " 
           << StoredOptions.at(detectorElement).StoreDigitizerDeadTime << " ns"  << G4endl;
  }
  if(StoredOptions.at(detectorElement).StoreDigitizerIntegrationWindow >= 0) {
    WCSimDigitize->SetDigitizerIntegrationWindow(StoredOptions.at(detectorElement).StoreDigitizerIntegrationWindow);
    G4cout << "\tDigitizer integration window set to " 
           << StoredOptions.at(detectorElement).StoreDigitizerIntegrationWindow << " ns" << G4endl;
  }
  WCSimDigitize->SetExtendIntegrationWindow(StoredOptions.at(detectorElement).StoreExtendDigitizerIntegrationWindow);
  G4cout<<"\tDigitizer integration window set to " << ((StoredOptions.at(detectorElement).StoreExtendDigitizerIntegrationWindow) ? "" : "not") << " extend when new hits arrive in an existing window" << G4endl;
  WCSimDigitize->SetDoPhotonIntegration(StoredOptions.at(detectorElement).StoreDoPhotonIntegration);
  G4cout<<"\tDigitizer set to ";
  if(StoredOptions.at(detectorElement).StoreExtendDigitizerIntegrationWindow){
    G4cout << "integrate photons on the same PMT close in time into digits" << G4endl;
  } else {
    G4cout << "convert each individual photon into its own digit" << G4endl;
  }
}

void WCSimWCDAQMessenger::AddDAQMessengerInstance(G4String detectorElementin){
  if(StoredOptions.count(detectorElementin)){
    G4cerr<<"Attempt to re-add existing DAQ messenger "<<detectorElementin<<G4endl;
    return;
  } else {
    //G4cout<<"AddDAQMessengerInstance for detectorElement "<<detectorElementin<<G4endl;
    Initialize(detectorElementin);
  }
}

void WCSimWCDAQMessenger::RemoveDAQMessengerInstance(G4String detectorElementin){
  if(StoredOptions.count(detectorElementin)){
    //G4cout<<"Removing DAQ Messenger instance " << detectorElementin << G4endl;
    StoredOptions.erase(detectorElementin);
  } else {
    G4cerr<<"Attempt to remove nonexistant element "<<detectorElementin<<" from WCDAQMessenger!"<<G4endl;
  }
}
