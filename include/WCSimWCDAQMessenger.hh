#ifndef WCSimWCDAQMessenger_h
#define WCSimWCDAQMessenger_h 1

#include "G4UImessenger.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithoutParameter.hh"

#include <map>

class G4UIcommand;
class WCSimEventAction;
class WCSimWCDigitizerBase;
class WCSimWCTriggerBase;

struct OptionsStore{
  // Trigger options
  G4int                 StoreNDigitsPostWindow;
  G4int                 StoreNDigitsPreWindow;
  G4bool                StoreNDigitsAdjustForNoise;
  G4int                 StoreNDigitsWindow;
  G4int                 StoreNDigitsThreshold;
  G4int                 StoreSaveFailuresPostWindow;
  G4int                 StoreSaveFailuresPreWindow;
  G4double              StoreSaveFailuresTime;
  G4int                 StoreSaveFailuresMode;
  G4String              StoreTriggerChoice;
  G4bool                StoreMultiDigitsPerTrigger;
  G4bool                MultiDigitsPerTriggerSet;
  // Digitizer options
  G4int                 StoreDigitizerIntegrationWindow;
  G4int                 StoreDigitizerDeadTime;
  G4String              StoreDigitizerChoice;
  G4bool                StoreExtendDigitizerIntegrationWindow;
};

class WCSimWCDAQMessenger: public G4UImessenger
{
public:
  static WCSimWCDAQMessenger* iInstance;
  WCSimWCDAQMessenger(WCSimEventAction*);
  
  static WCSimWCDAQMessenger* GetInstance();
  void AddDAQMessengerInstance(G4String detectorElement);
  void RemoveDAQMessengerInstance(G4String detectorElement);
  void SetDigitizerInstance(WCSimWCDigitizerBase* digitizerpoint, G4String detectorElement);
  void SetTriggerInstance(WCSimWCTriggerBase* triggerpoint, G4String detectorElement);
  void Initialize(G4String detectorElementin);
  G4String GetDetectorElement(){ return detectorElement; }

  ~WCSimWCDAQMessenger();

  void SetNewValue(G4UIcommand* command, G4String newValue);

  void SetTriggerOptions(G4String detectorElementin);
  void SetDigitizerOptions(G4String detectorElementin);

  void TellMeAboutTheDigitizer  (WCSimWCDigitizerBase* digitizer)   { WCSimDigitize = digitizer; }
  void TellMeAboutTheTrigger    (WCSimWCTriggerBase*   trigger)     { WCSimTrigger  = trigger; }

private:
  WCSimEventAction*     WCSimEvent;
  WCSimWCDigitizerBase* WCSimDigitize;
  WCSimWCTriggerBase*   WCSimTrigger;

  G4UIdirectory*        WCSimDAQDir;
  G4UIcmdWithAString*   DigitizerChoice;
  G4UIcmdWithAString*   TriggerChoice;
  G4UIcmdWithABool*     MultiDigitsPerTrigger;
  G4UIcmdWithABool*     ExtendDigitizerIntegrationWindow;
  G4bool                MultiDigitsPerTriggerSet;

  G4UIdirectory*        DigitizerDir;
  G4UIcmdWithAnInteger* DigitizerDeadTime;
  G4UIcmdWithAnInteger* DigitizerIntegrationWindow;

  G4UIdirectory*        SaveFailuresTriggerDir;
  G4UIcmdWithAnInteger* SaveFailuresTriggerMode;
  G4UIcmdWithADouble*   SaveFailuresTriggerTime;
  G4UIcmdWithAnInteger* SaveFailuresPreTriggerWindow;
  G4UIcmdWithAnInteger* SaveFailuresPostTriggerWindow;

  G4UIdirectory*        NDigitsTriggerDir;
  G4UIcmdWithAnInteger* NDigitsTriggerThreshold;
  G4UIcmdWithAnInteger* NDigitsTriggerWindow;
  G4UIcmdWithABool*     NDigitsTriggerAdjustForNoise;
  G4UIcmdWithAnInteger* NDigitsPreTriggerWindow;
  G4UIcmdWithAnInteger* NDigitsPostTriggerWindow;
  G4UIcmdWithAString*   SetDetectorElement;

  std::map<G4String, OptionsStore> StoredOptions;

  G4String initialiseString;
  G4bool   initialised;
  G4String detectorElement;
  
};

#endif
