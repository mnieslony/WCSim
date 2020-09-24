#ifndef WCSimPrimaryGeneratorMessenger_h
#define WCSimPrimaryGeneratorMessenger_h 1

class WCSimPrimaryGeneratorAction;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWith3Vector;
class G4UIcmdWith3VectorAndUnit;
class G4UIcmdWithADoubleAndUnit;

#include "G4UImessenger.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "globals.hh"
#include "G4Tokenizer.hh"

class WCSimPrimaryGeneratorMessenger: public G4UImessenger
{
 public:
  WCSimPrimaryGeneratorMessenger(WCSimPrimaryGeneratorAction* mpga);
  ~WCSimPrimaryGeneratorMessenger();
  
 public:
  void     SetNewValue(G4UIcommand* command, G4String newValues);
  G4String GetCurrentValue(G4UIcommand* command);
  
 private:
  WCSimPrimaryGeneratorAction* myAction;
  
 private: //commands
  G4UIdirectory*      mydetDirectory;
  G4UIcmdWithAString* genCmd;
  G4UIcmdWithAString* fileNameCmd;
  G4UIcmdWithAString* spectrumFileCmd;
  G4UIcmdWithAString* timeUnitCmd;
  G4UIcmdWithAString* isotopeCmd;
  G4UIcmdWithAString* radonScalingCmd;
  G4UIcmdWithADouble* radioactive_time_window_Cmd;
  G4UIcmdWithAnInteger* radonGeoSymCmd;
  G4UIcmdWithAString* geniefileDirectoryCmd;
  G4UIcmdWithAString* talysfileDirectoryCmd;
  G4UIcmdWithAnInteger* primariesStartEventCmd;

  G4UIcmdWith3VectorAndUnit* positionCmd;
  G4UIcmdWithADoubleAndUnit* radiusCmd;
  G4UIcmdWithADoubleAndUnit* heightCmd;
  G4UIcmdWith3Vector* rot1Cmd;
  G4UIcmdWith3Vector* rot2Cmd;
  
  void IsotopeCommand(G4String newValue);
  void RadonScalingCommand(G4String newValue);
};

#endif


