#include "WCSimPrimaryGeneratorMessenger.hh"
#include "WCSimPrimaryGeneratorAction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWith3Vector.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4ios.hh"

WCSimPrimaryGeneratorMessenger::WCSimPrimaryGeneratorMessenger(WCSimPrimaryGeneratorAction* pointerToAction)
:myAction(pointerToAction)
{
  mydetDirectory = new G4UIdirectory("/mygen/");
  mydetDirectory->SetGuidance("WCSim detector control commands.");

  genCmd = new G4UIcmdWithAString("/mygen/generator",this);
  genCmd->SetGuidance("Select primary generator.");
  //T. Akiri: Addition of laser, M.Nieslony: Addition of antinu, genie
  genCmd->SetGuidance(" Available generators : muline, gun, laser, gps, beam, antinu, genie");
  genCmd->SetParameterName("generator",true);
  genCmd->SetDefaultValue("beam");	// previously muline
  //T. Akiri: Addition of laser, M.Nieslony: Addition of antinu, genie
  genCmd->SetCandidates("muline gun laser gps beam antinu genie");

  fileNameCmd = new G4UIcmdWithAString("/mygen/vecfile",this);
  fileNameCmd->SetGuidance("Select the file of vectors.");
  fileNameCmd->SetGuidance(" Enter the file name of the vector file");
  fileNameCmd->SetParameterName("fileName",true);
  fileNameCmd->SetDefaultValue("inputvectorfile");

  spectrumFileCmd = new G4UIcmdWithAString("/mygen/Efile",this);
  spectrumFileCmd->SetGuidance("Select the file with energy spectrum.");
  spectrumFileCmd->SetGuidance("Enter the file name of the energy spectrum file.");
  spectrumFileCmd->SetParameterName("spectrumFile",true);
  spectrumFileCmd->SetDefaultValue("inputspectrum"); 
 
  positionCmd = new G4UIcmdWith3VectorAndUnit("/mygen/position",this);
  positionCmd->SetGuidance("Set gun Position");
  positionCmd->SetUnitCategory("Length");
  positionCmd->SetDefaultUnit("cm");
  positionCmd->SetUnitCandidates("mm cm m");

  radiusCmd = new G4UIcmdWithADoubleAndUnit("/mygen/cylradius",this);
  radiusCmd->SetGuidance("Set radius of the cylindrical distribution");
  radiusCmd->SetUnitCategory("Length");
  radiusCmd->SetDefaultUnit("cm");
  radiusCmd->SetUnitCandidates("mm cm m");
 
  heightCmd = new G4UIcmdWithADoubleAndUnit("/mygen/cylhalfz",this);
  heightCmd->SetGuidance("Set half the height of the cylindrical distribution");
  heightCmd->SetUnitCategory("Length");
  heightCmd->SetDefaultUnit("cm");
  heightCmd->SetUnitCandidates("mm cm m");

  rot1Cmd = new G4UIcmdWith3Vector("/mygen/rot1",this);
  rot1Cmd->SetGuidance("Set first rotation of coordinate system (x y z)");
  
  rot2Cmd = new G4UIcmdWith3Vector("/mygen/rot2",this);
  rot2Cmd->SetGuidance("Set second rotation of coordinate system (x y z)");

  primariesfileDirectoryCmd = new G4UIcmdWithAString("/mygen/primariesdirectory", this);
  primariesfileDirectoryCmd->SetGuidance("Specify the directory containing beam primary root files");
  primariesfileDirectoryCmd->SetParameterName("directoryName",true);
  primariesfileDirectoryCmd->SetDefaultValue("");
  
  neutrinosfileDirectoryCmd = new G4UIcmdWithAString("/mygen/neutrinosdirectory", this);
  neutrinosfileDirectoryCmd->SetGuidance("Specify the directory containing genie neutrino root files. Set this before setting the primariesDirectory. Both should be set at the same time.");
  neutrinosfileDirectoryCmd->SetParameterName("directoryName",true);
  neutrinosfileDirectoryCmd->SetDefaultValue("");

  primariesStartEventCmd = new G4UIcmdWithAnInteger("/mygen/primariesoffset", this);
  primariesStartEventCmd->SetGuidance("The starting entry number for reading primaries");
  primariesStartEventCmd->SetParameterName("primariesoffset",true);
  primariesStartEventCmd->SetDefaultValue(0);
  
  geniefileDirectoryCmd = new G4UIcmdWithAString("/mygen/geniedirectory", this);
  geniefileDirectoryCmd->SetGuidance("Specify the directory containing genie (non-beam) root files");
  geniefileDirectoryCmd->SetParameterName("directoryName",true);
  geniefileDirectoryCmd->SetDefaultValue("");
}

WCSimPrimaryGeneratorMessenger::~WCSimPrimaryGeneratorMessenger()
{
  delete genCmd;
  delete fileNameCmd;
  delete primariesfileDirectoryCmd;
  delete neutrinosfileDirectoryCmd;
  delete geniefileDirectoryCmd;
  delete mydetDirectory;
  delete primariesStartEventCmd;

  delete positionCmd;
  delete radiusCmd;
  delete heightCmd;
  delete rot1Cmd;
  delete rot2Cmd;
}

void WCSimPrimaryGeneratorMessenger::SetNewValue(G4UIcommand * command,G4String newValue)
{
  if( command==genCmd )
  {
    if (newValue == "muline")
    {
      G4cout<<"Setting generator source to muline"<<G4endl;
      myAction->SetMulineEvtGenerator(true);
      myAction->SetGunEvtGenerator(false);
      myAction->SetLaserEvtGenerator(false);
      myAction->SetBeamEvtGenerator(false);
      myAction->SetGPSEvtGenerator(false);
      myAction->SetAntiNuEvtGenerator(false);
      myAction->SetGenieEvtGenerator(false);
    }
    else if ( newValue == "gun")
    {
      G4cout<<"Setting generator source to gun"<<G4endl;
      myAction->SetMulineEvtGenerator(false);
      myAction->SetGunEvtGenerator(true);
      myAction->SetLaserEvtGenerator(false);
      myAction->SetBeamEvtGenerator(false);
      myAction->SetGPSEvtGenerator(false);
      myAction->SetAntiNuEvtGenerator(false);
      myAction->SetGenieEvtGenerator(false);
    }
    else if ( newValue == "laser")   //T. Akiri: Addition of laser
    {
      G4cout<<"Setting generator source to laser"<<G4endl;
      myAction->SetMulineEvtGenerator(false);
      myAction->SetGunEvtGenerator(false);
      myAction->SetLaserEvtGenerator(true);
      myAction->SetBeamEvtGenerator(false);
      myAction->SetGPSEvtGenerator(false);
      myAction->SetAntiNuEvtGenerator(false);
      myAction->SetGenieEvtGenerator(false);
    }
    else if ( newValue == "beam")
    {
      G4cout<<"Setting generator source to beam"<<G4endl;
      myAction->SetMulineEvtGenerator(false);
      myAction->SetGunEvtGenerator(false);
      myAction->SetLaserEvtGenerator(false);
      myAction->SetBeamEvtGenerator(true);
      myAction->SetGPSEvtGenerator(false);
      myAction->SetAntiNuEvtGenerator(false);
      myAction->SetGenieEvtGenerator(false);
    }
    else if ( newValue == "gps")
    {
      G4cout<<"Setting generator source to gps"<<G4endl;
      myAction->SetMulineEvtGenerator(false);
      myAction->SetGunEvtGenerator(false);
      myAction->SetLaserEvtGenerator(false);
      myAction->SetBeamEvtGenerator(false);
      myAction->SetGPSEvtGenerator(true);
      myAction->SetAntiNuEvtGenerator(false);
      myAction->SetGenieEvtGenerator(false);
    } 
    else if ( newValue == "antinu")
    {
      G4cout<<"Setting generator source to antinu"<<G4endl;
      myAction->SetMulineEvtGenerator(false);
      myAction->SetGunEvtGenerator(false);
      myAction->SetLaserEvtGenerator(false);
      myAction->SetBeamEvtGenerator(false);
      myAction->SetGPSEvtGenerator(false);
      myAction->SetAntiNuEvtGenerator(true);
      myAction->SetGenieEvtGenerator(false);
    }
    else if ( newValue == "genie")
    {
      G4cout<<"Setting generator source to genie"<<G4endl;
      myAction->SetMulineEvtGenerator(false);
      myAction->SetGunEvtGenerator(false);
      myAction->SetLaserEvtGenerator(false);
      myAction->SetBeamEvtGenerator(false);
      myAction->SetGPSEvtGenerator(false);
      myAction->SetAntiNuEvtGenerator(false);
      myAction->SetGenieEvtGenerator(true);
    }
  }

  if( command == spectrumFileCmd ){
    myAction->OpenSpectrumFile(newValue);
    G4cout << "Input energy spectrum file set to " << newValue << G4endl;
  }

  if( command == fileNameCmd )
  {
    myAction->OpenVectorFile(newValue);
    G4cout << "Input vector file set to " << newValue << G4endl;
  }
  
  if( command == primariesfileDirectoryCmd )
  {
    myAction->SetPrimaryFilesDirectory(newValue);
    myAction->SetNewPrimariesFlag(true);
    G4cout << "Input directory set to " << newValue << " (PrimariesFiles)" << G4endl;
  }
  
  if( command == neutrinosfileDirectoryCmd )
  {
    myAction->SetNeutrinoFilesDirectory(newValue);
    G4cout << "Input directory set to " << newValue << " (NeutrinoFiles)" << G4endl;
  }
  
  if( command == geniefileDirectoryCmd )
  {
    myAction->SetGenieFilesDirectory(newValue);
    G4cout << "Input directory set to " << newValue << " (GENIE)" << G4endl;
  }
  
  if( command == primariesStartEventCmd )
  {
    myAction->SetPrimariesOffset(primariesStartEventCmd->GetNewIntValue(newValue));
    G4cout << "Primary files will be read starting from entry "<<newValue << G4endl;
  }

  if ( command == positionCmd )
  {
    myAction->SetPosition(positionCmd->ConvertToDimensioned3Vector(newValue));
    G4cout << "Position for antinu event distribution set." << G4endl;
  }

  if ( command == radiusCmd )
  {
    myAction->SetRadius(radiusCmd->ConvertToDimensionedDouble(newValue));
    G4cout << "Radius for antinu event distribution set to "<<newValue << G4endl;
  }

  if ( command == heightCmd )
  {
    myAction->SetHalfZ(heightCmd->ConvertToDimensionedDouble(newValue));
    G4cout << "Height (half Z) for antinu event distribution set to " << newValue << G4endl;
  }

  if ( command == rot1Cmd )
  {
    myAction->SetRot1(rot1Cmd->ConvertTo3Vector(newValue));
    G4cout << "First rotation for antinu event distribution is set." << G4endl;
  }

  if ( command == rot2Cmd )
  {
    myAction->SetRot2(rot2Cmd->ConvertTo3Vector(newValue));
    G4cout << " Second rotation for antinu event distribution is set. " << G4endl;
  }

}

G4String WCSimPrimaryGeneratorMessenger::GetCurrentValue(G4UIcommand* command)
{
  G4String cv;
  
  if( command==genCmd )
  {
    if(myAction->IsUsingMulineEvtGenerator())
      { cv = "muline"; }
    else if(myAction->IsUsingGunEvtGenerator())
      { cv = "gun"; }
    else if(myAction->IsUsingLaserEvtGenerator())
      { cv = "laser"; }   //T. Akiri: Addition of laser
    else if(myAction->IsUsingBeamEvtGenerator())
      { cv = "beam"; }
    else if(myAction->IsUsingGPSEvtGenerator())
      { cv = "gps"; }
    else if(myAction->IsUsingAntiNuEvtGenerator())
      { cv = "antinu"; }
    else if(myAction->IsUsingGenieEvtGenerator())
      { cv = "genie"; }
  }
  
  return cv;
  G4cout<<"generator is currently "<<cv<<G4endl;
}

