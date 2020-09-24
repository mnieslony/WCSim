#include "WCSimPrimaryGeneratorMessenger.hh"
#include "WCSimPrimaryGeneratorAction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADouble.hh"
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
  //T. Akiri: Addition of laser
  genCmd->SetGuidance(" Available generators : muline, gun, laser, gps, radioactive, radon, antinu, genie");
  genCmd->SetParameterName("generator",true);
  genCmd->SetDefaultValue("muline");
  //T. Akiri: Addition of laser
  genCmd->SetCandidates("muline gun laser gps radioactive radon antinu genie");

  fileNameCmd = new G4UIcmdWithAString("/mygen/vecfile",this);
  fileNameCmd->SetGuidance("Select the file of vectors.");
  fileNameCmd->SetGuidance(" Enter the file name of the vector file");
  fileNameCmd->SetParameterName("fileName",true);
  fileNameCmd->SetDefaultValue("inputvectorfile");
  
  timeUnitCmd = new G4UIcmdWithAString("/mygen/time_unit",this);
  timeUnitCmd->SetGuidance("Define the units used for tme in the input file.");
  timeUnitCmd->SetGuidance("Can be picosecond, ps, ns, nanosecond, ms, millisecond, s, sec or second");
  timeUnitCmd->SetGuidance("Default if not set is nanosecond");
  timeUnitCmd->SetParameterName("unit",true);
  timeUnitCmd->SetDefaultValue("ns");
  
  radioactive_time_window_Cmd = new G4UIcmdWithADouble("/mygen/radioactive_time_window",this);
  radioactive_time_window_Cmd->SetGuidance("Select time window for radioactivity");
  radioactive_time_window_Cmd->SetParameterName("radioactive_time_window",true);
  radioactive_time_window_Cmd->SetDefaultValue(0.);
  
  isotopeCmd = new G4UIcmdWithAString("/mygen/isotope",this);
  isotopeCmd->SetGuidance("Select properties of radioactive isotope");
  isotopeCmd->SetGuidance("[usage] /mygen/isotope ISOTOPE LOCATION ACTIVITY");
  isotopeCmd->SetGuidance("     ISOTOPE : Tl208, Bi214, K40");
  isotopeCmd->SetGuidance("     LOCATION : water PMT");
  isotopeCmd->SetGuidance("     ACTIVITY : (int) activity of isotope (Bq) ");
  G4UIparameter* param;
  param = new G4UIparameter("ISOTOPE",'s',true);
  param->SetDefaultValue("Tl208");
  isotopeCmd->SetParameter(param);
  param = new G4UIparameter("LOCATION",'s',true);
  param->SetDefaultValue("water");
  isotopeCmd->SetParameter(param);
  param = new G4UIparameter("ACTIVITY",'d',true);
  param->SetDefaultValue("0");
  isotopeCmd->SetParameter(param);
  
  radonScalingCmd = new G4UIcmdWithAString("/mygen/radon_scaling",this);
  radonScalingCmd->SetGuidance("Select scalling scenario");
  radonScalingCmd->SetGuidance("[usage] /mygen/radon SCENARIO ");
  radonScalingCmd->SetGuidance("     SCENARIO : A, B, C");
  radonScalingCmd->SetCandidates("A B C");
  param = new G4UIparameter("SCENARIO",'s',true);
  param->SetDefaultValue("C");
  radonScalingCmd->SetParameter(param);
  
  radonGeoSymCmd = new G4UIcmdWithAnInteger("/mygen/radon_symmetry",this);
  radonGeoSymCmd->SetGuidance("Select scalling scenario");
  radonGeoSymCmd->SetGuidance("[usage] /mygen/radon SCENARIO ");
  radonGeoSymCmd->SetGuidance("     SYMMETRY : 1 ... ");
  param = new G4UIparameter("SYMMETRY",'d',true);
  param->SetDefaultValue("1");
  radonScalingCmd->SetParameter(param);

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

  primariesStartEventCmd = new G4UIcmdWithAnInteger("/mygen/primariesoffset", this);
  primariesStartEventCmd->SetGuidance("The starting entry number for reading primaries");
  primariesStartEventCmd->SetParameterName("primariesoffset",true);
  primariesStartEventCmd->SetDefaultValue(0);
  
  geniefileDirectoryCmd = new G4UIcmdWithAString("/mygen/geniedirectory", this);
  geniefileDirectoryCmd->SetGuidance("Specify the directory containing genie (non-beam) root files");
  geniefileDirectoryCmd->SetParameterName("directoryName",true);
  geniefileDirectoryCmd->SetDefaultValue("");

  talysfileDirectoryCmd = new G4UIcmdWithAString("/mygen/talysdirectory", this);
  talysfileDirectoryCmd->SetGuidance("Specify the directory containing talys root files");
  talysfileDirectoryCmd->SetParameterName("directoryName",true);
  talysfileDirectoryCmd->SetDefaultValue("");
}

WCSimPrimaryGeneratorMessenger::~WCSimPrimaryGeneratorMessenger()
{
  delete genCmd;
  delete mydetDirectory;
  delete talysfileDirectoryCmd;
  delete geniefileDirectoryCmd;
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
      myAction->SetMulineEvtGenerator(true);
      myAction->SetGunEvtGenerator(false);
      myAction->SetLaserEvtGenerator(false);
      myAction->SetGPSEvtGenerator(false);
      myAction->SetRadioactiveEvtGenerator(false);
      myAction->SetRadonEvtGenerator(false);
      myAction->SetAntiNuEvtGenerator(false);
      myAction->SetGenieEvtGenerator(false);
    }
    else if ( newValue == "gun")
    {
      myAction->SetMulineEvtGenerator(false);
      myAction->SetGunEvtGenerator(true);
      myAction->SetLaserEvtGenerator(false);
      myAction->SetGPSEvtGenerator(false);
      myAction->SetRadioactiveEvtGenerator(false);
      myAction->SetRadonEvtGenerator(false);
      myAction->SetAntiNuEvtGenerator(false);
      myAction->SetGenieEvtGenerator(false);
    }
    else if ( newValue == "laser")   //T. Akiri: Addition of laser
    {
      myAction->SetMulineEvtGenerator(false);
      myAction->SetGunEvtGenerator(false);
      myAction->SetLaserEvtGenerator(true);
      myAction->SetGPSEvtGenerator(false);
      myAction->SetRadioactiveEvtGenerator(false);
      myAction->SetRadonEvtGenerator(false);
      myAction->SetAntiNuEvtGenerator(false);
      myAction->SetGenieEvtGenerator(false);
    }
    else if ( newValue == "gps")
    {
      myAction->SetMulineEvtGenerator(false);
      myAction->SetGunEvtGenerator(false);
      myAction->SetLaserEvtGenerator(false);
      myAction->SetGPSEvtGenerator(true);
      myAction->SetRadioactiveEvtGenerator(false);
      myAction->SetRadonEvtGenerator(false);
      myAction->SetAntiNuEvtGenerator(false);
      myAction->SetGenieEvtGenerator(false);
    }
    else if ( newValue == "radioactive") //G. Pronost: Addition of Radioactivity (from F. Nova code)
    {
      myAction->SetMulineEvtGenerator(false);
      myAction->SetGunEvtGenerator(false);
      myAction->SetLaserEvtGenerator(false);
      myAction->SetGPSEvtGenerator(false);
      myAction->SetRadioactiveEvtGenerator(true);
      myAction->SetRadonEvtGenerator(false);
      myAction->SetAntiNuEvtGenerator(false);
      myAction->SetGenieEvtGenerator(false);
    }
    else if ( newValue == "radon" ) //G. Pronost: Addition of Radon generator (based on F. Nova's radioactive generator but dedicated to radioactive events in water)
    {
      myAction->SetMulineEvtGenerator(false);
      myAction->SetGunEvtGenerator(false);
      myAction->SetLaserEvtGenerator(false);
      myAction->SetGPSEvtGenerator(false);
      myAction->SetRadioactiveEvtGenerator(false);
      myAction->SetRadonEvtGenerator(true);
      myAction->SetAntiNuEvtGenerator(false);
      myAction->SetGenieEvtGenerator(false);
    }
    else if ( newValue == "antinu")
    {
      G4cout<<"Setting generator source to antinu"<<G4endl;
      myAction->SetMulineEvtGenerator(false);
      myAction->SetGunEvtGenerator(false);
      myAction->SetLaserEvtGenerator(false);
      myAction->SetGPSEvtGenerator(false);
      myAction->SetRadioactiveEvtGenerator(false);
      myAction->SetRadonEvtGenerator(false);
      myAction->SetAntiNuEvtGenerator(true);
      myAction->SetGenieEvtGenerator(false);
    }
    else if ( newValue == "genie")
    {
      G4cout<<"Setting generator source to genie"<<G4endl;
      myAction->SetMulineEvtGenerator(false);
      myAction->SetGunEvtGenerator(false);
      myAction->SetLaserEvtGenerator(false);
      myAction->SetGPSEvtGenerator(false);
      myAction->SetRadioactiveEvtGenerator(false);
      myAction->SetRadonEvtGenerator(false);
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
  
  if( command==isotopeCmd )
  {
    IsotopeCommand(newValue);
  }

  if( command==radioactive_time_window_Cmd )
  {
    myAction->SetRadioactiveTimeWindow(StoD(newValue));
  }
  
  if ( command==radonScalingCmd ) 
  {
    RadonScalingCommand(newValue);
  }
  
  if ( command==radonGeoSymCmd ) 
  {
    myAction->SetRadonSymmetry(radonGeoSymCmd->GetNewIntValue(newValue));
  }

  if ( command==timeUnitCmd)
  {
    myAction->SetTimeUnit(newValue);
    G4cout << "Time unit set to " << newValue << G4endl;
  }

  if( command == geniefileDirectoryCmd )
  {
    myAction->SetGenieFilesDirectory(newValue);
    G4cout << "Input directory set to " << newValue << " (GENIE)" << G4endl;
  }
  
  if( command == talysfileDirectoryCmd )
  {
    myAction->SetTalysFilesDirectory(newValue);
    G4cout << "Input directory set to " << newValue << " (Talys)" << G4endl;
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
    else if(myAction->IsUsingGPSEvtGenerator())
      { cv = "gps"; }
    else if(myAction->IsUsingRadonEvtGenerator())
      { cv = "radon"; }
    else if(myAction->IsUsingAntiNuEvtGenerator())
      { cv = "antinu"; }
    else if(myAction->IsUsingGenieEvtGenerator())
      { cv = "genie"; }
  }
  
  return cv;
}


void  WCSimPrimaryGeneratorMessenger::IsotopeCommand(G4String newValue)
{
  G4Tokenizer next( newValue );

  G4String isotope = next();
  G4String location = next();
  G4double activity = StoD(next());

  myAction->AddRadioactiveSource(isotope, location, activity);
}

void WCSimPrimaryGeneratorMessenger::RadonScalingCommand(G4String newValue)
{
  G4Tokenizer next( newValue );

  G4String scenario = next();
  G4int iScenario = 0;
   
  if ( scenario == "A" ) iScenario = 1; // Relative scaling with respect to full ID volume (Pessimistic)
  if ( scenario == "B" ) iScenario = 2; // Relative scaling with respect to fiducial volume
  if ( scenario == "C" ) iScenario = 3; // Absolute scaling with respect to ID border (Optimistic)
   
  myAction->SetRadonScenario(iScenario);
}
