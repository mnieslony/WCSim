#include "WCSimTuningMessenger.hh"
#include "WCSimTuningParameters.hh"

#include "G4UIdirectory.hh"
#include "G4UIcommand.hh"
#include "G4UIparameter.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithABool.hh" //jl145


WCSimTuningMessenger::WCSimTuningMessenger(WCSimTuningParameters* WCTuningPars):WCSimTuningParams(WCTuningPars) { 
  WCSimDir = new G4UIdirectory("/WCSim/tuning/");
  WCSimDir->SetGuidance("Commands to change tuning parameters");

  Rayff = new G4UIcmdWithADouble("/WCSim/tuning/rayff",this);
  Rayff->SetGuidance("Set the Rayleigh scattering parameter");
  Rayff->SetParameterName("Rayff",true);
  Rayff->SetDefaultValue(0.75);

  Bsrff = new G4UIcmdWithADouble("/WCSim/tuning/bsrff",this);
  Bsrff->SetGuidance("Set the Blacksheet reflection parameter");
  Bsrff->SetParameterName("Bsrff",true);
  Bsrff->SetDefaultValue(2.50);

  Abwff = new G4UIcmdWithADouble("/WCSim/tuning/abwff",this);
  Abwff->SetGuidance("Set the water attenuation parameter");
  Abwff->SetParameterName("Abwff",true);
  Abwff->SetDefaultValue(1.30);

  Rgcff = new G4UIcmdWithADouble("/WCSim/tuning/rgcff",this);
  Rgcff->SetGuidance("Set the cathode reflectivity parameter");
  Rgcff->SetParameterName("Rgcff",true);
  Rgcff->SetDefaultValue(0.32);

  Mieff = new G4UIcmdWithADouble("/WCSim/tuning/mieff",this);
  Mieff->SetGuidance("Set the Mie scattering parameter");
  Mieff->SetParameterName("Mieff",true);
  Mieff->SetDefaultValue(0.0);

  Teflonrff = new G4UIcmdWithADouble("/WCSim/tuning/teflonrff",this);
  Teflonrff->SetGuidance("Set teflon wrapping reflectivity parameter");
  Teflonrff->SetParameterName("Teflonrff",true);
  Teflonrff->SetDefaultValue(1.00);

  Holderrff = new G4UIcmdWithADouble("/WCSim/tuning/holderrff",this);
  Holderrff->SetGuidance("Set ANNIE holder reflectivity parameter");
  Holderrff->SetParameterName("Holderrff",true);
  Holderrff->SetDefaultValue(1.00);

  HolderrffLUX = new G4UIcmdWithADouble("/WCSim/tuning/holderrfflux",this);
  HolderrffLUX->SetGuidance("Set LUX/ETEL holder reflectivity parameter");
  HolderrffLUX->SetParameterName("HolderrffLUX",true);
  HolderrffLUX->SetDefaultValue(1.00);

  Linerrff = new G4UIcmdWithADouble("/WCSim/tuning/linerrff",this);
  Linerrff->SetGuidance("Set Liner reflectivity parameter");
  Linerrff->SetParameterName("Linerrff",true);
  Linerrff->SetDefaultValue(1.00);

  Holder = new G4UIcmdWithABool("/WCSim/tuning/holder",this);
  Holder->SetGuidance("Turn ANNIE holders on/off");
  Holder->SetParameterName("Holder",true);
  Holder->SetDefaultValue(0);

  //jl145 - for Top Veto
  TVSpacing = new G4UIcmdWithADouble("/WCSim/tuning/tvspacing",this);
  TVSpacing->SetGuidance("Set the Top Veto PMT Spacing, in cm.");
  TVSpacing->SetParameterName("TVSpacing",true);
  TVSpacing->SetDefaultValue(100.0);

  TopVeto = new G4UIcmdWithABool("/WCSim/tuning/topveto",this);
  TopVeto->SetGuidance("Turn Top Veto simulation on/off");
  TopVeto->SetParameterName("TopVeto",true);
  TopVeto->SetDefaultValue(0);

}

WCSimTuningMessenger::~WCSimTuningMessenger()
{
  delete Rayff;
  delete Bsrff;
  delete Abwff;
  delete Rgcff;
  delete Mieff;

  //ANNIE-specific variables
  delete Teflonrff;
  delete Holderrff;
  delete HolderrffLUX;
  delete Linerrff;
  delete Holder;

  //jl145 - for Top Veto
  delete TVSpacing;
  delete TopVeto;

  delete WCSimDir;
}

void WCSimTuningMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{    

  if(command == Rayff) {
	  // Set the Rayleigh scattering parameter
	  //	  printf("Input parameter %f\n",Rayff->GetNewDoubleValue(newValue));

  WCSimTuningParams->SetRayff(Rayff->GetNewDoubleValue(newValue));

  printf("Setting Rayleigh scattering parameter %f\n",Rayff->GetNewDoubleValue(newValue));

  }

 if(command == Bsrff) {
	  // Set the blacksheet reflection parameter
	  //	  printf("Input parameter %f\n",Bsrff->GetNewDoubleValue(newValue));

  WCSimTuningParams->SetBsrff(Bsrff->GetNewDoubleValue(newValue));

  printf("Setting blacksheet reflection parameter %f\n",Bsrff->GetNewDoubleValue(newValue));

  }

 if(command == Abwff) {
	  // Set the water attenuation  parameter
	  //	  printf("Input parameter %f\n",Abwff->GetNewDoubleValue(newValue));

  WCSimTuningParams->SetAbwff(Abwff->GetNewDoubleValue(newValue));

  printf("Setting water attenuation parameter %f\n",Abwff->GetNewDoubleValue(newValue));

  }

 if(command == Rgcff) {
	  // Set the cathode reflectivity parameter
	  //	  printf("Input parameter %f\n",Rgcff->GetNewDoubleValue(newValue));

  WCSimTuningParams->SetRgcff(Rgcff->GetNewDoubleValue(newValue));

  printf("Setting cathode reflectivity parameter %f\n",Rgcff->GetNewDoubleValue(newValue));

  }

 if(command == Mieff) {
	  // Set the Mie scattering parameter
	  //	  printf("Input parameter %f\n",Mieff->GetNewDoubleValue(newValue));

  WCSimTuningParams->SetMieff(Mieff->GetNewDoubleValue(newValue));

  printf("Setting Mie scattering parameter %f\n",Mieff->GetNewDoubleValue(newValue));

  }

  // ANNIE - Teflon reflectivity
  if (command == Teflonrff){
    // Set the Teflon reflectivity parameter
    //    printf("Input parameter %f\n",Teflonrff->GetNewDoubleValue(newValue));

  WCSimTuningParams->SetTeflonrff(Teflonrff->GetNewDoubleValue(newValue));

  printf("Setting ANNIE Teflon reflectivity parameter %f\n",Teflonrff->GetNewDoubleValue(newValue));

  }

  // ANNIE - Holder reflectivity
  if (command == Holderrff){
    // Set the ANNIE holder reflectivity parameter
    //    printf("Input parameter %f\n",Holderrff->GetNewDoubleValue(newValue));

  WCSimTuningParams->SetHolderrff(Holderrff->GetNewDoubleValue(newValue));

  printf("Setting ANNIE holder reflectivity parameter %f\n",Holderrff->GetNewDoubleValue(newValue));
  
  }

  // ANNIE - ETEL/LUX Holder reflectivity
  if (command == HolderrffLUX){
    // Set the LUX/ETEL holder reflectivity parameter
    //    printf("Input parameter %f\n",Holderrff->GetNewDoubleValue(newValue));

  WCSimTuningParams->SetHolderrffLUX(HolderrffLUX->GetNewDoubleValue(newValue));

  printf("Setting LUX/ETEL holder reflectivity parameter %f\n",HolderrffLUX->GetNewDoubleValue(newValue));
  
  }

  // ANNIE - Liner reflectivity
  if (command == Linerrff){
    // Set the ANNIE liner reflectivity parameter
    //    printf("Input parameter %f\n",Linerrff->GetNewDoubleValue(newValue));

  WCSimTuningParams->SetHolderrff(Linerrff->GetNewDoubleValue(newValue));

  printf("Setting ANNIE liner reflectivity parameter %f\n",Linerrff->GetNewDoubleValue(newValue));
  
  }

  else if(command == Holder) {
    // Turn the ANNIE PMT holders on/off in the simulation
    WCSimTuningParams->SetHolder(Holder->GetNewBoolValue(newValue));
    if(Holder->GetNewBoolValue(newValue))
      printf("Setting ANNIE PMT holders On\n");
    else
      printf("Setting ANNIE PMT holders Off\n");
  }

  //jl145 - For Top Veto

  else if(command == TVSpacing) {
    // Set the Top Veto PMT Spacing
    WCSimTuningParams->SetTVSpacing(TVSpacing->GetNewDoubleValue(newValue));
    printf("Setting Top Veto PMT Spacing %f\n",TVSpacing->GetNewDoubleValue(newValue));
  }

  else if(command == TopVeto) {
    // Set the Top Veto on/off
    WCSimTuningParams->SetTopVeto(TopVeto->GetNewBoolValue(newValue));
    if(TopVeto->GetNewBoolValue(newValue))
      printf("Setting Top Veto On\n");
    else
      printf("Setting Top Veto Off\n");
  }


}
