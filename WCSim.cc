#include "G4ios.hh"
#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"
#include "WCSimDetectorConstruction.hh"
#include "WCSimPhysicsList.hh"
#include "WCSimPhysicsMessenger.hh"
#include "WCSimPhysicsListFactory.hh"
#include "WCSimPhysicsListFactoryMessenger.hh"
#include "WCSimTuningParameters.hh"
#include "WCSimTuningMessenger.hh"
#include "WCSimPrimaryGeneratorAction.hh"
#include "WCSimEventAction.hh"
#include "WCSimRunAction.hh"
#include "WCSimStackingAction.hh"
#include "WCSimTrackingAction.hh"
#include "WCSimSteppingAction.hh"
#include "WCSimVisManager.hh"
#include "WCSimRandomParameters.hh"

void file_exists(const char * filename) {
  bool exists = access(filename, F_OK) != -1;
  if(!exists) {
    G4cerr << filename << " not found or inaccessible. Exiting" << G4endl;
    exit(-1);
  }
}

int main(int argc,char** argv)
{
  // Construct the default run manager
  G4RunManager* runManager = new G4RunManager;

  // get the pointer to the UI manager
  G4UImanager* UI = G4UImanager::GetUIpointer();

  // Set up the tuning parameters that need to be read before the detector
  //  construction is done
  WCSimTuningParameters* tuningpars = new WCSimTuningParameters();

  // Get the tuning parameters
  file_exists("tuning_parameters.mac");
  UI->ApplyCommand("/control/execute tuning_parameters.mac");

  // define random number generator parameters
  WCSimRandomParameters *randomparameters = new WCSimRandomParameters();

  // UserInitialization classes (mandatory)
  enum DetConfiguration {wfm=1,fwm=2};
  G4int WCSimConfiguration = fwm;
  G4cout<<"Making detector"<<G4endl;
  WCSimDetectorConstruction* WCSimdetector = new 
    WCSimDetectorConstruction(WCSimConfiguration,tuningpars);
  G4cout<<"Setting detector"<<G4endl;
  runManager->SetUserInitialization(WCSimdetector);
  G4cout<<"Making the physics"<<G4endl;
  // Added selectable physics lists 2010-07 by DMW
  // Set up the messenger hooks here, initialize the actual list after loading jobOptions.mac
  WCSimPhysicsListFactory *physFactory = new WCSimPhysicsListFactory();
  G4cout<<"Getting job options"<<G4endl;
  // Currently, default model is set to BINARY
  file_exists("jobOptions.mac");
  UI->ApplyCommand("/control/execute jobOptions.mac");
  G4cout<<"Initialising physics"<<G4endl;
  // Initialize the physics factory to register the selected physics.
  physFactory->InitializeList();
  runManager->SetUserInitialization(physFactory);

  // If the WCSim physics list was chosen in jobOptions.mac,
  // then it's hadronic model needs to be selected in jobOptions2.mac
  //=================================
  // Added by JLR 2005-07-05
  //=================================
  // Choice of hadronic interaction model for 
  // protons & neutrons. This file must be read in
  // by the program BEFORE the runManager is initialized.
  // If file does not exist, default model will be used.
  // Currently, default model is set to BINARY.
  G4cout<<"Getting further job options"<<G4endl;
  file_exists("jobOptions2.mac");
  UI->ApplyCommand("/control/execute jobOptions2.mac");
  G4cout<<"Creating visualisations"<<G4endl;
  // Visualization
  G4VisManager* visManager = new WCSimVisManager;
    G4cout<<"Initialising visualisations"<<G4endl;
  visManager->Initialize();
  G4cout<<"Loading torpedo bays"<<G4endl;
  // Set user action classes
  WCSimPrimaryGeneratorAction* myGeneratorAction = new 
    WCSimPrimaryGeneratorAction(WCSimdetector);
  G4cout<<"Acquiring target lock"<<G4endl;
  runManager->SetUserAction(myGeneratorAction);
  G4cout<<"Devising battle plan"<<G4endl;
  WCSimRunAction* myRunAction = new WCSimRunAction(WCSimdetector);
  G4cout<<"Submitting battle plan"<<G4endl;
  runManager->SetUserAction(myRunAction);
  G4cout<<"Setting event action"<<G4endl;
  runManager->SetUserAction(new WCSimEventAction(myRunAction, WCSimdetector,
						 myGeneratorAction));
  G4cout<<"Setting tracking action"<<G4endl;
  runManager->SetUserAction(new WCSimTrackingAction);
  G4cout<<"Setting stacking action"<<G4endl;
  runManager->SetUserAction(new WCSimStackingAction(WCSimdetector));
  G4cout<<"Setting stepping action"<<G4endl;
  runManager->SetUserAction(new WCSimSteppingAction);

  G4cout<<"Initialise"<<G4endl;
  // Initialize G4 kernel
  runManager->Initialize();

  if (argc==1)   // Define UI terminal for interactive mode  
  { 

    // Start UI Session
    G4UIsession* session =  new G4UIterminal(new G4UItcsh);

    G4cout<<"On screen!"<<G4endl;
    // Visualization Macro
    UI->ApplyCommand("/control/execute WCSim.mac");

    G4cout<<"Begin!"<<G4endl;
    // Start Interactive Mode
    session->SessionStart();

    delete session;
  }
  else           // Batch mode
  { 
    G4String command = "/control/execute ";
    G4String fileName = argv[1];

    UI->ApplyCommand(command+fileName);
  }

  delete visManager;

  delete runManager;
  return 0;
}


