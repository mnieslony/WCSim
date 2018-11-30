#include "WCSimPhysicsListFactory.hh"
#include "G4NeutronHPManager.hh"
#include "G4HadronicProcessStore.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

/* This code draws upon examples/extended/fields/field04 for inspiration */



WCSimPhysicsListFactory::WCSimPhysicsListFactory() :  G4VModularPhysicsList()
{
 defaultCutValue = 1.0*mm;
 SetVerboseLevel(0);
 
 PhysicsListName="NULL_LIST"; // default list is set in WCSimPhysicsListFactoryMessenger to FTFP_BERT
 factory = new G4PhysListFactory();
 // TODO create opticalPhyscics object?
 
 std::vector<G4String> ValidListsVector = factory->AvailablePhysLists();
 G4int nlists = ValidListsVector.size();
 G4cout << "There are " << nlists << " available physics lists, and they are: " << G4endl;
 for (G4int i=0; i<nlists; i++){
   G4cout << "  " << ValidListsVector[i] << G4endl;
   ValidListsString += ValidListsVector[i];
   ValidListsString += " ";
 }
 //G4cout << "ValidListsString=" << ValidListsString << G4endl;

 PhysicsMessenger = new WCSimPhysicsListFactoryMessenger(this, ValidListsString);

}

WCSimPhysicsListFactory::~WCSimPhysicsListFactory()
{
  delete PhysicsMessenger;
  PhysicsMessenger = NULL;

}

void WCSimPhysicsListFactory::ConstructParticle()
{
  G4VModularPhysicsList::ConstructParticle();
}

void WCSimPhysicsListFactory::ConstructProcess()
{
  G4VModularPhysicsList::ConstructProcess();
}

void WCSimPhysicsListFactory::SetCuts()
{
  if (verboseLevel >0){
      G4cout << "WCSimPhysicsListFactory::SetCuts:";
      G4cout << "CutLength : " << G4BestUnit(defaultCutValue,"Length") << G4endl;
  }

  // set cut values for gamma at first and for e- second and next for e+,
  // because some processes for e+/e- need cut values for gamma
  //
  SetCutValue(defaultCutValue, "gamma");
  SetCutValue(defaultCutValue, "e-");
  SetCutValue(defaultCutValue, "e+");

  if (verboseLevel>0) DumpCutValuesTable();

}

void WCSimPhysicsListFactory::SetList(G4String newvalue){
  G4cout << "Setting Physics list to " << newvalue << " and delaying initialization" << G4endl;
  PhysicsListName = newvalue;
}

void WCSimPhysicsListFactory::InitializeList(){
  G4cout << "Initializing physics list " << PhysicsListName << G4endl;

  G4VModularPhysicsList* phys = 0;

  if (factory->IsReferencePhysList(PhysicsListName)) {
    G4NeutronHPManager::GetInstance()->SetVerboseLevel(0);
    G4HadronicProcessStore::Instance()->SetVerbose(0);
    phys=factory->GetReferencePhysList(PhysicsListName);
    phys->SetVerboseLevel(0);
    for (G4int i = 0; ; ++i) {
      G4VPhysicsConstructor* elem =
        const_cast<G4VPhysicsConstructor*> (phys->GetPhysics(i));
      if (elem == NULL) break;
      G4cout << "RegisterPhysics: " << elem->GetPhysicsName() << G4endl;
      RegisterPhysics(elem);
    }
    G4cout << "RegisterPhysics: OpticalPhysics" << G4endl; 
    G4OpticalPhysics* opticalPhysics = new G4OpticalPhysics();
    RegisterPhysics(opticalPhysics);
    
    // can optionally turn off cerenkov/scintillation process like this:
    //opticalPhysics->Configure(kCerenkov, true);
    // or in G4.10.3+ can keep the process on to make *number of photons* available,
    // but do not put them on the tracking stack. This is a much faster equivalent to using:
    //   ` if(particleType==G4OpticalPhoton::OpticalPhotonDefinition()) return fKill; `
    // in the UserStackingAction.
    //opticalPhysics->SetCerenkovStackPhotons(true);
    
    // # photons produced is calculated from beta value at start of step, and has nonlinear dependence
    // limit beta change so that beta ~ constant and # photons generated is closer to correct
    opticalPhysics->SetMaxBetaChangePerStep(10.0);
    // similar to above for scintillation??
    opticalPhysics->SetMaxNumPhotonsPerStep(100);
    
    // prevent the stack getting too large by tracking photons first
    opticalPhysics->SetTrackSecondariesFirst(kScintillation,true);
    opticalPhysics->SetTrackSecondariesFirst(kCerenkov,true);
  } else {
    G4cout << "Physics list " << PhysicsListName << " is not understood" << G4endl;
  }
} 

void WCSimPhysicsListFactory::SaveOptionsToOutput(WCSimRootOptions * wcopt)
{
  wcopt->SetPhysicsListName(PhysicsListName);
}
