/* vim:set noexpandtab tabstop=4 wrap */
#include "WCSimPrimaryGeneratorAction.hh"
#include "WCSimDetectorConstruction.hh"
#include "WCSimPrimaryGeneratorMessenger.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4GeneralParticleSource.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleTypes.hh"
#include "G4IonTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4ThreeVector.hh"
#include "G4SPSAngDistribution.hh"
#include "G4SPSPosDistribution.hh"
#include "TFile.h"
#include "TLorentzVector.h"
#include "globals.hh"
#include "Randomize.hh"
#include <fstream>
#include <vector>
#include <string>
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include <math.h> 
#include <libgen.h>

#include "G4Navigator.hh"
#include "G4TransportationManager.hh"
#include "G4UImanager.hh"

// GENIE headers
#ifndef NO_GENIE
#include "GHEP/GHepParticle.h"
#include "GHEP/GHepUtils.h"
#include "Ntuple/NtpMCTreeHeader.h"
#include "Interaction/Interaction.h"
#endif

// when loading dirt primaries, skip entries that are from upstream rock interactions. 
#ifndef ONLY_TANK_EVENTS
//#define ONLY_TANK_EVENTS
#endif

// as help for reconstruction, a sample where light is only from muons, no other primary particles from the event.
#ifndef ONLY_MUONS
//#define ONLY_MUONS
#endif

using std::vector;
using std::string;
using std::fstream;

vector<string> tokenize( string separators, string input );

inline vector<string> readInLine(fstream& inFile, int lineSize, char* inBuf)
{
  // Read in line break it up into tokens
  inFile.getline(inBuf,lineSize);
  return tokenize(" $", inBuf);
}

inline void ReadInEnergySpectrum(fstream& EFile,std::vector<G4double>& vecE,std::vector<G4double>& vecProb)
{
	double temp_E, temp_Prob;
	while (EFile >> temp_E >> temp_Prob){
		if (EFile.eof()) break;
		std::cout <<"Reading in antineutrino energy "<<temp_E<<", and probability "<<temp_Prob<<std::endl;
		vecE.push_back(temp_E);
		vecProb.push_back(temp_Prob);
	}
}

inline float atof( const string& s ) {return std::atof( s.c_str() );}
inline int   atoi( const string& s ) {return std::atoi( s.c_str() );}

WCSimPrimaryGeneratorAction::WCSimPrimaryGeneratorAction(
					  WCSimDetectorConstruction* myDC)
  :myDetector(myDC), loadNewPrimaries(true), loadNewGenie(true), inputdata(0), geniedata(0), primariesDirectory(""), neutrinosDirectory(""), genieDirectory(""), vectorFileName("")
{
  //T. Akiri: Initialize GPS to allow for the laser use 
  MyGPS = new G4GeneralParticleSource();

  // Initialize to zero
  mode = 0;
  nvtxs = 0;
  for( Int_t u=0; u<MAX_N_PRIMARIES; u++){
    vtxsvol[u] = 0;
    vtxs[u] = G4ThreeVector(0.,0.,0.);
  }
  nuEnergy = 0.;
  _counterRock=0; // counter for generated in Rock
  _counterCublic=0; // counter generated
  
  //---Set defaults. Do once at beginning of session.
  
  G4int n_particle = 1;
  particleGun = new G4ParticleGun(n_particle);
  particleGun->SetParticleEnergy(1.0*CLHEP::GeV);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.0));
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4IonTable* ionTable = G4IonTable::GetIonTable();
  G4String particleName;
  particleGun->
    SetParticleDefinition(particleTable->FindParticle(particleName="mu+"));

  particleGun->
    SetParticlePosition(G4ThreeVector(0.*CLHEP::m,0.*CLHEP::m,0.*CLHEP::m));
    
  theSPSAng = new G4SPSAngDistribution;
  theSPSPos = new G4SPSPosDistribution;
  G4SPSRandomGenerator *RndGen = new G4SPSRandomGenerator;
  theSPSPos->SetBiasRndm(RndGen);
  theSPSAng->SetBiasRndm(RndGen);
  theSPSPos->SetPosDisType("Volume");
  theSPSPos->SetPosDisShape("Cylinder");
  thePosition = G4ThreeVector(0., 0., 0.);
      
  Espectrum.clear();
  ProbabilitySpec.clear();
  isFirstEvent = true;

  messenger = new WCSimPrimaryGeneratorMessenger(this);
  useMulineEvt = false;
  useGunEvt = false;
  useLaserEvt = false;
  useBeamEvt = true;
  useGPSEvt = false;
  useAntiNuEvt = false; 
  useGenieEvt = false; 
    
#ifndef NO_GENIE
  genierecordval = new genie::NtpMCEventRecord;
#endif
}

WCSimPrimaryGeneratorAction::~WCSimPrimaryGeneratorAction()
{
  if (IsGeneratingVertexInRock()){
    G4cout << "Fraction of Rock volume is : " << G4endl;
      G4cout << " Random number generated in Rock / in Cublic = " 
             << _counterRock << "/" << _counterCublic 
             << " = " << _counterRock/(G4double)_counterCublic << G4endl;
  }
  inputFile.close();
  inputSpecFile.close();
  delete particleGun;
  delete MyGPS;   //T. Akiri: Delete the GPS variable
  delete messenger;

  delete theSPSPos;
  delete theSPSAng;
  
  if(useBeamEvt){
    if(inputdata){
      inputdata->ResetBranchAddresses();
      metadata->ResetBranchAddresses();
      delete inputdata;
      delete metadata;
#ifndef NO_GENIE
      if(geniedata) geniedata->ResetBranchAddresses();
      if(geniedata) delete geniedata;
      if(genierecordval) delete genierecordval;
#endif
    }
  }

  if(useGenieEvt){
    if(geniedata){
      geniedata->ResetBranchAddresses();
      delete geniedata;
    }
   }

}

void WCSimPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{

  // We will need a particle table
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();

  // Temporary kludge to turn on/off vector text format 

  G4bool useNuanceTextFormat = true;


  // Do for every event

  if (useMulineEvt)
  {

    if ( !inputFile.is_open() )
    {
      G4cout << "Set a vector file using the command /mygen/vecfile name"
	     << G4endl;
      exit(-1);
    }

    //
    // Documentation describing the nuance text format can be found here: 
    // http://neutrino.phy.duke.edu/nuance-format/
    //
    // The format must be strictly adhered to for it to be processed correctly.
    // The lines and their meanings from begin through info are fixed, and then
    // a variable number of tracks may follow.
    //
    if (useNuanceTextFormat)
      {
	const int lineSize=100;
	char      inBuf[lineSize];
	vector<string> token(1);
	
	token = readInLine(inputFile, lineSize, inBuf);
	  
        if (token.size() == 0) 
	  {
	    G4cout << "end of nuance vector file!" << G4endl;
	  }
	else if (token[0] != "begin")
	  {
	    G4cout << "unexpected line begins with " << token[0] << G4endl;
	  }
	else   // normal parsing begins here
	  {
	    // Read the nuance line (ignore value now)

	    token = readInLine(inputFile, lineSize, inBuf);
	    mode = atoi(token[1]);

	    // Read the Vertex line
	    token = readInLine(inputFile, lineSize, inBuf);
	    vtxs[0] = G4ThreeVector(atof(token[1])*cm,
				    atof(token[2])*cm,
				    atof(token[3])*cm);
	    
            // true : Generate vertex in Rock , false : Generate vertex in WC tank
            SetGenerateVertexInRock(false);

	    // Next we read the incoming neutrino and target
	    
	    // First, the neutrino line

	    token=readInLine(inputFile, lineSize, inBuf);
	    beampdgs[0] = atoi(token[1]);
	    beamenergies[0] = atof(token[2])*MeV;
	    beamdirs[0] = G4ThreeVector(atof(token[3]),
					atof(token[4]),
					atof(token[5]));

	    // Now read the target line

	    token=readInLine(inputFile, lineSize, inBuf);
	    targetpdgs[0] = atoi(token[1]);
	    targetenergies[0] = atof(token[2])*MeV;
	    targetdirs[0] = G4ThreeVector(atof(token[3]),
					  atof(token[4]),
					  atof(token[5]));

	    // Read the info line, basically a dummy
	    token=readInLine(inputFile, lineSize, inBuf);
	    G4cout << "Vector File Record Number " << token[2] << G4endl;
            vecRecNumber = atoi(token[2]);
	    
	    // Now read the outgoing particles
	    // These we will simulate.


	    while ( token=readInLine(inputFile, lineSize, inBuf),
		    token[0] == "track" )
	      {
		// We are only interested in the particles
		// that leave the nucleus, tagged by "0"


		if ( token[6] == "0")
		  {
		    G4int pdgid = atoi(token[1]);
		    G4double energy = atof(token[2])*CLHEP::MeV;
		    G4ThreeVector dir = G4ThreeVector(atof(token[3]),
						      atof(token[4]),
						      atof(token[5]));
		    std::cout<<"PDGcode "<<pdgid<<"\n";
		    //must handle the case of an ion speratly from other particles
		    //check PDG code if we have an ion.
		    //PDG code format for ions ±10LZZZAAAI
		    char strPDG[11];
		    char strA[10]={0};
		    char strZ[10]={0};
		    

		    long int A=0,Z=0;
		    //		    A=strotl(strPDG,&str);
		    if(abs(pdgid) >= 1000000000)
		      {
			//ion
			sprintf(strPDG,"%i",abs(pdgid));
			strncpy(strZ, &strPDG[3], 3);
			strncpy(strA, &strPDG[6], 3);
			strA[3]='\0';
			strZ[3]='\0';
			A=atoi(strA);
			Z=atoi(strZ);
			G4ParticleDefinition* ion;
			ion =  G4IonTable::GetIonTable()->GetIon(Z, A, 0.);
			particleGun->SetParticleDefinition(ion);
			particleGun->SetParticleCharge(0);
		      }
		    else {
		      //not ion
		      particleGun->
			SetParticleDefinition(particleTable->
		      FindParticle(pdgid));
		    }
		    G4double mass = 
		      particleGun->GetParticleDefinition()->GetPDGMass();

		    G4double ekin = energy - mass;

		    particleGun->SetParticleEnergy(ekin);
		    //G4cout << "Particle: " << pdgid << " KE: " << ekin << G4endl;
		    particleGun->SetParticlePosition(vtxs[0]);
		    particleGun->SetParticleMomentumDirection(dir);
		    particleGun->GeneratePrimaryVertex(anEvent);
		  }
	      }
	  }
      }
    else 
      {    // old muline format  
	inputFile >> nuEnergy >> energy >> xPos >> yPos >> zPos 
		  >> xDir >> yDir >> zDir;
	
	G4double random_z = ((myDetector->GetWaterTubePosition())
			     - .5*(myDetector->GetWaterTubeLength()) 
			     + 1.*CLHEP::m + 15.0*CLHEP::m*G4UniformRand())/CLHEP::m;
	zPos = random_z;
	G4ThreeVector vtx = G4ThreeVector(xPos, yPos, random_z);
	G4ThreeVector dir = G4ThreeVector(xDir,yDir,zDir);

	particleGun->SetParticleEnergy(energy*CLHEP::MeV);
	particleGun->SetParticlePosition(vtx);
	particleGun->SetParticleMomentumDirection(dir);
	particleGun->GeneratePrimaryVertex(anEvent);
      }
  }

  else if (useGunEvt)
  {      // manual gun operation
    particleGun->GeneratePrimaryVertex(anEvent);

    //To prevent occasional seg fault from an un assigned targetpdg 
    targetpdgs[0] = 2212; //ie. proton

    G4ThreeVector P  =anEvent->GetPrimaryVertex()->GetPrimary()->GetMomentum();
    G4ThreeVector vtx=anEvent->GetPrimaryVertex()->GetPosition();
    G4double m       =anEvent->GetPrimaryVertex()->GetPrimary()->GetMass();
    G4int pdg        =anEvent->GetPrimaryVertex()->GetPrimary()->GetPDGcode();
   

    char strPDG[11];
    char strA[10]={0};
    char strZ[10]={0};
    
    
    long int A=0,Z=0;
    //		    A=strotl(strPDG,&str);
    if(abs(pdg) >= 1000000000)
      {
	//ion
	sprintf(strPDG,"%i",abs(pdg));
	strncpy(strZ, &strPDG[3], 3);
	strncpy(strA, &strPDG[6], 3);
	strA[3]='\0';
	strZ[3]='\0';
	A=atoi(strA);
	Z=atoi(strZ);

	G4ParticleDefinition* ion   = G4IonTable::GetIonTable()->GetIon(Z, A, 0);
	ion->SetPDGStable(false);
	ion->SetPDGLifeTime(0.);
	
	G4ParticleDefinition* ion2   = G4IonTable::GetIonTable()->GetIon(Z, A, 0);
	std::cout<<"ion2 "<<ion2->GetPDGLifeTime()<<"\n";
      }
    
    
    G4ThreeVector dir  = P.unit();
    G4double E         = std::sqrt((P.dot(P))+(m*m));

//     particleGun->SetParticleEnergy(E);
//     particleGun->SetParticlePosition(vtx);
//     particleGun->SetParticleMomentumDirection(dir);

    SetVtx(vtx);
    SetBeamEnergy(E);
    SetBeamDir(dir);
    SetBeamPDG(pdg);
  }
  else if (useLaserEvt)
    {
      targetpdgs[0] = 2212; //ie. proton 
      //T. Akiri: Create the GPS LASER event
      MyGPS->GeneratePrimaryVertex(anEvent);
      
      G4ThreeVector P   =anEvent->GetPrimaryVertex()->GetPrimary()->GetMomentum();
      G4double m        =anEvent->GetPrimaryVertex()->GetPrimary()->GetMass();      // this is rest mass
      G4ThreeVector vtx =anEvent->GetPrimaryVertex()->GetPosition();
      G4int pdg         =anEvent->GetPrimaryVertex()->GetPrimary()->GetPDGcode();
      
      G4ThreeVector dir  = P.unit();
      G4double E         = std::sqrt((P.dot(P))+(m*m));
      
      SetVtx(vtx);   // required to store the true vertex for Bonsai!
      SetBeamEnergy(E);
      SetBeamDir(dir);
      SetBeamPDG(pdg);
      const G4ParticleDefinition* parttype = anEvent->GetPrimaryVertex()->GetPrimary()->GetParticleDefinition();
      G4String particlename;
      if(parttype == G4OpticalPhoton::OpticalPhotonDefinition() ) {particlename="opticalphoton";}
      else { (parttype) ? particlename=parttype->GetParticleName() : particlename=std::to_string(pdg); }
      
      double tote = anEvent->GetPrimaryVertex()->GetPrimary()->GetTotalEnergy();
      double ke = anEvent->GetPrimaryVertex()->GetPrimary()->GetKineticEnergy();
      
      G4cout<<"Generating laser primary "<<particlename<<" with total energy "
            <<tote/MeV<<"MeV and kinetic energy "<<ke/MeV
            <<"MeV at ("<<vtx.x()/cm<<","<<vtx.y()/cm<<","<<vtx.z()/cm<<") in direction ("
            <<dir.x()<<", "<<dir.y()<<", "<<dir.z()<<") ";
      G4Navigator* theNavigator = G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking();
      G4VPhysicalVolume* primaryPV = theNavigator->LocateGlobalPointAndSetup(vtx);
      std::string vtxvol = ( (primaryPV) ? primaryPV->GetName() : "<<no-primaryPV>>" );
      G4cout<<" in "<<vtxvol<<G4endl;
    }
  
  else if (useBeamEvt)
    {
		// Load the next entry, with all required trees and files
		// ------------------------------------------------------
		loadbeamentry:
		if(loadNewPrimaries){ LoadNewPrimaries(); } // update TChain if a new file is loaded by messenger
		//inputdata has already had tree loaded at the end of last event's GeneratePrimaries call
		//localEntry will already be the value of the NEXT entry
		metadata->LoadTree(inputEntry);
		
		Int_t nextTreeNumber = inputdata->GetTreeNumber();
		if(treeNumber!=nextTreeNumber){
			G4cout<< "Reached end of Tree. Last entries' tree number was "
						<< treeNumber <<", this entries' tree number is "<< nextTreeNumber <<G4endl;
			dirtFileName = inputdata->GetCurrentFile()->GetName(); // new tree, new file
			char* dirtFileNameAsChar = strdup(dirtFileName.c_str());
			dirtFileName = basename(dirtFileNameAsChar);
			inputdata->SetBranchAddress("run",&runbranchval,&runBranch);
			inputdata->SetBranchAddress("ntank",&ntankbranchval,&nTankBranch);
			inputdata->SetBranchAddress("nupdg",&nupdgval,&nupdgBranch);
			inputdata->SetBranchAddress("nuvtxx",&nuvtxxval,&nuvtxxBranch);
			inputdata->SetBranchAddress("nuvtxy",&nuvtxyval,&nuvtxyBranch);
			inputdata->SetBranchAddress("nuvtxz",&nuvtxzval,&nuvtxzBranch);
			inputdata->SetBranchAddress("nuvtxt",&nuvtxtval,&nuvtxtBranch);
			inputdata->SetBranchAddress("vtxvol",&nupvval,&nuPVBranch);
			inputdata->SetBranchAddress("vtxmat",&numatval,&nuvtxmatBranch);
			inputdata->SetBranchAddress("entry",&genieentrybranchval,&genieentryBranch);
			metadata->SetBranchAddress("inputFluxName",&nufluxfilenameval,&nufluxfilenameBranch);
#ifndef NO_GENIE
			geniedata->SetBranchAddress("gmcrec",&genierecordval,&genierecordBranch);
			genierecordBranch->SetAutoDelete(kTRUE);
#else 
			genierecordBranch=(TBranch*)1;
#endif
			vtxxBranch=inputdata->GetBranch("vx");
			vtxyBranch=inputdata->GetBranch("vy");
			vtxzBranch=inputdata->GetBranch("vz");
			vtxtBranch=inputdata->GetBranch("vt");
			pxBranch=inputdata->GetBranch("px");
			pyBranch=inputdata->GetBranch("py");
			pzBranch=inputdata->GetBranch("pz");
			EBranch=inputdata->GetBranch("E");
			KEBranch=inputdata->GetBranch("kE");
			pdgBranch=inputdata->GetBranch("pdgtank");
			nuprimaryBranch=inputdata->GetBranch("primary");
		
			if(runBranch==0||nTankBranch==0||vtxxBranch==0||vtxyBranch==0||vtxzBranch==0||vtxtBranch==0||pxBranch==0||pyBranch==0||pzBranch==0||EBranch==0||KEBranch==0||pdgBranch==0||nupdgBranch==0||nuvtxxBranch==0||nuvtxyBranch==0||nuvtxzBranch==0||nuvtxtBranch==0||nuPVBranch==0||nuvtxmatBranch==0||nuprimaryBranch==0||nufluxfilenameBranch==0||genierecordBranch==0){
				G4cout<<"BRANCHES ARE ZOMBIES ARGH!"<<G4endl;
			} else { G4cout<<"entries in this tree: "<<vtxxBranch->GetEntries()<<G4endl; }
			
			entriesInThisTree = runBranch->GetEntries();
			treeNumber=nextTreeNumber;
		}
		
		G4cout<<"Loading primaries from entry "<<inputEntry<<", localentry "<<localEntry<<"/"<<entriesInThisTree<<G4endl;
		//runBranch->GetEntry(localEntry);
		//G4cout<<"Run is "<<runbranchval<<G4endl;
		
		nTankBranch->GetEntry(localEntry);
		nupdgBranch->GetEntry(localEntry);
		nuvtxxBranch->GetEntry(localEntry);
		nuvtxyBranch->GetEntry(localEntry);
		nuvtxzBranch->GetEntry(localEntry);
		nuvtxtBranch->GetEntry(localEntry);
		nuPVBranch->GetEntry(localEntry);
		nuvtxmatBranch->GetEntry(localEntry);
		genieentryBranch->GetEntry(localEntry);
		nufluxfilenameBranch->GetEntry(localEntry);
		
		// note info about this input event for recording into output file
		dirtEntryNum = localEntry;
		genieEntryNum = genieentrybranchval;
		genieFileName = nufluxfilenameval;
		
#ifdef ONLY_TANK_EVENTS
		if(strcmp(numatval,"TankWater")!=0){ // nu intx not in tank
		G4cout<<"---------------SKIPPING NON-TANK ENTRY----------------"<<G4endl;
		inputEntry++; 
		localEntry = inputdata->LoadTree(inputEntry);
		if(localEntry<0){
			// get the pointer to the UI manager
			G4UImanager* UI = G4UImanager::GetUIpointer();
			UI->ApplyCommand("/run/abort 0");	// abort without processing current event
			G4cout<<"@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#"<<G4endl;
			G4cout<<"@#@#@#@#@#@#@#@#@#@ REACHED END OF INPUT FILE! #@#@#@#@#@#@#@#@#@#@#@#"<<G4endl;
			G4cout<<"@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#"<<G4endl;
		} else { goto loadbeamentry; } // load the next entry
		}
#endif
		
#ifndef NO_GENIE
		Long64_t genielocalEntry = geniedata->LoadTree(genieentrybranchval);
		// load the appropriate genie entry. we assume 1:1 correspondance of genie:g4dirt files.
		// So the following should not be necessary, as the files should be loaded synchronously
		if(genielocalEntry<0){
			// get the pointer to the UI manager
			G4UImanager* UI = G4UImanager::GetUIpointer();
			UI->ApplyCommand("/run/abort 1");	// abort after processing current event
			G4cout<<"@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#"<<G4endl;
			G4cout<<"@#@#@#@#@#@#@#@#@#@ REACHED END OF GENIE FILE! #@#@#@#@#@#@#@#@#@#@#@#"<<G4endl;
			G4cout<<"@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#"<<G4endl;
		}
		genierecordBranch->GetEntry(genieentrybranchval);
#endif
		
		G4ParticleDefinition* parttype = particleTable->FindParticle(nupdgval);
		TLorentzVector neutrinovertex(nuvtxtval*CLHEP::s, nuvtxxval*CLHEP::m, nuvtxyval*CLHEP::m, nuvtxzval*CLHEP::m);	// position in m, times in s, convert to cm and ns
		G4cout<<"The origin interaction was a "<<(parttype->GetParticleName())<<" at ("<<nuvtxtval*1000000000.<<","<<nuvtxxval*100.<<","<<nuvtxyval*100.<<","<<nuvtxzval*100.<<")[ns, cm] in "<<nupvval<<" "<<numatval<<G4endl;
		G4cout<<"This entry has "<<ntankbranchval<<" primaries"<<G4endl;
		
		if(vtxxbranchval){delete[] vtxxbranchval;}
		if(vtxybranchval){delete[] vtxybranchval;}
		if(vtxzbranchval){delete[] vtxzbranchval;}
		if(vtxtbranchval){delete[] vtxtbranchval;}
		if(pxbranchval){delete[] pxbranchval;}
		if(pybranchval){delete[] pybranchval;}
		if(pzbranchval){delete[] pzbranchval;}
		if(ebranchval){delete[] ebranchval;}
		if(kebranchval){delete[] kebranchval;}
		if(pdgbranchval){delete[] pdgbranchval;}
		if(nuprimarybranchval){delete[] nuprimarybranchval;}
		
		vtxxbranchval = new Double_t[ntankbranchval];
		vtxybranchval = new Double_t[ntankbranchval];
		vtxzbranchval = new Double_t[ntankbranchval];
		vtxtbranchval = new Double_t[ntankbranchval];
		pxbranchval = new Double_t[ntankbranchval];
		pybranchval = new Double_t[ntankbranchval];
		pzbranchval = new Double_t[ntankbranchval];
		ebranchval = new Double_t[ntankbranchval];
		kebranchval = new Double_t[ntankbranchval];
		pdgbranchval = new Int_t[ntankbranchval];
		nuprimarybranchval = new Int_t[ntankbranchval];
		
		if(vtxxbranchval==0||vtxybranchval==0||vtxzbranchval==0||vtxtbranchval==0||pxbranchval==0||pybranchval==0||pzbranchval==0||ebranchval==0||kebranchval==0||pdgbranchval==0||nuprimarybranchval==0){
			G4cout<<"Arrays are zombies!"<<G4endl;
		}
		
		//G4cout<<"Setting branch addresses"<<G4endl;
		vtxxBranch->SetAddress(vtxxbranchval);
		vtxyBranch->SetAddress(vtxybranchval);
		vtxzBranch->SetAddress(vtxzbranchval);
		vtxtBranch->SetAddress(vtxtbranchval);
		pxBranch->SetAddress(pxbranchval);
		pyBranch->SetAddress(pybranchval);
		pzBranch->SetAddress(pzbranchval);
		EBranch->SetAddress(ebranchval);
		KEBranch->SetAddress(kebranchval);
		pdgBranch->SetAddress(pdgbranchval);
		nuprimaryBranch->SetAddress(nuprimarybranchval);
		
		//G4cout<<"Getting primary arrays"<<G4endl;
		vtxxBranch->GetEntry(localEntry);
		vtxyBranch->GetEntry(localEntry);
		vtxzBranch->GetEntry(localEntry);
		vtxtBranch->GetEntry(localEntry);
		pxBranch->GetEntry(localEntry);
		pyBranch->GetEntry(localEntry);
		pzBranch->GetEntry(localEntry);
		EBranch->GetEntry(localEntry);
		KEBranch->GetEntry(localEntry);
		pdgBranch->GetEntry(localEntry);
		nuprimaryBranch->GetEntry(localEntry);
		
		// OK, actually use the loaded information
		// ----------------------------------------
#ifdef ONLY_TANK_EVENTS
		// this *should* always be true, since tank nu events will generate primaries...
		Bool_t primariesinthisentry=false;
		for(int i=0;i<ntankbranchval;i++){
			if(nuprimarybranchval[i]==1){ primariesinthisentry=true; break; }
		}
		if(!primariesinthisentry){ // not genie primaries... (this shouldn't happen)
			G4cout<<"---------------SKIPPING ENTRY WITH NO GENIE PRIMARIES----------------"<<G4endl;
			inputEntry++;
			localEntry = inputdata->LoadTree(inputEntry);
			if(localEntry<0){
				// get the pointer to the UI manager
				G4UImanager* UI = G4UImanager::GetUIpointer();
				UI->ApplyCommand("/run/abort 0");	// abort without processing current event
				G4cout<<"@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#"<<G4endl;
				G4cout<<"@#@#@#@#@#@#@#@#@#@ REACHED END OF INPUT FILE! #@#@#@#@#@#@#@#@#@#@#@#"<<G4endl;
				G4cout<<"@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#"<<G4endl;
			} else { goto loadbeamentry; } // load the next entry
		}
#endif
#ifdef ONLY_MUONS
		// ensure the event has at least one muon and skip it if not
		Bool_t muonsinthisentry=false;
		for(int i=0;i<ntankbranchval;i++){
			if(pdgbranchval[i]==13){ muonsinthisentry=true; break; }
		}
		if(!muonsinthisentry){ // no muons in the event
			G4cout<<"---------------SKIPPING ENTRY WITH NO MUONS ----------------"<<G4endl;
			inputEntry++;
			localEntry = inputdata->LoadTree(inputEntry);
			if(localEntry<0){
				// get the pointer to the UI manager
				G4UImanager* UI = G4UImanager::GetUIpointer();
				UI->ApplyCommand("/run/abort 0");	// abort without processing current event
				G4cout<<"@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#"<<G4endl;
				G4cout<<"@#@#@#@#@#@#@#@#@#@ REACHED END OF INPUT FILE! #@#@#@#@#@#@#@#@#@#@#@#"<<G4endl;
				G4cout<<"@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#"<<G4endl;
			} else { goto loadbeamentry; } // load the next entry
		}
#endif
		// First the genie information (largely unused as not currently stored in wcsim output)
		// ===========================
#ifndef NO_GENIE
		// genie information available: (not all is currently stored)
		genie::EventRecord* gevtRec = genierecordval->event;
		genie::Interaction* genieint = gevtRec->Summary();
		// we write files with v2_8_6, but read with v2_12_0, so /*V*/ indicates reading is validated
	
		// process information:
//		TString procinfostring = genieint->ProcInfo().AsString();
//		TString scatteringtypestring = genieint->ScatteringTypeAsString();
//		TString interactiontypestring = genieint->InteractionTypeAsString();
//		Bool_t isQE = interaction.ProcInfo().IsQuasiElastic();
//		Double_t neutrinoq2 = genieint->Kine().Q2();
//		TLorentzVector& k1 = *(gevtRec->Probe()->P4());
//		TLorentzVector& k2 = *(gevtRec->FinalStatePrimaryLepton()->P4());
//		Double_t costhfsl = TMath::Cos( k2.Vect().Angle(k1.Vect()) );
//		TLorentzVector* genieVtx = gevtRec->Vertex();
//		G4double x = genieVtx->X() * m;         // same info as nuvtx in g4dirt file
//		G4double y = genieVtx->Y() * m;         // GENIE uses meters
//		G4double z = genieVtx->Z() * m;         // GENIE uses meters
//		G4double t = genieVtx->T() * second;    // GENIE uses seconds for time
		Int_t neutinteractioncode = genie::utils::ghep::NeutReactionCode(gevtRec); /*V*/
//		Int_t nuanceinteractioncode  = genie::utils::ghep::NuanceReactionCode(gevtRec);
		
		// neutrino information:
		Double_t probeenergy = genieint->InitState().ProbeE(genie::kRfLab) * GeV;  /*V*/
//		TSring probepartname = genieint->InitState().Probe()->GetName();
		Int_t probepdg = genieint->InitState().Probe()->PdgCode();                 /*V*/
		TLorentzVector* probemomentum = gevtRec->Probe()->P4();                    /*V*/
		TVector3 probethreemomentum = probemomentum->Vect();
		TVector3 probemomentumdir = probethreemomentum.Unit();
		// n.b.  genieint->InitState().Probe != gevtRec->Probe()
		
		// target nucleon:
//		int targetnucleonpdg = genieint->InitState().Tgt().HitNucPdg();
//		TString targetnucleonname;
//		if ( genie::pdg::IsNeutronOrProton(targetnucleonpdg) ) {
//			TParticlePDG * p = genie::PDGLibrary::Instance()->Find(targetnucleonpdg);
//			targetnucleonname = p->GetName();
//		} else {
//			targetnucleonname = targetnucleonpdg;
//		}
		TLorentzVector* targetnucleonmomentum=0;
		TVector3 targetnucleonthreemomentum(0.,0.,0.);
		Double_t targetnucleonenergy =0;
		if(gevtRec->HitNucleon()){
			targetnucleonmomentum = gevtRec->HitNucleon()->P4();               /*V*/
			targetnucleonthreemomentum = targetnucleonmomentum->Vect();
			targetnucleonenergy = targetnucleonmomentum->Energy() * GeV;
		}
		
		// target nucleus:
		Int_t targetnucleuspdg = genieint->InitState().Tgt().Pdg();                /*V*/
//		TParticlePDG * targetnucleus = genie::PDGLibrary::Instance()->Find( targetnucleuspdg );
//		TString targetnucleusname = "unknown";
//		if(targetnucleus){ targetnucleusname = nucleartarget->GetName(); }
//		Int_t targetnucleusZ = genieint->InitState().Tgt().Z();
//		Int_t targetnucleusA = genieint->InitState().Tgt().A();
	
		// remnant nucleus:
//		int remnucpos = gevtRec->RemnantNucleusPosition(); 
//		TString remnantnucleusname="n/a";
//		Double_t remnantnucleusenergy=-1. * GeV;
//		if(remnucpos>-1){
//			remnantnucleusname = gevtRec->Particle(remnucpos)->Name();
//			remnantnucleusenergy = gevtRec->Particle(remnucpos)->Energy(); //GeV
//		}
	
		// final state lepton:
//		int fsleppos = gevtRec->FinalStatePrimaryLeptonPosition();
//		TString fsleptonname="n/a";
//		Double_t fsleptonenergy=-1. * GeV;
//		if(fsleppos>-1){
//			fsleptonname = gevtRec->Particle(ipos)->Name();
//			fsleptonenergy = gevtRec->Particle(ipos)->Energy();
//		}
	
		// other remnants:
//		Int_t numfsprotons = genieint->ExclTag().NProtons();
//		Int_t numfsneutrons = genieint->ExclTag().NNeutrons();
//		Int_t numfspi0 = genieint->ExclTag().NPi0();
//		Int_t numfspiplus = genieint->ExclTag().NPiPlus();
//		Int_t numfspiminus = genieint->ExclTag().NPiMinus();
		
		//  The following information is retrieved from the PrimaryGeneratorAction in EndOfEventAction:
		//  For each neutrino vertex, the neutrino, target, + any nue, gamma and e daughters are stored.
		///////////////////////
		//	vecRecNumber           // entry number in input file. we should store filename too.
		//	mode;                  // neutrino interaction mode
		//	nvtxs;                 // number of neutrino(?) vertices.
		//	npar;                  // number of primary particles? not used.
		//	vtxsvol[nvtxs]         // neutrino vertex volume indices. not used, looked up in EndOfEventAction.
		//	vtxs[nvtxs]            // neutrino interaction vertices.
		//	beampdgs[nvtxs]        // pdgs of probe neutrinos
		//	beamenergies[nvtxs]    // energies of probe neutrinos, MeV (just as num)
		//	beamdirs[nvtxs]        // directions of probe neutrinos
		//	targetpdgs[nvtxs]      // pdgs of neutrino targets
		//	targetenergies[nvtxs]  // target nucleon energies, MeV (just as num)
		//	targetdirs[nvtxs]      // momentum direction unit vector of target
		///////////////////////
		
		vecRecNumber = genieentrybranchval;
		mode = neutinteractioncode;
		nvtxs = 1;
		npar = -1;                                              // ? not used.
		for(int i=0; i<nvtxs; i++){                             // we only ever have 1 neutrino intx
			vtxsvol[i] = -10;                                   // looked up in EndOfEventAction
			// neutrino vertices are stored in m not cm
			vtxs[i] = G4ThreeVector(nuvtxxval*CLHEP::m, nuvtxyval*CLHEP::m, nuvtxzval*CLHEP::m);
			beampdgs[i] = probepdg;
			beamenergies[i] = probeenergy;
			targetpdgs[i] = targetnucleuspdg;
			targetenergies[i] = targetnucleonenergy;
			G4ThreeVector probemomdir;                          // convert TVector3 to G4ThreeVector
			G4ThreeVector targetnucleonmomdir;
			for(int comp=0; comp<2; comp++){
				probemomdir[comp] = probemomentumdir[comp];
				targetnucleonmomdir[comp] = targetnucleonthreemomentum.Unit()[comp];
			}
			beamdirs[i] = probemomdir;
			targetdirs[i] = targetnucleonmomdir;
		}
#else
		// without genie we don't have the primary interaction information....
		vecRecNumber = genieentrybranchval;
		mode = -999;
		nvtxs = 1;
		npar = -1;                                              // ? not used.
		for(int i=0; i<nvtxs; i++){                             // we only ever have 1 neutrino intx
			vtxsvol[i] = -10;                               // looked up in EndOfEventAction
			vtxs[i] = G4ThreeVector(-999., -999., -999.);   // this is the NEUTRINO info
			beampdgs[i] = -999;                             // NOT THE PRIMARY PARTICLE INFO
			beamenergies[i] = -999.;
			beamdirs[i] = G4ThreeVector(-999., -999., -999.);
			targetpdgs[i] = -999;
			targetenergies[i] = -999.;
			targetdirs[i] = G4ThreeVector(-0., -0., -1.);
		}
#endif
		
		// Now read the outgoing particles: These we will simulate
		// =======================================================
		//G4cout<<"Looping over primaries"<<G4endl;
		for(int i=0;i<ntankbranchval;i++){
			//G4cout<<"Loading details of primary "<<i<<G4endl;
			vtxxval=vtxxbranchval[i]*CLHEP::cm;
			vtxyval=vtxybranchval[i]*CLHEP::cm;
			vtxzval=vtxzbranchval[i]*CLHEP::cm;
			vtxtval=vtxtbranchval[i]*CLHEP::ns;
			pxval=pxbranchval[i]*GeV;
			pyval=pybranchval[i]*GeV;
			pzval=pzbranchval[i]*GeV;
			eval=ebranchval[i]*GeV;
			keval=kebranchval[i]*GeV;
			pdgval=pdgbranchval[i];
			nuprimaryval=nuprimarybranchval[i];
			G4ThreeVector thevtx = G4ThreeVector(vtxxval, vtxyval, vtxzval);
			G4ThreeVector thepdir = G4ThreeVector(pxval, pyval, pzval);
			thepdir.unit();		// normalise to unit vector
			
			// correct for cases where keval == 0 by giving them a little boost
			if(keval==0){keval+=1.*eV; eval+=1.*eV; thepdir=G4ThreeVector(0.,0.,1.);}
			
			//must handle the case of an ion speratly from other particles
			//check PDG code if we have an ion: PDG code format for ions ±10LZZZAAAI
			char strPDG[11];
			char strA[10]={0};
			char strZ[10]={0};
			long int A=0,Z=0;
#ifdef ONLY_MUONS
			if(abs(pdgval)!=13){
				// skip non-muon primary
				continue;
			}
#endif
			if(abs(pdgval) >= 1000000000){
				//ion
				sprintf(strPDG,"%i",abs(pdgval));
				strncpy(strZ, &strPDG[3], 3);
				strncpy(strA, &strPDG[6], 3);
				strA[3]='\0';
				strZ[3]='\0';
				A=atoi(strA);
				Z=atoi(strZ);
				parttype =  G4IonTable::GetIonTable()->GetIon(Z, A, 0.);
				if(parttype){
					particleGun->SetParticleDefinition(parttype);
					particleGun->SetParticleCharge(0);
				} else {
					G4cerr << "skipping primary with PDG " << pdgval << G4endl;
					continue;
				}
			} else {
				//not ion
				parttype = particleTable->FindParticle(pdgval);
				if(parttype){
					particleGun->SetParticleDefinition(parttype);
				} else {
					G4cerr << "skipping primary with PDG " << pdgval << G4endl;
					continue;
				}
			}
			if(nuprimaryval==1){G4cout<<"genie primary ";} else {G4cout<<"genie secondary ";}
			G4cout<<keval/GeV<<" GeV ";
			if(parttype==0){G4cout<<"PDG: "<<pdgval;} else {G4cout<<parttype->GetParticleName();}
			G4Navigator* theNavigator = G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking();
			G4VPhysicalVolume* primaryPV = theNavigator->LocateGlobalPointAndSetup(thevtx); 
			G4cout<<" in "<< primaryPV->GetName()<<G4endl;
			

			particleGun->SetParticleEnergy(keval);       // !!!kinetic!!! energy
			particleGun->SetParticlePosition(thevtx);
			//particleGun->SetParticleTime(vtxtval);     // set event time t=0 for prompt trigger.
			particleGun->SetParticleMomentumDirection(thepdir);
			particleGun->GeneratePrimaryVertex(anEvent); //anEvent provided by G4 when invoking the method
			//G4cout<<"Vertex set"<<G4endl;
		}
		
		// check if this is to be the last event and soft abort the run
		// (abort at the end of this event) if so
		inputEntry++;
		localEntry = inputdata->LoadTree(inputEntry);
		if(localEntry<0){
			// get the pointer to the UI manager
			G4UImanager* UI = G4UImanager::GetUIpointer();
			UI->ApplyCommand("/run/abort 1");	// abort after processing current event
			G4cout<<"@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#"<<G4endl;
			G4cout<<"@#@#@#@#@#@#@#@#@#@ REACHED END OF INPUT FILE! #@#@#@#@#@#@#@#@#@#@#@#"<<G4endl;
			G4cout<<"@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#"<<G4endl;
		}
		
	}

  else if (useGPSEvt)
    {
      MyGPS->GeneratePrimaryVertex(anEvent);

      G4ThreeVector P   =anEvent->GetPrimaryVertex()->GetPrimary()->GetMomentum();
      G4ThreeVector vtx =anEvent->GetPrimaryVertex()->GetPosition();
      G4double m        =anEvent->GetPrimaryVertex()->GetPrimary()->GetMass(); // this is rest mass
      G4int pdg         =anEvent->GetPrimaryVertex()->GetPrimary()->GetPDGcode();
      
      G4ThreeVector dir  = P.unit();
      G4double E         = std::sqrt((P.dot(P))+(m*m));
      
      SetVtx(vtx);
      SetBeamEnergy(E);
      SetBeamDir(dir);
      SetBeamPDG(pdg);
            
      double tote = anEvent->GetPrimaryVertex()->GetPrimary()->GetTotalEnergy();
      double ke = anEvent->GetPrimaryVertex()->GetPrimary()->GetKineticEnergy();
      
      int nprimaryvertices = anEvent->GetNumberOfPrimaryVertex();
      G4cout<<"Generating event with "<<nprimaryvertices<<" primary vertices"<<G4endl;
      for(int evtvtxi=0; evtvtxi<nprimaryvertices; evtvtxi++){
        G4PrimaryVertex* thevertex = anEvent->GetPrimaryVertex(evtvtxi);
        int nparticles = thevertex->GetNumberOfParticle();
        G4cout<<"Primary vertex "<<evtvtxi<<" has "<<nparticles<<" particles"<<G4endl;
        vtx  = thevertex->GetPosition();
        
        for(int parti=0; parti<nparticles; parti++){
          G4PrimaryParticle* theprimary = thevertex->GetPrimary(parti);
          pdg  = theprimary->GetPDGcode();
          tote = theprimary->GetTotalEnergy();
          ke   = theprimary->GetKineticEnergy();
          dir  = theprimary->GetMomentum().unit();
          
          G4ParticleDefinition* parttype = particleTable->FindParticle(pdg);
          G4String particlename;
          particlename = (parttype!=0) ? (std::string(parttype->GetParticleName())) : (std::to_string(pdg));
          G4cout<<"Generating primary "<<particlename<<" with total energy "
                <<tote/MeV<<"MeV and kinetic energy "<<ke/MeV
                <<"MeV at ("<<vtx.x()<<", "<<vtx.y()<<", "<<vtx.z()<<") in direction ("
                <<dir.x()<<", "<<dir.y()<<", "<<dir.z()<<") "<<G4endl;
        }
      }
    } else if (useGenieEvt){

      G4cout <<"Genie + TALYS event generator: Initialising..."<<G4endl;

      if (localEntry!=0) loadNewGenie = false;
      else {
        LoadTalysFiles();    //Load Talys root files once at the beginning!
      }
      G4cout <<"loadNewGenie: "<<loadNewGenie<<G4endl;      

      // Load the next entry, with all required trees and files
      // -----------------------------------------------------
      
      loadgenieentry:
      if(loadNewGenie){
        G4cout <<"Loading new GENIE File"<<G4endl;
        LoadNewGENIEFile(); 
      } // update TChain if a new file is loaded by messenger
      
      G4cout <<"loadNewGenie (after loading GENIE File): "<<loadNewGenie<<G4endl;
     
      //inputdata has already had tree loaded at the end of last event's GeneratePrimaries call
      //localEntry will already be the value of the NEXT entry
      
      G4cout <<"Get tree number..."<<G4endl;
      Int_t nextTreeNumber = geniedata->GetTreeNumber();
      G4cout <<"treeNumber: "<<treeNumber<<", nextTreeNumber: "<<nextTreeNumber<<G4endl;
      if(treeNumber!=nextTreeNumber){
  
        G4cout<< "Reached end of Tree. Last entries' tree number was "
                                                << treeNumber <<", this entries' tree number is "<< nextTreeNumber <<G4endl;
        genieFileName = inputdata->GetCurrentFile()->GetName(); // new tree, new file
        char* genieFileNameAsChar = strdup(genieFileName.c_str());
        genieFileName = basename(genieFileNameAsChar);

        geniedata->SetBranchAddress("iev",&genie_entry,&genie_entry_branch);
        geniedata->SetBranchAddress("neu",&genie_neutrinoflav,&genie_neutrinoflav_branch);
        geniedata->SetBranchAddress("Ev",&genie_neutrinoE,&genie_neutrinoE_branch);
        geniedata->SetBranchAddress("pxv",&genie_neutrinopx,&genie_neutrinopx_branch);
        geniedata->SetBranchAddress("pyv",&genie_neutrinopy,&genie_neutrinopy_branch);
        geniedata->SetBranchAddress("pzv",&genie_neutrinopz,&genie_neutrinopz_branch);
        geniedata->SetBranchAddress("cc",&genie_cc,&genie_cc_branch);
        geniedata->SetBranchAddress("nc",&genie_nc,&genie_nc_branch);
        geniedata->SetBranchAddress("Z",&genie_Z,&genie_Z_branch);
        geniedata->SetBranchAddress("A",&genie_A,&genie_A_branch);
        geniedata->SetBranchAddress("nfn",&genie_final_n,&genie_final_n_branch);
        geniedata->SetBranchAddress("nfp",&genie_final_p,&genie_final_p_branch);
        geniedata->SetBranchAddress("nfpip",&genie_final_pip,&genie_final_pip_branch);
        geniedata->SetBranchAddress("nfpim",&genie_final_pim,&genie_final_pim_branch);
        geniedata->SetBranchAddress("nfpi0",&genie_final_pi0,&genie_final_pi0_branch);
        geniedata->SetBranchAddress("nfkp",&genie_final_kp,&genie_final_kp_branch);
        geniedata->SetBranchAddress("nfkm",&genie_final_km,&genie_final_km_branch);
        geniedata->SetBranchAddress("nfk0",&genie_final_k0,&genie_final_k0_branch);
        geniedata->SetBranchAddress("nf",&genie_num_final,&genie_num_final_branch);
        geniedata->SetBranchAddress("pdgf",&genie_pdg_final,&genie_pdg_final_branch);
        geniedata->SetBranchAddress("Ef",&genie_E_final,&genie_E_final_branch);
        geniedata->SetBranchAddress("pxf",&genie_px_final,&genie_px_final_branch);
        geniedata->SetBranchAddress("pyf",&genie_py_final,&genie_py_final_branch);
        geniedata->SetBranchAddress("pzf",&genie_pz_final,&genie_pz_final_branch);
        geniedata->SetBranchAddress("vtxx",&genie_vtxx,&genie_vtxx_branch);
        geniedata->SetBranchAddress("vtxy",&genie_vtxy,&genie_vtxy_branch);
        geniedata->SetBranchAddress("vtxz",&genie_vtxz,&genie_vtxz_branch);
        geniedata->SetBranchAddress("vtxt",&genie_vtxt,&genie_vtxt_branch);

        //Load GENIE gst branches to check whether they are zombies
        
        if(genie_entry_branch==0 || genie_neutrinoflav_branch==0 || genie_neutrinoE_branch==0 || genie_neutrinopx_branch==0 || genie_neutrinopy_branch==0 || genie_neutrinopz_branch==0 || genie_cc_branch==0 || genie_nc_branch==0 || genie_Z_branch==0 || genie_A_branch==0 || genie_final_n_branch==0 || genie_final_p_branch==0 || genie_final_pip_branch==0 || genie_final_pim_branch==0 || genie_final_pi0_branch==0 || genie_final_kp_branch==0 || genie_final_km_branch==0 || genie_final_k0_branch==0 || genie_pdg_final_branch==0 || genie_E_final_branch==0 || genie_px_final_branch==0 || genie_py_final_branch==0 || genie_pz_final_branch==0 || genie_vtxx_branch==0 || genie_vtxy_branch==0 || genie_vtxz_branch==0){
          G4cout<<"LoadNewGENIEFile: BRANCHES ARE ZOMBIES ARGH!"<<G4endl;
        } else {
	  G4cout <<"Current genie file, num entries: "<<genie_vtxx_branch->GetEntries()<<G4endl;
        }
        entriesInThisTree = genie_vtxx_branch->GetEntries();
        treeNumber=nextTreeNumber;
      }

      G4cout<<"Loading primaries from entry "<<inputEntry<<", localentry "<<localEntry<<"/"<<entriesInThisTree<<G4endl;

      genie_entry_branch->GetEntry(localEntry);
      genie_neutrinoflav_branch->GetEntry(localEntry);
      genie_neutrinoE_branch->GetEntry(localEntry);
      genie_neutrinopx_branch->GetEntry(localEntry);
      genie_neutrinopy_branch->GetEntry(localEntry);
      genie_neutrinopz_branch->GetEntry(localEntry);
      genie_cc_branch->GetEntry(localEntry);
      genie_nc_branch->GetEntry(localEntry);
      genie_Z_branch->GetEntry(localEntry);
      genie_A_branch->GetEntry(localEntry);
      genie_final_n_branch->GetEntry(localEntry);
      genie_final_p_branch->GetEntry(localEntry);
      genie_final_pip_branch->GetEntry(localEntry);
      genie_final_pim_branch->GetEntry(localEntry);
      genie_final_pi0_branch->GetEntry(localEntry);
      genie_final_kp_branch->GetEntry(localEntry);
      genie_final_km_branch->GetEntry(localEntry);
      genie_final_k0_branch->GetEntry(localEntry);
      genie_num_final_branch->GetEntry(localEntry);
      genie_vtxx_branch->GetEntry(localEntry);
      genie_vtxy_branch->GetEntry(localEntry);
      genie_vtxz_branch->GetEntry(localEntry);
      genie_vtxt_branch->GetEntry(localEntry);

      genie_vtxy -= 0.1441;
      genie_vtxz += 1.681;


      G4cout <<"genie_entry: "<<genie_entry<<", genie_neutrinoflav: "<<genie_neutrinoflav<<", genie_neutrinoE: "<<genie_neutrinoE<<", genie_neutrinopx: "<<genie_neutrinopx<<", genie_neutrinopy: "<<genie_neutrinopy<<", genie_neutrinopz: "<<genie_neutrinopz<<", genie_cc: "<<genie_cc<<", genie_nc: "<<genie_nc<<", genie_Z: "<<genie_Z<<", genie_A: "<<genie_A<<", genie_final_n: "<<genie_final_n<<", genie_final_p: "<<genie_final_p<<G4endl;

      //Define neutrino particle properties
      G4ParticleDefinition* parttype = particleTable->FindParticle(genie_neutrinoflav);
      TLorentzVector neutrinovertex(genie_vtxt*CLHEP::s, genie_vtxx*CLHEP::m, genie_vtxy*CLHEP::m, genie_vtxz*CLHEP::m);  // position in m, times in s, convert to cm and ns
      G4cout<<"The origin interaction was a "<<(parttype->GetParticleName())<<" at ("<<genie_vtxt*1000000000.<<","<<genie_vtxx*100.<<","<<genie_vtxy*100.<<","<<genie_vtxz*100.<<")[ns, cm] in (Z,A) = ("<<genie_Z<<","<<genie_A<<")"<<G4endl;
      G4cout<<"This entry has "<<genie_num_final<<" primaries"<<G4endl;

      if (genie_Z == 8 && genie_A == 16){
        G4cout <<"Oxygen nucleus detected! Read in talys-information."<<G4endl;
        G4cout <<"Number of final n: "<<genie_final_n<<", number of final p: "<<genie_final_p<<G4endl;
        G4cout <<"TALYS file name: "<<CalculateResNucleus(genie_final_n,genie_final_p);
      } else {
        G4cout <<"No oxygen nucleus, don't read in TALYS information. (Z = "<<genie_Z<<", A = "<<genie_A<<")"<<G4endl;
      }


      if(pxbranchval){delete[] pxbranchval;}
      if(pybranchval){delete[] pybranchval;}
      if(pzbranchval){delete[] pzbranchval;}
      if(ebranchval){delete[] ebranchval;}
      if(pdgbranchval){delete[] pdgbranchval;}

      pxbranchval = new Double_t[genie_num_final];
      pybranchval = new Double_t[genie_num_final];
      pzbranchval = new Double_t[genie_num_final];
      ebranchval = new Double_t[genie_num_final];
      pdgbranchval = new Int_t[genie_num_final];

      if(pxbranchval==0||pybranchval==0||pzbranchval==0||ebranchval==0||pdgbranchval==0){
        G4cout<<"Arrays are zombies!"<<G4endl;
      }

      G4cout <<"Setting branch addresses"<<G4endl;
      genie_px_final_branch->SetAddress(pxbranchval);
      genie_py_final_branch->SetAddress(pybranchval);
      genie_pz_final_branch->SetAddress(pzbranchval);
      genie_E_final_branch->SetAddress(ebranchval);
      genie_pdg_final_branch->SetAddress(pdgbranchval);

      G4cout <<"Getting entries"<<G4endl;

      genie_px_final_branch->GetEntry(localEntry);
      genie_py_final_branch->GetEntry(localEntry);
      genie_pz_final_branch->GetEntry(localEntry);
      genie_E_final_branch->GetEntry(localEntry);
      genie_pdg_final_branch->GetEntry(localEntry);
      
      double genie_vtxx_temp, genie_vtxy_temp, genie_vtxz_temp, genie_vtxt_temp;

      for (int i_final=0; i_final<genie_num_final; i_final++){

        G4cout <<"Loading details of GENIE primary" << G4endl;
        genie_vtxx_temp=genie_vtxx*100*CLHEP::cm;
        genie_vtxy_temp=genie_vtxy*100*CLHEP::cm;
        genie_vtxz_temp=genie_vtxz*100*CLHEP::cm;
        genie_vtxt_temp=genie_vtxt*CLHEP::ns;
        pxval=pxbranchval[i_final]*GeV;
        pyval=pybranchval[i_final]*GeV;
        pzval=pzbranchval[i_final]*GeV;
        eval=ebranchval[i_final]*GeV;
        pdgval=pdgbranchval[i_final];
        G4ThreeVector thevtx = G4ThreeVector(genie_vtxx_temp, genie_vtxy_temp, genie_vtxz_temp);
        G4ThreeVector thepdir = G4ThreeVector(pxval, pyval, pzval);
        thepdir.unit();         // normalise to unit vector
     
        G4cout <<"thevtx: "<<genie_vtxx<<","<<genie_vtxy<<","<<genie_vtxz<<G4endl;
        G4cout <<"thepdir: "<<pxval<<","<<pyval<<","<<pzval<<G4endl;
        G4cout <<"thepdir.mag(): "<<thepdir.mag()<<G4endl;
        if (thepdir.mag()<1) {
          thepdir.setX(0);
          thepdir.setY(0);
          thepdir.setZ(1);
        }

        //Get PDG number & containers to store possible ion information
        char strPDG[11];
        char strA[10]={0};
        char strZ[10]={0};
        long int A=0,Z=0;
        G4cout <<"PDG value is "<<pdgval<<G4endl;

	//Check if primary particle is ion or not
        if(abs(pdgval) >= 1000000000){
          
          G4cout <<"Ion case!"<<G4endl;
          //ION case
          sprintf(strPDG,"%i",abs(pdgval));
          strncpy(strZ, &strPDG[3], 3);
          strncpy(strA, &strPDG[6], 3);
          strA[3]='\0';
          strZ[3]='\0';
          A=atoi(strA);
          Z=atoi(strZ);
          parttype =  G4IonTable::GetIonTable()->GetIon(Z, A, 0.);
          if(parttype){
            particleGun->SetParticleDefinition(parttype);
            particleGun->SetParticleCharge(0);
          } else {
            G4cerr << "skipping primary with PDG " << pdgval << G4endl;
            continue;
          }
        } else {

          G4cout <<"Non-ION case1"<<G4endl;
          //NON-ION case
          parttype = particleTable->FindParticle(pdgval);
          if(parttype){
            particleGun->SetParticleDefinition(parttype);
          } else {
            G4cerr << "skipping primary with PDG " << pdgval << G4endl;
            continue;
          }
        }

	G4double mass = parttype->GetPDGMass()*0.001*GeV;
        keval = eval - mass;
        G4cout<<keval<<" ("<<keval/GeV<<" GeV) ";
        if(parttype==0){G4cout<<"PDG: "<<pdgval;} else {G4cout<<parttype->GetParticleName();}

        //Get navigator for transportation
        G4Navigator* theNavigator = G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking();
        
        //Get interaction volume
        G4VPhysicalVolume* primaryPV = theNavigator->LocateGlobalPointAndSetup(thevtx);
        G4cout<<" in "<< primaryPV->GetName()<<G4endl;

        if(keval==0){keval+=1.*eV; eval+=1.*eV; thepdir=G4ThreeVector(0.,0.,1.);}
        
	//Set properties of the primary particle
	particleGun->SetParticleEnergy(keval);
        particleGun->SetParticlePosition(thevtx);
        particleGun->SetParticleMomentumDirection(thepdir);

        //Fire the particle away
        particleGun->GeneratePrimaryVertex(anEvent); 
      }      

      G4cout<<"inputEntry: "<<inputEntry<<G4endl;
      inputEntry++;
      localEntry = geniedata->LoadTree(inputEntry);
      G4cout <<"localEntry: "<<localEntry<<G4endl;
      if(localEntry<0){
        // get the pointer to the UI manager
	G4UImanager* UI = G4UImanager::GetUIpointer();
	UI->ApplyCommand("/run/abort 1");	// abort after processing current event
	G4cout<<"@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#"<<G4endl;
	G4cout<<"@#@#@#@#@#@#@#@#@#@ REACHED END OF INPUT FILE! #@#@#@#@#@#@#@#@#@#@#@#"<<G4endl;
	G4cout<<"@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#"<<G4endl;
      }
  

    }
    else if (useAntiNuEvt){

      if ( !inputSpecFile.is_open() )
      {
        G4cout << "Set an energy spectrum file using the command /mygen/Efile name"
             << G4endl;
        exit(-1);
      }

       if (isFirstEvent){
	ReadInEnergySpectrum(inputSpecFile,Espectrum,ProbabilitySpec);

  	G4double me = 0.511;
  	Espectrum_positron.clear();

  	//convert antineutrino energies to positron energies
  	for (unsigned int iE = 0; iE < Espectrum.size(); iE++){
    		double positron_energy = Espectrum.at(iE) - 2*me - 0.782;
    		Espectrum_positron.push_back(positron_energy);
  	}
	isFirstEvent = false;
	}


	theSPSAng->SetAngDistType("iso");
        theSPSAng->SetPosDistribution(theSPSPos);
        thePosition = theSPSPos->GenerateOne();

        G4ParticleDefinition* positron = particleTable->FindParticle(-11);
        G4double KinE = 0.;
        KinE = ShootEnergyPositronCustom()/MeV;
        G4double energy = KinE + positron->GetPDGMass();

        theDirection = theSPSAng->GenerateOne();
        G4double pmom = std::sqrt(pow(energy,2.) - pow(positron->GetPDGMass(),2.));
        G4double px = pmom*theDirection.x();
        G4double py = pmom*theDirection.y();
        G4double pz = pmom*theDirection.z();

        G4PrimaryVertex*   vertex   = new G4PrimaryVertex(thePosition,0);
        G4PrimaryParticle* particle = new G4PrimaryParticle(positron,px,py,pz);     
	vertex->SetPrimary( particle );
	anEvent->AddPrimaryVertex( vertex );

	SetNvtxs(3);	//antineutrino + positron + neutron
      
        //anti-electron neutrino
	SetVtxs(0,thePosition);
        SetBeamEnergy(0.,0);	//dummy value for neutrino
        SetBeamDir(G4ThreeVector(0,0,1),0);  //dummy value for neutrino
        SetBeamPDG(-12,0);      
	
	//positron
	SetVtxs(1,thePosition);
	SetBeamEnergy(energy,1);
	SetBeamDir(theDirection,1);
	SetBeamPDG(-11,1);

	G4cout <<"Anti-Neutrino event: Generating positron at position ("<<thePosition.x()<<","<<thePosition.y()<<","<<thePosition.z()<<", energy "<<energy<<", dir = ("<<theDirection.x()<<","<<theDirection.y()<<","<<theDirection.z()<<")"<<G4endl;

        G4ParticleDefinition* neutron = particleTable->FindParticle(2112);
	KinE = ShootEnergyNeutron()/MeV;
	energy = KinE + neutron->GetPDGMass();
	theDirection = theSPSAng->GenerateOne();
	
	pmom = std::sqrt(pow(energy,2.) - pow(neutron->GetPDGMass(),2.));
	px = pmom*theDirection.x();
	py = pmom*theDirection.y();
	pz = pmom*theDirection.z();

	vertex = new G4PrimaryVertex(thePosition,0);
	particle = new G4PrimaryParticle(neutron,px,py,pz);
	vertex->SetPrimary( particle );
	anEvent->AddPrimaryVertex( vertex );

	//neutron
	SetVtxs(2,thePosition);
	SetBeamEnergy(energy,2);
	SetBeamDir(theDirection,2);
	SetBeamPDG(2112,2);

	G4cout <<"Anti-Neutrino event: Generating neutron at position ("<<thePosition.x()<<","<<thePosition.y()<<","<<thePosition.z()<<", energy "<<energy<<", dir = ("<<theDirection.x()<<","<<theDirection.y()<<","<<theDirection.z()<<")"<<G4endl;
	
	//print properties of generated particles
	double tote = anEvent->GetPrimaryVertex()->GetPrimary()->GetTotalEnergy();
      	double ke = anEvent->GetPrimaryVertex()->GetPrimary()->GetKineticEnergy();

      	int nprimaryvertices = anEvent->GetNumberOfPrimaryVertex();
      	G4cout<<"Generating event with "<<nprimaryvertices<<" primary vertices"<<G4endl;
      	for(int evtvtxi=0; evtvtxi<nprimaryvertices; evtvtxi++){
        	G4PrimaryVertex* thevertex = anEvent->GetPrimaryVertex(evtvtxi);
        	int nparticles = thevertex->GetNumberOfParticle();
        	G4cout<<"Primary vertex "<<evtvtxi<<" has "<<nparticles<<" particles"<<G4endl;
        	G4ThreeVector vtx  = thevertex->GetPosition();

        	for(int parti=0; parti<nparticles; parti++){
          		G4PrimaryParticle* theprimary = thevertex->GetPrimary(parti);
          		G4int pdg  = theprimary->GetPDGcode();
          		tote = theprimary->GetTotalEnergy();
         		ke   = theprimary->GetKineticEnergy();
          		G4ThreeVector dir  = theprimary->GetMomentum().unit();

          		G4ParticleDefinition* parttype = particleTable->FindParticle(pdg);
          		G4String particlename;
          		particlename = (parttype!=0) ? (std::string(parttype->GetParticleName())) : (std::to_string(pdg));
          		G4cout<<"Generating primary "<<particlename<<" with total energy "
                		<<tote/MeV<<"MeV and kinetic energy "<<ke/MeV
                		<<"MeV at ("<<vtx.x()<<", "<<vtx.y()<<", "<<vtx.z()<<") in direction ("
                		<<dir.x()<<", "<<dir.y()<<", "<<dir.z()<<") "<<G4endl;
        	}
      }
    }   
}



void WCSimPrimaryGeneratorAction::SaveOptionsToOutput(WCSimRootOptions * wcopt)
{
  if(useMulineEvt)
    wcopt->SetVectorFileName(vectorFileName);
  else
    wcopt->SetVectorFileName("");
  wcopt->SetGeneratorType(GetGeneratorTypeString());
}

G4String WCSimPrimaryGeneratorAction::GetGeneratorTypeString()
{
  if(useMulineEvt)
    return "muline";
  else if(useGunEvt)
    return "gun";
  else if(useGPSEvt)
    return "gps";
  else if(useLaserEvt)
    return "laser";
  else if(useBeamEvt)
    return "beam";
  else if (useAntiNuEvt)
    return "antinu";
  else if (useGenieEvt)
    return "genie";
  return "";
}

// Returns a vector with the tokens
vector<string> tokenize( string separators, string input ) 
{
  std::size_t startToken = 0, endToken; // Pointers to the token pos
  vector<string> tokens;  // Vector to keep the tokens
  
  if( separators.size() > 0 && input.size() > 0 ) 
    {
    
      while( startToken < input.size() )
	{
	  // Find the start of token
	  startToken = input.find_first_not_of( separators, startToken );
      
	  // If found...
	  if( startToken != input.npos ) 
	    {
	      // Find end of token
	      endToken = input.find_first_of( separators, startToken );
	      if( endToken == input.npos )
		// If there was no end of token, assign it to the end of string
		endToken = input.size();
        
	      // Extract token
	      tokens.push_back( input.substr( startToken, endToken - startToken ) );
        
	      // Update startToken
	      startToken = endToken;
	    }
	}
    }
  
  return tokens;
}

void WCSimPrimaryGeneratorAction::LoadNewPrimaries(){
	if(primariesDirectory==""){
		G4cout<<"No primary files specified! Cannot generate beam events!"<<G4endl;
		assert(false);
	}
	//TODO: should also figure out how to check if loading the tchain is successful
	G4cout<<"loading new primary TChain from: "<<primariesDirectory<<G4endl;
	if(inputdata){ inputdata->ResetBranchAddresses(); delete inputdata; }
	inputdata = new TChain("tankflux");	// input is name of tree in contributing files
	inputdata->Add(primariesDirectory);
	inputdata->LoadTree(0);
	if(metadata){ metadata->ResetBranchAddresses(); delete metadata; }
	metadata = new TChain("tankmeta");
	metadata->Add(primariesDirectory);
	metadata->LoadTree(0);
#ifndef NO_GENIE
	if(geniedata){ geniedata->ResetBranchAddresses(); delete geniedata; }
	geniedata = new TChain("gtree");
	geniedata->Add(neutrinosDirectory);
	geniedata->LoadTree(0);
#endif
	
	inputdata->SetBranchAddress("run",&runbranchval,&runBranch);
	inputdata->SetBranchAddress("ntank",&ntankbranchval,&nTankBranch);
	inputdata->SetBranchAddress("nupdg",&nupdgval,&nupdgBranch);
	inputdata->SetBranchAddress("nuvtxx",&nuvtxxval,&nuvtxxBranch);
	inputdata->SetBranchAddress("nuvtxy",&nuvtxyval,&nuvtxyBranch);
	inputdata->SetBranchAddress("nuvtxz",&nuvtxzval,&nuvtxzBranch);
	inputdata->SetBranchAddress("nuvtxt",&nuvtxtval,&nuvtxtBranch);
	inputdata->SetBranchAddress("vtxvol",&nupvval,&nuPVBranch);
	inputdata->SetBranchAddress("vtxmat",&numatval,&nuvtxmatBranch);
	inputdata->SetBranchAddress("entry",&genieentrybranchval,&genieentryBranch);
	metadata->SetBranchAddress("inputFluxName",&nufluxfilenameval,&nufluxfilenameBranch);
#ifndef NO_GENIE
	geniedata->SetBranchAddress("gmcrec",&genierecordval,&genierecordBranch);
#else 
	genierecordBranch=(TBranch*)1;
#endif
	
	vtxxBranch=inputdata->GetBranch("vx");
	vtxyBranch=inputdata->GetBranch("vy");
	vtxzBranch=inputdata->GetBranch("vz");
	vtxtBranch=inputdata->GetBranch("vt");
	pxBranch=inputdata->GetBranch("px");
	pyBranch=inputdata->GetBranch("py");
	pzBranch=inputdata->GetBranch("pz");
	EBranch=inputdata->GetBranch("E");
	KEBranch=inputdata->GetBranch("kE");
	pdgBranch=inputdata->GetBranch("pdgtank");
	nuprimaryBranch=inputdata->GetBranch("primary");
	
	if(runBranch==0||nTankBranch==0||vtxxBranch==0||vtxyBranch==0||vtxzBranch==0||vtxtBranch==0||pxBranch==0||pyBranch==0||pzBranch==0||EBranch==0||KEBranch==0||pdgBranch==0||nupdgBranch==0||nuvtxxBranch==0||nuvtxyBranch==0||nuvtxzBranch==0||nuvtxtBranch==0||nuPVBranch==0||nuvtxmatBranch==0||nuprimaryBranch==0 || nufluxfilenameBranch==0||genierecordBranch==0){
		G4cout<<"BRANCHES ARE ZOMBIES ARGH!"<<G4endl;
	}
	
	entriesInThisTree = runBranch->GetEntries();
	inputEntry=primariesoffset;
	primariesoffset=0;
	runBranch->GetEntry(inputEntry);
	G4cout<<"first run: "<<runbranchval<<G4endl;
	treeNumber=inputdata->GetTreeNumber();
	localEntry = inputdata->LoadTree(inputEntry);
	dirtFileName = inputdata->GetCurrentFile()->GetName(); // new tree, new file
	char* dirtFileNameAsChar = strdup(dirtFileName.c_str());
	dirtFileName = basename(dirtFileNameAsChar);
	
	loadNewPrimaries=false;
}

void WCSimPrimaryGeneratorAction::LoadNewGENIEFile(){

	//Alternative implementation for loading GENIE files when not using an additional primary file (only information from the ghep.root file)
        //Use gst-converted file version

        G4cout <<"genieDirectory is "<<genieDirectory<<G4endl;

	if(genieDirectory==""){
                G4cout<<"No genie files specified! Cannot generate genie events!"<<G4endl;
                assert(false);
        }

        if(geniedata){ geniedata->ResetBranchAddresses(); delete geniedata; }
        geniedata = new TChain("gst");
        geniedata->Add(genieDirectory);
        geniedata->LoadTree(0);
	

	G4cout <<"Setting branch addresses for geniedata tree"<<G4endl;
	//Set all branch addresses for important variables
	geniedata->SetBranchAddress("iev",&genie_entry,&genie_entry_branch);
	geniedata->SetBranchAddress("neu",&genie_neutrinoflav,&genie_neutrinoflav_branch);
	geniedata->SetBranchAddress("Ev",&genie_neutrinoE,&genie_neutrinoE_branch);
	geniedata->SetBranchAddress("pxv",&genie_neutrinopx,&genie_neutrinopx_branch);
	geniedata->SetBranchAddress("pyv",&genie_neutrinopy,&genie_neutrinopy_branch);
	geniedata->SetBranchAddress("pzv",&genie_neutrinopz,&genie_neutrinopz_branch);
	geniedata->SetBranchAddress("cc",&genie_cc,&genie_cc_branch);
	geniedata->SetBranchAddress("nc",&genie_nc,&genie_nc_branch);
	geniedata->SetBranchAddress("Z",&genie_Z,&genie_Z_branch);
	geniedata->SetBranchAddress("A",&genie_A,&genie_A_branch);
	geniedata->SetBranchAddress("nfn",&genie_final_n,&genie_final_n_branch);
	geniedata->SetBranchAddress("nfp",&genie_final_p,&genie_final_p_branch);
	geniedata->SetBranchAddress("nfpip",&genie_final_pip,&genie_final_pip_branch);
	geniedata->SetBranchAddress("nfpim",&genie_final_pim,&genie_final_pim_branch);
	geniedata->SetBranchAddress("nfpi0",&genie_final_pi0,&genie_final_pi0_branch);
	geniedata->SetBranchAddress("nfkp",&genie_final_kp,&genie_final_kp_branch);
	geniedata->SetBranchAddress("nfkm",&genie_final_km,&genie_final_km_branch);
	geniedata->SetBranchAddress("nfk0",&genie_final_k0,&genie_final_k0_branch);
	geniedata->SetBranchAddress("nf",&genie_num_final,&genie_num_final_branch);
	geniedata->SetBranchAddress("pdgf",&genie_pdg_final,&genie_pdg_final_branch);
	geniedata->SetBranchAddress("Ef",&genie_E_final,&genie_E_final_branch);
	geniedata->SetBranchAddress("pxf",&genie_px_final,&genie_px_final_branch);
	geniedata->SetBranchAddress("pyf",&genie_py_final,&genie_py_final_branch);
	geniedata->SetBranchAddress("pzf",&genie_pz_final,&genie_pz_final_branch);
	geniedata->SetBranchAddress("vtxx",&genie_vtxx,&genie_vtxx_branch);
	geniedata->SetBranchAddress("vtxy",&genie_vtxy,&genie_vtxy_branch);
	geniedata->SetBranchAddress("vtxz",&genie_vtxz,&genie_vtxz_branch);
	geniedata->SetBranchAddress("vtxt",&genie_vtxt,&genie_vtxt_branch);

	//Load GENIE gst branches to check whether they are zombies

	if(genie_entry_branch==0 || genie_neutrinoflav_branch==0 || genie_neutrinoE_branch==0 || genie_neutrinopx_branch==0 || genie_neutrinopy_branch==0 || genie_neutrinopz_branch==0 || genie_cc_branch==0 || genie_nc_branch==0 || genie_Z_branch==0 || genie_A_branch==0 || genie_final_n_branch==0 || genie_final_p_branch==0 || genie_final_pip_branch==0 || genie_final_pim_branch==0 || genie_final_pi0_branch==0 || genie_final_kp_branch==0 || genie_final_km_branch==0 || genie_final_k0_branch==0 || genie_pdg_final_branch==0 || genie_E_final_branch==0 || genie_px_final_branch==0 || genie_py_final_branch==0 || genie_pz_final_branch==0 || genie_vtxx_branch==0 || genie_vtxy_branch==0 || genie_vtxz_branch==0){
                G4cout<<"LoadNewGENIEFile: BRANCHES ARE ZOMBIES ARGH!"<<G4endl;
        }

        G4cout <<"Getting entries in tree..."<<G4endl;
	entriesInThisTree=geniedata->GetEntries();
        G4cout <<"Entries in this tree: "<<entriesInThisTree<<G4endl;
	inputEntry = primariesoffset;
	primariesoffset=0;
	//geniedata->GetEntry(inputEntry);
        
        G4cout<<"first run: "<<runbranchval<<G4endl;
        treeNumber=geniedata->GetTreeNumber();
        localEntry = geniedata->LoadTree(inputEntry);
        genieFileName = geniedata->GetCurrentFile()->GetName(); // new tree, new file
        std::cout <<"Current genieFileName: "<<genieFileName<<std::endl;
        char* genieFileNameAsChar = strdup(genieFileName.c_str());
        genieFileName = basename(genieFileNameAsChar);	

	loadNewGenie=false;

}

void WCSimPrimaryGeneratorAction::LoadTalysFiles(){

	//Load nuclear de-excitation information from TALYS-generated root-files
	//Version of talys file depends on how many neutrons & protons have been knocked out of the nucleus
	//
	
        if(talysDirectory==""){
		G4cout<<"No talys files specified! Cannot look-up nuclear de-excitation modes!"<<G4endl;
		assert(false);
	}

        std::string filename_O15 = "O15.root";
        std::string filename_N15 = "N15.root";
        std::string filename_N14 = "N14.root";
        std::string filename_Li9 = "Li9.root";
        std::string filename_Li7 = "Li7.root";
        std::string filename_C14 = "C14.root";
        std::string filename_C13 = "C13.root";
        std::string filename_C11 = "C11.root";
        std::string filename_C10 = "C10.root";
        std::string filename_Be10 = "Be10.root";
        std::string filename_Be9 = "Be9.root";
        std::string filename_B11 = "B11.root";
        std::string filename_B10 = "B10.root";
        std::string filename_B9 = "B9.root";
        
        std::string filename_O15gamma = "O15gamma.root";
        std::string filename_N15gamma = "N15gamma.root";
        std::string filename_N14gamma = "N14gamma.root";
        std::string filename_Li9gamma = "Li9gamma.root";
        std::string filename_Li7gamma = "Li7gamma.root";
        std::string filename_C14gamma = "C14gamma.root";
        std::string filename_C13gamma = "C13gamma.root";
        std::string filename_C11gamma = "C11gamma.root";
        std::string filename_C10gamma = "C10gamma.root";
        std::string filename_Be10gamma = "Be10gamma.root";
        std::string filename_Be9gamma = "Be9gamma.root";
        std::string filename_B11gamma = "B11gamma.root";
        std::string filename_B10gamma = "B10gamma.root";
        std::string filename_B9gamma = "B9gamma.root";
	
        std::string filepath_O15 =talysDirectory + filename_O15;
        std::string filepath_N15 =talysDirectory + filename_N15;
        std::string filepath_N14 =talysDirectory + filename_N14;
        std::string filepath_Li9 =talysDirectory + filename_Li9;
        std::string filepath_Li7 =talysDirectory + filename_Li7;
        std::string filepath_C14 =talysDirectory + filename_C14;
        std::string filepath_C13 =talysDirectory + filename_C13;
        std::string filepath_C11 =talysDirectory + filename_C11;
        std::string filepath_C10 =talysDirectory + filename_C10;
        std::string filepath_Be10 =talysDirectory + filename_Be10;
        std::string filepath_Be9 =talysDirectory + filename_Be9;
        std::string filepath_B11 =talysDirectory + filename_B11;
        std::string filepath_B10 =talysDirectory + filename_B10;
        std::string filepath_B9 =talysDirectory + filename_B9;
        
        std::string filepath_O15gamma =talysDirectory + filename_O15gamma;
        std::string filepath_N15gamma =talysDirectory + filename_N15gamma;
        std::string filepath_N14gamma =talysDirectory + filename_N14gamma;
        std::string filepath_Li9gamma =talysDirectory + filename_Li9gamma;
        std::string filepath_Li7gamma =talysDirectory + filename_Li7gamma;
        std::string filepath_C14gamma =talysDirectory + filename_C14gamma;
        std::string filepath_C13gamma =talysDirectory + filename_C13gamma;
        std::string filepath_C11gamma =talysDirectory + filename_C11gamma;
        std::string filepath_C10gamma =talysDirectory + filename_C10gamma;
        std::string filepath_Be10gamma =talysDirectory + filename_Be10gamma;
        std::string filepath_Be9gamma =talysDirectory + filename_Be9gamma;
        std::string filepath_B11gamma =talysDirectory + filename_B11gamma;
        std::string filepath_B10gamma =talysDirectory + filename_B10gamma;
        std::string filepath_B9gamma =talysDirectory + filename_B9gamma;

        f_O15 = new TFile(filepath_O15.c_str(),"READ");
        f_N15 = new TFile(filepath_N15.c_str(),"READ");
        f_N14 = new TFile(filepath_N14.c_str(),"READ");
        f_Li9 = new TFile(filepath_Li9.c_str(),"READ");
        f_Li7 = new TFile(filepath_Li7.c_str(),"READ");
        f_C14 = new TFile(filepath_C14.c_str(),"READ");
        f_C13 = new TFile(filepath_C13.c_str(),"READ");
        f_C11 = new TFile(filepath_C11.c_str(),"READ");
        f_C10 = new TFile(filepath_C10.c_str(),"READ");
        f_Be10 = new TFile(filepath_Be10.c_str(),"READ");
        f_Be9 = new TFile(filepath_Be9.c_str(),"READ");
        f_B11 = new TFile(filepath_B11.c_str(),"READ");
        f_B10 = new TFile(filepath_B10.c_str(),"READ");
        f_B9 = new TFile(filepath_B9.c_str(),"READ");
        
        f_O15gamma = new TFile(filepath_O15gamma.c_str(),"READ");
        f_N15gamma = new TFile(filepath_N15gamma.c_str(),"READ");
        f_N14gamma = new TFile(filepath_N14gamma.c_str(),"READ");
        f_Li9gamma = new TFile(filepath_Li9gamma.c_str(),"READ");
        f_Li7gamma = new TFile(filepath_Li7gamma.c_str(),"READ");
        f_C14gamma = new TFile(filepath_C14gamma.c_str(),"READ");
        f_C13gamma = new TFile(filepath_C13gamma.c_str(),"READ");
        f_C11gamma = new TFile(filepath_C11gamma.c_str(),"READ");
        f_C10gamma = new TFile(filepath_C10gamma.c_str(),"READ");
        f_Be10gamma = new TFile(filepath_Be10gamma.c_str(),"READ");
        f_Be9gamma = new TFile(filepath_Be9gamma.c_str(),"READ");
        f_B11gamma = new TFile(filepath_B11gamma.c_str(),"READ");
        f_B10gamma = new TFile(filepath_B10gamma.c_str(),"READ");
        f_B9gamma = new TFile(filepath_B9gamma.c_str(),"READ");
        

	talys_O15 = (TTree*) f_O15->Get("TreeNucldeex");
        talys_O15->SetName("talys_O15");
	talys_N15 = (TTree*) f_N15->Get("TreeNucldeex");
        talys_N15->SetName("talys_N15");
	talys_N14 = (TTree*) f_N14->Get("TreeNucldeex");
        talys_N14->SetName("talys_N14");
	talys_Li9 = (TTree*) f_Li9->Get("TreeNucldeex");
        talys_Li9->SetName("talys_Li9");
	talys_Li7 = (TTree*) f_Li7->Get("TreeNucldeex");
        talys_Li7->SetName("talys_Li7");
	talys_C14 = (TTree*) f_C14->Get("TreeNucldeex");
        talys_C14->SetName("talys_C14");
	talys_C13 = (TTree*) f_C13->Get("TreeNucldeex");
        talys_C13->SetName("talys_C13");
	talys_C11 = (TTree*) f_C11->Get("TreeNucldeex");
        talys_C11->SetName("talys_C11");
	talys_C10 = (TTree*) f_C10->Get("TreeNucldeex");
        talys_C10->SetName("talys_C10");
        talys_Be10 = (TTree*) f_Be10->Get("TreeNucldeex");
        talys_Be10->SetName("talys_Be10");
        talys_Be9 = (TTree*) f_Be9->Get("TreeNucldeex");
        talys_Be9->SetName("talys_Be9");
        talys_B11 = (TTree*) f_B11->Get("TreeNucldeex");
        talys_B11->SetName("talys_B11");
        talys_B10 = (TTree*) f_B10->Get("TreeNucldeex");
        talys_B10->SetName("talys_B10");
        talys_B9 = (TTree*) f_B9->Get("TreeNucldeex");
        talys_B9->SetName("talys_B9");
	
        talys_O15gamma = (TTree*) f_O15gamma->Get("TreeNucldeex");
        talys_O15gamma->SetName("talys_O15gamma");
	talys_N15gamma = (TTree*) f_N15gamma->Get("TreeNucldeex");
        talys_N15gamma->SetName("talys_N15gamma");
	talys_N14gamma = (TTree*) f_N14gamma->Get("TreeNucldeex");
        talys_N14gamma->SetName("talys_N14gamma");
	talys_Li9gamma = (TTree*) f_Li9gamma->Get("TreeNucldeex");
        talys_Li9gamma->SetName("talys_Li9gamma");
	talys_Li7gamma = (TTree*) f_Li7gamma->Get("TreeNucldeex");
        talys_Li7gamma->SetName("talys_Li7gamma");
	talys_C14gamma = (TTree*) f_C14gamma->Get("TreeNucldeex");
        talys_C14gamma->SetName("talys_C14gamma");
	talys_C13gamma = (TTree*) f_C13gamma->Get("TreeNucldeex");
        talys_C13gamma->SetName("talys_C13gamma");
	talys_C11gamma = (TTree*) f_C11gamma->Get("TreeNucldeex");
        talys_C11gamma->SetName("talys_C11gamma");
	talys_C10gamma = (TTree*) f_C10gamma->Get("TreeNucldeex");
        talys_C10gamma->SetName("talys_C10gamma");
        talys_Be10gamma = (TTree*) f_Be10gamma->Get("TreeNucldeex");
        talys_Be10gamma->SetName("talys_Be10gamma");
        talys_Be9gamma = (TTree*) f_Be9gamma->Get("TreeNucldeex");
        talys_Be9gamma->SetName("talys_Be9gamma");
        talys_B11gamma = (TTree*) f_B11gamma->Get("TreeNucldeex");
        talys_B11gamma->SetName("talys_B11gamma");
        talys_B10gamma = (TTree*) f_B10gamma->Get("TreeNucldeex");
        talys_B10gamma->SetName("talys_B10gamma");
        talys_B9gamma = (TTree*) f_B9gamma->Get("TreeNucldeex");
        talys_B9gamma->SetName("talys_B9gamma");

        G4cout <<"Loaded complete de-excitation information for relevenat nuclei from talys *gamma-files."<<G4endl;

}

G4double WCSimPrimaryGeneratorAction::ShootEnergyNeutron() {

	G4double sum = 0.;
	G4double norm = 0.;
	G4double Probability [51];
	G4double EnergyBin [51];
	G4double deltaX, x, y;
	
	G4double ANeutronEnergySp [51] =    
 	{ 0., 274.4119, 342.3826, 257.9844, 187.8962, 141.1938,
   	107.2338,  78.0118,  62.4729,  45.4234,  37.3674,
   	28.1552,  27.4770,  17.0363,  17.4592,   9.9503,
   	9.9544,   7.9960,   4.7464,   5.1454,   3.6949,
   	2.1576,   1.7950,   1.1628,   1.1344,   0.9440,
   	0.5982,   0.4545,   0.3037,   0.3071,   0.2996,
   	0.1526,   0.1976,   0.0576,   0.1063,   0.0565,
   	0.0597,   0.0334,   0.0270,   0.0326,   0.0113,
   	0.0121,   0.0025,   0.0063,   0.0022,   0.0011,
   	0.0000,   0.0000,   0.0000,   0.0000,   0.0000};

	norm= 0;

     	for(G4int i = 0; i< 51 ; i++) {
      		norm += ANeutronEnergySp[i];
     	}

    	for(G4int i = 0; i< 51 ; i++) {
      		sum += ANeutronEnergySp[i]/norm;
      		Probability[i]=sum;
      		EnergyBin[i]=G4float(i)*4;
		//energy range of emitted neutrons will be 0MeV - 0.2 MeV (4*50MeV/1000)
    	}

  	G4double val = G4UniformRand();

  	for(G4int i=0;i<51;i++) {
    		if(Probability[i] >= val) {
	    		if(i == 0) {  
		    		return EnergyBin[0]/1000;
	    		}
			deltaX = val - Probability[i] ;
      			y = EnergyBin[i] - EnergyBin[i-1] ;
      			x = Probability[i] - Probability[i-1] ;
      			return ((deltaX*y/x + EnergyBin[i])/1000) ;
    		}
  	}

	return 0;
}


G4double WCSimPrimaryGeneratorAction::ShootEnergyPositronCustom(){


	G4double deltaX, x, y;	
	G4double val = G4UniformRand();

	for (G4int i = 0; i < Espectrum_positron.size(); i++){
		if (ProbabilitySpec.at(i) >= val){
			if (i == 0) return Espectrum_positron.at(0);
			deltaX = val - ProbabilitySpec.at(i);
			y = Espectrum_positron.at(i) - Espectrum_positron.at(i-1);
			x = ProbabilitySpec.at(i) - ProbabilitySpec.at(i-1);
			return (deltaX*y/x + Espectrum_positron.at(i));
		}
	}

}

std::string WCSimPrimaryGeneratorAction::CalculateResNucleus(int num_n, int num_p){

	std::string res_nucleus = "none";

        if (num_n == 1 && num_p == 0) res_nucleus = "O15";
	if (num_n == 0 && num_p == 1) res_nucleus = "N15";
	if (num_n == 1 && num_p == 1) res_nucleus = "N14";
	if (num_n == 0 && num_p == 2) res_nucleus = "C14";
	if (num_n == 1 && num_p == 2) res_nucleus = "C13";
	if (num_n == 3 && num_p == 2) res_nucleus = "C11";
	if (num_n == 4 && num_p == 2) res_nucleus = "C10";
	if (num_n == 2 && num_p == 3) res_nucleus = "B11";
	if (num_n == 3 && num_p == 3) res_nucleus = "B10";
	if (num_n == 4 && num_p == 3) res_nucleus = "B9";
	if (num_n == 2 && num_p == 4) res_nucleus = "Be10";
	if (num_n == 3 && num_p == 4) res_nucleus = "Be9";
	if (num_n == 2 && num_p == 5) res_nucleus = "Li9";
	if (num_n == 4 && num_p == 5) res_nucleus = "Li7";

	return res_nucleus;

}
