#include "WCSimPrimaryGeneratorAction.hh"
#include "WCSimDetectorConstruction.hh"
#include "WCSimPrimaryGeneratorMessenger.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4GeneralParticleSource.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4ThreeVector.hh"
#include "TLorentzVector.h"
#include "globals.hh"
#include "Randomize.hh"
#include <fstream>
#include <vector>
#include <string>
#include "G4SystemOfUnits.hh"
#include <math.h> 

#include "G4Navigator.hh"
#include "G4TransportationManager.hh"

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

inline float atof( const string& s ) {return std::atof( s.c_str() );}
inline int   atoi( const string& s ) {return std::atoi( s.c_str() );}

WCSimPrimaryGeneratorAction::WCSimPrimaryGeneratorAction(
					  WCSimDetectorConstruction* myDC)
  :myDetector(myDC)
{
  //T. Akiri: Initialize GPS to allow for the laser use 
  MyGPS = new G4GeneralParticleSource();

  // Initialize to zero
  mode = 0;
  vtxvol = 0;
  vtx = G4ThreeVector(0.,0.,0.);
  nuEnergy = 0.;
  _counterRock=0; // counter for generated in Rock
  _counterCublic=0; // counter generated
  
  //---Set defaults. Do once at beginning of session.
  
  G4int n_particle = 1;
  particleGun = new G4ParticleGun(n_particle);
  particleGun->SetParticleEnergy(1.0*CLHEP::GeV);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.0));

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  particleGun->
    SetParticleDefinition(particleTable->FindParticle(particleName="mu+"));

  particleGun->
    SetParticlePosition(G4ThreeVector(0.*CLHEP::m,0.*CLHEP::m,0.*CLHEP::m));
    
  messenger = new WCSimPrimaryGeneratorMessenger(this);
  useMulineEvt = false;
  useNormalEvt = false;
  useBeamEvt = true;
  
  if(useBeamEvt){
	inputdata = new TChain("tankflux");	// input is name of tree in contributing files
	//inputdata->Add("/pnfs/annie/persistent/users/moflaher/g4dirt/annie_tank_flux.*.root");
	inputdata->Add("/annie/app/users/edrakopo/annie_tank_flux.2000.root");
	inputdata->LoadTree(0);
	
	inputdata->SetBranchAddress("run",&runbranchval,&runBranch);
	inputdata->SetBranchAddress("ntank",&ntankbranchval,&nTankBranch);
	inputdata->SetBranchAddress("nupdg",&nupdgval,&nupdgBranch);
	inputdata->SetBranchAddress("nuvtxx",&nuvtxxval,&nuvtxxBranch);
	inputdata->SetBranchAddress("nuvtxy",&nuvtxyval,&nuvtxyBranch);
	inputdata->SetBranchAddress("nuvtxz",&nuvtxzval,&nuvtxzBranch);
	inputdata->SetBranchAddress("nuvtxt",&nuvtxtval,&nuvtxtBranch);
	inputdata->SetBranchAddress("vtxvol",&nupvval,&nuPVBranch);
	inputdata->SetBranchAddress("vtxmat",&numatval,&nuvtxmatBranch);
	entriesInThisTree = runBranch->GetEntries();
	
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
	geniePrimaryBranch=inputdata->GetBranch("primary");
	
	
	if(runBranch==0||nTankBranch==0||vtxxBranch==0||vtxyBranch==0||vtxzBranch==0||vtxtBranch==0||pxBranch==0||pyBranch==0||pzBranch==0||EBranch==0||KEBranch==0||pdgBranch==0||nupdgBranch==0||nuvtxxBranch==0||nuvtxyBranch==0||nuvtxzBranch==0||nuvtxtBranch==0||nuPVBranch==0||nuvtxmatBranch==0||geniePrimaryBranch==0){
		G4cout<<"BRANCHES ARE ZOMBIES ARGH!"<<G4endl;
	}
	
	inputEntry=0;
	runBranch->GetEntry(inputEntry);
	G4cout<<"first run: "<<runbranchval<<G4endl;
	treeNumber=inputdata->GetTreeNumber();
	
//	inputfilereader = new TTreeReader(inputdata);
//	Int_t numinputevents = inputfilereader->GetEntries(true);
//	G4cout<<numinputevents<<" events in the input file."<<G4endl;
//	TTreeReader::EEntryStatus entrystats = inputfilereader->SetEntry(0);
//	G4cout<<"entrystat: "<<entrystats<<G4endl;
//	Bool_t nextsuccess = inputfilereader->Next();
//	if(nextsuccess){G4cout<<"success"<<G4endl;} else {G4cout<<"failure"<<G4endl;}
//	
//	TTreeReaderValue<Int_t> rdrrun((*inputfilereader), "run");
//	TTreeReaderValue<Int_t> rdrentry((*inputfilereader), "entry");
//	TTreeReaderValue<Int_t> rdriter((*inputfilereader), "iter");
//	TTreeReaderValue<Int_t> rdrniter((*inputfilereader), "niter");
//	TTreeReaderValue<Int_t> rdrnupdg((*inputfilereader), "nupdg");
//	TTreeReaderValue<Double_t> rdrnuvtxx((*inputfilereader), "nuvtxx");
//	TTreeReaderValue<Double_t> rdrnuvtxy((*inputfilereader), "nuvtxy");
//	TTreeReaderValue<Double_t> rdrnuvtxz((*inputfilereader), "nuvtxz");
//	TTreeReaderValue<Double_t> rdrnuvtxt((*inputfilereader), "nuvtxt");
//	TTreeReaderValue<Int_t> rdrintank((*inputfilereader), "intank");
//	TTreeReaderValue<Int_t> rdrinhall((*inputfilereader), "inhall");
//	TTreeReaderValue<Char_t> rdrvtxvol((*inputfilereader), "vtxvol");
//	TTreeReaderValue<Char_t> rdrvtxmat((*inputfilereader), "vtxmat");
//	TTreeReaderValue<Int_t> rdrntank((*inputfilereader), "ntank");
//	TTreeReaderArray<Int_t> rdrpdgtank((*inputfilereader), "pdgtank");
//	TTreeReaderArray<Int_t> rdrprimary((*inputfilereader), "primary");
//	TTreeReaderArray<Double_t> rdrvx((*inputfilereader), "vx");
//	TTreeReaderArray<Double_t> rdrvy((*inputfilereader), "vy");
//	TTreeReaderArray<Double_t> rdrvz((*inputfilereader), "vz");
//	TTreeReaderArray<Double_t> rdrvt((*inputfilereader), "vt");
//	TTreeReaderArray<Double_t> rdrpx((*inputfilereader), "px");
//	TTreeReaderArray<Double_t> rdrpy((*inputfilereader), "py");
//	TTreeReaderArray<Double_t> rdrpz((*inputfilereader), "pz");
//	TTreeReaderArray<Double_t> rdrE((*inputfilereader), "E");
//	TTreeReaderArray<Double_t> rdrkE((*inputfilereader), "kE");
//	
//	entrystats = inputfilereader->SetEntry(0);
//	G4cout<<"entrystat: "<<entrystats<<G4endl;
//	nextsuccess = inputfilereader->Next();
//	if(nextsuccess){G4cout<<"success"<<G4endl;} else {G4cout<<"failure"<<G4endl;}
//	
//	if(rdrrun.GetSetupStatus()<0) G4cout<<rdrrun.GetSetupStatus()<<G4endl;
//	if(rdrentry.GetSetupStatus()<0) G4cout<<rdrentry.GetSetupStatus()<<G4endl;
//	if(rdriter.GetSetupStatus()<0) G4cout<<rdriter.GetSetupStatus()<<G4endl;
//	if(rdrniter.GetSetupStatus()<0) G4cout<<rdrniter.GetSetupStatus()<<G4endl;
//	if(rdrnupdg.GetSetupStatus()<0) G4cout<<rdrnupdg.GetSetupStatus()<<G4endl;
//	if(rdrnuvtxx.GetSetupStatus()<0) G4cout<<rdrnuvtxx.GetSetupStatus()<<G4endl;
//	if(rdrnuvtxy.GetSetupStatus()<0) G4cout<<rdrnuvtxy.GetSetupStatus()<<G4endl;
//	if(rdrnuvtxz.GetSetupStatus()<0) G4cout<<rdrnuvtxz.GetSetupStatus()<<G4endl;
//	if(rdrnuvtxt.GetSetupStatus()<0) G4cout<<rdrnuvtxt.GetSetupStatus()<<G4endl;
//	if(rdrintank.GetSetupStatus()<0) G4cout<<rdrintank.GetSetupStatus()<<G4endl;
//	if(rdrinhall.GetSetupStatus()<0) G4cout<<rdrinhall.GetSetupStatus()<<G4endl;
//	if(rdrvtxvol.GetSetupStatus()<0) G4cout<<rdrvtxvol.GetSetupStatus()<<G4endl;
//	if(rdrvtxmat.GetSetupStatus()<0) G4cout<<rdrvtxmat.GetSetupStatus()<<G4endl;
//	if(rdrntank.GetSetupStatus()<0) G4cout<<rdrntank.GetSetupStatus()<<G4endl;
//	if(rdrpdgtank.GetSetupStatus()<0) G4cout<<rdrpdgtank.GetSetupStatus()<<G4endl;
//	if(rdrprimary.GetSetupStatus()<0) G4cout<<rdrprimary.GetSetupStatus()<<G4endl;
//	if(rdrvx.GetSetupStatus()<0) G4cout<<rdrvx.GetSetupStatus()<<G4endl;
//	if(rdrvy.GetSetupStatus()<0) G4cout<<rdrvy.GetSetupStatus()<<G4endl;
//	if(rdrvz.GetSetupStatus()<0) G4cout<<rdrvz.GetSetupStatus()<<G4endl;
//	if(rdrvt.GetSetupStatus()<0) G4cout<<rdrvt.GetSetupStatus()<<G4endl;
//	if(rdrpx.GetSetupStatus()<0) G4cout<<rdrpx.GetSetupStatus()<<G4endl;
//	if(rdrpy.GetSetupStatus()<0) G4cout<<rdrpy.GetSetupStatus()<<G4endl;
//	if(rdrpz.GetSetupStatus()<0) G4cout<<rdrpz.GetSetupStatus()<<G4endl;
//	if(rdrE.GetSetupStatus()<0) G4cout<<rdrE.GetSetupStatus()<<G4endl;
//	if(rdrkE.GetSetupStatus()<0) G4cout<<rdrkE.GetSetupStatus()<<G4endl;
//	
//	rdrrunp = &rdrrun;
//	rdrentryp = &rdrentry;
//	rdriterp = &rdriter;
//	rdrniterp = &rdrniter;
//	rdrnupdgp = &rdrnupdg;
//	rdrnuvtxxp = &rdrnuvtxx;
//	rdrnuvtxyp = &rdrnuvtxy;
//	rdrnuvtxtp = &rdrnuvtxt;
//	rdrintankp = &rdrintank;
//	rdrinhallp = &rdrinhall;
//	rdrvtxvolp = &rdrvtxvol;
//	rdrntankp = &rdrntank;
//	rdrpdgtankp = &rdrpdgtank;
//	rdrprimaryp = &rdrprimary;
//	rdrvxp = &rdrvx;
//	rdrvyp = &rdrvy;
//	rdrvzp = &rdrvz;
//	rdrvtp = &rdrvt;
//	rdrpxp = &rdrpx;
//	rdrpyp = &rdrpy;
//	rdrpzp = &rdrpz;
//	rdrEp = &rdrE;
//	rdrkEp = &rdrkE;
	
  }
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
  delete particleGun;
  delete MyGPS;   //T. Akiri: Delete the GPS variable
  delete messenger;
  
  if(useBeamEvt){
    inputdata->ResetBranchAddresses();
    delete inputdata;
    delete inputfilereader;
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
	    vtx = G4ThreeVector(atof(token[1])*CLHEP::cm,
				atof(token[2])*CLHEP::cm,
				atof(token[3])*CLHEP::cm);
	    
            // true : Generate vertex in Rock , false : Generate vertex in WC tank
            SetGenerateVertexInRock(false);

	    // Next we read the incoming neutrino and target
	    
	    // First, the neutrino line

	    token=readInLine(inputFile, lineSize, inBuf);
	    beampdg = atoi(token[1]);
	    beamenergy = atof(token[2])*CLHEP::MeV;
	    beamdir = G4ThreeVector(atof(token[3]),
				    atof(token[4]),
				    atof(token[5]));

	    // Now read the target line

	    token=readInLine(inputFile, lineSize, inBuf);
	    targetpdg = atoi(token[1]);
	    targetenergy = atof(token[2])*CLHEP::MeV;
	    targetdir = G4ThreeVector(atof(token[3]),
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
		    particleGun->
		      SetParticleDefinition(particleTable->
					    FindParticle(pdgid));
		    G4double mass = 
		      particleGun->GetParticleDefinition()->GetPDGMass();

		    G4double ekin = energy - mass;

		    particleGun->SetParticleEnergy(ekin);
		    //G4cout << "Particle: " << pdgid << " KE: " << ekin << G4endl;
		    particleGun->SetParticlePosition(vtx);
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

  else if (useNormalEvt)
  {      // manual gun operation
    particleGun->GeneratePrimaryVertex(anEvent);

    G4ThreeVector P  =anEvent->GetPrimaryVertex()->GetPrimary()->GetMomentum();
    G4ThreeVector vtx=anEvent->GetPrimaryVertex()->GetPosition();
    G4double m       =anEvent->GetPrimaryVertex()->GetPrimary()->GetMass();
    G4int pdg        =anEvent->GetPrimaryVertex()->GetPrimary()->GetPDGcode();

    G4ThreeVector dir  = P.unit();
    G4double E         = std::sqrt((P.dot(P))+(CLHEP::m*CLHEP::m));

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
      //T. Akiri: Create the GPS LASER event
      MyGPS->GeneratePrimaryVertex(anEvent);
      
      G4ThreeVector P   =anEvent->GetPrimaryVertex()->GetPrimary()->GetMomentum();
      G4ThreeVector vtx =anEvent->GetPrimaryVertex()->GetPosition();
      G4int pdg         =anEvent->GetPrimaryVertex()->GetPrimary()->GetPDGcode();
      
      G4ThreeVector dir  = P.unit();
      G4double E         = std::sqrt((P.dot(P)));
      
      SetVtx(vtx);
      SetBeamEnergy(E);
      SetBeamDir(dir);
      SetBeamPDG(pdg);
    }
  
  else if (useBeamEvt)
    {   // read in from ROOT file output of rob hatcher's annie dirt simulations
//  	if(!inputfilereader->Next()){
//  		G4cout<<" END OF INPUT FILE!!!!!!"<<G4endl;
//  		inputfilereader->SetEntry(0);
//  	}
//  
//      if (inputfilereader->GetEntryStatus() == kEntryValid) {
//         G4cout << "Loaded entry " << inputfilereader->GetCurrentEntry() << '\n';
//      } else {
//         switch (inputfilereader->GetEntryStatus()) {
//         kEntryValid:
//            // Handled above.
//            break;
//         kEntryNotLoaded:
//            G4cerr << "Error: TTreeReader has not loaded any data yet!\n";
//            break;
//         kEntryNoTree:
//            G4cerr << "Error: TTreeReader cannot find a tree names \"MyTree\"!\n";
//            break;
//         kEntryNotFound:
//            // Can't really happen as TTreeReader::Next() knows when to stop.
//            G4cerr << "Error: The entry number doe not exist\n";
//            break;
//         kEntryChainSetupError:
//            G4cerr << "Error: TTreeReader cannot access a chain element, e.g. file without the tree\n";
//            break;
//         kEntryChainFileError:
//            G4cerr << "Error: TTreeReader cannot open a chain element, e.g. missing file\n";
//            break;
//         kEntryDictionaryError:
//            G4cerr << "Error: TTreeReader cannot find the dictionary for some data\n";
//            break;
//         }
//         G4Exception("WCLitePrimaryGeneratorAction::GeneratePrimaries","Expl01",FatalException,"Error reading input file");
//         return;
//      }
//  	G4cout<<"getting entry"<<G4endl;
//  	for(int i=0;i<(**rdrntankp);i++){
//  	// generate each of the primaries from this event
//  		G4cout<<"entry got"<<G4endl;
//  		
//  		G4ThreeVector rdrvtx = G4ThreeVector((*rdrvxp)[i], (*rdrvyp)[i], (*rdrvzp)[i]);
//  		G4ThreeVector rdrdir = G4ThreeVector((*rdrpxp)[i], (*rdrpyp)[i], (*rdrpzp)[i]);
//  		rdrdir.unit();		// normalise to unit vector
//  		G4cout<<"setting particle gun"<<G4endl;
//  		particleGun->SetParticleDefinition( particleTable->FindParticle( (*rdrpdgtankp)[i] ) );	//FindParticle finds by either name or pdg.
//  		particleGun->SetParticleEnergy((*rdrkEp)[i]);
//  		particleGun->SetParticlePosition(rdrvtx);
//  		particleGun->SetParticleTime((*rdrvtp)[i]);
//  		particleGun->SetParticleMomentumDirection(rdrdir);
//  		particleGun->GeneratePrimaryVertex(anEvent);	//anEvent is provided by geant4 when invoking the method
//  		
//  	}
		
		Long64_t localEntry = inputdata->LoadTree(inputEntry);
		if(localEntry<0){
			G4cout<<"@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#"<<G4endl;
			G4cout<<"@#@#@#@#@#@#@#@#@#@ REACHED END OF INPUT FILE! #@#@#@#@#@#@#@#@#@#@#@#"<<G4endl;
			G4cout<<"@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#"<<G4endl;
			inputEntry=0;
			localEntry = inputdata->LoadTree(inputEntry);
			// should be done below since this opens a new tree
//			inputdata->SetBranchAddress("run",&runbranchval,&runBranch);
//			inputdata->SetBranchAddress("ntank",&ntankbranchval,&nTankBranch);
//			inputdata->SetBranchAddress("nupdg",&nupdgval,&nupdgBranch);
//			inputdata->SetBranchAddress("nuvtxx",&nuvtxxval,&nuvtxxBranch);
//			inputdata->SetBranchAddress("nuvtxy",&nuvtxyval,&nuvtxyBranch);
//			inputdata->SetBranchAddress("nuvtxz",&nuvtxzval,&nuvtxzBranch);
//			inputdata->SetBranchAddress("nuvtxt",&nuvtxtval,&nuvtxtBranch);
//			inputdata->SetBranchAddress("vtxvol",&nupvval,&nuPVBranch);
//			inputdata->SetBranchAddress("vtxmat",&numatval,&nuvtxmatBranch);
//			vtxxBranch=inputdata->GetBranch("vx");
//			vtxyBranch=inputdata->GetBranch("vy");
//			vtxzBranch=inputdata->GetBranch("vz");
//			vtxtBranch=inputdata->GetBranch("vt");
//			pxBranch=inputdata->GetBranch("px");
//			pyBranch=inputdata->GetBranch("py");
//			pzBranch=inputdata->GetBranch("pz");
//			EBranch=inputdata->GetBranch("E");
//			KEBranch=inputdata->GetBranch("kE");
//			pdgBranch=inputdata->GetBranch("pdgtank");
//			geniePrimaryBranch=inputdata->GetBranch("primary");
//			entriesInThisTree = runBranch->GetEntries();
		}
		
		Int_t nextTreeNumber = inputdata->GetTreeNumber();
		if(treeNumber!=nextTreeNumber){
			G4cout<< "Reached end of Tree. Last entries' tree number was "
						<< treeNumber <<", this entries' tree number is "<< nextTreeNumber <<G4endl;
			G4cout<<"Getting new tree branches"<<G4endl;
			// maybe all this isn't necessary??
			inputdata->SetBranchAddress("run",&runbranchval,&runBranch);
			inputdata->SetBranchAddress("ntank",&ntankbranchval,&nTankBranch);
			inputdata->SetBranchAddress("nupdg",&nupdgval,&nupdgBranch);
			inputdata->SetBranchAddress("nuvtxx",&nuvtxxval,&nuvtxxBranch);
			inputdata->SetBranchAddress("nuvtxy",&nuvtxyval,&nuvtxyBranch);
			inputdata->SetBranchAddress("nuvtxz",&nuvtxzval,&nuvtxzBranch);
			inputdata->SetBranchAddress("nuvtxt",&nuvtxtval,&nuvtxtBranch);
			inputdata->SetBranchAddress("vtxvol",&nupvval,&nuPVBranch);
			inputdata->SetBranchAddress("vtxmat",&numatval,&nuvtxmatBranch);
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
			geniePrimaryBranch=inputdata->GetBranch("primary");
			entriesInThisTree = runBranch->GetEntries();
		
			if(runBranch==0||nTankBranch==0||vtxxBranch==0||vtxyBranch==0||vtxzBranch==0||vtxtBranch==0||pxBranch==0||pyBranch==0||pzBranch==0||EBranch==0||KEBranch==0||pdgBranch==0||nupdgBranch==0||nuvtxxBranch==0||nuvtxyBranch==0||nuvtxzBranch==0||nuvtxtBranch==0||nuPVBranch==0||nuvtxmatBranch==0||geniePrimaryBranch==0){
				G4cout<<"BRANCHES ARE ZOMBIES ARGH!"<<G4endl;
			} else { G4cout<<"entries in this tree: "<<vtxxBranch->GetEntries()<<G4endl; }
			
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
		
		G4ParticleDefinition* parttype = particleTable->FindParticle(nupdgval);
		TLorentzVector neutrinovertex(nuvtxtval, nuvtxxval, nuvtxyval, nuvtxzval);
		G4cout<<"The origin interaction was a "<<(parttype->GetParticleName())<<" at ("<<nuvtxtval<<","<<nuvtxxval<<","<<nuvtxyval<<","<<nuvtxzval<<") in "<<nupvval<<" "<<numatval<<G4endl;
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
		if(genieprimarybranchval){delete[] genieprimarybranchval;}
		
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
		genieprimarybranchval = new Int_t[ntankbranchval];
		
		if(vtxxbranchval==0||vtxybranchval==0||vtxzbranchval==0||vtxtbranchval==0||pxbranchval==0||pybranchval==0||pzbranchval==0||ebranchval==0||kebranchval==0||pdgbranchval==0||genieprimarybranchval==0){
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
		geniePrimaryBranch->SetAddress(genieprimarybranchval);
		
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
		geniePrimaryBranch->GetEntry(localEntry);
		
		//G4cout<<"Looping over primaries"<<G4endl;
		for(int i=0;i<ntankbranchval;i++){
			//G4cout<<"Loading details of primary "<<i<<G4endl;
			vtxxval=vtxxbranchval[i]*cm;
			vtxyval=vtxybranchval[i]*cm;
			vtxzval=vtxzbranchval[i]*cm;
			vtxtval=vtxtbranchval[i]*ns;
			pxval=pxbranchval[i]*GeV;
			pyval=pybranchval[i]*GeV;
			pzval=pzbranchval[i]*GeV;
			eval=ebranchval[i]*GeV;
			keval=kebranchval[i]*GeV;
			pdgval=pdgbranchval[i];
			genieprimaryval=genieprimarybranchval[i];
			G4ThreeVector thevtx = G4ThreeVector(vtxxval, vtxyval, vtxzval);
			G4ThreeVector thepdir = G4ThreeVector(pxval, pyval, pzval);
			thepdir.unit();		// normalise to unit vector
			
			if(keval==0){keval+=(pow(1.,-6.))*eV; eval+=(pow(1.,-6.))*eV; thepdir=G4ThreeVector(0.,0.,1.);}
			parttype = particleTable->FindParticle(pdgval);
			if(genieprimaryval==1){G4cout<<"genie primary ";} else {G4cout<<"genie secondary ";}
			G4cout<<eval/GeV<<" GeV ";
			if(parttype==0){G4cout<<"PDG: "<<pdgval;} else {G4cout<<parttype->GetParticleName();}
			G4Navigator* theNavigator = G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking();
			G4VPhysicalVolume* primaryPV = theNavigator->LocateGlobalPointAndSetup(thevtx); 
			G4cout<<" in "<< primaryPV->GetName()<<G4endl;
			
			particleGun->SetParticleDefinition( particleTable->FindParticle(pdgval) );	//FindParticle finds by either name or pdg.
			particleGun->SetParticleEnergy(eval*MeV);
			particleGun->SetParticlePosition(thevtx);
			particleGun->SetParticleTime(vtxtval);
			particleGun->SetParticleMomentumDirection(thepdir);
			particleGun->GeneratePrimaryVertex(anEvent);	//anEvent is provided by geant4 when invoking the method
			//G4cout<<"Vertex set"<<G4endl;
		}
		
		inputEntry++;
		
	}
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

