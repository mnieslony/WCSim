#include "WCSimPrimaryGeneratorAction.hh"
#include "WCSimDetectorConstruction.hh"
#include "WCSimPrimaryGeneratorMessenger.hh"
#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4GeneralParticleSource.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4ThreeVector.hh"
#include "G4SPSAngDistribution.hh"
#include "G4SPSPosDistribution.hh"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TRandom3.h"
#include "globals.hh"
#include "Randomize.hh"
#include <fstream>
#include <vector>
#include <string>
#include <math.h> 
#include <libgen.h>

#include "G4Navigator.hh"
#include "G4TransportationManager.hh"
#include "G4UImanager.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

using std::vector;
using std::string;
using std::fstream;

// GENIE headers
#ifndef NO_GENIE
#include "GHEP/GHepParticle.h"
#include "GHEP/GHepUtils.h"
#include "Ntuple/NtpMCTreeHeader.h"
#include "Interaction/Interaction.h"
#endif

vector<string> tokenize( string separators, string input );

inline vector<string> readInLine(fstream& inFile, int lineSize, char* inBuf)
{
  // Read in line break it up into tokens
  // Any line starting with # is ignored
	while(true)                                               
	{  
		if (inFile.getline(inBuf,lineSize))
		{
			if(inBuf[0]!='#')                                  
				return tokenize(" $\r", inBuf);
		}
		else
		{
			vector<string> nullLine;                               
			return nullLine; 
		}
	}
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

inline double atof( const string& S ) {return std::atof( S.c_str() );}
inline int    atoi( const string& S ) {return std::atoi( S.c_str() );}

WCSimPrimaryGeneratorAction::WCSimPrimaryGeneratorAction(
					  WCSimDetectorConstruction* myDC)
  :myDetector(myDC), vectorFileName("")
{
  //T. Akiri: Initialize GPS to allow for the laser use 
  MyGPS = new G4GeneralParticleSource();

  // Initialize to zero
  mode[0] = 0;
  nvtxs = 0;
  for( Int_t u=0; u<MAX_N_VERTICES; u++){
    vtxsvol[u] = 0;
    vtxs[u] = G4ThreeVector(0.,0.,0.);
  }
  nuEnergy = 0.;
  _counterRock=0; // counter for generated in Rock
  _counterCublic=0; // counter generated
  
  //---Set defaults. Do once at beginning of session.
  
  G4int n_particle = 1;
  particleGun = new G4ParticleGun(n_particle);
  particleGun->SetParticleEnergy(1.0*GeV);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.0));

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4IonTable* ionTable = G4IonTable::GetIonTable();
  G4String particleName;
  particleGun->
    SetParticleDefinition(particleTable->FindParticle(particleName="mu+"));

  particleGun->
    SetParticlePosition(G4ThreeVector(0.*m,0.*m,0.*m));
   
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
  useMulineEvt 		= true;
  useGunEvt    		= false;
  useLaserEvt  		= false;
  useGPSEvt    		= false;
  useRadioactiveEvt  	= false;
  useRadonEvt        	= false;
  useAntiNuEvt		= false;
  useGenieEvt		= false; 
 
  // Radioactive and Radon generator variables:
  radioactive_sources.clear();
  myRn222Generator	= 0;
  fRnScenario		= 0;
  fRnSymmetry		= 1;
  // Time units for vertices   
  fTimeUnit=CLHEP::nanosecond;

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
  delete particleGun;
  delete MyGPS;   //T. Akiri: Delete the GPS variable
  delete messenger;

  delete theSPSPos;
  delete theSPSAng;

  if (useGenieEvt){
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
  G4IonTable* ionTable = G4IonTable::GetIonTable();
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
	// The original documentation describing the nuance text format can be found here: 
    // http://neutrino.phy.duke.edu/nuance-format/
    //
    // Information specific to WCSim can be found in the file Nuance_MC_Format.txt in
    // the doc directory.
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

    	if (token.size() == 0 || token[0] == "stop" ) 
    	{
    		G4cout << "End of nuance vector file - run terminated..."<< G4endl;
    		G4RunManager::GetRunManager()-> AbortRun();
    	}
    	else if (token[0] != "begin" )
    	{
    		G4cout << "unexpected line begins with " << token[0] << " we were expecting \" begin \" "<<G4endl;
    	}
		else   // normal parsing begins here
		{
		// Read the nuance line 
        // should be nuance <value>
        // but could be just  
        // nuance 
		// if value is given set mode to equal it.

			token = readInLine(inputFile, lineSize, inBuf);
			int iVertex=0;
			while(token[0]=="nuance" && iVertex < MAX_N_VERTICES)
			{
				if(token.size()>1)
					mode[iVertex] = atoi(token[1]);
	            // Read the Vertex line
				token = readInLine(inputFile, lineSize, inBuf);
				vtxs[iVertex] = G4ThreeVector(atof(token[1])*cm,
					atof(token[2])*cm,
					atof(token[3])*cm);
				G4double VertexTime=atof(token[4])*fTimeUnit; 
				vertexTimes[iVertex]=VertexTime;
                // true : Generate vertex in Rock , false : Generate vertex in WC tank
				SetGenerateVertexInRock(false);
	            // Next we read the incoming neutrino and target
	            // First, the neutrino line
				token=readInLine(inputFile, lineSize, inBuf);
				beampdgs[iVertex] = atoi(token[1]);
				beamenergies[iVertex] = atof(token[2])*MeV;
				beamdirs[iVertex] = G4ThreeVector(atof(token[3]),
					atof(token[4]),
					atof(token[5]));
				SetBeamEnergy(beamenergies[iVertex]);         
                SetBeamDir(beamdirs[iVertex]);
	            // Now read the target line
				token=readInLine(inputFile, lineSize, inBuf);
				targetpdgs[iVertex] = atoi(token[1]);
				targetenergies[iVertex] = atof(token[2])*MeV;
				targetdirs[iVertex] = G4ThreeVector(atof(token[3]),
					atof(token[4]),
					atof(token[5]));
	            // Read the info line, basically a dummy
				token=readInLine(inputFile, lineSize, inBuf);
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
						G4double tempEnergy = atof(token[2])*MeV;
						G4ThreeVector dir = G4ThreeVector(atof(token[3]),
							atof(token[4]),
							atof(token[5]));
		                //must handle the case of an ion seperatly from other particles
		                //check PDG code if we have an ion.
		                //PDG code format for ions ±10LZZZAAAI
						char strPDG[11];
						char strA[10]={0};
						char strZ[10]={0};
						long int A=0,Z=0;
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
							ion =  ionTable->GetIon(Z, A, 0.);
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

						G4double ekin = tempEnergy - mass;

						particleGun->SetParticleEnergy(ekin);
						particleGun->SetParticlePosition(vtxs[iVertex]);
						particleGun->SetParticleMomentumDirection(dir);
						particleGun->SetParticleTime(VertexTime);
						particleGun->GeneratePrimaryVertex(anEvent);

					}
				}
				iVertex++;
				if(iVertex > MAX_N_VERTICES)
					G4cout<<" CAN NOT DEAL WITH MORE THAN "<<MAX_N_VERTICES<<" VERTICES - TRUNCATING EVENT HERE "<<G4endl;
			}
			nvtxs=iVertex;
			SetNvtxs(nvtxs);

		}
	}
    else 
    {    // old muline format  
    	inputFile >> nuEnergy >> energy >> xPos >> yPos >> zPos 
    	>> xDir >> yDir >> zDir;

    	G4double random_z = ((myDetector->GetWaterTubePosition())
    		- .5*(myDetector->GetWaterTubeLength()) 
    		+ 1.*m + 15.0*m*G4UniformRand())/m;
    	zPos = random_z;
    	G4ThreeVector vtx = G4ThreeVector(xPos, yPos, random_z);
    	G4ThreeVector dir = G4ThreeVector(xDir,yDir,zDir);

    	particleGun->SetParticleEnergy(energy*MeV);
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
    G4double mass    =anEvent->GetPrimaryVertex()->GetPrimary()->GetMass();
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
    G4double E         = std::sqrt((P.dot(P))+(mass*mass));

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
      G4ThreeVector vtx =anEvent->GetPrimaryVertex()->GetPosition();
      G4int pdg         =anEvent->GetPrimaryVertex()->GetPrimary()->GetPDGcode();
      
      //     G4ThreeVector dir  = P.unit();
      G4double E         = std::sqrt((P.dot(P)));
      
      //SetVtx(vtx);
      SetBeamEnergy(E);
      //SetBeamDir(dir);
      SetBeamPDG(pdg);
    }
  else if (useGPSEvt)
    {
      MyGPS->GeneratePrimaryVertex(anEvent);
      
      G4ThreeVector P   =anEvent->GetPrimaryVertex()->GetPrimary()->GetMomentum();
      G4ThreeVector vtx =anEvent->GetPrimaryVertex()->GetPosition();
      G4double mass     =anEvent->GetPrimaryVertex()->GetPrimary()->GetMass();
      G4int pdg         =anEvent->GetPrimaryVertex()->GetPrimary()->GetPDGcode();
      
      G4ThreeVector dir  = P.unit();
      G4double E         = std::sqrt((P.dot(P))+(mass*mass));
      
      SetVtx(vtx);
      SetBeamEnergy(E);
      SetBeamDir(dir);
      SetBeamPDG(pdg);
    }
  else if (useRadioactiveEvt)
    {
      
      // initialize GPS properties
      MyGPS->ClearAll();
      
      MyGPS->SetMultipleVertex(true);
      
      std::vector<WCSimPmtInfo*> *pmts=NULL;
      
      std::vector<struct radioactive_source>::iterator it;
      
      for ( it = radioactive_sources.begin(); it != radioactive_sources.end(); it++ ){
      	G4String IsotopeName = it->IsotopeName;
      	G4String IsotopeLocation = it->IsotopeLocation;
      	G4double IsotopeActivity = it->IsotopeActivity;

      	double average= IsotopeActivity * GetRadioactiveTimeWindow();
      	if (IsotopeLocation.compareTo("PMT") == 0){
      		pmts = myDetector->Get_Pmts();
      		average *= pmts->size();
      	}
	  
	// random poisson number of vertices based on average
	int n_vertices = CLHEP::RandPoisson::shoot(average);

	//	n_vertices = 1; // qqq

	for(int u=0; u<n_vertices; u++){
	    
	  MyGPS->AddaSource(1.);
	    
	  MyGPS->SetCurrentSourceto(MyGPS->GetNumberofSource() - 1);

	  if (IsotopeName.compareTo("Tl208") == 0){
	    MyGPS->SetParticleDefinition(G4IonTable::GetIonTable()->GetIon( 81, 208, 0));
	  }else if (IsotopeName.compareTo("Bi214") == 0){
	    MyGPS->SetParticleDefinition(G4IonTable::GetIonTable()->GetIon( 83, 214, 0));
	  }else if (IsotopeName.compareTo("K40") == 0){
	    MyGPS->SetParticleDefinition(G4IonTable::GetIonTable()->GetIon( 19, 40, 0));
	  }else if (IsotopeName.compareTo("Rn220") == 0){
	    MyGPS->SetParticleDefinition(G4IonTable::GetIonTable()->GetIon( 86, 220, 0));
	  }else if (IsotopeName.compareTo("Po216") == 0){
	    MyGPS->SetParticleDefinition(G4IonTable::GetIonTable()->GetIon( 84, 216, 0));
	  }else if (IsotopeName.compareTo("Pb212") == 0){
	    MyGPS->SetParticleDefinition(G4IonTable::GetIonTable()->GetIon( 82, 212, 0));
	  }else if (IsotopeName.compareTo("Bi212") == 0){
	    MyGPS->SetParticleDefinition(G4IonTable::GetIonTable()->GetIon( 83, 212, 0));
	  }else if (IsotopeName.compareTo("Po212") == 0){
	    MyGPS->SetParticleDefinition(G4IonTable::GetIonTable()->GetIon( 84, 212, 0));
	  }else if (IsotopeName.compareTo("Rn222") == 0){
	    MyGPS->SetParticleDefinition(G4IonTable::GetIonTable()->GetIon( 86, 222, 0));
	  }else if (IsotopeName.compareTo("Po218") == 0){
	    MyGPS->SetParticleDefinition(G4IonTable::GetIonTable()->GetIon( 84, 218, 0));
	  }else if (IsotopeName.compareTo("At218") == 0){
	    MyGPS->SetParticleDefinition(G4IonTable::GetIonTable()->GetIon( 85, 218, 0));
	  }else if (IsotopeName.compareTo("Pb214") == 0){
	    MyGPS->SetParticleDefinition(G4IonTable::GetIonTable()->GetIon( 82, 214, 0));
	  }else if (IsotopeName.compareTo("Po214") == 0){
	    MyGPS->SetParticleDefinition(G4IonTable::GetIonTable()->GetIon( 84, 214, 0));
	  }else if (IsotopeName.compareTo("Tl210") == 0){
	    MyGPS->SetParticleDefinition(G4IonTable::GetIonTable()->GetIon( 81, 210, 0));
	  }else if (IsotopeName.compareTo("Pb210") == 0){
	    MyGPS->SetParticleDefinition(G4IonTable::GetIonTable()->GetIon( 82, 210, 0));
	  }else if (IsotopeName.compareTo("Bi210") == 0){
	    MyGPS->SetParticleDefinition(G4IonTable::GetIonTable()->GetIon( 83, 210, 0));
	  }else if (IsotopeName.compareTo("Po210") == 0){
	    MyGPS->SetParticleDefinition(G4IonTable::GetIonTable()->GetIon( 84, 210, 0));
	  }else if (IsotopeName.compareTo("Hg206") == 0){
	    MyGPS->SetParticleDefinition(G4IonTable::GetIonTable()->GetIon( 80, 206, 0));
	  }else if (IsotopeName.compareTo("Tl206") == 0){
	    MyGPS->SetParticleDefinition(G4IonTable::GetIonTable()->GetIon( 81, 206, 0));
	  }else if (IsotopeName.compareTo("Rn219") == 0){
	    MyGPS->SetParticleDefinition(G4IonTable::GetIonTable()->GetIon( 86, 219, 0));
	  }else if (IsotopeName.compareTo("Po215") == 0){
	    MyGPS->SetParticleDefinition(G4IonTable::GetIonTable()->GetIon( 84, 215, 0));
	  }else if (IsotopeName.compareTo("At215") == 0){
	    MyGPS->SetParticleDefinition(G4IonTable::GetIonTable()->GetIon( 85, 215, 0));
	  }else if (IsotopeName.compareTo("Pb211") == 0){
	    MyGPS->SetParticleDefinition(G4IonTable::GetIonTable()->GetIon( 82, 211, 0));
	  }else if (IsotopeName.compareTo("Bi211") == 0){
	    MyGPS->SetParticleDefinition(G4IonTable::GetIonTable()->GetIon( 83, 211, 0));
	  }else if (IsotopeName.compareTo("Po211") == 0){
	    MyGPS->SetParticleDefinition(G4IonTable::GetIonTable()->GetIon( 84, 211, 0));
	  }else if (IsotopeName.compareTo("Tl207") == 0){
	    MyGPS->SetParticleDefinition(G4IonTable::GetIonTable()->GetIon( 81, 207, 0));
	  }else if (IsotopeName.compareTo("Th232") == 0){
	    MyGPS->SetParticleDefinition(G4IonTable::GetIonTable()->GetIon( 90, 232, 0));
	  }else if (IsotopeName.compareTo("Ra228") == 0){
	    MyGPS->SetParticleDefinition(G4IonTable::GetIonTable()->GetIon( 88, 228, 0));
	  }else if (IsotopeName.compareTo("Ac228") == 0){
	    MyGPS->SetParticleDefinition(G4IonTable::GetIonTable()->GetIon( 89, 228, 0));
	  }else if (IsotopeName.compareTo("Th228") == 0){
	    MyGPS->SetParticleDefinition(G4IonTable::GetIonTable()->GetIon( 90, 228, 0));
	  }else if (IsotopeName.compareTo("Ra224") == 0){
	    MyGPS->SetParticleDefinition(G4IonTable::GetIonTable()->GetIon( 88, 224, 0));
	  }else if (IsotopeName.compareTo("U238") == 0){
	    MyGPS->SetParticleDefinition(G4IonTable::GetIonTable()->GetIon( 92, 238, 0));
	  }else if (IsotopeName.compareTo("Th234") == 0){
	    MyGPS->SetParticleDefinition(G4IonTable::GetIonTable()->GetIon( 90, 234, 0));
	  }else if (IsotopeName.compareTo("Pa234") == 0){
	    MyGPS->SetParticleDefinition(G4IonTable::GetIonTable()->GetIon( 91, 234, 0));
	  }else if (IsotopeName.compareTo("U234") == 0){
	    MyGPS->SetParticleDefinition(G4IonTable::GetIonTable()->GetIon( 92, 234, 0));
	  }else if (IsotopeName.compareTo("Th230") == 0){
	    MyGPS->SetParticleDefinition(G4IonTable::GetIonTable()->GetIon( 90, 230, 0));
	  }else if (IsotopeName.compareTo("Ra226") == 0){
	    MyGPS->SetParticleDefinition(G4IonTable::GetIonTable()->GetIon( 88, 226, 0));
	  }else if (IsotopeName.compareTo("U235") == 0){
	    MyGPS->SetParticleDefinition(G4IonTable::GetIonTable()->GetIon( 92, 235, 0));
	  }else if (IsotopeName.compareTo("Th231") == 0){
	    MyGPS->SetParticleDefinition(G4IonTable::GetIonTable()->GetIon( 90, 231, 0));
	  }else if (IsotopeName.compareTo("Pa231") == 0){
	    MyGPS->SetParticleDefinition(G4IonTable::GetIonTable()->GetIon( 91, 231, 0));
	  }else if (IsotopeName.compareTo("Ac227") == 0){
	    MyGPS->SetParticleDefinition(G4IonTable::GetIonTable()->GetIon( 89, 227, 0));
	  }else if (IsotopeName.compareTo("Th227") == 0){
	    MyGPS->SetParticleDefinition(G4IonTable::GetIonTable()->GetIon( 90, 227, 0));
	  }else if (IsotopeName.compareTo("Fr223") == 0){
	    MyGPS->SetParticleDefinition(G4IonTable::GetIonTable()->GetIon( 87, 223, 0));
	  }else if (IsotopeName.compareTo("Ra223") == 0){
	    MyGPS->SetParticleDefinition(G4IonTable::GetIonTable()->GetIon( 88, 223, 0));
	  }else if (IsotopeName.compareTo("At219") == 0){
	    MyGPS->SetParticleDefinition(G4IonTable::GetIonTable()->GetIon( 85, 219, 0));
	  }else if (IsotopeName.compareTo("Bi215") == 0){
	    MyGPS->SetParticleDefinition(G4IonTable::GetIonTable()->GetIon( 83, 215, 0));
	  }

	  if (IsotopeLocation.compareTo("water") == 0){
	    MyGPS->GetCurrentSource()->GetEneDist()->SetEnergyDisType("Mono");
	    MyGPS->GetCurrentSource()->GetEneDist()->SetMonoEnergy(0.);
	    MyGPS->GetCurrentSource()->GetPosDist()->SetPosDisType("Point");
	    MyGPS->GetCurrentSource()->GetPosDist()->SetCentreCoords(G4ThreeVector(0, 0, 0));
	    MyGPS->GetCurrentSource()->GetPosDist()->SetPosDisType("Volume");
	    MyGPS->GetCurrentSource()->GetPosDist()->SetPosDisShape("Cylinder");
	    G4String WCIDCollectionName = myDetector->GetIDCollectionName();
	    WCSimPMTObject *PMT = myDetector->GetPMTPointer(WCIDCollectionName);
	    MyGPS->GetCurrentSource()->GetPosDist()->SetRadius(myDetector->GetGeo_Dm(3)*CLHEP::cm - 2.*PMT->GetRadius());
	    MyGPS->GetCurrentSource()->GetPosDist()->SetHalfZ(myDetector->GetGeo_Dm(2)*CLHEP::cm/2. - 2.*PMT->GetRadius());
	    MyGPS->GetCurrentSource()->GetPosDist()->SetPosRot1(G4ThreeVector(1, 0, 0));
	    MyGPS->GetCurrentSource()->GetPosDist()->SetPosRot2(G4ThreeVector(0, 1, 0));

	  }
	  else if (IsotopeLocation.compareTo("PMT") == 0){
	    int npmts = pmts->size();
	    int random_pmt_id = CLHEP::RandFlat::shootInt(1,npmts);
	    WCSimPmtInfo* pmtinfo = (WCSimPmtInfo*)pmts->at( random_pmt_id - 1 );
	    G4ThreeVector random_pmt_center(pmtinfo->Get_transx()*CLHEP::cm, pmtinfo->Get_transy()*CLHEP::cm, pmtinfo->Get_transz()*CLHEP::cm);
	    double random_cos_theta = CLHEP::RandFlat::shoot(0., 1.);
	    double random_sin_theta = sqrt(1. - pow(random_cos_theta,2));
	    random_sin_theta *= (CLHEP::RandFlat::shootBit() == 0 ? -1 : 1);
	    double random_phi = CLHEP::RandFlat::shoot(0., 2.*CLHEP::pi*CLHEP::rad);
	    G4String WCIDCollectionName = myDetector->GetIDCollectionName();
	    WCSimPMTObject *PMT = myDetector->GetPMTPointer(WCIDCollectionName);
	    double PMT_radius = PMT->GetRadius();
	    double glassThickness = PMT->GetPMTGlassThickness();
	    double expose = PMT->GetExposeHeight();
	    double sphereRadius = (expose*expose+ PMT_radius*PMT_radius)/(2*expose);
	    double Rmin = sphereRadius-glassThickness;
	    double Rmax = sphereRadius;
	    double random_R = CLHEP::RandFlat::shoot(Rmin, Rmax);
	    G4ThreeVector orientation(pmtinfo->Get_orienx(), pmtinfo->Get_orieny(), pmtinfo->Get_orienz());
	    G4ThreeVector axis_1 = orientation.orthogonal();
	    G4ThreeVector axis_2 = orientation.cross(axis_1);
	    G4ThreeVector position = random_pmt_center + random_R*(orientation*random_cos_theta + axis_1*random_sin_theta*cos(random_phi) + axis_2*random_sin_theta*sin(random_phi));
	      
	    //G4cout << " random id " << random_pmt_id << " of " << npmts << " costheta " << random_cos_theta << " sintheta " << random_sin_theta << " phi " << random_phi << " WCIDCollectionName " << WCIDCollectionName << " PMT_radius " << PMT_radius << " expose " << expose << " sphereRadius " << sphereRadius << " Rmin " << Rmin << " Rmax " << Rmax << " random_R " << random_R << " orientation (" << orientation.x() << ", " << orientation.y() << ", " << orientation.z() << ") center (" << random_pmt_center.x() << ", " << random_pmt_center.y() << ", " << random_pmt_center.z() << ") position (" << position.x() << ", " << position.y() << ", " << position.z() << ") " << G4endl;
	      
	    MyGPS->GetCurrentSource()->GetEneDist()->SetEnergyDisType("Mono");
	    MyGPS->GetCurrentSource()->GetEneDist()->SetMonoEnergy(0.);
	    MyGPS->GetCurrentSource()->GetPosDist()->SetPosDisType("Point");
	    MyGPS->GetCurrentSource()->GetPosDist()->SetCentreCoords(position);
	  }
	    
	}

//	G4cout << " is " << IsotopeName << " of " << radioactive_sources.size() << " loc " << IsotopeLocation << " a " << IsotopeActivity << " nv " << n_vertices << G4endl;

      }

      G4int number_of_sources = MyGPS->GetNumberofSource();

      // this will generate several primary vertices
      MyGPS->GeneratePrimaryVertex(anEvent);

      SetNvtxs(number_of_sources);
      for( G4int u=0; u<number_of_sources; u++){
	targetpdgs[u] = 2212; //ie. proton 

      	G4ThreeVector P   =anEvent->GetPrimaryVertex(u)->GetPrimary()->GetMomentum();
      	G4ThreeVector vtx =anEvent->GetPrimaryVertex(u)->GetPosition();
      	G4int pdg         =anEvent->GetPrimaryVertex(u)->GetPrimary()->GetPDGcode();
      
      	//       G4ThreeVector dir  = P.unit();
      	G4double E         = std::sqrt((P.dot(P)));

//	G4cout << " vertex " << u << " of " << number_of_sources << " (" << vtx.x() << ", " << vtx.y() << ", " << vtx.z() << ")" << G4endl;

      	SetVtxs(u,vtx);
      	SetBeamEnergy(E,u);
      	//       SetBeamDir(dir);
      	SetBeamPDG(pdg,u);
      }

    }
  else if (useRadonEvt)
    { //G. Pronost: Add Radon (adaptation of Radioactive event)
    
      // Currently only one generator is possible
      // In order to have several, we need to find a solution for the fitting graphes (which are static currently)
      // Idea: array of fitting graphes? (each new generators having a specific ID)
      if ( !myRn222Generator ) {
      	myRn222Generator = new WCSimGenerator_Radioactivity(myDetector);
      	myRn222Generator->Configuration(fRnScenario);
      }
      
      //G4cout << " Generate radon events " << G4endl;
      // initialize GPS properties
      MyGPS->ClearAll();
      
      MyGPS->SetMultipleVertex(true);
      
      
      std::vector<struct radioactive_source>::iterator it;
      
      G4String IsotopeName = "Rn222";
      G4double IsotopeActivity = myRn222Generator->GetMeanActivity() * 1e-3; // mBq to Bq
      G4double iEventAvg = IsotopeActivity * GetRadioactiveTimeWindow();

      //G4cout << " Average " << iEventAvg << G4endl;
      // random poisson number of vertices based on average
      int n_vertices = CLHEP::RandPoisson::shoot(iEventAvg);

      if ( n_vertices < 1 ) {
      	 n_vertices = 1;
      }
      
      for(int u=0; u<n_vertices; u++){
	
	MyGPS->AddaSource(1.);	
	MyGPS->SetCurrentSourceto(MyGPS->GetNumberofSource() - 1);
	
	// Bi214 (source of electron in Rn222 decay chain, assumed to be in equilibrium)
	MyGPS->SetParticleDefinition(G4IonTable::GetIonTable()->GetIon( 83, 214, 0));
	
	// Get position (first position take few seconds to be produced, there after there is no trouble)
	//G4cout << "GetRandomVertex" << G4endl;
	G4ThreeVector position = myRn222Generator->GetRandomVertex(fRnSymmetry);
	//G4cout << "Done: " << position << G4endl;
	// energy 
	MyGPS->GetCurrentSource()->GetEneDist()->SetEnergyDisType("Mono");
	MyGPS->GetCurrentSource()->GetEneDist()->SetMonoEnergy(0.);
	    
	// position 
	MyGPS->GetCurrentSource()->GetPosDist()->SetPosDisType("Point");
	MyGPS->GetCurrentSource()->GetPosDist()->SetCentreCoords(position);

	//G4cout << u << " is " << IsotopeName << " loc " << position  << G4endl;

      }
      G4int number_of_sources = MyGPS->GetNumberofSource();

      // this will generate several primary vertices
      MyGPS->GeneratePrimaryVertex(anEvent);

      SetNvtxs(number_of_sources);
      for( G4int u=0; u<number_of_sources; u++){
	targetpdgs[u] = 2212; //ie. proton 

      	G4ThreeVector P   =anEvent->GetPrimaryVertex(u)->GetPrimary()->GetMomentum();
      	G4ThreeVector vtx =anEvent->GetPrimaryVertex(u)->GetPosition();
      	G4int pdg         =anEvent->GetPrimaryVertex(u)->GetPrimary()->GetPDGcode();
      
      	//       G4ThreeVector dir  = P.unit();
      	G4double E         = std::sqrt((P.dot(P)));
	
	//G4cout << " vertex " << u << " of " << number_of_sources << " (" << vtx.x() << ", " << vtx.y() << ", " << vtx.z() << ") with pdg: " << pdg << G4endl;

      	SetVtxs(u,vtx);
      	SetBeamEnergy(E,u);
      	//       SetBeamDir(dir);
      	SetBeamPDG(pdg,u);
      }

    } else if (useGenieEvt){

      G4cout <<"Genie + TALYS event generator: Initialising..."<<G4endl;

      if (localEntry!=0) loadNewGenie = false;
      else {
        LoadTalysFiles();    //Load Talys root files once at the beginning!
      }
      G4cout <<"loadNewGenie: "<<loadNewGenie<<G4endl; 

      // Load the next entry, with all required trees and files
      //-----------------------------------------------------

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
	//neutrino
	geniedata->SetBranchAddress("neu",&genie_neutrinoflav,&genie_neutrinoflav_branch);
        geniedata->SetBranchAddress("Ev",&genie_neutrinoE,&genie_neutrinoE_branch);
        geniedata->SetBranchAddress("pxv",&genie_neutrinopx,&genie_neutrinopx_branch);
        geniedata->SetBranchAddress("pyv",&genie_neutrinopy,&genie_neutrinopy_branch);
        geniedata->SetBranchAddress("pzv",&genie_neutrinopz,&genie_neutrinopz_branch);

	//target
	geniedata->SetBranchAddress("cc",&genie_cc,&genie_cc_branch);
        geniedata->SetBranchAddress("nc",&genie_nc,&genie_nc_branch);
        geniedata->SetBranchAddress("Z",&genie_Z,&genie_Z_branch);
        geniedata->SetBranchAddress("A",&genie_A,&genie_A_branch);
        geniedata->SetBranchAddress("hitnuc",&genie_hitnuc,&genie_hitnuc_branch);
        geniedata->SetBranchAddress("qel",&genie_qel,&genie_qel_branch);
        geniedata->SetBranchAddress("res",&genie_res,&genie_res_branch);
        geniedata->SetBranchAddress("dis",&genie_dis,&genie_dis_branch);
        geniedata->SetBranchAddress("coh",&genie_coh,&genie_coh_branch);
        geniedata->SetBranchAddress("imd",&genie_imd,&genie_imd_branch);

	//lepton
	geniedata->SetBranchAddress("fspl",&genie_fspl,&genie_fspl_branch);
	geniedata->SetBranchAddress("El",&genie_El,&genie_El_branch);
	geniedata->SetBranchAddress("pxl",&genie_pxl,&genie_pxl_branch);
	geniedata->SetBranchAddress("pyl",&genie_pyl,&genie_pyl_branch);
	geniedata->SetBranchAddress("pzl",&genie_pzl,&genie_pzl_branch);

	//Counts of final state particles
	geniedata->SetBranchAddress("nfn",&genie_final_n,&genie_final_n_branch);
        geniedata->SetBranchAddress("nfp",&genie_final_p,&genie_final_p_branch);
        geniedata->SetBranchAddress("nfpip",&genie_final_pip,&genie_final_pip_branch);
        geniedata->SetBranchAddress("nfpim",&genie_final_pim,&genie_final_pim_branch);
        geniedata->SetBranchAddress("nfpi0",&genie_final_pi0,&genie_final_pi0_branch);
        geniedata->SetBranchAddress("nfkp",&genie_final_kp,&genie_final_kp_branch);
        geniedata->SetBranchAddress("nfkm",&genie_final_km,&genie_final_km_branch);
        geniedata->SetBranchAddress("nfk0",&genie_final_k0,&genie_final_k0_branch);

	//Detailed information about final state particles
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

        if(genie_entry_branch==0 || genie_neutrinoflav_branch==0 || genie_neutrinoE_branch==0 || genie_neutrinopx_branch==0 || genie_neutrinopy_branch==0 || genie_neutrinopz_branch==0 || genie_cc_branch==0 || genie_nc_branch==0 || genie_Z_branch==0 || genie_A_branch==0 || genie_hitnuc_branch == 0 || genie_qel_branch == 0 || genie_res_branch == 0 || genie_dis_branch == 0 || genie_coh_branch == 0 || genie_imd_branch == 0 || genie_fspl_branch == 0 || genie_El_branch == 0 || genie_pxl_branch == 0 || genie_pyl_branch == 0 || genie_pzl_branch == 0 || genie_final_n_branch==0 || genie_final_p_branch==0 || genie_final_pip_branch==0 || genie_final_pim_branch==0 || genie_final_pi0_branch==0 || genie_final_kp_branch==0 || genie_final_km_branch==0 || genie_final_k0_branch==0 || genie_pdg_final_branch==0 || genie_E_final_branch==0 || genie_px_final_branch==0 || genie_py_final_branch==0 || genie_pz_final_branch==0 || genie_vtxx_branch==0 || genie_vtxy_branch==0 || genie_vtxz_branch==0){
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
      genie_hitnuc_branch->GetEntry(localEntry);
      genie_qel_branch->GetEntry(localEntry);
      genie_res_branch->GetEntry(localEntry);
      genie_dis_branch->GetEntry(localEntry);
      genie_coh_branch->GetEntry(localEntry);
      genie_imd_branch->GetEntry(localEntry);
      genie_fspl_branch->GetEntry(localEntry);
      genie_El_branch->GetEntry(localEntry);
      genie_pxl_branch->GetEntry(localEntry);
      genie_pyl_branch->GetEntry(localEntry);
      genie_pzl_branch->GetEntry(localEntry);
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

      //genie_vtxy -= 0.1441;	//just needed in ANNIE because of coordinate system shift
      //genie_vtxz += 1.681;

      //Generate random position within specified volume (easy scale up from ANNIE to Super-K) 
      G4cout <<"Position before randomizing of position: ("<<genie_vtxx<<","<<genie_vtxy<<","<<genie_vtxz<<")"<<G4endl;
      theSPSAng->SetAngDistType("iso");
      theSPSAng->SetPosDistribution(theSPSPos);
      thePosition = theSPSPos->GenerateOne();
      genie_vtxx = thePosition.x()/1000.;
      genie_vtxy = thePosition.y()/1000.;
      genie_vtxz = thePosition.z()/1000.;
      G4cout <<"Position after randomizing: ("<<genie_vtxx<<","<<genie_vtxy<<","<<genie_vtxz<<")"<<G4endl;

      G4cout <<"genie_entry: "<<genie_entry<<", genie_neutrinoflav: "<<genie_neutrinoflav<<", genie_neutrinoE: "<<genie_neutrinoE<<", genie_neutrinopx: "<<genie_neutrinopx<<", genie_neutrinopy: "<<genie_neutrinopy<<", genie_neutrinopz: "<<genie_neutrinopz<<", genie_cc: "<<genie_cc<<", genie_nc: "<<genie_nc<<", genie_Z: "<<genie_Z<<", genie_A: "<<genie_A<<", genie_hitnuc: "<<genie_hitnuc<<", genie_final_n: "<<genie_final_n<<", genie_final_p: "<<genie_final_p<<G4endl;

      //Define neutrino particle properties
      G4ParticleDefinition* parttype = particleTable->FindParticle(genie_neutrinoflav);
      TLorentzVector neutrinovertex(genie_vtxt*CLHEP::s, genie_vtxx*CLHEP::m, genie_vtxy*CLHEP::m, genie_vtxz*CLHEP::m);  // position in m, times in s, convert to cm and ns
      G4cout<<"The origin interaction was a "<<(parttype->GetParticleName())<<" at ("<<genie_vtxt*1000000000.<<","<<genie_vtxx*100.<<","<<genie_vtxy*100.<<","<<genie_vtxz*100.<<")[ns, cm] in (Z,A) = ("<<genie_Z<<","<<genie_A<<")"<<G4endl;
      G4cout<<"This entry has "<<genie_num_final<<" primaries"<<G4endl;


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

      //Get TALYS de-excitation information

      G4cout <<"Clearing vectors"<<G4endl;
      talys_pdg.clear();	//additional de-excitation final state pdgs
      talys_momdir.clear();     //additional de-excitation final state momenta
      talys_energy.clear();     //de-excitation final state energies

      G4cout <<"Check de excitation prob"<<G4endl;
      bool do_deexcitation=false;
      double random_deex = G4UniformRand();
      G4cout <<"random_deex = "<<random_deex<<", ExcitationProb = "<<ExcitationProb<<G4endl;
      if (random_deex < ExcitationProb){
        do_deexcitation=true;
      }

      G4cout <<"do_deexcitation: "<<do_deexcitation<<G4endl;
      if (do_deexcitation){
        int resNuclA=0, resNuclZ=0;
        double charge=0;
        for (int i_final = 0; i_final < genie_num_final; i_final++){
          charge+=particleTable->FindParticle(pdgbranchval[i_final])->GetPDGCharge();
        }
        resNuclZ=genie_Z-charge;
        resNuclA=genie_A-genie_final_n-genie_final_p;
        G4cout <<"resNuclZ = "<<resNuclZ<<", resNuclA = "<<resNuclA<<G4endl;

        if (genie_hitnuc == 2112 || genie_hitnuc == 2212) //neutron or proton was knocked out
        {
          if (genie_Z == 8 && genie_A == 16){
            G4cout <<"Oxygen nucleus detected! Read in talys-information."<<G4endl;
            G4cout <<"Number of final n: "<<genie_final_n<<", number of final p: "<<genie_final_p<<G4endl;
            std::string res_nucleus_name = CalculateResNucleus(resNuclA,resNuclZ);
	    G4cout <<"TALYS file name: "<<res_nucleus_name<<G4endl;
            if (res_nucleus_name != "none"){
		LoadDeexcitationProb();
                GenerateDeexcitation(&talys_pdg,&talys_momdir,&talys_energy,resNuclA,resNuclZ);
            }
            else {
              G4cout <<"TALYS de-excitation information could not be loaded for residual nucleus (A,Z) = ("<<resNuclA<<","<<resNuclZ<<"). Use no de-excitation in this event"<<G4endl;
            }
          } else {
            G4cout <<"No oxygen nucleus, don't read in TALYS information. (Z = "<<genie_Z<<", A = "<<genie_A<<")"<<G4endl;
          }
        }
        G4cout <<"Summary of TALYS-generated de-excitation particles:"<<G4endl;
        G4cout <<"Number of particles: "<<talys_pdg.size()<<G4endl;
        for (unsigned int i_par =0; i_par < talys_pdg.size(); i_par ++){
          G4cout <<"PDG: "<<talys_pdg.at(i_par)<<", Energy: "<<talys_energy.at(i_par)<<", Momentum direction: ("<<talys_momdir.at(i_par).getX()<<","<<talys_momdir.at(i_par).getY()<<","<<talys_momdir.at(i_par).getZ()<<")"<<G4endl;
        }
      }

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
        
          //ION case  
          G4cout <<"Ion case!"<<G4endl;
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
        //G4VPhysicalVolume* primaryPV = theNavigator->LocateGlobalPointAndSetup(thevtx);
        //G4cout<<" in "<< primaryPV->GetName()<<G4endl;

        if(keval==0){keval+=1.*eV; eval+=1.*eV; thepdir=G4ThreeVector(0.,0.,1.);}
        
	G4cout <<"keval after correction: "<<keval<<G4endl;
	//Set properties of the primary particle
	particleGun->SetParticleEnergy(keval);
        particleGun->SetParticlePosition(thevtx);
        particleGun->SetParticleMomentumDirection(thepdir);

        //Fire the particle away
        particleGun->GeneratePrimaryVertex(anEvent); 
      } 

      //Simulate primary lepton
      int pdg_lepton = genie_fspl;
      G4double kinE_lepton = genie_El*GeV;
      genie_vtxx_temp=genie_vtxx*100*CLHEP::cm;
      genie_vtxy_temp=genie_vtxy*100*CLHEP::cm;
      genie_vtxz_temp=genie_vtxz*100*CLHEP::cm;      
      G4ThreeVector lepton_vertex = G4ThreeVector(genie_vtxx_temp,genie_vtxy_temp,genie_vtxz_temp);     
      G4double px_lepton = genie_pxl*GeV;
      G4double py_lepton = genie_pyl*GeV;
      G4double pz_lepton = genie_pzl*GeV;
      G4ThreeVector lepton_dir = G4ThreeVector(px_lepton,py_lepton,pz_lepton);

      G4ParticleDefinition* primary_lepton = particleTable->FindParticle(pdg_lepton);
      particleGun->SetParticleDefinition(primary_lepton);
      particleGun->SetParticleEnergy(kinE_lepton);
      particleGun->SetParticlePosition(lepton_vertex);
      particleGun->SetParticleMomentumDirection(lepton_dir);

      G4cout <<"Add primary lepton particle with pdg "<<pdg_lepton<<", energy "<<kinE_lepton<<", vtx = ("<<genie_vtxx_temp<<","<<genie_vtxy_temp<<","<<genie_vtxz_temp<<"), dir = ("<<px_lepton<<","<<py_lepton<<","<<pz_lepton<<")"<<G4endl;

      //Fire the primary lepton away
      particleGun->GeneratePrimaryVertex(anEvent);

      //Simulate additional TALYS-deexcitation particles:	
     for (unsigned int i_deexc = 0; i_deexc < talys_pdg.size(); i_deexc++){

        G4double momentum,mass;
        genie_vtxx_temp=genie_vtxx*100*CLHEP::cm;
        genie_vtxy_temp=genie_vtxy*100*CLHEP::cm;
        genie_vtxz_temp=genie_vtxz*100*CLHEP::cm;
        G4double kinE = talys_energy[i_deexc]*MeV;

       /*if (talys_pdg.at(i_deexc)==2112) n_deexNeutron++;
           else if (talys_pdg.at(i_deexc)==2212) n_deexProton++;
           else if (talys_pdg.at(i_deexc)==1000010020) n_deexDeuterium++;
 	   else if (talys_pdg.at(i_deexc)==1000010030) n_deexTritium++;
           else if (talys_pdg.at(i_deexc)==1000020030) n_deexHelium++;
           else if (talys_pdg.at(i_deexc)==1000020040) n_deexAlpha++;
       */

        G4ParticleDefinition* deexc_particle = particleTable->FindParticle(talys_pdg.at(i_deexc));
        if (deexc_particle) //does particle type exist?

        {
          mass=deexc_particle->GetPDGMass();
          if (fabs(mass)<1E-9) momentum = kinE;
          else momentum = sqrt(2*mass*kinE);//use non-relativistic formula since the procued nuclei are pretty heavy and won't be at relativistic speeds
	  
          G4ThreeVector deexc_vtx = G4ThreeVector(genie_vtxx_temp, genie_vtxy_temp, genie_vtxz_temp);
          //Set properties of this de-excitation particle
          particleGun->SetParticleDefinition(deexc_particle);
          particleGun->SetParticleEnergy(kinE);
          particleGun->SetParticlePosition(deexc_vtx);
          particleGun->SetParticleMomentumDirection(talys_momdir[i_deexc]);

          G4cout <<"Add de-excitation particle with pdg "<<talys_pdg[i_deexc]<<", energy "<<talys_energy[i_deexc]<<", vtx = ("<<genie_vtxx_temp<<","<<genie_vtxy_temp<<","<<genie_vtxz_temp<<"), dir = ("<<talys_momdir[i_deexc].getX()<<","<<talys_momdir[i_deexc].getY()<<","<<talys_momdir[i_deexc].getZ()<<")"<<G4endl;

          //Fire the particle away
          particleGun->GeneratePrimaryVertex(anEvent);
        }
        else {
          G4cerr <<"Could not find de-excitation particle with pdg "<<talys_pdg[i_deexc]<<", energy "<<talys_energy[i_deexc]<<G4endl;
        }
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

    } else if (useAntiNuEvt){

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
        geniedata->SetBranchAddress("hitnuc",&genie_hitnuc,&genie_hitnuc_branch);
        geniedata->SetBranchAddress("qel",&genie_qel,&genie_qel_branch);
        geniedata->SetBranchAddress("res",&genie_res,&genie_res_branch);
	geniedata->SetBranchAddress("dis",&genie_dis,&genie_dis_branch);
	geniedata->SetBranchAddress("coh",&genie_coh,&genie_coh_branch);
	geniedata->SetBranchAddress("imd",&genie_imd,&genie_imd_branch);
        geniedata->SetBranchAddress("fspl",&genie_fspl,&genie_fspl_branch);
	geniedata->SetBranchAddress("El",&genie_El,&genie_El_branch);
	geniedata->SetBranchAddress("pxl",&genie_pxl,&genie_pxl_branch);
	geniedata->SetBranchAddress("pyl",&genie_pyl,&genie_pyl_branch);
	geniedata->SetBranchAddress("pzl",&genie_pzl,&genie_pzl_branch);
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
	if(genie_entry_branch==0 || genie_neutrinoflav_branch==0 || genie_neutrinoE_branch==0 || genie_neutrinopx_branch==0 || genie_neutrinopy_branch==0 || genie_neutrinopz_branch==0 || genie_cc_branch==0 || genie_nc_branch==0 || genie_Z_branch==0 || genie_A_branch==0 || genie_hitnuc_branch==0 || genie_qel_branch==0 || genie_res_branch==0 || genie_dis_branch==0 || genie_coh_branch==0 || genie_imd_branch==0 || genie_fspl_branch==0 || genie_El_branch==0 || genie_pxl_branch==0 || genie_pyl_branch==0 || genie_pzl_branch==0 || genie_final_n_branch==0 || genie_final_p_branch==0 || genie_final_pip_branch==0 || genie_final_pim_branch==0 || genie_final_pi0_branch==0 || genie_final_kp_branch==0 || genie_final_km_branch==0 || genie_final_k0_branch==0 || genie_pdg_final_branch==0 || genie_E_final_branch==0 || genie_px_final_branch==0 || genie_py_final_branch==0 || genie_pz_final_branch==0 || genie_vtxx_branch==0 || genie_vtxy_branch==0 || genie_vtxz_branch==0){
                G4cout<<"LoadNewGENIEFile: BRANCHES ARE ZOMBIES ARGH!"<<G4endl;
        }

        G4cout <<"Getting entries in tree..."<<G4endl;
	entriesInThisTree=geniedata->GetEntries();
        G4cout <<"Entries in this tree: "<<entriesInThisTree<<G4endl;
	inputEntry = primariesoffset;
	primariesoffset=0;
	//geniedata->GetEntry(inputEntry);
	
       G4cout<<"first run: "<<runbranchval<<G4endl;
        localEntry = geniedata->LoadTree(inputEntry);
        treeNumber=geniedata->GetTreeNumber();
        genieFileName = geniedata->GetCurrentFile()->GetName(); // new tree, new file
        std::cout <<"Current genieFileName: "<<genieFileName<<std::endl;
        char* genieFileNameAsChar = strdup(genieFileName.c_str());
        genieFileName = basename(genieFileNameAsChar);	

	loadNewGenie=false;

}	



void WCSimPrimaryGeneratorAction::LoadTalysFiles(){

	//Load nuclear de-excitation information from TALYS-generated root-files
	//Version of talys file depends on how many neutrons & protons have been knocked out of the nucleus
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

	G4cout <<"Read first file (015)"<<G4endl;
        f_O15 = new TFile(filepath_O15.c_str(),"READ");
	G4cout <<"Read second file (N15)"<<G4endl;
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
        talys_treemap.emplace("O15",talys_O15);
	talys_N15 = (TTree*) f_N15->Get("TreeNucldeex");
        talys_N15->SetName("talys_N15");
        talys_treemap.emplace("N15",talys_N15);
	talys_N14 = (TTree*) f_N14->Get("TreeNucldeex");
        talys_N14->SetName("talys_N14");
        talys_treemap.emplace("N14",talys_N14);
	talys_Li9 = (TTree*) f_Li9->Get("TreeNucldeex");
        talys_Li9->SetName("talys_Li9");
        talys_treemap.emplace("Li9",talys_Li9);
	talys_Li7 = (TTree*) f_Li7->Get("TreeNucldeex");
        talys_Li7->SetName("talys_Li7");
        talys_treemap.emplace("Li7",talys_Li7);
	talys_C14 = (TTree*) f_C14->Get("TreeNucldeex");
        talys_C14->SetName("talys_C14");
        talys_treemap.emplace("C14",talys_C14);
	talys_C13 = (TTree*) f_C13->Get("TreeNucldeex");
        talys_C13->SetName("talys_C13");
        talys_treemap.emplace("C13",talys_C13);
	talys_C11 = (TTree*) f_C11->Get("TreeNucldeex");
        talys_C11->SetName("talys_C11");
        talys_treemap.emplace("C11",talys_C11);
	talys_C10 = (TTree*) f_C10->Get("TreeNucldeex");
        talys_C10->SetName("talys_C10");
        talys_treemap.emplace("C10",talys_C10);
        talys_Be10 = (TTree*) f_Be10->Get("TreeNucldeex");
        talys_Be10->SetName("talys_Be10");
        talys_treemap.emplace("Be10",talys_Be10);
        talys_Be9 = (TTree*) f_Be9->Get("TreeNucldeex");
        talys_Be9->SetName("talys_Be9");
        talys_treemap.emplace("Be9",talys_Be9);
        talys_B11 = (TTree*) f_B11->Get("TreeNucldeex");
        talys_B11->SetName("talys_B11");
        talys_treemap.emplace("B11",talys_B11);
        talys_B10 = (TTree*) f_B10->Get("TreeNucldeex");
        talys_B10->SetName("talys_B10");
        talys_treemap.emplace("B10",talys_B10);
        talys_B9 = (TTree*) f_B9->Get("TreeNucldeex");
        talys_B9->SetName("talys_B9");
        talys_treemap.emplace("B9",talys_B9);
	
        talys_O15gamma = (TTree*) f_O15gamma->Get("TreeNucldeex");
        talys_O15gamma->SetName("talys_O15gamma");
        talys_gammatreemap.emplace("O15",talys_O15gamma);
	talys_N15gamma = (TTree*) f_N15gamma->Get("TreeNucldeex");
        talys_N15gamma->SetName("talys_N15gamma");
        talys_gammatreemap.emplace("N15",talys_N15gamma);
	talys_N14gamma = (TTree*) f_N14gamma->Get("TreeNucldeex");
        talys_N14gamma->SetName("talys_N14gamma");
        talys_gammatreemap.emplace("N14",talys_N14gamma);
	talys_Li9gamma = (TTree*) f_Li9gamma->Get("TreeNucldeex");
        talys_Li9gamma->SetName("talys_Li9gamma");
        talys_gammatreemap.emplace("Li9",talys_Li9gamma);
	talys_Li7gamma = (TTree*) f_Li7gamma->Get("TreeNucldeex");
        talys_Li7gamma->SetName("talys_Li7gamma");
        talys_gammatreemap.emplace("Li7",talys_Li7gamma);
	talys_C14gamma = (TTree*) f_C14gamma->Get("TreeNucldeex");
        talys_C14gamma->SetName("talys_C14gamma");
        talys_gammatreemap.emplace("C14",talys_C14gamma);
	talys_C13gamma = (TTree*) f_C13gamma->Get("TreeNucldeex");
        talys_C13gamma->SetName("talys_C13gamma");
        talys_gammatreemap.emplace("C13",talys_C13gamma);
	talys_C11gamma = (TTree*) f_C11gamma->Get("TreeNucldeex");
        talys_C11gamma->SetName("talys_C11gamma");
        talys_gammatreemap.emplace("C11",talys_C11gamma);
	talys_C10gamma = (TTree*) f_C10gamma->Get("TreeNucldeex");
        talys_C10gamma->SetName("talys_C10gamma");
        talys_gammatreemap.emplace("C10",talys_C10gamma);
        talys_Be10gamma = (TTree*) f_Be10gamma->Get("TreeNucldeex");
        talys_Be10gamma->SetName("talys_Be10gamma");
        talys_gammatreemap.emplace("Be10",talys_Be10gamma);
        talys_Be9gamma = (TTree*) f_Be9gamma->Get("TreeNucldeex");
        talys_Be9gamma->SetName("talys_Be9gamma");
        talys_gammatreemap.emplace("Be9",talys_Be9gamma);
        talys_B11gamma = (TTree*) f_B11gamma->Get("TreeNucldeex");
        talys_B11gamma->SetName("talys_B11gamma");
        talys_gammatreemap.emplace("B11",talys_B11gamma);
        talys_B10gamma = (TTree*) f_B10gamma->Get("TreeNucldeex");
        talys_B10gamma->SetName("talys_B10gamma");
        talys_gammatreemap.emplace("B10",talys_B10gamma);
        talys_B9gamma = (TTree*) f_B9gamma->Get("TreeNucldeex");
        talys_B9gamma->SetName("talys_B9gamma");
        talys_gammatreemap.emplace("B9",talys_B9gamma);

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

std::string WCSimPrimaryGeneratorAction::CalculateResNucleus(int resNuclA, int resNuclZ){

	std::string res_nucleus = "none";

        if (resNuclA == 15 && resNuclZ == 8) res_nucleus = "O15";
	if (resNuclA == 15 && resNuclZ == 7) res_nucleus = "N15";
	if (resNuclA == 14 && resNuclZ == 7) res_nucleus = "N14";
	if (resNuclA == 14 && resNuclZ == 6) res_nucleus = "C14";
	if (resNuclA == 13 && resNuclZ == 6) res_nucleus = "C13";
	if (resNuclA == 11 && resNuclZ == 6) res_nucleus = "C11";
	if (resNuclA == 10 && resNuclZ == 6) res_nucleus = "C10";
	if (resNuclA == 11 && resNuclZ == 5) res_nucleus = "B11";
	if (resNuclA == 10 && resNuclZ == 5) res_nucleus = "B10";
	if (resNuclA == 9 && resNuclZ == 5) res_nucleus = "B9";
	if (resNuclA == 10 && resNuclZ == 4) res_nucleus = "Be10";
	if (resNuclA == 9 && resNuclZ == 4) res_nucleus = "Be9";
	if (resNuclA == 9 && resNuclZ == 3) res_nucleus = "Li9";
	if (resNuclA == 7 && resNuclZ == 3) res_nucleus = "Li7";

	if (res_nucleus != "none"){
	  talys_current = talys_treemap[res_nucleus];
	  talys_currentgamma = talys_gammatreemap[res_nucleus];
	}

	return res_nucleus;

}

void WCSimPrimaryGeneratorAction::LoadDeexcitationProb(){

  //vectors storing information about nucleusXX.root
  
  deex_channel.clear();
  deex_part_pdg.clear();
  deex_part_energy.clear();

  //temporary storing vectors
    int channel;
  std::vector<double> *gamma_energy=nullptr;
  std::vector<double> *neutron_energy=nullptr;
  std::vector<double> *proton_energy=nullptr;
  std::vector<double> *deuteron_energy=nullptr;
  std::vector<double> *tritium_energy=nullptr;
  std::vector<double> *helium_energy=nullptr;
  std::vector<double> *alpha_energy=nullptr;
  std::vector<int> pdg;
  std::vector<double> energy;
  int n_entries=0;

  if (talys_current){ //make sure tree exists and was read in correctly
 talys_current->SetBranchAddress("Channel",&channel);
    talys_current->SetBranchAddress("GammaEnergy",&gamma_energy);
    talys_current->SetBranchAddress("NeutronEnergy",&neutron_energy);
    talys_current->SetBranchAddress("ProtonEnergy",&proton_energy);
    talys_current->SetBranchAddress("DeuteronEnergy",&deuteron_energy);
    talys_current->SetBranchAddress("TritiumEnergy",&tritium_energy);
    talys_current->SetBranchAddress("Helium3Energy",&helium_energy);
    talys_current->SetBranchAddress("AlphaEnergy",&alpha_energy);
    n_entries = (int) talys_current->GetEntries();
  }

  if (n_entries == 0){
    G4cerr << "Nuclear de-excitation data could not be loaded"<<G4endl;
    return;
  }

  for (int i_entry=0; i_entry < n_entries; i_entry++){

    talys_current->GetEntry(i_entry);
    pdg.clear();
    energy.clear();
    deex_channel.push_back(channel);
    for (unsigned int i=0; i< gamma_energy->size(); i++){
      pdg.push_back(22);
      energy.push_back(gamma_energy->at(i));
    }
    for (unsigned int i=0; i< neutron_energy->size(); i++){
      pdg.push_back(2112);
      energy.push_back(neutron_energy->at(i));
    }
    for (unsigned int i=0; i< proton_energy->size(); i++){
      pdg.push_back(2212);
      energy.push_back(proton_energy->at(i));
    }
    for (unsigned int i=0; i< deuteron_energy->size(); i++){
      pdg.push_back(1000010020);
      energy.push_back(deuteron_energy->at(i));
    }
    for (unsigned int i=0; i< tritium_energy->size(); i++){
      pdg.push_back(1000010030);
      energy.push_back(tritium_energy->at(i));
    }
    for (unsigned int i=0; i< helium_energy->size(); i++){
      pdg.push_back(1000020030);
      energy.push_back(helium_energy->at(i));
    }
    for (unsigned int i=0; i< alpha_energy->size(); i++){
      pdg.push_back(1000020040);
      energy.push_back(alpha_energy->at(i));
    }

    deex_part_pdg.push_back(pdg);
    deex_part_energy.push_back(energy);

  }
  
  //vectors storing information about nucleusgammaXX.root
  deex_resnuclZ.clear();
  deex_resnuclA.clear();
  deex_resnuclLevel.clear();
  deex_resnuclEnergy.clear();
  deex_resnuclPopulation.clear();

  //temp vectors for tree information
  int Z,A;
  std::vector<int> tempLevel;
  std::vector<double> tempEnergy;
  std::vector<double> tempPopulation;
  std::vector<int>* pLevel=nullptr;
  std::vector<double>* pEnergy=nullptr;
  std::vector<double>* pPopulation=nullptr;
  n_entries = 0;

  //Loop through talys gamma tree and store the information
    if (talys_currentgamma){
    talys_currentgamma->SetBranchAddress("ResNuclZ",&Z);
    talys_currentgamma->SetBranchAddress("ResNuclA",&A);
    talys_currentgamma->SetBranchAddress("ResNuclLevel",&pLevel);
    talys_currentgamma->SetBranchAddress("ResNuclEnergy",&pEnergy);   
    talys_currentgamma->SetBranchAddress("ResNuclPopulation",&pPopulation);   
    n_entries = (int) talys_currentgamma->GetEntries();
  } 

  if (n_entries==0){
    G4cerr <<"Nuclear de-excitation data could not be loaded."<<G4endl;
    return;
  }

  for (int i_entry=0; i_entry < n_entries; i_entry++){
    talys_currentgamma->GetEntry(i_entry);
    deex_resnuclZ.push_back(Z);
    deex_resnuclA.push_back(A);
    tempLevel.clear();
    tempEnergy.clear();
    tempPopulation.clear();
    for (unsigned int i = 0; i< pLevel->size(); i++){
      tempLevel.push_back(pLevel->at(i));
      tempEnergy.push_back(pEnergy->at(i));
      tempPopulation.push_back(pPopulation->at(i));
    }
    deex_resnuclLevel.push_back(tempLevel);
    deex_resnuclEnergy.push_back(tempEnergy);
    deex_resnuclPopulation.push_back(tempPopulation);
  }

}

void WCSimPrimaryGeneratorAction::GenerateDeexcitation(std::vector<int> *Talys_pdg, std::vector<G4ThreeVector> *Talys_momdir, std::vector<double> *Talys_energy, int resNuclA, int resNuclZ){

  //GENERAL Deexcitation Tree
  
  int rand,Z,A,channel;
  A = resNuclA;
  Z = resNuclZ; 
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  rand = (int)floor(G4UniformRand()*deex_part_pdg.size());

  Talys_pdg->insert(Talys_pdg->end(),deex_part_pdg[rand].begin(),deex_part_pdg[rand].end());
  Talys_energy->insert(Talys_energy->end(),deex_part_energy[rand].begin(),deex_part_energy[rand].end());

  channel = deex_channel[rand];

  for (unsigned int i=0; i< deex_part_pdg[rand].size(); i++){
  
    Talys_momdir->push_back(G4RandomDirection());
    G4ParticleDefinition *particle = particleTable->FindParticle(deex_part_pdg[rand][i]);
    Z-=particle->GetPDGCharge();
    A-=particle->GetBaryonNumber();

  }

  int nuclDeexZ = Z;
  int nuclDeexA = A;

  //GAMMA Deexcitation Tree
  
  double sum_branch=0,x=0;
  double energy_temp;
  std::vector<double> branch_ratio;
  int n=0;

  //Get population of discrete energy levels of residual nucleus
  for (unsigned int i_Z=0; i_Z < deex_resnuclZ.size(); i_Z++){
    if (deex_resnuclZ.at(i_Z)==nuclDeexZ && deex_resnuclA.at(i_Z)==nuclDeexA){
      n=i_Z;
      break;
    }
  }

  //Calculate branching rtio of energy levels from populations of levels
    for (unsigned int i_level=0; i_level < deex_resnuclLevel[n].size(); i_level++){
    sum_branch+=deex_resnuclPopulation[n][i_level];
    branch_ratio.push_back(sum_branch);
  } 

  x = G4UniformRand()*sum_branch;
  for (unsigned int i = 0; i< branch_ratio.size(); i++){
    if (x < branch_ratio[i]){
      energy_temp = deex_resnuclEnergy[n][i];
    }
  }

  //Add Gamma with random momentum to de-excitation particle list
  if (energy_temp){
    Talys_pdg->push_back(22);
    Talys_momdir->push_back(G4RandomDirection());
    Talys_energy->push_back(energy_temp);
  }
}

