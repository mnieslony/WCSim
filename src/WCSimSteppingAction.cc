#include <stdlib.h>
#include <stdio.h>

#include "WCSimSteppingAction.hh"

#include "G4Track.hh"
#include "G4VProcess.hh"
#include "G4VParticleChange.hh"
#include "G4SteppingVerbose.hh"
#include "G4SteppingManager.hh"
#include "G4PVParameterised.hh"
#include "G4PVReplica.hh"
#include "G4SDManager.hh"
#include "G4RunManager.hh"

#include "mrdPMTSD.hh"
#include "faccPMTSD.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4TrackStatus.hh"

void WCSimSteppingAction::UserSteppingAction(const G4Step* aStep)
{
  //DISTORTION must be used ONLY if INNERTUBE or INNERTUBEBIG has been defined in BidoneDetectorConstruction.cc
  
  const G4Event* evt = G4RunManager::GetRunManager()->GetCurrentEvent();

  const G4Track* track       = aStep->GetTrack();
  G4VPhysicalVolume* volume  = track->GetVolume();
  // G4cout<<"SteppingAction: physical volume is: "<<volume->GetName()<<G4endl;

  G4SDManager* SDman   = G4SDManager::GetSDMpointer();
  G4HCofThisEvent* HCE = evt->GetHCofThisEvent();

  //debugging 
//  G4Track* theTrack = aStep->GetTrack();
//  const G4DynamicParticle* aParticle = theTrack->GetDynamicParticle();
//  G4ThreeVector aMomentum = aParticle->GetMomentumDirection();
//  G4double vx = aMomentum.x();
//  G4int ix = std::isnan(vx);
//  if(ix != 0){
//    G4cout << " PROBLEM! " << theTrack->GetCreatorProcess()->GetProcessName() <<
//  std::flush << G4endl;
//  }

  // record hits manually for photons crossing boundaries into 'pmt' surfaces in the veto / mrd sandbox style
  // ========================================================================================================
  
  // retrieve the pointer to the optical boundary process: this code only gets executed once
  // use 'static' variable - will retain it's value between calls ...
  static G4OpBoundaryProcess* boundary=NULL;
  // ... so that once boundary has a value, the following will be skipped:
  // Retrieve pointer to the 'optical boundary process' so that we may retrive it's status
  // (invoked, what type of action etc) at the end of each step
  if(!boundary){
    G4ProcessManager* pm = aStep->GetTrack()->GetDefinition()->GetProcessManager();
    G4int nprocesses = pm->GetProcessListLength();
    G4ProcessVector* pv = pm->GetProcessList();
    G4int i;
    for( i=0;i<nprocesses;i++){
      if((*pv)[i]->GetProcessName()=="OpBoundary"){
        boundary = (G4OpBoundaryProcess*)(*pv)[i];
        //G4cout<<"optical boundary process: "<<boundary<<G4endl;
        break;
      }
    }
  }

  G4StepPoint* thePostPoint = aStep->GetPostStepPoint();
  G4VPhysicalVolume* thePostPV = thePostPoint->GetPhysicalVolume();
  // NOTE: if the step is limited by hitting a boundary, thePostPV will point to the volume it is ENTERING, 
  // NOT the volume this STEP WAS THROUGH.
  if(!thePostPV){
    //out of world
    fExpectedNextStatus=Undefined;
    return;
  }
  
  
  G4ParticleDefinition* particleType = track->GetDefinition();
  if(particleType==G4OpticalPhoton::OpticalPhotonDefinition()){ 	//Optical photon only
   // for photons entering the 'sensitive' volumes
    if((thePostPV->GetName()=="mrdsdsurf_phys")||(thePostPV->GetName()=="vetoSurface_phys") ){
      //|| (thePostPV->GetName()=="WCLAPPD") ){	// &&(thePostPoint->GetProcessDefinedStep()->GetProcessName() == "OpAbsorption")
      aStep->GetTrack()->SetTrackStatus(fStopAndKill);	//needs a non-const track pointer
      G4OpBoundaryProcessStatus boundaryStatus=Undefined;
      boundaryStatus=boundary->GetStatus();
      // double check to ensure particle at boundary (may not be valid for vers <Geant4.6.0-p1? omittable >?)
      if(thePostPoint->GetStepStatus()==fGeomBoundary){
        if(fExpectedNextStatus==StepTooSmall){
          if(boundaryStatus!=StepTooSmall){
            G4cout<<"expected status step too small, but boundary status isn't!"<<G4endl;
            /*G4ExceptionDescription ed;     // exceptiondescription not in geant4.9.4
            ed << "SteppingAction::UserSteppingAction(): "
               << "No reallocation step after reflection!"
               << G4endl;
            G4Exception("SteppingAction::UserSteppingAction()", "Expl01", FatalException,ed,
            "Something is wrong with the surface normal or geometry"); */
            G4Exception("SteppingAction::UserSteppingAction()","Expl01",FatalException,"No reallocation step after reflection!");
            // G4Exception(const char* issure, const char* errorCode, G4ExceptionSeverity severity, const char* comments);
          }
        }
        fExpectedNextStatus=Undefined;
        // perform actions based on boundary status... 
        switch(boundaryStatus){
          case Absorption:
              //G4cout<<"Absorption of optical photon on boundary of "<<thePostPV->GetName()<<G4endl;
              //info->SetMRDdetected(boundaryAbsorbed);
              break;
          case Detection: {
              //G4cout<<"Detection of optical photon on boundary of "<<thePostPV->GetName()<<G4endl;
              G4SDManager* SDman = G4SDManager::GetSDMpointer();
              if(thePostPV->GetName()=="mrdsdsurf_phys"){
                G4String sdName="MRDPMTSD";
                mrdPMTSD* pmtSD = (mrdPMTSD*)SDman->FindSensitiveDetector(sdName);
                if(pmtSD){
                  pmtSD->ProcessHits_PhotonHit(aStep, NULL);
                  //G4cout<<"called process photon hit!"<<G4endl;
                  //info->SetMRDdetected(hitPMT);
                }
              } else if(thePostPV->GetName()=="vetoSurface_phys"){
                G4String sdName="FACCPMTSD";
                faccPMTSD* pmtSD = (faccPMTSD*)SDman->FindSensitiveDetector(sdName);
                if(pmtSD){
                  pmtSD->ProcessHits_PhotonHit(aStep, NULL);
                  //G4cout<<"called process photon hit!"<<G4endl;
                  //info->SetFACCdetected(hitPMT);
                }
              }
              break;
          }
          case FresnelReflection:
          case TotalInternalReflection:
          case LambertianReflection:
          case LobeReflection:
          case SpikeReflection:
             // trackInformation->IncReflections();
             break;
          case BackScattering: {
               // trackInformation->IncReflections();
               fExpectedNextStatus=StepTooSmall;
               //G4cout<<"boundary status backscattering - ignore the next step (too small)"<<G4endl;
               break;
          }
          default:
              break;
        }	//close switch statement
     } // post-step status not on boundary
   } // post-step volume was not in a 'sensitive' volume
  } // not an optical photon
}


G4int WCSimSteppingAction::G4ThreeVectorToWireTime(G4ThreeVector *pos3d,
						    G4ThreeVector lArPos,
						    G4ThreeVector start,
						    G4int i)
{
  G4double x0 = start[0]-lArPos[0];
  G4double y0 = start[1]-lArPos[1];
//   G4double y0 = 2121.3;//mm
  G4double z0 = start[2]-lArPos[2];

  G4double dt=0.8;//mm
//   G4double midt = 2651.625;
  G4double pitch = 3;//mm
//   G4double midwir = 1207.10;
  G4double c45 = 0.707106781;
  G4double s45 = 0.707106781;

  G4double w1;
  G4double w2;
  G4double t;

//   G4double xField(0.);
//   G4double yField(0.);
//   G4double zField(0.);

//   if(detector->getElectricFieldDistortion())
//     {
//       Distortion(pos3d->getX(),
// 		 pos3d->getY());
      
//       xField = ret[0];
//       yField = ret[1];
//       zField = pos3d->getZ();
      
//       w1 = (int)(((zField+z0)*c45 + (x0-xField)*s45)/pitch); 
//       w2 = (int)(((zField+z0)*c45 + (x0+xField)*s45)/pitch); 
//       t = (int)(yField+1);

//       //G4cout<<" x orig "<<pos3d->getX()<<" y orig "<<(pos3d->getY()+y0)/dt<<G4endl;
//       //G4cout<<" x new "<<xField<<" y new "<<yField<<G4endl;
//     }
//   else 
//     {
      
  w1 = (int) (((pos3d->getZ()+z0)*c45 + (x0-pos3d->getX())*s45)/pitch); 
  w2 = (int)(((pos3d->getZ()+z0)*c45 + (x0+pos3d->getX())*s45)/pitch); 
  t  = (int)((pos3d->getY()+y0)/dt +1);
//     }

  if (i==0)
    return (int)w1;
  else if (i==1)
    return (int)w2;
  else if (i==2)
    return (int)t;
  else return 0;
} 


void WCSimSteppingAction::Distortion(G4double /*x*/,G4double /*y*/)
{
 
//   G4double theta,steps,yy,y0,EvGx,EvGy,EField,velocity,tSample,dt;
//   y0=2121.3;//mm
//   steps=0;//1 mm steps
//   tSample=0.4; //micros
//   dt=0.8;//mm
//   LiquidArgonMedium medium;  
//   yy=y;
//   while(y<y0 && y>-y0 )
//     {
//       EvGx=FieldLines(x,y,1);
//       EvGy=FieldLines(x,y,2);
//       theta=atan(EvGx/EvGy);
//       if(EvGy>0)
// 	{
// 	  x+=sin(theta);
// 	  y+=cos(theta);
// 	}
//       else
// 	{
// 	  y-=cos(theta);
// 	  x-=sin(theta);
// 	}
//       EField=sqrt(EvGx*EvGx+EvGy*EvGy);//kV/mm
//       velocity=medium.DriftVelocity(EField*10);// mm/microsec
//       steps+=1/(tSample*velocity);

//       //G4cout<<" step "<<steps<<" x "<<x<<" y "<<y<<" theta "<<theta<<" Gx "<<eventaction->Gx->Eval(x,y)<<" Gy "<<eventaction->Gy->Eval(x,y)<<" EField "<<EField<<" velocity "<<velocity<<G4endl;
//     }

//   //numbers
//   //EvGx=FieldLines(0,1000,1);
//   //EvGy=FieldLines(0,1000,2);
//   //EField=sqrt(EvGx*EvGx+EvGy*EvGy);//kV/mm
//   //velocity=medium.DriftVelocity(EField*10);// mm/microsec
//   //G4double quenching;
//   //quenching=medium.QuenchingFactor(2.1,EField*10);
//   //G4cout<<" Gx "<<EvGx<<" Gy "<<EvGy<<" EField "<<EField<<" velocity "<<velocity<<" quenching "<<quenching<<G4endl;


//   ret[0]=x;
//   if(yy>0)
//     ret[1]=2*y0/dt -steps;
//   else
//     ret[1]=steps; 
}


double WCSimSteppingAction::FieldLines(G4double /*x*/,G4double /*y*/,G4int /*coord*/)
{ //0.1 kV/mm = field
  //G4double Radius=302;//mm
//   G4double Radius=602;//mm
//   if(coord==1) //x coordinate
//     return  (0.1*(2*Radius*Radius*x*abs(y)/((x*x+y*y)*(x*x+y*y))));
//   else //y coordinate
//     return 0.1*((abs(y)/y)*(1-Radius*Radius/((x*x+y*y)*(x*x+y*y))) + abs(y)*(2*Radius*Radius*y/((x*x+y*y)*(x*x+y*y))));
  return 0;
}

G4String ToName(G4OpBoundaryProcessStatus boundaryStatus){
G4String processnames[14]={"Undefined","FresnelRefraction","FresnelReflection","TotalInternalReflection","LambertianReflection","LobeReflection","SpikeReflection","BackScattering","Absorption","Detection","NotAtBoundary","SameMaterial","StepTooSmall","NoRINDEX"};
return processnames[boundaryStatus];}

G4String ToName2(G4StepStatus stepStatus){
G4String processnames[8]={"fWorldBoundary","fGeomBoundary","fAtRestDoItProc","fAlongStepDoItProc","fPostStepDoItProc","fUserDefinedLimit","fExclusivelyForcedProc","fUndefined"};
return processnames[stepStatus];}
