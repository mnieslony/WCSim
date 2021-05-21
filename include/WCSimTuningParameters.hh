#ifndef WCSimTuningParameters_h
#define WCSimTuningParameters_h 1
#include "WCSimTuningMessenger.hh"
#include "WCSimRootOptions.hh"
#include "globals.hh"

class WCSimTuningParameters
{
public:
  WCSimTuningParameters();
  ~WCSimTuningParameters();


  // Setters and getters
  G4double GetRayff() {return rayff;}
  void SetRayff(G4double rparam) {rayff=rparam;}

  G4double GetBsrff() {return bsrff;}
  void SetBsrff(G4double rparam) {bsrff=rparam;}

  G4double GetAbwff() {return abwff;}
  void SetAbwff(G4double rparam) {abwff=rparam;}

  G4double GetRgcff() {return rgcff;}
  void SetRgcff(G4double rparam) {rgcff=rparam;}

  G4double GetMieff() {return mieff;}
  void SetMieff(G4double rparam) {mieff=rparam;}

  //ANNIE-specific setters and getters
  G4double GetTeflonrff() {return teflonrff;}
  void SetTeflonrff(G4double rparam) {teflonrff=rparam;}

  G4double GetHolderrff() {return holderrff;}
  void SetHolderrff(G4double rparam) {holderrff=rparam;}

  G4double GetHolderrffLUX() {return holderrfflux;}
  void SetHolderrffLUX(G4double rparam) {holderrfflux=rparam;}

  G4double GetLinerrff() {return linerrff;}
  void SetLinerrff(G4double rparam) {linerrff=rparam;}

  G4bool GetHolder() {return holder;}
  void SetHolder(G4double hparam) {holder=hparam;}

  //For Top Veto - jl145
  G4double GetTVSpacing() {return tvspacing;}
  void SetTVSpacing(G4double tparam) {tvspacing=tparam;}

  G4bool GetTopVeto() {return topveto;}
  void SetTopVeto(G4double tparam) {topveto=tparam;}

  void SaveOptionsToOutput(WCSimRootOptions * wcopt);

private:

  // The messenger
  WCSimTuningMessenger* TuningMessenger;

  // The parameters that need to be set before WCSimDetectorConstruction
  // is created

  G4double rayff;
  G4double bsrff;
  G4double abwff;
  G4double rgcff;
  G4double mieff;

  // ANNIE-specfic tuning parameters
  G4double teflonrff;   //Teflon-wrapped Inner Structure --> Reflectivity tuning factor
  G4double holderrff;   //ANNIE holders --> Reflectivity tuning factor
  G4double holderrfflux;   //LUX/ETEL holders --> Reflectivity tuning factor
  G4double linerrff;    //Liner --> Reflectivity tuning factor
  G4bool holder;        //Should ANNIE PMT holders be implemented in the simulation?

  //For Top Veto - jl145
  G4double tvspacing;
  G4bool topveto;

};

#endif







