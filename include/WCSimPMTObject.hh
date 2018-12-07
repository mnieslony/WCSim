#ifndef WCSimWCPMTObject_h
#define WCSimWCPMTObject_h 1

#include "WCSimDetectorConstruction.hh"
#include "globals.hh"
#include "Randomize.hh"
#include <map>
#include <vector>


class WCSimPMTObject
{

public:
  virtual G4String GetPMTName()=0;
  virtual G4double GetExposeHeight()=0;
  virtual G4double GetRadius()=0;
  virtual G4float* Getqpe()=0;
  virtual G4float* GetQE()=0;
  virtual G4float* GetQEWavelength()=0;
  virtual G4float  GetmaxQE()=0;
  virtual G4float  GetCollectionEfficiency(float);
  virtual float    HitTimeSmearing(float)=0;
  virtual G4double GetPMTGlassThickness()=0;
  virtual G4float  GetDarkRate()=0;
  virtual G4float  GetDarkRateConversionFactor()=0;
  virtual G4double GetGelThickness(){ return 0.;};
  virtual G4double GetShamferRadius(){ return 0.;};
protected:
  virtual G4float* GetCollectionEfficiencyArray();
  virtual G4float* GetCollectionEfficiencyAngle();
  G4float Interpolate_func(G4float, G4int, G4float*, G4float*);
};

class PMT20inch : public WCSimPMTObject
{

public:
  
  PMT20inch();
  ~PMT20inch();
 
public:
  G4String GetPMTName() ;
  G4double GetExposeHeight();
  G4double GetRadius();
  G4float* Getqpe();
  G4float* GetQE();
  G4float* GetQEWavelength();
  G4float  GetmaxQE();
  float    HitTimeSmearing(float);
  G4double GetPMTGlassThickness();
  G4float  GetDarkRate();
  G4float  GetDarkRateConversionFactor();


};

class PMT8inch : public WCSimPMTObject
{

public:
  
  PMT8inch();
  ~PMT8inch();
 
public:
  G4String GetPMTName(); 
  G4double GetExposeHeight(); 
  G4double GetRadius();
  G4float* Getqpe();
  G4float* GetQE();
  G4float* GetQEWavelength();
  G4float  GetmaxQE();
  float    HitTimeSmearing(float);
  G4double GetPMTGlassThickness();
  G4float  GetDarkRate();
  G4float  GetDarkRateConversionFactor();
};

 class PMT10inch : public WCSimPMTObject
{

public: 
  PMT10inch();
  ~PMT10inch();
 
public:
  G4String GetPMTName(); 
  G4double GetExposeHeight(); 
  G4double GetRadius(); 
  G4float* Getqpe();
  G4float* GetQE();
  G4float* GetQEWavelength();
  G4float  GetmaxQE();
  float    HitTimeSmearing(float);
  G4double GetPMTGlassThickness();
  G4float  GetDarkRate();
  G4float  GetDarkRateConversionFactor();
 };

 class PMT10inchHQE : public WCSimPMTObject
{

public: 
  PMT10inchHQE();
  ~PMT10inchHQE();
 
public:
  G4String GetPMTName(); 
  G4double GetExposeHeight(); 
  G4double GetRadius(); 
  G4float* Getqpe();
  G4float* GetQE();
  G4float* GetQEWavelength();
  G4float  GetmaxQE();
  float    HitTimeSmearing(float);
  G4double GetPMTGlassThickness();
  G4float  GetDarkRate();
  G4float  GetDarkRateConversionFactor();
 };

 class PMT12inchHQE : public WCSimPMTObject
{

public: 
  PMT12inchHQE();
  ~PMT12inchHQE();
 
public:
  G4String GetPMTName(); 
  G4double GetExposeHeight(); 
  G4double GetRadius(); 
  G4float* Getqpe();
  G4float* GetQE();
  G4float* GetQEWavelength();
  G4float  GetmaxQE();
  float    HitTimeSmearing(float);
  G4double GetPMTGlassThickness();
  G4float  GetDarkRate();
  G4float  GetDarkRateConversionFactor();
 };

class HPD20inchHQE : public WCSimPMTObject
{

public:
  
  HPD20inchHQE();
  ~HPD20inchHQE();
 
public:
  G4String GetPMTName() ;
  G4double GetExposeHeight();
  G4double GetRadius();
  G4float* Getqpe();
  G4float* GetQE();
  G4float* GetQEWavelength();
  G4float  GetmaxQE();
  float    HitTimeSmearing(float);
  G4double GetPMTGlassThickness();
  G4float  GetDarkRate();
  G4float  GetDarkRateConversionFactor();
protected:
  G4float* GetCollectionEfficiencyArray();
};

class HPD12inchHQE : public WCSimPMTObject
{

public:
  
  HPD12inchHQE();
  ~HPD12inchHQE();
 
public:
  G4String GetPMTName() ;
  G4double GetExposeHeight();
  G4double GetRadius();
  G4float* Getqpe();
  G4float* GetQE();
  G4float* GetQEWavelength();
  G4float  GetmaxQE();
  float    HitTimeSmearing(float);
  G4double GetPMTGlassThickness();
  G4float  GetDarkRate();
  G4float  GetDarkRateConversionFactor();
protected:
  G4float* GetCollectionEfficiencyArray();
};

class BoxandLine20inchHQE : public WCSimPMTObject
{

public:
  
  BoxandLine20inchHQE();
  ~BoxandLine20inchHQE();
 
public:
  G4String GetPMTName() ;
  G4double GetExposeHeight();
  G4double GetRadius();
  G4float* Getqpe();
  G4float* GetQE();
  G4float* GetQEWavelength();
  G4float  GetmaxQE();
  float    HitTimeSmearing(float);
  G4double GetPMTGlassThickness();
  G4float  GetDarkRate();
  G4float  GetDarkRateConversionFactor();
protected:
  G4float* GetCollectionEfficiencyArray();
};

class BoxandLine12inchHQE : public WCSimPMTObject
{

public:
  
  BoxandLine12inchHQE();
  ~BoxandLine12inchHQE();
 
public:
  G4String GetPMTName() ;
  G4double GetExposeHeight();
  G4double GetRadius();
  G4float* Getqpe();
  G4float* GetQE();
  G4float* GetQEWavelength();
  G4float  GetmaxQE();
  float    HitTimeSmearing(float);
  G4double GetPMTGlassThickness();
  G4float  GetDarkRate();
  G4float  GetDarkRateConversionFactor();
protected:
  G4float* GetCollectionEfficiencyArray();
};

class FlatFacedPMT2inch : public WCSimPMTObject
{
public:
FlatFacedPMT2inch();
~FlatFacedPMT2inch();

public:
  G4String GetPMTName(); 
  G4double GetExposeHeight(); 
  G4double GetRadius();
  G4float* Getqpe();
  G4float* GetQE();
  G4float* GetQEWavelength();
  G4float  GetmaxQE();
  float    HitTimeSmearing(float);
  G4double GetPMTGlassThickness();
  G4float  GetDarkRate();
  G4float  GetDarkRateConversionFactor();
  G4double GetGelThickness();
  G4double GetShamferRadius();
};

class FlatFacedPMT4inch : public WCSimPMTObject
{
public:
FlatFacedPMT4inch();
~FlatFacedPMT4inch();

public:
  G4String GetPMTName(); 
  G4double GetExposeHeight(); 
  G4double GetRadius();
  G4float* Getqpe();
  G4float* GetQE();
  G4float* GetQEWavelength();
  G4float  GetmaxQE();
  float    HitTimeSmearing(float);
  G4double GetPMTGlassThickness();
  G4float  GetDarkRate();
  G4float  GetDarkRateConversionFactor();
  G4double GetGelThickness();
  G4double GetShamferRadius();
};

class PMT1cm : public WCSimPMTObject
{
public:
PMT1cm();
~PMT1cm();

public:
  G4String GetPMTName(); 
  G4double GetExposeHeight(); 
  G4double GetRadius();
  G4float* Getqpe();
  G4float* GetQE();
  G4float* GetQEWavelength();
  G4float  GetmaxQE();
  float    HitTimeSmearing(float);
  G4double GetPMTGlassThickness();
  G4float  GetDarkRate();
  G4float  GetDarkRateConversionFactor();
};

class PMT_R5912HQE : public WCSimPMTObject
{
public:
PMT_R5912HQE();
~PMT_R5912HQE();

public:
  G4String GetPMTName(); 
  G4double GetExposeHeight(); 
  G4double GetRadius();
  G4float* Getqpe();
  G4float* GetQE();
  G4float* GetQEWavelength();
  G4float  GetmaxQE();
  float    HitTimeSmearing(float);
  G4double GetPMTGlassThickness();
  G4float  GetDarkRate();
  G4float  GetDarkRateConversionFactor();
protected:
  G4float* GetCollectionEfficiencyArray();
};

class PMT_D784KFLB : public WCSimPMTObject
{
public:
PMT_D784KFLB();
~PMT_D784KFLB();

public:
  G4String GetPMTName(); 
  G4double GetExposeHeight(); 
  G4double GetRadius();
  G4float* Getqpe();
  G4float* GetQE();
  G4float* GetQEWavelength();
  G4float  GetmaxQE();
  float    HitTimeSmearing(float);
  G4double GetPMTGlassThickness();
  G4float  GetDarkRate();
  G4float  GetDarkRateConversionFactor();
protected:
  G4float* GetCollectionEfficiencyArray();
};

class PMT_R7081 : public WCSimPMTObject
{
public:
PMT_R7081();
~PMT_R7081();

public:
  G4String GetPMTName(); 
  G4double GetExposeHeight(); 
  G4double GetRadius();
  G4float* Getqpe();
  G4float* GetQE();
  G4float* GetQEWavelength();
  G4float  GetmaxQE();
  float    HitTimeSmearing(float);
  G4double GetPMTGlassThickness();
  G4float  GetDarkRate();
  G4float  GetDarkRateConversionFactor();
protected:
  G4float* GetCollectionEfficiencyArray();
};

 class PMT_R7081HQE : public WCSimPMTObject  // actually just the same as PMT10inchHQE
{

public: 
  PMT_R7081HQE();
  ~PMT_R7081HQE();
 
public:
  G4String GetPMTName(); 
  G4double GetExposeHeight(); 
  G4double GetRadius(); 
  G4float* Getqpe();
  G4float* GetQE();
  G4float* GetQEWavelength();
  G4float  GetmaxQE();
  float    HitTimeSmearing(float);
  G4double GetPMTGlassThickness();
  G4float  GetDarkRate();
  G4float  GetDarkRateConversionFactor();
protected:
  G4float* GetCollectionEfficiencyArray();
 };

#endif
