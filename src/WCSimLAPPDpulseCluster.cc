#include "WCSimLAPPDpulse.hh"
#include "WCSimLAPPDpulseCluster.hh"
#include "TObject.h"
#include "TString.h"
#include "TFile.h"
#include "TH1.h"
#include "TF1.h"


#include <iostream>
#include <cmath>


ClassImp(WCSimLAPPDpulse)

WCSimLAPPDpulseCluster::WCSimLAPPDpulseCluster()
{
  this->Reset();

  for(int i=0; i<30; i++){
    fLAPPDstripPulseCount.push_back(0);	
    std::vector<int> myinput;
    fLAPPDstripPulseCoordinate.push_back(myinput);
  }
//std::cout<<"the vector array SIZE: "<<fLAPPDstripPulseCoordinate.size()<<std::endl;
}

WCSimLAPPDpulseCluster::~WCSimLAPPDpulseCluster()
{
  this->Reset();
}

void WCSimLAPPDpulseCluster::Reset()
{
  for(auto apulse : fLAPPDpulseList){
    delete apulse;
    apulse=0;
  }
  fLAPPDpulseList.clear();
  fLAPPDstripPulseCount.clear();
  fLAPPDstripPulseCoordinate.clear();
  _nchannels=0;

}

/* 
static bool CompareStripnums(LAPPDpulse *sn1, LAPPDpulse *sn2)
{
  return ( sn1->Getstripnum() > sn2->Getstripnum() );
}


void LAPPDpulseCluster::SortCluster()
{
  std::sort(fLAPPDpulseList.begin(), fLAPPDpulseList.end(), CompareStripnums);
}
*/

void WCSimLAPPDpulseCluster::AddPulse(WCSimLAPPDpulse* pulse)
{
  int CHnum = pulse->Getstripnum();
  fLAPPDpulseList.push_back(pulse);
  if(fLAPPDstripPulseCount.at(CHnum-1)==0) _nchannels++;
  fLAPPDstripPulseCount.at(CHnum-1)++;
  fLAPPDstripPulseCoordinate[(CHnum-1)].push_back(fLAPPDpulseList.size()-1);	
}

WCSimLAPPDpulse* WCSimLAPPDpulseCluster::GetPulse(Int_t n)
{
  return (WCSimLAPPDpulse*)(fLAPPDpulseList.at(n));
}


Int_t WCSimLAPPDpulseCluster::GetNPulsesStrip(int stripnum)
{
  return fLAPPDstripPulseCoordinate[(stripnum-1)].size();
}


Int_t WCSimLAPPDpulseCluster::GetPulseNum(int stripnum, int pulsenum)
{
  if( (pulsenum+1) > fLAPPDstripPulseCoordinate[(stripnum-1)].size() ){
     std::cout<<"STRIPNUM TOO HIGH"<<std::endl;
     return -5555;
  } else{ 
     return fLAPPDstripPulseCoordinate[(stripnum-1)].at(pulsenum);
  }
}


Int_t WCSimLAPPDpulseCluster::GetNPulses()
{
  return fLAPPDpulseList.size();
}


