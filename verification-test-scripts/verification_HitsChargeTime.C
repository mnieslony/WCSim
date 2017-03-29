#include <iostream>
#include <TH1F.h>
#include <stdio.h>     
#include <stdlib.h>    

// Simple example of reading a generated Root file
int verification_HitsChargeTime(std::string filename, std::string directoryname, int &numtubeshit, int &numdigits, int &numpes, int &numevents, TH1F *hits, TH1F *time, TH1F *charge, bool verbose=false)
{
  
  //for(int file=0; file<2; file++){
    
//    // Clear global scope
//    cout<<"clearing global scope"<<endl;
//    gROOT->Reset();
    
    // Set style preferences 
    gStyle->SetOptStat(0);
    gStyle->SetCanvasColor(0);
    gStyle->SetTitleColor(1);
    gStyle->SetStatColor(0);
    gStyle->SetFrameFillColor(0);
    gStyle->SetPadColor(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetTitleSize(0.04);
    gStyle->SetCanvasBorderMode(0);
    gStyle->SetFrameBorderMode(0);
    gStyle->SetFrameLineWidth(2);
    gStyle->SetPadBorderMode(0);
    gStyle->SetPalette(1);
    gStyle->SetTitleAlign(23);
    gStyle->SetTitleX(.5);
    gStyle->SetTitleY(0.99);
    gStyle->SetTitleBorderSize(0);
    gStyle->SetTitleFillColor(0);
    gStyle->SetHatchesLineWidth(2);
    gStyle->SetLineWidth(1.5);
    gStyle->SetTitleFontSize(0.07);
    gStyle->SetLabelSize(0.05,"X");
    gStyle->SetLabelSize(0.05,"Y");
    gStyle->SetTitleSize(0.04,"X");
    gStyle->SetTitleSize(0.04,"Y");
    gStyle->SetTitleBorderSize(0);
    gStyle->SetCanvasBorderMode(0);
    
    cout<<"opening file"<<endl;
    TFile *f = new TFile(filename.c_str(),"read");
    if (!f->IsOpen()){
      cout << "Error, could not open input file: " << filename << endl;
      return -1;
    }
    
    cout<<"getting tree"<<endl;
    TTree  *wcsimT = (TTree*)f->Get("wcsimT");
    int nevent = wcsimT->GetEntries();
    cout<<nevent<<" events in the tree"<<endl;
    
    // Create a WCSimRootEvent to put stuff from the tree in and set the branch address for reading from the tree
    cout<<"making new WCSimRootEvent"<<endl;
    WCSimRootEvent *wcsimrootsuperevent = new WCSimRootEvent();
    cout<<"setting branch address"<<endl;
    wcsimT->SetBranchAddress("wcsimrootevent",&wcsimrootsuperevent);
    
    // Force deletion to prevent memory leak when issuing multiple calls to GetEvent()
    cout<<"setting autodelete"<<endl;
    wcsimT->GetBranch("wcsimrootevent")->SetAutoDelete(kTRUE);
    
    // Print the first event from the modified WCSim version
    cout<<"getting event 0"<<endl;
    wcsimT->GetEvent(0);
    WCSimRootTrigger *wcsimrootevent = wcsimrootsuperevent->GetTrigger(0);
    cout << "Stats for the first event in file " << filename << endl;
    cout << "Number of tube hits " << wcsimrootevent->GetNumTubesHit() << endl;
    cout << "Number of digitized tube hits " << wcsimrootevent->GetNumDigiTubesHit() << endl;
    cout << "Number of photoelectron hit times " << wcsimrootevent->GetCherenkovHitTimes()->GetEntries() << endl;
    
    numtubeshit = wcsimrootevent->GetNumTubesHit();
    numdigits = wcsimrootevent->GetNumDigiTubesHit();
    numpes = wcsimrootevent->GetCherenkovHitTimes()->GetEntries();
    numevents = nevent;
    
    cout<<"filling histograms"<<endl;
    // Now loop over events and fill the histograms
    for (int ev=0; ev<nevent; ev++){
      // Read the event from the tree into the WCSimRootEvent instance
      wcsimT->GetEvent(ev);
      wcsimrootevent = wcsimrootsuperevent->GetTrigger(0);
      if(verbose){
        printf("********************************************************");
        printf("Evt, date %d %d\n", wcsimrootevent->GetHeader()->GetEvtNum(),
               wcsimrootevent->GetHeader()->GetDate());
        printf("Mode %d\n", wcsimrootevent->GetMode());
        printf("Number of subevents %d\n",
               wcsimrootsuperevent->GetNumberOfSubEvents());
        
        printf("Vtxvol %d\n", wcsimrootevent->GetVtxvol());
        printf("Vtx %f %f %f\n", wcsimrootevent->GetVtx(0),
               wcsimrootevent->GetVtx(1),wcsimrootevent->GetVtx(2));
      }
      
      for (int index = 0 ; index < wcsimrootsuperevent->GetNumberOfEvents(); index++){
          wcsimrootevent = wcsimrootsuperevent->GetTrigger(index);
          int ncherenkovdigihits = wcsimrootevent->GetNcherenkovdigihits();
          hits->Fill(ncherenkovdigihits);
          
          float totalq = 0.;
          float totalt = 0.;
          // Loop through elements in the TClonesArray of WCSimRootCherenkovHits
          for (int i=0; i< ncherenkovdigihits; i++){
              TObject *Digi = (wcsimrootevent->GetCherenkovDigiHits())->At(i);
              WCSimRootCherenkovDigiHit *wcsimrootcherenkovdigihit = 
                dynamic_cast<WCSimRootCherenkovDigiHit*>(Digi);
              
              float q = wcsimrootcherenkovdigihit->GetQ();
              float t = wcsimrootcherenkovdigihit->GetT();
              totalq+=q;
              totalt+=t;
          }
          float av_time = (ncherenkovdigihits > 0) ? totalt/ncherenkovdigihits : 0;
          float av_q = (ncherenkovdigihits > 0) ? totalq/ncherenkovdigihits : 0;
          charge->Fill(av_q);  
          time->Fill(av_time);
      }
      
      // reinitialize super event between loops.
      wcsimrootsuperevent->ReInitialize();
    }// End of loop over events
    
    f->Close();
    delete f;
    f=0;
    
  //}
  
  cout<<"returning"<<endl;
  cout << "***********************************************************" << endl;
  return 1;
}
