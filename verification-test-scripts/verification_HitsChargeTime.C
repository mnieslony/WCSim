#include <iostream>
#include <TH1F.h>
#include <stdio.h>     
#include <stdlib.h>    

// Simple example of reading a generated Root file
int verification_HitsChargeTime(const char *filename1="wcsimtest.root", const char *filename2="../../WCSim_clean/verification-test-scripts/wcsimtest.root", bool verbose=false)
{
  
  // Load the library with class dictionary info
  // (create with "gmake shared")
  std::string wcsimdirenv = gSystem->Getenv("WCSIMDIR");
  std::string validationscriptpath;
  if(wcsimdirenv !=  NULL){
    cout<<"wcsimdirenv="<<wcsimdirenv<<endl;
    validationscriptpath = wcsimdirenv+"/verification-test-scripts/";
    cout<<"validationscriptpath="<<validationscriptpath<<endl;
  } else {
    cout<<"WCSIMDIR not defined: please define it."<<endl;
    exit(-1);
  }
//  if(wcsimdirenv !=  NULL){
//    gSystem->Load("${WCSIMDIR}/libWCSimRoot.so");
//  }else{
//    gSystem->Load("../libWCSimRoot.so");
//  }
  
  // vectors to store the stats for event comparisons
  std::vector<int> numtubeshit;
  std::vector<int> numdigits;
  std::vector<int> numpes;
  std::vector<int> numevents;
  
  // Histograms for the modified WCSim version
  TH1F *hits1 = new TH1F("PMT Hits", "# Digitized Hits", 500, 0, 3000);
  TH1F *time1 = new TH1F("Average Time", "Average Time", 600, 900, 2000);
  TH1F *charge1 = new TH1F("Q/# Digitized PMT", "Average Charge", 200, 0, 5);
  // ... and for the "clean" version
  TH1F *hits2 = new TH1F("PMT Hits 2", "Digitized Hits", 500, 0, 3000);
  TH1F *time2 = new TH1F("Average Time 2", "Average Time", 600, 900, 2000);
  TH1F *charge2 = new TH1F("Q/# Digitized PMT 2", "Average Charge", 200, 0, 5);
  
  for(int file=0; file<2; file++){
  
    // Clear global scope
    gROOT->Reset();
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
    
    //gInterpreter->AddIncludePath(theincludepath); ??
    std::string directoryname, filename;
    TH1F *hits, *time, *charge;
    if(file==0){
      directoryname="upstream/";
      filename=directoryname+filename1;
      hits=hits1;
      time=time1;
      charge=charge1;
    } else {
      directoryname="annie/";
      filename=directoryname+filename2;
      hits=hits2;
      time=time2;
      charge=charge2;
    }
    
    std::vector<std::string> includefiles { "WCSimRootEvent.hh", "WCSimRootGeom.hh", "WCSimPmtInfo.hh",
    "WCSimEnumerations.hh", "WCSimRootLinkDef.hh", "WCSimRootOptions.hh"};
    //if(file==1) includefiles.pop_back();
    for(auto includefile : includefiles){
      std::string includecommand="#include \""+validationscriptpath+directoryname+includefile+"\"";
      gROOT->ProcessLine(includecommand.c_str());
    }
    std::string librarycommand=validationscriptpath+directoryname+"libWCSimRoot.so";
    gSystem->Load(librarycommand.c_str());
    
    std::string fullfilepath = validationscriptpath+filename;
    TFile *f = new TFile(fullfilepath.c_str(),"read");
    if (!f->IsOpen()){
      cout << "Error, could not open input file: " << filename << endl;
      return -1;
    }
    
    TTree  *wcsimT = (TTree*)f->Get("wcsimT");
    int nevent = wcsimT->GetEntries();
    
    // Create a WCSimRootEvent to put stuff from the tree in and set the branch address for reading from the tree
    WCSimRootEvent *wcsimrootsuperevent = new WCSimRootEvent();
    wcsimT->SetBranchAddress("wcsimrootevent",&wcsimrootsuperevent);
    
    // Force deletion to prevent memory leak when issuing multiple calls to GetEvent()
    wcsimT->GetBranch("wcsimrootevent")->SetAutoDelete(kTRUE);
    
    // Print the first event from the modified WCSim version
    wcsimT->GetEvent(0);
    WCSimRootTrigger *wcsimrootevent = wcsimrootsuperevent->GetTrigger(0);
    cout << "Stats for the first event in file " << filename << endl;
    cout << "Number of tube hits " << wcsimrootevent->GetNumTubesHit() << endl;
    cout << "Number of digitized tube hits " << wcsimrootevent->GetNumDigiTubesHit() << endl;
    cout << "Number of photoelectron hit times " << wcsimrootevent->GetCherenkovHitTimes()->GetEntries() << endl;
    cout << "***********************************************************" << endl;
    
    numtubeshit.push_back(wcsimrootevent->GetNumTubesHit());
    numdigits.push_back(wcsimrootevent->GetNumDigiTubesHit());
    numpes.push_back(wcsimrootevent->GetCherenkovHitTimes()->GetEntries());
    numevents.push_back(nevent);
    
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
    
  }
  
  // Compare these two events
  cout <<  "***********************************************************" << endl;
  if (abs(numtubeshit.at(0) - numtubeshit.at(1))>1.0e-6){cout << "FIRST EVENT TEST FAILED: Number of hit tubes do not match" << endl;}
  else {cout << "FIRST EVENT TEST PASSED: Number of hit tubes matches" << endl;}
  if (abs(numdigits.at(0) - numdigits.at(1))>1.0e-6){cout << "FIRST EVENT TEST FAILED: Number of digitized tubes do not match" << endl; }
  else {cout << "FIRST EVENT TEST PASSED: Number of digitized tubes matches" << endl; }
  if (abs(numpes.at(0) - numpes.at(1))> 1.0e-6){cout << "FIRST EVENT TEST FAILED: Number of hit times do not match" << endl;}
  else {cout << "FIRST EVENT TEST PASSED: Number of hit times matches" << endl;}

  if (numevents.at(0) != numevents.at(1)) {
    cout <<  "***********************************************************" << endl;
    cout << "The input files don’t contain the same number of events. Only the first events were compared. To see histograms of the number of hits, deposited charge and hit time, please choose two input files which contain the same number of events." << endl;
    return -1;
  }
  
  Double_t ks_hits = hits1->KolmogorovTest(hits2);
  Double_t ks_charge = charge1->KolmogorovTest(charge2);
  Double_t ks_time = time1->KolmogorovTest(time2);
  cout << "***********************************************************" << endl;
  cout << "ks test for # of digitized hits: " << ks_hits << endl;
  cout << "ks test for average charge: " << ks_charge << endl;
  cout << "ks test for average time: " << ks_time << endl;

  //  TCanvas c1("c1"); 
  float win_scale = 0.75;
  int n_wide(2);
  int n_high(2);
  TCanvas* c1 = new TCanvas("c1", "Test Plots", 500*n_wide*win_scale, 500*n_high*win_scale);
  c1->Draw();
  c1->Divide(2,2);
  c1->cd(1); 
  hits2->SetLineColor(kRed);
  hits1->Draw();
  c1->cd(1); hits2->Draw("SAME");

 TLegend *leg = new TLegend(0.2,0.7,0.55,0.85, "");
 leg->SetFillColor(0);
 leg->SetBorderSize(0);
 leg->AddEntry(hits1,filename1, "l");
 leg->AddEntry(hits2,filename2, "l");
 leg->Draw();
 
 c1->cd(2);
 charge1->GetXaxis()->SetTitle("Total Charge / # digitized hits");
 charge1->Draw();
 charge2->SetLineColor(kRed);
 c1->cd(2);
 charge2->Draw("SAME");
 c1->cd(3);
 time1->GetXaxis()->SetTitle("Total Time / # digitized hits (ns)");  
 time1->Draw();
 time2->SetLineColor(kRed);
 c1->cd(3);
 time2->Draw("SAME");
  
 return 1;
}
