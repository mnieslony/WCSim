{
  std::string comparisondirectoryname="annie"; // this is the directory name of the folder in validationscriptpath.
  
  
  cout<<"getting WCSIMDIR"<<endl;
  std::string wcsimdirenv = gSystem->Getenv("WCSIMDIR");
  std::string validationscriptpath;
  if(wcsimdirenv !=  NULL){
    validationscriptpath = wcsimdirenv+"/verification-test-scripts/";
  } else {
    cout<<"WCSIMDIR not defined: please define it."<<endl;
    exit(-1);
  }
  
  cout<<"defining pointers"<<endl;
  #include <unistd.h>
  std::string comparisonfilename="comparisonfile.root";
  TFile* comparisonfile=0;
  std::string file;
  TTree* tree=0;
  TBranch *numtubeshitb=0;
  TBranch *numdigitsb=0;
  TBranch *numpesb=0;
  TBranch *numeventsb=0;
  TBranch *hitsb=0;
  TBranch *timesb=0;
  TBranch *chargeb=0;
  
  cout<<"making variables"<<endl;
  int numtubeshit;
  int numdigits;
  int numpes;
  int numevents;
  TH1F* hits   = new TH1F("PMT Hits", "# Digitized Hits", 500, 0, 3000);
  TH1F* times   = new TH1F("Average Time", "Average Time", 600, 900, 2000);
  TH1F* charge = new TH1F("Q/# Digitized PMT", "Average Charge", 200, 0, 5);
  
  cout<<"searching for comparison file "<<(validationscriptpath+comparisonfilename).c_str()<<endl;
  if(access((validationscriptpath+comparisonfilename).c_str(), F_OK)!=-1){
    cout<<"comparison file found: retrieving it"<<endl;
    comparisonfile = TFile::Open((validationscriptpath+comparisonfilename).c_str(),"UPDATE");
    tree=(TTree*)comparisonfile->Get("comparisontree");
    
    cout<<"setting branch addresses for first event stats"<<endl;
    tree->SetBranchAddress("numtubeshit", &numtubeshit, &numtubeshitb);
    tree->SetBranchAddress("numdigits", &numdigits, &numdigitsb);
    tree->SetBranchAddress("numpes", &numpes, &numpesb);
    tree->SetBranchAddress("numevents", &numevents, &numeventsb);
    
    cout<<"setting branch addresses for histograms"<<endl;
    tree->SetBranchAddress("hits_hist", &hits, &hitsb);
    tree->SetBranchAddress("time_hist", &times, &timesb);
    tree->SetBranchAddress("charge_hist", &charge, &chargeb);
    
    file=comparisondirectoryname;
    
  } else {
    cout<<"comparison file not found: making it"<<endl;
    comparisonfile = new TFile((validationscriptpath+comparisonfilename).c_str(),"RECREATE","File comparing upstream with another");
    tree = new TTree("comparisontree","Comparison Tree");
    
    cout<<"creating branches for first event stats"<<endl;
    numtubeshitb = tree->Branch("numtubeshit", &numtubeshit);
    numdigitsb = tree->Branch("numdigits", &numdigits);
    numpesb = tree->Branch("numpes", &numpes);
    numeventsb = tree->Branch("numevents", &numevents);
    
    cout<<"creating branches for histograms"<<endl;
    hitsb   = tree->Branch("hits_hist", &hits);
    timesb   = tree->Branch("time_hist", &times);
    chargeb = tree->Branch("charge_hist", &charge);
    
    file="upstream";
  }
  
  //for(int file=0; file<2; file++){
    std::string filename, directoryname;
    if(file.compare("upstream")==0){
      cout<<"reading stats from upstream file"<<endl;
      directoryname="upstream/";
      filename=validationscriptpath+directoryname+"wcsimtest.root";
    } else {
      cout<<"reading stats from alternate file"<<endl;
      directoryname=file+"/";
      filename=validationscriptpath+directoryname+"wcsimtest.root";
    }
    
    cout<<"including header files before reading file "<<filename<<endl;
    std::string theincludepath=validationscriptpath+directoryname;
    gInterpreter->AddIncludePath(theincludepath.c_str());
    std::vector<std::string> includefiles { "WCSimEnumerations.hh", "WCSimRootEvent.hh", "WCSimRootGeom.hh", "WCSimPmtInfo.hh", "WCSimRootLinkDef.hh", "WCSimRootOptions.hh"};
    if(file.compare("annie")==0) includefiles.pop_back(); // annie doesn't yet use WCSimRootOptions.hh
    for(auto includefile : includefiles){
      std::string includecommand="#include \""+theincludepath+includefile+"\"";
      cout<<"processing \'"<<includecommand<<"\'"<<endl;
      gROOT->ProcessLine(includecommand.c_str());
    }
    // must load headers before the library
    std::string librarycommand=validationscriptpath+directoryname+"libWCSimRoot.so";
    cout<<"loading library \'"<<librarycommand<<"\'"<<endl;
    gSystem->Load(librarycommand.c_str());
    
    gROOT->ProcessLine(".L /annie/app/users/moflaher/wcsim/wcsim/verification-test-scripts/verification_HitsChargeTime.C");
    verification_HitsChargeTime(filename, directoryname, numtubeshit, numdigits, numpes, numevents, hits, times, charge);
    
    numtubeshitb->Fill();
    numdigitsb->Fill();
    numpesb->Fill();
    numeventsb->Fill();
    hitsb->Fill();
    timesb->Fill();
    chargeb->Fill();
    
    //assert(false);
  //}
  
  if(numtubeshitb->GetEntries()==2){
  
    tree->SetEntries(2);
    comparisonfile->cd();
    tree->Write("",TObject::kOverwrite);
    
    cout<<"comparing histograms"<<endl;
    // pull both events into memory for comparison
    if(file.compare("upstream")==0){
      tree->GetEntry(1);
    } else {
      tree->GetEntry(0);
    }
    int numtubeshit1=numtubeshit;
    int numdigits1=numdigits;
    int numpes1=numpes;
    int numevents1=numevents;
    TH1F* hits1 = new TH1F(*hits);
    TH1F* times1 = new TH1F(*times);
    TH1F* charge1 = new TH1F(*charge);
    if(file.compare("upstream")==0){
      tree->GetEntry(0);
    } else {
      tree->GetEntry(1);
    }
    int numtubeshit2=numtubeshit;
    int numdigits2=numdigits;
    int numpes2=numpes;
    int numevents2=numevents;
    TH1F* hits2 = new TH1F(*hits);
    TH1F* times2 = new TH1F(*times);
    TH1F* charge2 = new TH1F(*charge);
    
    // Compare these two events
    //cout <<  "***********************************************************" << endl;
    if (abs(numtubeshit1 - numtubeshit2)>1.0e-6){cout << "FIRST EVENT TEST FAILED: Number of hit tubes do not match" << endl;}
    else {cout << "FIRST EVENT TEST PASSED: Number of hit tubes matches" << endl;}
    if (abs(numdigits1 - numdigits2)>1.0e-6){cout << "FIRST EVENT TEST FAILED: Number of digitized tubes do not match" << endl; }
    else {cout << "FIRST EVENT TEST PASSED: Number of digitized tubes matches" << endl; }
    if (abs(numpes1 - numpes2)> 1.0e-6){cout << "FIRST EVENT TEST FAILED: Number of hit times do not match" << endl;}
    else {cout << "FIRST EVENT TEST PASSED: Number of hit times matches" << endl;}

    if (numevents1 != numevents2) {
      cout <<  "***********************************************************" << endl;
      cout << "The input files don’t contain the same number of events. Only the first events were compared. To see histograms of the number of hits, deposited charge and hit time, please choose two input files which contain the same number of events." << endl;
      return -1;
    }
    
    Double_t ks_hits = hits1->KolmogorovTest(hits2);
    Double_t ks_charge = charge1->KolmogorovTest(charge2);
    Double_t ks_times = times1->KolmogorovTest(times2);
    cout << "***********************************************************" << endl;
    cout << "ks test for # of digitized hits: " << ks_hits << endl;
    cout << "ks test for average charge: " << ks_charge << endl;
    cout << "ks test for average time: " << ks_times << endl;

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
    leg->AddEntry(hits1,"upstream", "l");
    leg->AddEntry(hits2,"comparison", "l");
    leg->Draw();
    
    c1->cd(2);
    charge1->GetXaxis()->SetTitle("Total Charge / # digitized hits");
    charge1->Draw();
    charge2->SetLineColor(kRed);
    c1->cd(2);
    charge2->Draw("SAME");
    c1->cd(3);
    times1->GetXaxis()->SetTitle("Total Time / # digitized hits (ns)");  
    times1->Draw();
    times2->SetLineColor(kRed);
    c1->cd(3);
    times2->Draw("SAME");
    
    comparisonfile->cd();
    c1->cd();
    cd->SaveAs("comparisonplots.png");
    hits1->Write();
    hits2->Write();
    charge1->Write();
    charge2->Write();
    time1->Write();
    time2->Write();
    
  } else {
    
    tree->SetEntries(1);
    comparisonfile->cd();
    tree->Write("",TObject::kOverwrite);
    comparisonfile->Close();
    gApplication->Terminate();
    
  }
  
}
