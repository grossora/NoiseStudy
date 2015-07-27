//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Jun 10 18:16:19 2015 by ROOT version 6.02/05
// from TTree _wf_tree/Waveform Tree
// found on file: run_253_tree.root
//////////////////////////////////////////////////////////

#ifndef NoiseAnalyzer_h
#define NoiseAnalyzer_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH2.h>

// Header file for the classes stored in the TTree if any.
#include <vector>
#include <iostream>

using namespace std;

class NoiseAnalyzer {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain
   Int_t           fRunNumber;

   TH1F           *avgMeanHist;
   TH1F           *avgRMSHist;
   TH1F           *avgMaxHist;
   TH1F           *chirpFracHist;
   TH1F           *zigzagFracHist;
   TH1F           *waveFracHist;

   TH1F           *avgMeanElecHist;
   TH1F           *avgRMSElecHist;
   TH1F           *avgMaxElecHist;
   TH1F           *chirpFracElecHist;
   TH1F           *zigzagFracElecHist;
   TH1F           *waveFracElecHist;

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           _event;
   Int_t           _crate;
   Int_t           _fem;
   Int_t           _chan;
   vector<short>   *adc_v;
   Double_t        _mean;
   Double_t        _rms;

   // List of branches
   TBranch        *b_event;   //!
   TBranch        *b_crate;   //!
   TBranch        *b_fem;   //!
   TBranch        *b_chan;   //!
   TBranch        *b_adc_v;   //!
   TBranch        *b_mean;   //!
   TBranch        *b_rms;   //!

   NoiseAnalyzer(Int_t runNumber = 253, TTree *tree=0);
   virtual ~NoiseAnalyzer();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);

   // New functions
   virtual Double_t   GetMaxVal();
   virtual Double_t   ChirpAlg();
   virtual Double_t   ZigzagAlg();
   virtual Double_t   WaveAlg();
   virtual vector<TH1F*>  GetHists();
};

#endif

#ifdef NoiseAnalyzer_cxx
NoiseAnalyzer::NoiseAnalyzer(Int_t runNumber, TTree *tree) : fChain(0) 
{
   fRunNumber = runNumber;

   TFile *f = new TFile(Form("data/run_%d_tree.root",runNumber)); 
   f->GetObject("_wf_tree",tree);
   if(f)
   {
     Init(tree);
   }
   else
     cout << "FILE DOES NOT EXIST." << endl;
}

NoiseAnalyzer::~NoiseAnalyzer()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t NoiseAnalyzer::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t NoiseAnalyzer::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void NoiseAnalyzer::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   adc_v = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("_event", &_event, &b_event);
   fChain->SetBranchAddress("_crate", &_crate, &b_crate);
   fChain->SetBranchAddress("_fem", &_fem, &b_fem);
   fChain->SetBranchAddress("_chan", &_chan, &b_chan);
   fChain->SetBranchAddress("adc_v", &adc_v, &b_adc_v);
   fChain->SetBranchAddress("_mean", &_mean, &b_mean);
   fChain->SetBranchAddress("_rms", &_rms, &b_rms);
   Notify();
}

Bool_t NoiseAnalyzer::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void NoiseAnalyzer::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t NoiseAnalyzer::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef NoiseAnalyzer_cxx
