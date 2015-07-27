#include "NoiseAnalyzer.h"
#include "TH2.h"
#include "TFile.h"
#include "TStyle.h"

vector<TH1F*> doOneRun(Int_t runNum);

int main(int argc, char **argv)
{
  const Int_t maxChannels = 8300;
  const Int_t maxChannelsElec = 8640;
  //const Int_t minRunNum = 253;
  //const Int_t maxRunNum = 315;
  //const Int_t minRunNum = 369;
  //const Int_t maxRunNum = 392;
  //const Int_t minRunNum = 470;
  //const Int_t maxRunNum = 523;
  //const Int_t minRunNum = 532;
  //const Int_t maxRunNum = 593;
  //const Int_t minRunNum = 594;
  //const Int_t maxRunNum = 636;

  const Int_t minRunNum = 900;
  const Int_t maxRunNum = 900; 

  //const Int_t minRunNum = 713;
  //const Int_t maxRunNum = 771;

  TFile* outputFile2D = new TFile("noiseStudyHists.root","RECREATE");

  TH2F *meanHist2D = new TH2F("meanHist2D","ADC Value Mean (Averaged Over Events);Channel;Run #", maxChannels,-0.5,maxChannels-0.5,maxRunNum-minRunNum+1,minRunNum-0.5,maxRunNum+0.5);
  TH2F *rmsHist2D = new TH2F("rmsHist2D", "ADC Value RMS (Averaged Over Events);Channel;Run #", maxChannels,-0.5,maxChannels-0.5,maxRunNum-minRunNum+1,minRunNum-0.5,maxRunNum+0.5);
  TH2F *maxHist2D = new TH2F("maxHist2D", "ADC Value Max-Mean (Averaged Over Events);Channel;Run #", maxChannels,-0.5,maxChannels-0.5,maxRunNum-minRunNum+1,minRunNum-0.5,maxRunNum+0.5);
  TH2F *chirpHist2D = new TH2F("chirpHist2D", "Chirpiness (Averaged Over Events);Channel;Run #", maxChannels,-0.5,maxChannels-0.5,maxRunNum-minRunNum+1,minRunNum-0.5,maxRunNum+0.5);
  TH2F *zigzagHist2D = new TH2F("zigzagHist2D", "Zigzagginess (Averaged Over Events);Channel;Run #", maxChannels,-0.5,maxChannels-0.5,maxRunNum-minRunNum+1,minRunNum-0.5,maxRunNum+0.5);
  TH2F *waveHist2D = new TH2F("waveHist2D", "Waviness (Averaged Over Events);Channel;Run #", maxChannels,-0.5,maxChannels-0.5,maxRunNum-minRunNum+1,minRunNum-0.5,maxRunNum+0.5);

  TH2F *meanElecHist2D = new TH2F("meanElecHist2D","ADC Value Mean (Averaged Over Events);ElecChannel;Run #", maxChannelsElec,-0.5,maxChannelsElec-0.5,maxRunNum-minRunNum+1,minRunNum-0.5,maxRunNum+0.5);
  TH2F *rmsElecHist2D = new TH2F("rmsElecHist2D", "ADC Value RMS (Averaged Over Events);ElecChannel;Run #", maxChannelsElec,-0.5,maxChannelsElec-0.5,maxRunNum-minRunNum+1,minRunNum-0.5,maxRunNum+0.5);
  TH2F *maxElecHist2D = new TH2F("maxElecHist2D", "ADC Value Max-Mean (Averaged Over Events);ElecChannel;Run #", maxChannelsElec,-0.5,maxChannelsElec-0.5,maxRunNum-minRunNum+1,minRunNum-0.5,maxRunNum+0.5);
  TH2F *chirpElecHist2D = new TH2F("chirpElecHist2D", "Chirpiness (Averaged Over Events);ElecChannel;Run #", maxChannelsElec,-0.5,maxChannelsElec-0.5,maxRunNum-minRunNum+1,minRunNum-0.5,maxRunNum+0.5);
  TH2F *zigzagElecHist2D = new TH2F("zigzagElecHist2D", "Zigzagginess (Averaged Over Events);ElecChannel;Run #", maxChannelsElec,-0.5,maxChannelsElec-0.5,maxRunNum-minRunNum+1,minRunNum-0.5,maxRunNum+0.5);
  TH2F *waveElecHist2D = new TH2F("waveElecHist2D", "Waviness (Averaged Over Events);ElecChannel;Run #", maxChannelsElec,-0.5,maxChannelsElec-0.5,maxRunNum-minRunNum+1,minRunNum-0.5,maxRunNum+0.5);

  vector<TH1F*> histSet;
  for(Int_t i = minRunNum; i <= maxRunNum; i++)
  {
    cout << "Analyzing Run #" << i << endl;
    histSet = doOneRun(i);

    for(Int_t j = 0; j < maxChannels; j++)
    {
      meanHist2D->Fill(j,i,histSet[0]->GetBinContent(j+1));
      rmsHist2D->Fill(j,i,histSet[1]->GetBinContent(j+1));
      maxHist2D->Fill(j,i,histSet[2]->GetBinContent(j+1));
      chirpHist2D->Fill(j,i,histSet[3]->GetBinContent(j+1));
      zigzagHist2D->Fill(j,i,histSet[4]->GetBinContent(j+1));
      waveHist2D->Fill(j,i,histSet[5]->GetBinContent(j+1));

      meanElecHist2D->Fill(j,i,histSet[6]->GetBinContent(j+1));
      rmsElecHist2D->Fill(j,i,histSet[7]->GetBinContent(j+1));
      maxElecHist2D->Fill(j,i,histSet[8]->GetBinContent(j+1));
      chirpElecHist2D->Fill(j,i,histSet[9]->GetBinContent(j+1));
      zigzagElecHist2D->Fill(j,i,histSet[10]->GetBinContent(j+1));
      waveElecHist2D->Fill(j,i,histSet[11]->GetBinContent(j+1));
    }
  }

  outputFile2D->cd();

  meanHist2D->Write();
  rmsHist2D->Write();
  maxHist2D->Write();
  chirpHist2D->Write();
  zigzagHist2D->Write();
  waveHist2D->Write();

  meanElecHist2D->Write();
  rmsElecHist2D->Write();
  maxElecHist2D->Write();
  chirpElecHist2D->Write();
  zigzagElecHist2D->Write();
  waveElecHist2D->Write();

  outputFile2D->Close();

  return 0;
}

vector<TH1F*> doOneRun(Int_t runNum)
{
  NoiseAnalyzer *noiseAna = new NoiseAnalyzer(runNum);
  noiseAna->Loop(); 

  vector<TH1F*> hists = noiseAna->GetHists();
  
  delete noiseAna;
  return hists;
}
