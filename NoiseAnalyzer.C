#define NoiseAnalyzer_cxx
#include "NoiseAnalyzer.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>

#include <iostream>
#include <fstream>

using namespace std;

const Int_t maxChannels = 8300;
const Int_t maxChannelsElec = 8640;
const Int_t maxTicks = 9594;

const Bool_t saveWaveforms = false;
const Bool_t runChirpAlg = false;
const Bool_t runZigzagAlg = false;
const Bool_t runWaveAlg = false;
const Int_t maxEvents = 1;
//const Bool_t saveWaveforms = true;
//const Bool_t runChirpAlg = true;
//const Bool_t runZigzagAlg = true;
//const Bool_t runWaveAlg = true;

void NoiseAnalyzer::Loop()
{
  if (fChain == 0) return;

  ifstream infile;
  infile.open("chanInfo.txt");

  avgMeanHist = new TH1F("avgMeanHist",";Channel;ADC Value Mean (Averaged Over Events)", maxChannels,-0.5,maxChannels-0.5);
  avgRMSHist = new TH1F("avgRMSHist", ";Channel;ADC Value RMS (Averaged Over Events)", maxChannels,-0.5,maxChannels-0.5);
  avgMaxHist = new TH1F("avgMaxHist", ";Channel;ADC Value Max-Mean (Averaged Over Events)", maxChannels,-0.5,maxChannels-0.5);
  chirpFracHist = new TH1F("chirpFracHist", ";Channel;Chirpiness (Averaged Over Events)", maxChannels,-0.5,maxChannels-0.5);
  zigzagFracHist = new TH1F("zigzagFracHist", ";Channel;Zigzagginess (Averaged Over Events)", maxChannels,-0.5,maxChannels-0.5);
  waveFracHist = new TH1F("waveFracHist", ";Channel;Waviness (Averaged Over Events)", maxChannels,-0.5,maxChannels-0.5);

  avgMeanElecHist = new TH1F("avgMeanElecHist",";ElecChannel;ADC Value Mean (Averaged Over Events)", maxChannelsElec,-0.5,maxChannelsElec-0.5);
  avgRMSElecHist = new TH1F("avgRMSElecHist", ";ElecChannel;ADC Value RMS (Averaged Over Events)", maxChannelsElec,-0.5,maxChannelsElec-0.5);
  avgMaxElecHist = new TH1F("avgMaxElecHist", ";ElecChannel;ADC Value Max-Mean (Averaged Over Events)", maxChannelsElec,-0.5,maxChannelsElec-0.5);
  chirpFracElecHist = new TH1F("chirpFracElecHist", ";ElecChannel;Chirpiness (Averaged Over Events)", maxChannelsElec,-0.5,maxChannelsElec-0.5);
  zigzagFracElecHist = new TH1F("zigzagFracElecHist", ";ElecChannel;Zigzagginess (Averaged Over Events)", maxChannelsElec,-0.5,maxChannelsElec-0.5);
  waveFracElecHist = new TH1F("waveFracElecHist", ";ElecChannel;Waviness (Averaged Over Events)", maxChannelsElec,-0.5,maxChannelsElec-0.5);

  std::string string_crate;
  std::string string_slot;
  std::string string_femch;
  std::string string_larch;
  std::string string_larwire;
  std::string string_wirenum;
  std::string string_plane;
  std::string string_length;

  Int_t varCrate[10000];
  Int_t varSlot[10000];
  Int_t varFemch[10000];
  Int_t varLarch[10000];
  Int_t varLarwire[10000];
  Int_t varWirenum[10000];
  Int_t varPlane[10000];
  Double_t varLength[10000];

  Int_t numChan = 0;
  while(infile >> string_crate)
  {
    numChan++;

    infile >> string_slot;
    infile >> string_femch;
    infile >> string_larch;
    infile >> string_larwire;
    infile >> string_wirenum;
    infile >> string_plane;
    infile >> string_length;

    varCrate[numChan-1] = atoi(string_crate.c_str());
    varSlot[numChan-1] = atoi(string_slot.c_str());
    varFemch[numChan-1] = atoi(string_femch.c_str());
    varLarch[numChan-1] = atoi(string_larch.c_str());
    varLarwire[numChan-1] = atoi(string_larwire.c_str());
    varWirenum[numChan-1] = atoi(string_wirenum.c_str());
    varPlane[numChan-1] = atoi(string_plane.c_str());
    varLength[numChan-1] = atof(string_length.c_str());
  }
  infile.close();

  TFile* outputFile = new TFile(Form("output_run%d.root",fRunNumber),"RECREATE");

  TH1F** meanHist = new TH1F * [maxEvents];
  TH1F** rmsHist = new TH1F * [maxEvents];
  TH1F** maxHist = new TH1F * [maxEvents];
  TH1F** chirpHist = new TH1F * [maxEvents];
  TH1F** zigzagHist = new TH1F * [maxEvents];
  TH1F** waveHist = new TH1F * [maxEvents];
  TH1F** meanElecHist = new TH1F * [maxEvents];
  TH1F** rmsElecHist = new TH1F * [maxEvents];
  TH1F** maxElecHist = new TH1F * [maxEvents];
  TH1F** chirpElecHist = new TH1F * [maxEvents];
  TH1F** zigzagElecHist = new TH1F * [maxEvents];
  TH1F** waveElecHist = new TH1F * [maxEvents];
  for(Int_t i = 0; i < maxEvents; i++)
  {
    meanHist[i] = new TH1F(Form("meanHist_event%d",i),";Channel;ADC Value Mean", maxChannels,-0.5,maxChannels-0.5);
    rmsHist[i] = new TH1F(Form("rmsHist_event%d",i), ";Channel;ADC Value RMS", maxChannels,-0.5,maxChannels-0.5);
    maxHist[i] = new TH1F(Form("maxHist_event%d",i), ";Channel;ADC Value Max-Mean", maxChannels,-0.5,maxChannels-0.5);
    chirpHist[i] = new TH1F(Form("chirpHist_event%d",i), ";Channel;Chirpiness", maxChannels,-0.5,maxChannels-0.5);
    zigzagHist[i] = new TH1F(Form("zigzagHist_event%d",i), ";Channel;Zigzagginess", maxChannels,-0.5,maxChannels-0.5);
    waveHist[i] = new TH1F(Form("waveHist_event%d",i), ";Channel;Waviness", maxChannels,-0.5,maxChannels-0.5);

    meanElecHist[i] = new TH1F(Form("meanElecHist_event%d",i),";ElecChannel;ADC Value Mean", maxChannelsElec,-0.5,maxChannelsElec-0.5);
    rmsElecHist[i] = new TH1F(Form("rmsElecHist_event%d",i), ";ElecChannel;ADC Value RMS", maxChannelsElec,-0.5,maxChannelsElec-0.5);
    maxElecHist[i] = new TH1F(Form("maxElecHist_event%d",i), ";ElecChannel;ADC Value Max-Mean", maxChannelsElec,-0.5,maxChannelsElec-0.5);
    chirpElecHist[i] = new TH1F(Form("chirpElecHist_event%d",i), ";ElecChannel;Chirpiness", maxChannelsElec,-0.5,maxChannelsElec-0.5);
    zigzagElecHist[i] = new TH1F(Form("zigzagElecHist_event%d",i), ";ElecChannel;Zigzagginess", maxChannelsElec,-0.5,maxChannelsElec-0.5);
    waveElecHist[i] = new TH1F(Form("waveElecHist_event%d",i), ";ElecChannel;Waviness", maxChannelsElec,-0.5,maxChannelsElec-0.5);
  }

  TH1F** noiseHist_event0 = new TH1F * [maxChannels];
  for(Int_t i = 0; i < maxChannels; i++)
  {
    noiseHist_event0[i] = new TH1F(Form("noiseHist_event0_channel%d",i),";Tick;ADC Value", maxTicks,-0.5,maxTicks-0.5);
  }

  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  Int_t channelNum;
  Int_t channelNumElec;
  int FEMslotOffset[9] = {4,0,0,0,0,0,0,0,1};
  for (Long64_t jentry=0; jentry<nentries;jentry++)
  {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);
    nbytes += nb;

    channelNum = 15*64*_crate+64*_fem+_chan;
    channelNumElec = 15*64*_crate+64*(14-_fem-FEMslotOffset[_crate])+_chan;

    //if((channelNum > 4190) && (channelNum < 4415))
    //{
    //  cout << _crate+1 << " " << _fem+4 << " " << _chan << endl;
    //}

    if(varLarch[channelNum] != -1)
    {
      Double_t maxVal = GetMaxVal();

      meanHist[_event]->Fill(varLarch[channelNum],_mean);
      rmsHist[_event]->Fill(varLarch[channelNum],_rms);
      maxHist[_event]->Fill(varLarch[channelNum],maxVal-_mean);
      meanElecHist[_event]->Fill(channelNumElec,_mean);
      rmsElecHist[_event]->Fill(channelNumElec,_rms);
      maxElecHist[_event]->Fill(channelNumElec,maxVal-_mean);

      avgMeanHist->Fill(varLarch[channelNum],_mean);
      avgRMSHist->Fill(varLarch[channelNum],_rms);
      avgMaxHist->Fill(varLarch[channelNum],maxVal-_mean);
      avgMeanElecHist->Fill(channelNumElec,_mean);
      avgRMSElecHist->Fill(channelNumElec,_rms);
      avgMaxElecHist->Fill(channelNumElec,maxVal-_mean);
      
      if(saveWaveforms == true)
      {
        for(Int_t i = 0; i < adc_v->size(); i++)
	{
          if(_event == 0)
	  {
            noiseHist_event0[varLarch[channelNum]]->Fill(i,(*adc_v)[i]);
	  }
	}
      }

      if(runChirpAlg == true)
      {
        Double_t chirpVal = ChirpAlg();

        if(chirpVal > 0.0)
	{
          chirpHist[_event]->Fill(varLarch[channelNum],chirpVal);
          chirpElecHist[_event]->Fill(channelNumElec,chirpVal);

          chirpFracHist->Fill(varLarch[channelNum],chirpVal);
          chirpFracElecHist->Fill(channelNumElec,chirpVal);
	}
      }

      if(runZigzagAlg == true)
      {
        Double_t zigzagVal = ZigzagAlg();

        if(zigzagVal > 0.0)
	{
          zigzagHist[_event]->Fill(varLarch[channelNum],zigzagVal);
          zigzagElecHist[_event]->Fill(channelNumElec,zigzagVal);

          zigzagFracHist->Fill(varLarch[channelNum],zigzagVal);
          zigzagFracElecHist->Fill(channelNumElec,zigzagVal);
	}

        Double_t waveVal = WaveAlg();

        if(waveVal > 0.0)
	{
          waveHist[_event]->Fill(varLarch[channelNum],waveVal);
          waveElecHist[_event]->Fill(channelNumElec,waveVal);

          waveFracHist->Fill(varLarch[channelNum],waveVal);
          waveFracElecHist->Fill(channelNumElec,waveVal);
	}
      }
    }
  }
  avgMeanHist->Scale(1.0/((Double_t)maxEvents));
  avgRMSHist->Scale(1.0/((Double_t)maxEvents));
  avgMaxHist->Scale(1.0/((Double_t)maxEvents));
  chirpFracHist->Scale(1.0/((Double_t)maxEvents));
  zigzagFracHist->Scale(1.0/((Double_t)maxEvents));
  waveFracHist->Scale(1.0/((Double_t)maxEvents));
  avgMeanElecHist->Scale(1.0/((Double_t)maxEvents));
  avgRMSElecHist->Scale(1.0/((Double_t)maxEvents));
  avgMaxElecHist->Scale(1.0/((Double_t)maxEvents));
  chirpFracElecHist->Scale(1.0/((Double_t)maxEvents));
  zigzagFracElecHist->Scale(1.0/((Double_t)maxEvents));
  waveFracElecHist->Scale(1.0/((Double_t)maxEvents));

  gStyle->SetOptStat(0);

  for(Int_t i = 0; i < maxEvents; i++)
  {
    meanHist[i]->SetOption("HIST");
    rmsHist[i]->SetOption("HIST");
    maxHist[i]->SetOption("HIST");
    meanHist[i]->Write();
    rmsHist[i]->Write();
    maxHist[i]->Write();

    meanElecHist[i]->SetOption("HIST");
    rmsElecHist[i]->SetOption("HIST");
    maxElecHist[i]->SetOption("HIST");
    meanElecHist[i]->Write();
    rmsElecHist[i]->Write();
    maxElecHist[i]->Write();
    
    if(runChirpAlg == true)
    {
      chirpHist[i]->SetOption("HIST");
      chirpHist[i]->Write();

      chirpElecHist[i]->SetOption("HIST");
      chirpElecHist[i]->Write();
    }

    if(runZigzagAlg == true)
    {
      zigzagHist[i]->SetOption("HIST");
      zigzagHist[i]->Write();

      waveHist[i]->SetOption("HIST");
      waveHist[i]->Write();

      zigzagElecHist[i]->SetOption("HIST");
      zigzagElecHist[i]->Write();

      waveElecHist[i]->SetOption("HIST");
      waveElecHist[i]->Write();
    }
  }

  if(saveWaveforms == true)
  {
    for(Int_t i = 0; i < maxChannels; i++)
    {
      noiseHist_event0[i]->SetOption("HIST");
      noiseHist_event0[i]->Write();
    }
  }

  outputFile->Close();
}

Double_t NoiseAnalyzer::GetMaxVal()
{
  Double_t ADCval;
  Double_t maxADC = -1.0;
  for(Int_t i = 0; i < adc_v->size(); i++)
  {
    ADCval = (*adc_v)[i];

    if(ADCval > maxADC)
      maxADC = ADCval;
  }

  return maxADC;
}

Double_t NoiseAnalyzer::ChirpAlg()
{
  const Int_t windowSize = 20;
  const Double_t chirpAmp = 0.0;
  const Double_t chirpMinRMS = 0.50;
  const Double_t chirpMaxRMS = 4.0;
  const Double_t chirpRatioMaxMin = 0.0;
  const Double_t chirpMinRMS_forNumLowCalc = 1.0;

  Int_t counter = 0;
  Double_t ADCval;
  Double_t meanAmp = 0.0;
  Double_t runningAmpMean = 0.0;
  Double_t runningAmpRMS = 0.0;
  Double_t maxAmp = -1.0;
  Double_t maxRMS = -1.0;
  Double_t minRMS = 9999999.0;
  Int_t minRMSBin = -1;
  Int_t maxRMSBin = -1;
  Int_t numLowRMS = 0;
  Int_t firstLowRMSBin = -1;
  Int_t lastLowRMSBin = -1;
  Bool_t lowRMSFlag = false;
  for(Int_t i = 0; i < adc_v->size(); i++)
  {
    ADCval = (*adc_v)[i];
    meanAmp += ADCval;
    runningAmpMean += ADCval;
    runningAmpRMS += TMath::Power(ADCval,2.0);

    if(ADCval > maxAmp)
      maxAmp = ADCval;

    counter++;
    if(counter == windowSize)
    {
      runningAmpMean /= (Double_t)windowSize;
      runningAmpRMS /= (Double_t)windowSize;
      runningAmpRMS = TMath::Sqrt(runningAmpRMS-TMath::Power(runningAmpMean,2.0));

      if(runningAmpRMS > maxRMS)
      {
        maxRMS = runningAmpRMS;
        maxRMSBin = i-windowSize+1;
      }

      if(runningAmpRMS < minRMS)
      {
        minRMS = runningAmpRMS;
        minRMSBin = i-windowSize+1;
      }

      if(runningAmpRMS < chirpMinRMS_forNumLowCalc)
      {
        numLowRMS++;
        if(lowRMSFlag == false)
 	{
          lowRMSFlag = true;
          firstLowRMSBin = i-windowSize+1;
          lastLowRMSBin = i-windowSize+1;
	}
	else
	{
          lastLowRMSBin = i-windowSize+1;
	}
      }

      counter = 0;
      runningAmpMean = 0.0;
      runningAmpRMS = 0.0;
    }
  }
  meanAmp /= (Double_t)adc_v->size();

  Double_t ratioMaxMin = maxRMS/minRMS;

  if((maxAmp-meanAmp > chirpAmp) && (minRMS < chirpMinRMS) && (maxRMS > chirpMaxRMS) && (ratioMaxMin > chirpRatioMaxMin))
  {
    //cout << minRMS << " " << maxRMS << endl;
    //cout << "MIN RMS START BIN:  " << minRMSBin << endl;
    //cout << "MAX RMS START BIN:  " << maxRMSBin << endl;
    //cout << "NUM LOW RMS SECTIONS:  " << numLowRMS << endl;
    //cout << "FIRST LOW RMS BIN:  " << firstLowRMSBin << endl;
    //cout << "LAST LOW RMS BIN:  " << lastLowRMSBin << endl;

    return ((Double_t) numLowRMS)/(((Double_t) maxTicks)/((Double_t) windowSize));
  }
  else
    return 0.0;
}

Double_t NoiseAnalyzer::ZigzagAlg()
{
  const Int_t startSidebandBin = 3700;
  const Int_t endSidebandBin = 4000;
  const Int_t startPeakBin = 4000;
  const Int_t endPeakBin = 4744;

  Double_t sidebandAvg;
  Double_t peakAvg;

  TH1F *currentHist = new TH1F("","",maxTicks,-0.5,maxTicks-0.5);
  for(Int_t i = 0; i < adc_v->size(); i++)
  {
    currentHist->SetBinContent(i+1,(*adc_v)[i]);
  }

  TH1F *currentFFTHist = new TH1F("","",maxTicks,-0.5,maxTicks-0.5);
  currentFFTHist = (TH1F*)currentHist->FFT(currentFFTHist,"MAG");

  Double_t integral = currentFFTHist->Integral(startPeakBin,startPeakBin+50);
  Double_t maxIntegral = integral;
  for(int i = startPeakBin+1; i <= endPeakBin; i++)
  {
    integral -= currentFFTHist->GetBinContent(i-1);
    integral += currentFFTHist->GetBinContent(i+50);

    if(integral > maxIntegral)
      maxIntegral = integral;
  }
  peakAvg = maxIntegral/50.0;
  sidebandAvg = TMath::Max(0.0001,currentFFTHist->Integral(startSidebandBin,endSidebandBin)/((Double_t) endSidebandBin-startSidebandBin+1));

  delete currentHist;
  delete currentFFTHist;

  if(peakAvg/sidebandAvg > 0.0)
    return peakAvg/sidebandAvg;
  else
    return 0.0;
}

Double_t NoiseAnalyzer::WaveAlg()
{
  const Int_t startSidebandBin = 300;
  const Int_t endSidebandBin = 500;
  const Int_t startPeakBin = 2;
  const Int_t endPeakBin = 248;

  Double_t sidebandAvg;
  Double_t peakAvg;

  TH1F *currentHist = new TH1F("","",maxTicks,-0.5,maxTicks-0.5);
  for(Int_t i = 0; i < adc_v->size(); i++)
  {
    currentHist->SetBinContent(i+1,(*adc_v)[i]);
  }

  TH1F *currentFFTHist = new TH1F("","",maxTicks,-0.5,maxTicks-0.5);
  currentFFTHist = (TH1F*)currentHist->FFT(currentFFTHist,"MAG");

  Double_t integral = currentFFTHist->Integral(startPeakBin,startPeakBin+50);
  Double_t maxIntegral = integral;
  for(int i = startPeakBin+1; i <= endPeakBin; i++)
  {
    integral -= currentFFTHist->GetBinContent(i-1);
    integral += currentFFTHist->GetBinContent(i+50);

    if(integral > maxIntegral)
      maxIntegral = integral;
  }
  peakAvg = maxIntegral/50.0;
  sidebandAvg = TMath::Max(0.0001,currentFFTHist->Integral(startSidebandBin,endSidebandBin)/((Double_t) endSidebandBin-startSidebandBin+1));

  delete currentHist;
  delete currentFFTHist;

  if(peakAvg/sidebandAvg > 0.0)
    return peakAvg/sidebandAvg;
  else
    return 0.0;
}

vector<TH1F*> NoiseAnalyzer::GetHists()
{
  vector<TH1F*> histVec;

  histVec.push_back((TH1F*)avgMeanHist->Clone());
  histVec.push_back((TH1F*)avgRMSHist->Clone());
  histVec.push_back((TH1F*)avgMaxHist->Clone());
  histVec.push_back((TH1F*)chirpFracHist->Clone());
  histVec.push_back((TH1F*)zigzagFracHist->Clone());
  histVec.push_back((TH1F*)waveFracHist->Clone());

  histVec.push_back((TH1F*)avgMeanElecHist->Clone());
  histVec.push_back((TH1F*)avgRMSElecHist->Clone());
  histVec.push_back((TH1F*)avgMaxElecHist->Clone());
  histVec.push_back((TH1F*)chirpFracElecHist->Clone());
  histVec.push_back((TH1F*)zigzagFracElecHist->Clone());
  histVec.push_back((TH1F*)waveFracElecHist->Clone());

  return histVec;
}
