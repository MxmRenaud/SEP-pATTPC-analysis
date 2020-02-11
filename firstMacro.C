/*
* ================================================================
*
*	Filename:    firstMacro.C
*	
*	Description: Program to read the .root files and make them comprehensible to a human mind. 
*
*	Created:     7/11/2019
*               
*	Author:      Maxime Renaud, NSH, UND - KUL
*	Co-author:   Shilun Jin, NSH, UND
*
* ================================================================
* 
* TODO look into remanence for plotted histograms so you can delete your objects
* NOTE: because I'm an uncivilised Barbarian, this is supposed to be invoked as a root macro, not a compiled program. 
* 
* 
* NOTE HOWTO:at the beginning of the firstMacro(), a few self-explanatory boolians
* can be activated/deactivated depending on what you want the program to do.
*/

#include <stdlib.h>
#include <stdio.h>
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TBranch.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TSpectrum.h"
#include "TF1.h"
#include "TMath.h"
#include "TGraph.h"
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <time.h>
#include "TStyle.h"
#include "TLine.h"
#include "TSpectrum.h"
#include <cstdio>
#include "Fit/FitResult.h"
#include "TFitResult.h"
#include "TVirtualFitter.h"
#include "TMatrixT.h"
#include "TLegend.h"
#include "TApplication.h"
#include <cstring>
#include <TClass.h>

#include "fourier-fcts.c"


using namespace std;

Int_t globalNumberErrors;
Int_t globalNumberOfPUevents;
Int_t globalNumberOfLongEvents;

Int_t whichAuxChanIsTAC(TH1I* histo){ //is the TAC in the auxillary channel 1 or 2 ?
    //not valid, some MCP signals have very low plateaus and the similar averages
//     bool boo = false;
//     TF1 *fa0 = new TF1("fa0", "pol0(0)",0,511);
//     histo->Fit(fa0,"0");
//     if (fa0->GetParameter(0) > 2000){boo = true;}
//     delete fa0;
//     if (boo == true) return 2;
//     else return 1;
    //locate the max aqnd it's position: if in the back and over 3000, then it's MCP not TAC
    bool boo = false;
    Int_t maxReached=0;
    Int_t savePosition;
    for (int i =0;i<512;i++){
        if (histo->GetBinContent(i) > maxReached){
            maxReached = histo->GetBinContent(i);
            savePosition = i;
        }
    }
    if (maxReached >= 3000 && savePosition > 480){return 2;}
    else if ((maxReached < 3000 && savePosition > 480) || (maxReached >= 3000 && savePosition < 480)){
        cout<<"\n\n================ WARNING ! ================\nRevise your TAC-or-MCP conditions !\n================ WARNING ! ================\n\n";
        globalNumberErrors++;
        return 2;
    }
    else {return 1;}
}


Int_t pileUp(TH1I* histo){
    TF1 *fa0 = new TF1("fa0", "pol0(0)",0,511);
    histo->Fit(fa0,"0Q");
    Float_t barem = fa0->GetParameter(0)*1.33;
//     cout<<"TAC test-value : "<<barem<<endl;
    delete fa0;
    Int_t nbrOfEvent = 0;
    for (int i=0; i<511;i++){
       if (histo->GetBinContent(i) <= barem && histo->GetBinContent(i+1) > barem){nbrOfEvent++;}
    }
//     cout<<nbrOfEvent<<endl;
    return nbrOfEvent;
}


Double_t* locateTACsignal(TH1F* histo_input, Int_t howMany){
    Float_t rAtiO = 0.4;
    Float_t wIdTh = 3;
    Int_t rangeMin = 2;
    Int_t rangeMax = 450;
    TSpectrum *s = new TSpectrum(howMany);
    histo_input->GetXaxis()->SetRangeUser(rangeMin,rangeMax);
    s->Search(histo_input, wIdTh, "R",rAtiO);
    
    Double_t* XPics = s->GetPositionX();
    Double_t v;
    Int_t i,j;
    //--Sorting by Decreasing order
    for (i=0;i<howMany;i++){ 
        for (j=i+1;j<howMany;j++){
            if (XPics[i]<XPics[j]){
                v=XPics[j];
                XPics[j]=XPics[i];
                XPics[i]=v;
            }
        }
    }
    return XPics;
}
   

Float_t findDerivative(TH1F* histo_input, Double_t startingPoint, Int_t numberOfBins){
    if (numberOfBins > 0 && numberOfBins < 511/*histo_input->GetMaximumStored()*/){
        Float_t deri;
//         histo_input->GetXaxis()->SetRangeUser((int)startingPoint,(int)startingPoint+numberOfBins);
        TF1 *fa1 = new TF1("fa1", "pol1(0)",(int)startingPoint,(int)startingPoint+numberOfBins);
        histo_input->Fit(fa1,"QR");
        deri = fa1->GetParameter(1);
        
        delete fa1;
        return deri;
    }
    else {
        cout<<"\nIn findDerivative, numberOfBins is "<<numberOfBins<<", and that shit ain't gonna fly fam'.\n\n";
        return 0.;
        globalNumberErrors++;
    }
}


// Function that return average of an array. 
float AveragE(Float_t a[], int n, int start){ 
    if (start-n < 0){cout<<"\nInput error for fct 'AveragE(Float_t a[], int n, int start)'.\n";return 0;globalNumberErrors++;}
    else{
        float sum = 0; 
        for (int i=start-n; i<start; i++){sum = sum + a[i];}
        return sum/n;
    }
} 

//functions for TF1's
double the_gausppar(double* var, double* param){
  return param[0]*TMath::Gaus(var[0],param[1],param[2]);
}
double the_steppar(double* var, double* param){
    return param[0]*((var[0]<param[1])?1:0);
}



//=================================================================================================================================


Int_t firstMacro(){

    
//END---------WHAT DO Y0U WISH TO DO HERE, MORTAL ? -------------------
    
    bool runOnTheCRC = false; //deactivates mandatory user input & adapt file paths NOTE not optimal for the latter part.
    bool multiRunAnalysis = true; //select between two file, 'list-preAnal.txt' and 'list-analysis.txt', for convenience in swtching between debugging/actual run
    
    bool fillIndividualChannels = false;
    bool realignTracks = true; //BEST TRACK LENGTH FINDING ALGORYTHM IS LOCKED IN HERE, as well
    bool acceptPeaksBelowBeam = false; // do you want to keep events that peak much further than the primary beam ? see 'pointOfPeakRejection' Will exclude most events that are just noise and some disintegration events.
    bool useFusionGates = false; //attempt to only keep events fulfilling certain conditions [hard-coded atm]
    bool substractBckgrd = true;
    bool beQuite = true; //reduced terminal output
    bool addSignalCharacOnPlots = true;
    bool deconvoluteThis = true;
    
    bool drawWhichPadIsIt = false;
    bool drawIndividualChannels = false;
    bool drawSummedChannels = false; //sums all channels of a single event into a track and draw.
    bool drawAuxillaryChannels = false;
    bool drawTBFeventRemanent = true; //drawTimeBucketFormatEventRemanent = draw a summed profile of all tracks for all events
    bool drawDerivatives = false; //only one event at a time
    bool drawLengthVpeakHeight = true;
    bool drawLengthVtrackCharge = true;
    bool drawLengthVbackgrd = true;
    bool drawLengthVtac = true;
    bool drawChargeVtac = true;
    bool drawPeakHeightVtac = true;
    bool drawChargeVpeakHeight = true;
    bool draw3DTrack = false; //WARNING please only one event at a time
    bool drawDoubleProjectTrack = false; //TODO. TODO. TODO TODO TODO TODOOOOOOOOOOOOOOOOOOOOOO !
    
    bool saveHistograms = true;
    
//END---THE CONTRACT HAS BEEN SEALED, MORTAL. I SHALL NOW PROCEED...-----
    
    
    
  const Int_t derivativeCalcWindowSize = 5;
  const Int_t estimatedWidth = 9.; //estimated width of the track's beginning rising edge at high TimeBicket
  const Int_t pointOfPeakRejection = 75;
  const Int_t EVENT = 0;//run_0090/!\ 11602;  // 14030->60; /!\ 14030 & 14041 & 14048 & 14053(noDipBelow); /!\ 14031 & 14032 & 14052 (1TAC2Peak) 14078 (1TAC3peaks); /!\ 4041 & 14056 (weirdMCP) /!\ 14030 & 14053 & 14059 (singleShort)
  char eRrOr;
  const Float_t meanBeamPeakHeight = 256.5;//EDIT value on 02/04, prev.: 253.5;
  const Float_t meanBeamCharge = 12095.5;//EDIT value on 02/04, prev.: 11259.5;
  const Float_t driftV = 20.; //mm/us
  const Float_t timeBucketSize = 0.160; //us
  const Float_t Li7ChargePercentage = 0.2537; //charge dep by Li7 as percentage Li8 beam charge
  const Float_t Li7peakHeight = 0.9998; //bragg peak height of Li7 as percentage Li8 bragg peak height
  const Float_t Li7avTrLenght = 57.5; //Li7 bragg peak position in Time Buckets


  if (fillIndividualChannels == false && drawIndividualChannels == true){cout<<"\n~~~~~~~~~~~~~~~~~~~~~WARNING !~~~~~~~~~~~~~~~~~~~\n\nY'need to fill the channels to draw them, asshat.\nGo activate 'fillIndividualChannels' and try again.\n\n\tCan you manage to get THAT right, at least ?\n\n";return 0;}
  if (realignTracks == true && drawTBFeventRemanent == true && runOnTheCRC == false){
      cout<<"\nWARNING ! You asked for both the tracks to be re-aligned before drawing (realignTracks) AND for the tracks to be drawn without being re-aligned (drawTBFeventRemanent).\nDo you wish to proceed ? (y/n) : ";
      cin>>eRrOr;
      if (eRrOr == 'n') {return 0;}
      else if (eRrOr == 'y') {cout<<"\nProceeding..."<<endl;}
      else {cout<<"\n... you think you're being funny ?\n\n\n\tYou. Aren't.\n\n";return 0;}
  }
  if (useFusionGates == true && (drawLengthVpeakHeight == false || drawLengthVtrackCharge == false || realignTracks == false)){
      cout<<"\nWARNING\nFor now  please activate drawLengthVtrackCharge, drawLengthVpeakHeight & realignTracks to run useFusionGates. will patch later.\n\n";return 0;}
  if (drawChargeVtac == true && drawLengthVtrackCharge == false){
      cout<<"\nWARNING ! 'drawChargeVtac == true && drawLengthVtrackCharge == false' is not a supported combo.\nPlease correct.\n\nThank you kindly.\n";return 0;}
  if (drawPeakHeightVtac == true && drawLengthVpeakHeight == false){
      cout<<"\nWARNING ! 'drawPeakHeightVtac == true && drawLengthVpeakHeight == false' ain't no good. Sorry.\n";return 0;}
  if (draw3DTrack == true && realignTracks == false){cout<<"\nWell met, stranger. I hereby inform you that the option 'draw3DTrack' requires option 'realignTracks' to be enabled, in order to be executed.\nThis message will now terminate your run.\n\nOur deepest apologies.\n";return 0;}




  //--CALIBRATION ELEMENTS--------
  
  //Go fetch the mapping matrix
  const int numberOfPads = 2021;
  float mapp[numberOfPads][7];
  const char *fMapStreamChar;
  string fStreamName;
  int line = 0;
  if (runOnTheCRC == true){fStreamName = "/afs/crc.nd.edu/user/m/mrenaud1/Public/SEP-pATTPC-analysis/TPC-SEP2019-fullMap.dat";}
  else{fStreamName = "TPC-SEP2019-fullMap.dat";}
  fMapStreamChar = fStreamName.c_str(); //WARNING http://www.cplusplus.com/forum/general/100714/
  fstream mapStream(fMapStreamChar);
  if(mapStream.is_open()){//structure Cobo<< "\t" << Asad<< "\t" << Aget<< "\t" << Channel<< "\t" << PadNum << "\t" << xPos(mm) << "\t" << yPos(mm)
      while(mapStream>>mapp[line][0]>>mapp[line][1]>>mapp[line][2]>>mapp[line][3]>>mapp[line][4]>>mapp[line][5]>>mapp[line][6]){
          if (line == numberOfPads){
          cout<<"\n\nWARNING !\nMore lines in the 'mapStream' than space in the 'mapp' matrix, there are.\nFucked up, you have.\nFix it, you should.\n";
          return 0;
          }
          line++;
      }
  }
  else{
   cout<<"\nWarning ! Failed to open TPC-SEP2019-fullMap.dat, please start panicking at your earliest convenience.\n";
   return 0;
  }
  
  const int numberOfEnergies = 16;
  float fusionCharac[numberOfEnergies][4];
  const char *fCharacStreamChar;
  if (runOnTheCRC == true){fStreamName = "/afs/crc.nd.edu/user/m/mrenaud1/Public/SEP-pATTPC-analysis/signal-Li8Ar40to48Sc.dat";}
  else{fStreamName = "../ExpPrep/signal-Li8Ar40to48Sc.dat";}
  fCharacStreamChar = fStreamName.c_str(); //WARNING http://www.cplusplus.com/forum/general/100714/
  fstream characStream(fCharacStreamChar);
  line = 0;
  if ((useFusionGates == true || addSignalCharacOnPlots == true) && characStream.is_open()){//energy[kev], charge[% of BeamE], peakHeight[fraction of beam peak], dist travelled [estim., mm]
      while (characStream>>fusionCharac[line][0]>>fusionCharac[line][1]>>fusionCharac[line][2]>>fusionCharac[line][3]){
          if (line == numberOfEnergies){
              cout<<"\n\nWARNING !\nMore lines in the 'characStream' than space in the 'fusionCharac' matrix, there are.\nFucked up, you have.\nFix it, you should.\n";
              return 0;
          }
          line++;
      }
  }
  else if (useFusionGates == true && characStream.is_open() == false){
      cout<<"\nWarning ! Failed to open signal-Li8Ar40To48Sc.dat, please start panicking at your earliest convenience.\n";
      return 0;
  }
  


  
  //Definition of objects
  const char *fListChar;
  if (runOnTheCRC == true){
	if (multiRunAnalysis == false) {fStreamName = "/afs/crc.nd.edu/user/m/mrenaud1/Public/SEP-pATTPC-analysis/list-preAnal.txt";}
	else {fStreamName = "/afs/crc.nd.edu/user/m/mrenaud1/Public/SEP-pATTPC-analysis/list-analysis.txt";}
  }
  else{
	if (multiRunAnalysis) {fStreamName = "list-preAnal.txt";}
	else {fStreamName = "list-analysis.txt";}
  }
  fListChar = fStreamName.c_str(); //WARNING http://www.cplusplus.com/forum/general/100714/
  ifstream List(fListChar);


  string WhichFile;
  const char *FileName;
  
  TFile *file;
  TTree *tree;
  
  TH1I* whichChannel = new TH1I("channel", "channels read;entry;channel",301,-1,299);
  
  TH1I* timeBucketForm[100];
  TH1F* timeBucketFormAveraged;
  TH1I* timeBucketFormAuxillary[6];
  TH1F* TBFeventRemanent;
  TH1F* TBFrealigned;
  TH1F* preTBF;
  TH1F* whatsTheDerivative;
  TH2F* lengthVSheight;
  TH2F* lengthVScharge;
  TH2F* lengthVbackground;
  TH2F* lengthVtac;
  TH2F* chargeVtac;
  TH2F* peakHeightVtac;
  TH2F* chargeVpeakHeight;
  TH3F* track3D;
  
  timeBucketFormAveraged = new TH1F("TBFAv","Averaged channels;time bucket;Energy[arb.]",512,0,511);
  timeBucketFormAveraged->Fill(0.);
  if (realignTracks){
      TBFrealigned = new TH1F("TBrealigned","Realigned tracks;altered time bucket;Energy[arb]",512,0,511);
      TBFrealigned->Fill(0.);
  }
  if (fillIndividualChannels){
      for (unsigned int i=0; i<100;i++){
        timeBucketForm[i] = new TH1I(Form("TBF%i",i), Form("channel %i;time bucket;Energy[arb.]",i),512,0,511);
    }
  }
  for (unsigned int i=0; i<6;i++){
        timeBucketFormAuxillary[i] = new TH1I(Form("TBFA%i",i), Form("Aux. chan %itime bucket;Energy[arb.]",i),512,0,511);
  }
  TBFeventRemanent = new TH1F("TBFrem","Added events;time bucket;Energy[arb]",512,0,511);
  preTBF = new TH1F("preTBF","Added events;time bucket;Energy[arb]",512,0,511);
  whatsTheDerivative = new TH1F("derivative","derivative;time bucket;derivative",512,0,511);
  
  if (drawLengthVpeakHeight == true && substractBckgrd == false) lengthVSheight = new TH2F("LVH","Track length V. peak height; Tr. Length [time buckets];max peak height [arb]",512,0,511,1100,0,1100);
  if (drawLengthVpeakHeight == true && substractBckgrd == true) lengthVSheight = new TH2F("LVH","Track length V. peak height; Tr. Length [time buckets];max peak height [arb]",512,0,511,700,0,700);
  if (drawLengthVtrackCharge == true && substractBckgrd == false) lengthVScharge = new TH2F("LVC","Track length V. charge; Tr. Length [time buckets];dep. charge [arb]",512,0,511,4500,0,90000);//2500,0,25000);
  if (drawLengthVtrackCharge == true && substractBckgrd == true) lengthVScharge = new TH2F("LVC","Track length V. charge; Tr. Length [time buckets];dep. charge [arb]",512,0,511,2500,0,25000);
  if (drawLengthVbackgrd) lengthVbackground = new TH2F("LVbck","Track length V. backgrd; Tr. Length [time buckets];sub. bckgrd [arb]",512,0,511,500,200,1200);
  if (drawLengthVtac) lengthVtac = new TH2F("LVtac","Track Length V. TAC;Tr. Length[TB];tac [arb]",512,0,511,850,0,850);
  if (drawPeakHeightVtac) peakHeightVtac = new TH2F("phVtac","Peak height V. TAC;max peak height [arb];tac [arb]",1100,0,1100,850,0,850);
  if (drawChargeVtac) chargeVtac = new TH2F("QVtac","Charge V. TAC;dep. charge [arb];tac [arb]",2500,0,25000,850,0,850);
  if (drawChargeVpeakHeight == true && substractBckgrd == true) chargeVpeakHeight = new TH2F("CVPH","Charge V. peak height; dep. charge [arb];max peak height [arb]",2500,0,25000,700,0,700);
  if (drawChargeVpeakHeight == true && substractBckgrd == false) chargeVpeakHeight = new TH2F("CVPH","Charge V. peak height; dep. charge [arb];max peak height [arb]",4500,0,90000,1100,0,1100);
  if (draw3DTrack) track3D = new TH3F("track3D","track representation;Time buckets;x [mm];y [mm]",512,0,511,244,-121,121,244,-121,121);

//   TH1F* timeBucketForm = new TH1F("TBF", "channel 1;time bucket;Energy[arb.]",512,0,511);
  

/* Unused histograms
 * 
 * 
 * TH2F* EadcCh;
  TH2F* adcPulse;
  if (needADC == 1){
      EadcCh = new TH2F("adc/ch", "adc w.r.t channel;Channel; E[arb.]",96,1,97,4100,0,4100);
      adcPulse = new TH2F("adcPulse", "pulser-data plot;pulse;adc;#Count",3100,0,6200,2050,0,4100);
  }
  TH2F* SscalCh;
  TH1F* SPulser;
  if (needScaler == 1){
      SscalCh = new TH2F("scaler/ch","scaler use profile;Channel; ?[very arb.]",16,1,17,5000,0,68800);
      SPulser = new TH1F("Pulser","pulser profile;[arb];counts",13000,0,13000);
  }
  TH2F* TtdcCh;
  if (needTDC == 1){TtdcCh = new TH2F("tdc/ch","tdc use profile;Channel; T[arb.]",96,1,97,5000,0,23000);}
  TH1F* PSDstrip[24];
  TH2F* stripCalibHist[24];   //for calibration, need to create 1 point per event with hit that goes E-back v. E-strip
  TH1F* EspectrDet[96];
  
  for (int stripNumber = 0;stripNumber<24; stripNumber++){
	  PSDstrip[stripNumber] = new TH1F(Form("PSDstrip_%d",stripNumber+1),Form("Pos. Sens. Det.%d E spectr",stripNumber+1),2050,0,4100);
	  stripCalibHist[stripNumber] = new TH2F(Form("calibStripHisto_%d", stripNumber+1),Form("E-E spectr., Det.%d;E back [arb.];E PSD [arb.]", stripNumber+1),4100,0,4100,450,0,450);//300,0,300,450,0,450);
  }
  for (int detNumber = 0;detNumber<96;detNumber++){
	  EspectrDet[detNumber] = new TH1F(Form("EspectrDet_%d",detNumber+1),Form("E spectrum of det %d",detNumber+1),4100,0,4100);
  }
  
  TH1F* adcCh = new TH1F("adc", "adc channel use;Channel;#Count",96,1,97);
  TH1F* tdcCh = new TH1F("tdc", "tdc channel use;Channel;#Count", 96,1,97);
  
  TH1F* PSDtdc[30];
  
  PSDtdc[stripNumber] = new TH1F(Form("PSDtdc_%d",stripNumber),Form("Pos. Sens. Det.%d t response",stripNumber),1000,0,65536);
 
  TH2F* Eadc2Ch = new TH2F("adc2/ch", "sec. adc w.r.t channel;Channel; E[arb.]",96,0,96,4200,0,4200);
  
  TH1F* scalCh = new TH1F("scaler", "scaler thing use;Channel num.;#Count",16,0,16);
  			 scalCh->Fill(ch[j]);
  example of 3d hist: TH3F* S_QE = new TH3F("S_QE", "Signal avec QE;Y [m];Z [m];#hits", 8,-0.026,0.026, 8,-0.026,0.026,100,0,0.25);
*/

  //def of (de)convolute
  float anArray1[16],anArray2[16],anArray3[32];
  float dataTest[16];// {0.,0.00150362,0.018439,0.0835921,0.142119,0.0933979,0.00707567,0.000207203,2.26662e-06,1.07102e-08,1.95578e-08,2.6077e-08,2.98023e-08,0.0183157,0.00300689,0.00021422};
  float responseTest[16] = {0.};
  float answerTest[32] = {0.};
  for (int i = 0;i<16;i++){
//       if (i<8) responseTest[i] = TMath::Gaus(i,0,1);
//       if (i>=8) responseTest[i] = 0.;
//       if (i>=8 && i<16) responseTest[i] = TMath::Gaus(-16+i,0,1);
      responseTest[i] = TMath::Gaus(i,7,2);
      dataTest[i] = 0.;
      anArray1[i] = i;
      if (i<16) anArray2[i] = i;
  }
//   dataTest[2] = dataTest[4] = 1;
  dataTest[7] = 4;
  for (int i=0;i<32;i++){
      cout<<endl<<"answerTest["<<i<<"] = "<<answerTest[i];
      anArray3[i] = i;
      if(i<16){cout<<",\tdataTest["<<i<<"] = "<<dataTest[i];/*dataTest2[i] = dataTest[i];*/}
      if(i<16){cout<<",\tresponseTest["<<i<<"] = "<<responseTest[i];/*responseTest2[i] = responseTest[i];*/}
  }
  const float *theArray1 = anArray1;
  const float *theArray2 = anArray2;
  const float *theArray3 = anArray3;
  const float *dataTest2 = dataTest;
  const float *responseTest2 = responseTest;
  TGraph* grD = new TGraph(16,theArray1,dataTest);
  TGraph* grR = new TGraph(16,theArray2,responseTest2);
  
  cout<<endl<<endl<<endl;
  convlv(dataTest,16,responseTest,16,1,answerTest);
  for (int i=0;i<32;i++){
      cout<<endl<<"answerTest["<<i<<"] = "<<answerTest[i];
      if(i<16){cout<<",\tdataTest["<<i<<"] = "<<dataTest[i];}
      if(i<16){cout<<",\tresponseTest["<<i<<"] = "<<responseTest[i];}
  }
  cout<<endl;
  const float *answerTest2 = answerTest;
  TGraph* grA = new TGraph(32,theArray3,answerTest2);
  
  TCanvas* Ctestc = new TCanvas("Ctestc","Ctestc",600,600);
  grA->SetLineColor(kGreen);
  grA->Draw("AC*");
  grD->SetLineColor(kBlue);
  grD->Draw("same");
  grR->SetLineColor(kRed);
  grR->Draw("same");
  
  return 0;
  
  
  
  
  
  //Def of variables /!\ mind that variable have to match for root extraction, storing Int_t->Double_t won't fly
  Int_t num_hits;
  Int_t whatChan;
  Int_t data[2055][517];
  globalNumberErrors = 0;
  globalNumberOfPUevents = 0;
  globalNumberOfLongEvents = 0;
  bool done = 0;
  Int_t auxChanDone = 0;
  bool checkCleared = 0; //WARNING when 'noPileUp condition relaxed, this check will be a problem
  Int_t doubleCheckPU = 0;
  Double_t* daPeak;
  Int_t tacPeak;
  Double_t TACval = 0;
  Int_t daPeakElement = 0;
  Float_t storedDerivative[100], derivativeAverage[100];
  Int_t howManyInDerivArray = 0;
  Float_t maxDerivativePosition[2];
  Float_t temp[3]; //variables for random use
  bool fusionTestValidated = 0;
  Float_t padsPosition[100][512][2];
  int checkTrack3D;

  

  while (List >> WhichFile){
	
	FileName = WhichFile.c_str(); //TFile needs conversion to char
	if (FileName[0] != '*'){//filenames commented by adding a '*' in front won't be read
   cout<<"\nOpening next file, "<<FileName<<", processing...\n";
   
        
	//Reading & Initialising
	file = new TFile(FileName);
	if (file->IsOpen()){
    tree = (TTree*)(file->Get("T"));  
    tree->SetBranchAddress("br_len", &num_hits);
    tree->SetBranchAddress("myints", data);
    
    //IMPORTANT! structure: Cobo<< "\t" << Asad<< "\t" << Aget<< "\t" << Channel<< "\t" << PadNum
    
    
    //Recursively go through tree, getting "Fill" after "Fill"findDerivative
    
    for (Int_t i=EVENT; i</*EVENT+1000*/tree->GetEntries(); i++){
        tree->GetEntry(i);
        if (beQuite != 1) cout<<"\nEvent #"<<i<<", num_hits : "<<num_hits<<endl<<endl;
        if (num_hits > 100){cout<<"\nWARNING, event "<<i<<"has too many entries.";globalNumberErrors++;}//TODO add error handeling
        
        //resetting variables/histograms
        auxChanDone = 0;
        for (int i=0;i<512;i++){
            timeBucketFormAveraged->SetBinContent(i,0.);
            whatsTheDerivative->SetBinContent(i,0.);
            preTBF->SetBinContent(i,0.);
            if (i<100) {storedDerivative[i] = 0; derivativeAverage[i] = 0;}
            if (i<3) {temp[i] = 0;}
            if (i<2) {maxDerivativePosition[i] = 0;}
        }
//         fill(storedDerivative,storedDerivative+100,0);
//         fill(derivativeAverage,derivativeAverage+100,0);
        checkCleared = 0;
        howManyInDerivArray = 0; //how many non-zero entries in the array containing the derivative of the track
        doubleCheckPU = 0;
        daPeakElement = 0;
        fusionTestValidated = 0;
        TACval = 0;

        
        //loop through channels in event
        for (Int_t j=0; j<num_hits; j++) {
            
            done = 0;
//             for (int i=0;i<512;i++){padPositions[i][0] = 0.;padPositions[i][1] = 0.;}
            
            for (Int_t k=0; k<517; k++){
                /*temp cout data    if (k%100==0){cout<<"\nnum_hits : "<<num_hits<<", data["<<j<<"]["<<k<<"] = "<<data[j][k];}*/
                
                if (done==0 && k == 4){//ID what the channel is 
                    if (data[j][0] < 2){
                        for (Int_t check = 0; check<numberOfPads; check++){
                            if ((data[j][0] == mapp[check][0]) && (data[j][1] == mapp[check][1]) && (data[j][2] == mapp[check][2]) && (data[j][3] == mapp[check][3])){
                                whatChan = mapp[check][4];
                                checkTrack3D = check;
                            }
                        }
                    }  
                    else if(data[j][0] == 2){whatChan =3000+data[j][3];} // ID the auxillary channel
                    else {cout<<"\n*Posh british accent* I say, this is rather unexpected indeed old chap'...\n\tEntry num. : "<<i<<", j : "<<j<<", k : "<<k<<endl<<endl;return 0;}
                    whichChannel->SetBinContent(j,whatChan);
                    done=1;
//                     cout<<"\n\twhichChannel here is : "<<data[j][k];
                }// \end ID what the channel is
                if (done==1 && k>4 && data[j][0] != 2){//fill in time bucket format for non-auxillary channels
                    if (fillIndividualChannels == true && j<100){timeBucketForm[j]->SetBinContent(k-5,data[j][k]);}
                    timeBucketFormAveraged->AddBinContent(k-5,(float)(data[j][k])/max(num_hits-6,1)); //WARNING diff gains inner-outer
                }  
                if (done==1 && k>4 && data[j][0] == 2 && auxChanDone<6){ //fill in time bucket format for auxillary channels
                    timeBucketFormAuxillary[auxChanDone]->SetBinContent(k-5,data[j][k]);
                    if(k==516){auxChanDone++;} //move to next histogram
                }
                else if(auxChanDone >= 6){cout<<"\nError, which is sad.\nAnd 'sad' backwards is 'das', and das not good...\nLook at the Auxillaries.\n";return 0;}
                if (k>5 && draw3DTrack == true && num_hits <100){
                    padsPosition[j][k-4][0] = mapp[checkTrack3D][5];
                    padsPosition[j][k-4][1] = mapp[checkTrack3D][6];
                }
                

            }
        }
        
        //section to see if event has pile-up & attempt clean-up/renorm
        temp[0] = whichAuxChanIsTAC(timeBucketFormAuxillary[1]);
        Int_t pileUpInEvent = pileUp(timeBucketFormAuxillary[(int)temp[0]]);
        if (beQuite != 1)cout<<"TAC is in aux. chan. "<<temp[0]<<endl;
        if (beQuite != 1)cout<<"#ofEvents : "<<pileUpInEvent;
        if (pileUpInEvent == 1){ if (beQuite != 1) cout<<"\t-> Accepted"<<endl;}
        else {if (beQuite != 1){cout<<"\t-> Rejected"<<endl<<endl;} globalNumberOfPUevents++;}
        
        if (done == 1 && pileUpInEvent == 1){//if at end of event+event was clean && no pile-up, perform last check; then add the average to the signal
            preTBF->Add(timeBucketFormAveraged);
            daPeak = locateTACsignal(timeBucketFormAveraged, pileUpInEvent+2);
            if (beQuite != 1) cout<<"Watch : "<<daPeak[0]<<"\t"<<daPeak[1]<<"\t"<<daPeak[2]<<endl;
            for (int j=0;j<pileUpInEvent+2;j++){
                if (daPeak[j] >=1 && daPeak[j] <= 513){doubleCheckPU++;}
                else if (daPeak[j] >= 513){daPeakElement++;}//if you're here the real peak you wantnot the first, so go fetch next entry, not the first "peak" -> daPeak[daPeakElement] instead of daPeak[0]
            }
            
            //Double-check pileUp
            if (pileUpInEvent == doubleCheckPU){
                if (acceptPeaksBelowBeam == true){checkCleared = 1;}
                else if (acceptPeaksBelowBeam == false && daPeak[daPeakElement] > pointOfPeakRejection){checkCleared =1;}
                else {if (beQuite != 1) {cout<<"\n\tWARNING ! Peak beyond rejection point.\t-> Rejected"<<endl<<endl;} globalNumberErrors++;globalNumberOfLongEvents++;}
            }
            else if (doubleCheckPU == 0){if (beQuite != 1) {cout<<"\n\tWARNING ! No signal !\t-> Rejected"<<endl<<endl;} globalNumberErrors++;}
            else{
                if (beQuite != 1){cout<<"\n\tWARNING ! Multiple peaks for one TAC signal !\t-> Rejected"<<endl<<endl;}
                globalNumberErrors++;
                globalNumberOfPUevents++;
            }
            
            //Get TAC value
            tacPeak = timeBucketFormAuxillary[(int)temp[0]]->GetMaximumBin();
            TACval = timeBucketFormAuxillary[(int)temp[0]]->GetBinContent(tacPeak) - timeBucketFormAuxillary[(int)temp[0]]->GetBinContent(tacPeak-6);
            if (beQuite != 1) cout<<"TAC position : "<<tacPeak<<"\tTAC value : "<<TACval<<endl;
            
            
            
            if (checkCleared == 1){//actual clean & treatement
                temp[0]=0;
                //fit noise before track and substract it.
                if (beQuite != 1) cout<<"location of peaks : "<<daPeak[daPeakElement]<<endl;
                TF1 *fa0 = new TF1("fa0", "pol0(0)",5,daPeak[daPeakElement]-20);
                timeBucketFormAveraged->Fit(fa0,"RQ"); //fit options: no draw, limited range, quiet
                if (drawLengthVbackgrd) temp[2] = fa0->GetParameter(0);
                if (beQuite != 1) cout<<"fa0->GetParameter(0) = "<<fa0->GetParameter(0)<<endl;
                fa0->SetRange(0,511);
                if (substractBckgrd){preTBF->Add(fa0,-1,"");}
                preTBF->SetBinContent(0,0);
                
                //locate start of track and substract after track noise
                for (int bin=0;bin<(511.-daPeak[daPeakElement])/5.;bin++){//use the derivative of the track signal to help ID the end (sharp drop)
//                     cout<<daPeak[daPeakElement]+bin*5<<"\t"<<findDerivative(preTBF,daPeak[daPeakElement]+bin*5,derivativeCalcWindowSize)<<endl;
                    storedDerivative[bin] = findDerivative(preTBF,daPeak[daPeakElement]+bin*5,derivativeCalcWindowSize);
                    whatsTheDerivative->SetBinContent(daPeak[daPeakElement]+bin*5,storedDerivative[bin]);
                    howManyInDerivArray++;
                }
                fa0->SetRange(daPeak[daPeakElement]+10,505);
                preTBF->Fit(fa0,"RQ");
                if (drawLengthVbackgrd == true && substractBckgrd == true) temp[2] += fa0->GetParameter(0);
                else if (drawLengthVbackgrd == true && substractBckgrd == false) temp[2] = fa0->GetParameter(0);
                if (beQuite != 1) cout<<"fa0->GetParameter(0), v2 = "<<fa0->GetParameter(0)<<endl;
                for (int j=daPeak[daPeakElement];j<511;j++){
                    if (preTBF->GetBinContent(j) <= max(fa0->GetParameter(0)*1.5,7.)){
//                         cout<<"\nj : "<<j;
                        for (int u=0; u<howManyInDerivArray-3; u++){
//                             cout<<u<<"\t"<<storedDerivative[u]<<"\t"<<AveragE(storedDerivative,3,3+u)<<endl; 
                            derivativeAverage[u]=AveragE(storedDerivative,3,3+u); //TODO optimise this section ! too many loops !
                        }
                        break;
                    }
                }
                for (int p = howManyInDerivArray-3;p>0;p--){
                    if (derivativeAverage[p] < -1.01){
                        if (beQuite != 1) cout<<"tentative : "<<p<<"\t"<<daPeak[daPeakElement]+(p+1)*5<<endl; //p +1(for stepping) 
                        fa0->SetRange(daPeak[daPeakElement]+(p+1)*5+10,505);
                        preTBF->Fit(fa0,"RQ");
//                         cout<<"fa0->GetParameter(0), v3 = "<<fa0->GetParameter(0)<<endl;
                        fa0->SetRange(0,511);
                        if (substractBckgrd){preTBF->Add(fa0,-1.01,"");}
                        
                        if (realignTracks){
                            temp[0] = 0;
                            temp[1] = 5.;
                            for (int bin=0;bin<40;bin++){
                                if (temp[0] >= findDerivative(preTBF,daPeak[daPeakElement]+(p+1)*5-20+bin,3)){
                                    maxDerivativePosition[0]=daPeak[daPeakElement]+(p+1)*5-20+bin; //(peak position)+(to reconstruct position of beginning track)+(-20+bin: scan that passage. Thus, here, find position of max derivative -> start track)
                                    temp[0] = findDerivative(preTBF,daPeak[daPeakElement]+(p+1)*5-20+bin,3); //store derivative
                                }
                                if (temp[1] < findDerivative(preTBF,(int)daPeak[daPeakElement]-40+bin,3)){
                                    maxDerivativePosition[1]=(int)daPeak[daPeakElement]-40+bin; //find end of track
                                    temp[1] = 1000000.; //break
                                }
                            }
                            if (beQuite != 1) cout<<"\nmaxDerivativePosition[0] = "<<maxDerivativePosition[0]<<", maxDerivativePosition[1] = "<<maxDerivativePosition[1]<<endl; //respectively position of start and end of track
                            
//                             cout<<"\n\t\tNotice me Senpai !\n"; 
                            
                            if (drawLengthVbackgrd) lengthVbackground->Fill(maxDerivativePosition[0]+estimatedWidth-maxDerivativePosition[1],temp[2]); //temp[2] holds the background/baseline
                            if (drawLengthVpeakHeight){
                                if (useFusionGates == false){lengthVSheight->Fill((maxDerivativePosition[0]+estimatedWidth-maxDerivativePosition[1]),preTBF->GetMaximum());}
                                else {
                                    for (line=15;line>=0;line--){
                                        if ((maxDerivativePosition[0]+estimatedWidth-maxDerivativePosition[1]) <= fusionCharac[line][3]/(driftV*timeBucketSize)){
                                            fusionTestValidated = 1;
                                            temp[0] = preTBF->Integral(maxDerivativePosition[1],maxDerivativePosition[0]+estimatedWidth,"");
                                            break;
                                        }
                                        globalNumberOfLongEvents++;
                                    }
                                    //current condition : >=90% of calculated peak height, and <=110% of calculated charge
                                    if (fusionTestValidated == true && temp[0]<92.545455*(maxDerivativePosition[0]+estimatedWidth-maxDerivativePosition[1])-755.81819){lengthVSheight->Fill(maxDerivativePosition[0]+estimatedWidth-maxDerivativePosition[1],preTBF->GetMaximum());}
                                }
                                if (drawPeakHeightVtac){peakHeightVtac->Fill(preTBF->GetMaximum(),TACval);}
                            }
                            if (drawLengthVtrackCharge){
                                if (useFusionGates == false){lengthVScharge->Fill((maxDerivativePosition[0]+estimatedWidth-maxDerivativePosition[1]),preTBF->Integral/*AndError*/(maxDerivativePosition[1],maxDerivativePosition[0]+estimatedWidth,""));}
                                else {
                                    if (fusionTestValidated == true && temp[0]<92.545455*(maxDerivativePosition[0]+estimatedWidth-maxDerivativePosition[1])-755.81819){
                                        lengthVScharge->Fill(maxDerivativePosition[0]+estimatedWidth-maxDerivativePosition[1],preTBF->Integral/*AndError*/(maxDerivativePosition[1],maxDerivativePosition[0]+estimatedWidth,""));
                                    }
                                }
                                if (drawChargeVtac){chargeVtac->Fill(preTBF->Integral/*AndError*/(maxDerivativePosition[1],maxDerivativePosition[0]+estimatedWidth,""),TACval);}
                            }
                            if (drawLengthVtac){lengthVtac->Fill((maxDerivativePosition[0]+estimatedWidth-maxDerivativePosition[1]),TACval);}
                            if (drawChargeVpeakHeight){chargeVpeakHeight->Fill(preTBF->Integral/*AndError*/(maxDerivativePosition[1],maxDerivativePosition[0]+estimatedWidth,""),preTBF->GetMaximum());}
                            if (draw3DTrack == true && num_hits < 100){
                                for (int i=0;i<num_hits;i++){
                                    for (int j=maxDerivativePosition[1];j<maxDerivativePosition[0]+estimatedWidth;j++){
                                        track3D->Fill(j+(510-maxDerivativePosition[0]+estimatedWidth),padsPosition[i][j][0],padsPosition[i][j][1]);
                                    }
                                }
                            }
                            else if (draw3DTrack == true && num_hits >= 100){
                                if (!beQuite) cout<<"\nWARNING ! Event incompatible with 'draw3DTrack' option. Dropping option.\n";
                                draw3DTrack = false;
                            }
                        }// END 'realignTracks'
                        else if (drawLengthVpeakHeight == true && realignTracks == false){
//                             cout<<"On fill LVPH avec : "<<(daPeak[daPeakElement]+(p+1)*5-(daPeak[daPeakElement]-10))<<"\t"<<preTBF->GetMaximum()<<endl;
                            lengthVSheight->Fill((daPeak[daPeakElement]+(p+1)*5-(daPeak[daPeakElement]-10)),preTBF->GetMaximum());
                        }
                        if (drawLengthVtrackCharge == true && realignTracks == false){
//                             cout<<"On fill LVC avec : "<<(daPeak[daPeakElement]+(p+1)*5-(daPeak[daPeakElement]-10))<<"\t"<<preTBF->Integral/*AndError*/(daPeak[daPeakElement]-10,daPeak[daPeakElement]+(p+1)*5,"")<<endl;
                            lengthVScharge->Fill((daPeak[daPeakElement]+(p+1)*5-(daPeak[daPeakElement]-10)),preTBF->Integral/*AndError*/(daPeak[daPeakElement]-10,daPeak[daPeakElement]+(p+1)*5,""));
                        }
                        if (realignTracks == false && drawLengthVbackgrd == true) lengthVbackground->Fill((daPeak[daPeakElement]+(p+1)*5-(daPeak[daPeakElement]-10)),temp[2]);
                        if (realignTracks == false && drawLengthVtac == true) lengthVtac->Fill((daPeak[daPeakElement]+(p+1)*5-(daPeak[daPeakElement]-10)),TACval);
                        if (realignTracks == false && drawChargeVpeakHeight == true){chargeVpeakHeight->Fill(preTBF->Integral/*AndError*/(daPeak[daPeakElement]-10,daPeak[daPeakElement]+(p+1)*5,""),preTBF->GetMaximum());}
                        
                        p = 1; //break;
                    }
                }
               
                
                //clean up negatives && filling 
                if (drawTBFeventRemanent){
                    if (useFusionGates == false){
                        TBFeventRemanent->Add(preTBF);
                        for (int i=0;i<513;i++){if (TBFeventRemanent->GetBinContent(i)<0){TBFeventRemanent->SetBinContent(i,0.);}}
                    }
                    else if (fusionTestValidated == true  && temp[0]<92.545455*(maxDerivativePosition[0]+estimatedWidth-maxDerivativePosition[1])-755.81819){
                        TBFeventRemanent->Add(preTBF);
                        for (int i=0;i<513;i++){if (TBFeventRemanent->GetBinContent(i)<0){TBFeventRemanent->SetBinContent(i,0.);}}
                    }
                }
                if (realignTracks){
                    if (useFusionGates == false){
                        for (int b=0;b<((int)maxDerivativePosition[0]+estimatedWidth);b++){
                            TBFrealigned->AddBinContent(b+(510-(maxDerivativePosition[0]+estimatedWidth)),max(preTBF->GetBinContent(b),0.));
                        }
                    }
                    else if (fusionTestValidated == true  && temp[0]<92.545455*(maxDerivativePosition[0]+estimatedWidth-maxDerivativePosition[1])-755.81819){
                        for (int b=0;b<((int)maxDerivativePosition[0]+estimatedWidth);b++){
                            TBFrealigned->AddBinContent(b+(510-(maxDerivativePosition[0]+estimatedWidth)),max(preTBF->GetBinContent(b),0.));
                        }
                    }
                }
                delete fa0;
               
            }// END checkCleared
        } // END noPileUp
        
    } // END of loop through events
   } // END if("file.is_open()")
   else {cout << "\n\nWARNING ! FAILED TO OPEN FILE "<<WhichFile<<", ABORTING.\nNow go and find out where you messed up. You prick.\n\n";break;return 0;}
   
   
   file->Close();
   delete file;
   
   
   }//FileName not *-ed
  }//while loop
  
  //Def of projections
  

  
  
  //Def fit fonction
/*  //Step-function convolution
  Float_t ampli = 50.;
  Float_t posi = 500.;
  TF1 *steppy = new TF1("steppy",the_steppar,0,512,2);
  steppy->SetParameter(0,ampli);
  steppy->SetParameter(1,posi);
  TF1 *detResp = new TF1("detResp",the_gausppar,0,512,3);
  detResp->SetParameter(2,1./timeBucketSize);
  TF1Convolution *step_conv = new TF1Convolution(steppy,detResp,485,512,true); 
  TF1* f_conv = new TF1("f_conv",*step_conv,485,512,step_conv->GetNpar());
  // END step-function*/
  
  
  
  
// -----------------------------DEF OF CANVAS AND PLOT---------------------------------------

  //Defining style of graph
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(1111111);//now showing integral under-and-overflow ! //old: 1001110
  gStyle->SetOptFit(1);
  gStyle->ToggleEventStatus();
  
  //save path definition
  string fStreamNameS;
  if (runOnTheCRC == false){fStreamName = "/home/mrenaud1/Documents/Li8Analysis/requestedGraphs/autoSaved-";}
  else {fStreamName = "/afs/crc.nd.edu/user/m/mrenaud1/Public/stockageGraphAnal/autoSaved-";}
  //prepare for point additions
  float_t ptx[numberOfEnergies], pty[numberOfEnergies], ptz[numberOfEnergies];
  for (int i=0;i<numberOfEnergies;i++){ptz[i] = i;}
  
  
  
  
  if (drawWhichPadIsIt){
      TCanvas* Cchan = new TCanvas("Cchan","Cchan",400,400);
        Cchan->ToggleEventStatus();
        whichChannel->Draw();
  }
  else {delete whichChannel;}
  
  if (drawIndividualChannels){//draws the first 100 used individual channels in that event
      
    TCanvas* CIC1 = new TCanvas("cic1","cic1",1000,1000);
        CIC1->Divide(4,5); 
        CIC1->ToggleEventStatus();
        for (unsigned int divNumb=1;divNumb<21;divNumb++){
            CIC1->cd(divNumb);
            timeBucketForm[divNumb-1]->Draw();
        }
    TCanvas* CIC2 = new TCanvas("cic2","cic2",1000,1000);
        CIC2->Divide(4,5); 
        CIC2->ToggleEventStatus();
        for (unsigned int divNumb=1;divNumb<21;divNumb++){
            CIC2->cd(divNumb);
            timeBucketForm[divNumb-1+20]->Draw();
        }
    TCanvas* CIC3 = new TCanvas("cic3","cic3",1000,1000);
        CIC3->Divide(4,5); 
        CIC3->ToggleEventStatus();
        for (unsigned int divNumb=1;divNumb<21;divNumb++){
            CIC3->cd(divNumb);
            timeBucketForm[divNumb-1+40]->Draw();
        }
    TCanvas* CIC4 = new TCanvas("cic4","cic4",1000,1000);
        CIC4->Divide(4,5); 
        CIC4->ToggleEventStatus();
        for (unsigned int divNumb=1;divNumb<21;divNumb++){
            CIC4->cd(divNumb);
            timeBucketForm[divNumb-1+60]->Draw();
        }
    TCanvas* CIC5 = new TCanvas("cic5","cic5",1000,1000);
        CIC5->Divide(4,5); 
        CIC5->ToggleEventStatus();
        for (unsigned int divNumb=1;divNumb<21;divNumb++){
            CIC5->cd(divNumb);
            timeBucketForm[divNumb-1+80]->Draw();
        }
  }
  
  if (drawSummedChannels){ //single event
      TCanvas* Call = new TCanvas("call","call",500,500);
        Call->ToggleEventStatus();
        timeBucketFormAveraged->Draw();
  }
  else {delete timeBucketFormAveraged;}
  
  if (drawAuxillaryChannels){
    TCanvas* CauxChan = new TCanvas("CauxChan","CauxChan",800,800);
        CauxChan->Divide(2,3);
        CauxChan->ToggleEventStatus();
        for (unsigned int divNumb=1;divNumb<7;divNumb++){
            CauxChan->cd(divNumb);
            timeBucketFormAuxillary[divNumb-1]->Draw();
        }
  }
  else {delete *timeBucketFormAuxillary;}
  
  if (drawTBFeventRemanent){ //summed up over all events
      TCanvas* CTBFevRem = new TCanvas("CTBFevRem","CTBFevRem",600,600);
        CTBFevRem->ToggleEventStatus();
        TBFeventRemanent->Draw();
        if (saveHistograms){
            fStreamNameS = fStreamName + "unalignedProfile.root";
            const char *nameCTBFsvRem = fStreamNameS.c_str(); //WARNING http://www.cplusplus.com/forum/general/100714/
            CTBFevRem->SaveAs(nameCTBFsvRem);
            delete CTBFevRem;
        }
  }
  else {delete TBFeventRemanent;}
  
  if (realignTracks){
      TCanvas* CTBFR = new TCanvas("CTBFR","CTBFR",600,600);
        CTBFR->ToggleEventStatus();
        TBFrealigned->Draw();
        if (saveHistograms){
            fStreamNameS = fStreamName + "realignedProfile.root";
            const char *nameCTBFR = fStreamNameS.c_str(); //WARNING http://www.cplusplus.com/forum/general/100714/
            CTBFR->SaveAs(nameCTBFR);
            delete CTBFR;
        }
  }
  
  if (drawDerivatives){
      TCanvas* Cderi = new TCanvas("deri","deri",500,500);
        Cderi->ToggleEventStatus();
        whatsTheDerivative->GetXaxis()->SetRangeUser(0,505);
        whatsTheDerivative->Draw();
  }
  else {delete whatsTheDerivative;}
  
  
  //---------2D plots------------------------------------------
  if (drawLengthVpeakHeight){
      TCanvas* CLVPH = new TCanvas("CLVPH","CLVPH",600,600);
        CLVPH->SetLogz();
        lengthVSheight->Draw("colz");
        if (addSignalCharacOnPlots){
            TMarker *mLi = new TMarker(Li7avTrLenght,Li7peakHeight*meanBeamPeakHeight,23);
            mLi->Draw("same");
            for (int i=0;i<numberOfEnergies;i++){
                ptx[i] = fusionCharac[i][3]/(driftV*timeBucketSize);
                pty[i] = fusionCharac[i][2]*meanBeamPeakHeight;
//                 cout<<"\nptx["<<i<<"] = "<<ptx[i]<<"\t, pty["<<i<<"] = "<<pty[i]<<endl;
                TMarker *m = new TMarker(ptx[i],pty[i],29);
                m->SetMarkerColor(4);
                m->Draw("same");
            }
        }
        if (saveHistograms){
            fStreamNameS = fStreamName + "TrLengthVpeakHeight.root";
            const char *nameCLVPH = fStreamNameS.c_str(); //WARNING http://www.cplusplus.com/forum/general/100714/
            CLVPH->SaveAs(nameCLVPH);
            delete CLVPH;
        }
  }
  
  if (drawLengthVtrackCharge){
      TCanvas* CLVC = new TCanvas("CLVC","CLVC",600,600);
        CLVC->SetLogz();
        lengthVScharge->Draw("colz");
        if (addSignalCharacOnPlots){
            TMarker *nLi = new TMarker(Li7avTrLenght,Li7ChargePercentage*meanBeamCharge,23);
            nLi->Draw("same");
            for (int i=0;i<numberOfEnergies;i++){
                ptx[i] = fusionCharac[i][3]/(driftV*timeBucketSize);
                pty[i] = fusionCharac[i][1]*meanBeamCharge/100.;
//                 cout<<"\nptx["<<i<<"] = "<<ptx[i]<<"\t, pty["<<i<<"] = "<<pty[i]<<endl;
                TMarker *n = new TMarker(ptx[i],pty[i],29);
                n->SetMarkerColor(2);
                n->Draw("same");
            }
        }
        if (saveHistograms){
            fStreamNameS = fStreamName + "TrLengthVcharge.root";
            const char *nameCLVC = fStreamNameS.c_str(); //WARNING http://www.cplusplus.com/forum/general/100714/
            CLVC->SaveAs(nameCLVC);
            delete CLVC;
        }
  }
  
  if (drawLengthVbackgrd){
      TCanvas* CLVbck = new TCanvas("CLVbck","CLVbck",600,600);
        CLVbck->SetLogz();
        lengthVbackground->Draw("colz");
        if (saveHistograms){
            fStreamNameS = fStreamName + "TrLengthVbackground.root";
            const char *nameCLVbck = fStreamNameS.c_str(); //WARNING http://www.cplusplus.com/forum/general/100714/
            CLVbck->SaveAs(nameCLVbck);
            delete CLVbck;
        }
  }
  
  if (drawLengthVtac){
      TCanvas* CLVTAC = new TCanvas("CLVTAC","CLVTAC",600,600);
        CLVTAC->SetLogz();
        lengthVtac->Draw("colz");
        if (saveHistograms){
            fStreamNameS = fStreamName + "TrLengthVtac.root";
            const char *nameCLVTAC = fStreamNameS.c_str(); //WARNING http://www.cplusplus.com/forum/general/100714/
            CLVTAC->SaveAs(nameCLVTAC);
            delete CLVTAC;
        }
  }
  
  if (drawChargeVtac){
      TCanvas* CQVTAC = new TCanvas("CQVTAC","CQVTAC",600,600);
        CQVTAC->SetLogz();
        chargeVtac->Draw("colz");
        if (saveHistograms){
            fStreamNameS = fStreamName + "chargeVtac.root";
            const char *nameCQVTAC = fStreamNameS.c_str(); //WARNING http://www.cplusplus.com/forum/general/100714/
            CQVTAC->SaveAs(nameCQVTAC);
            delete CQVTAC;
        }
  }
  
  if (drawPeakHeightVtac){
      TCanvas* CPHVTAC = new TCanvas("CPHVTAC","CPHVTAC",600,600);
        CPHVTAC->SetLogz();
        peakHeightVtac->Draw("colz");
        if (saveHistograms){
            fStreamNameS = fStreamName + "peakHeightVtac.root";
            const char *nameCPHVTAC = fStreamNameS.c_str(); //WARNING http://www.cplusplus.com/forum/general/100714/
            CPHVTAC->SaveAs(nameCPHVTAC);
            delete CPHVTAC;
        }
  }
  
  if (drawChargeVpeakHeight){
      TCanvas* CQVPH = new TCanvas("CQVPH","CQVPH",600,600);
        CQVPH->SetLogz();
        chargeVpeakHeight->Draw("colz");
        if (addSignalCharacOnPlots){
            TMarker *oLi = new TMarker(Li7ChargePercentage*meanBeamCharge,Li7peakHeight*meanBeamPeakHeight,23);
            oLi->Draw("same");
            for (int i=0;i<numberOfEnergies;i++){
                ptx[i] = fusionCharac[i][1]*meanBeamCharge/100.;
                pty[i] = fusionCharac[i][2]*meanBeamPeakHeight;
//                 cout<<"\nptx["<<i<<"] = "<<ptx[i]<<"\t, pty["<<i<<"] = "<<pty[i]<<endl;
                TMarker *o = new TMarker(ptx[i],pty[i],29);
                o->SetMarkerColor(3);
                o->Draw("same");
            }
        }
        if (saveHistograms){
            fStreamNameS = fStreamName + "chargeVpeakHeight.root";
            const char *nameCQVPH = fStreamNameS.c_str(); //WARNING http://www.cplusplus.com/forum/general/100714/
            CQVPH->SaveAs(nameCQVPH);
            delete CQVPH;
        }
  }
  
  if (draw3DTrack){
      TCanvas* C3DT = new TCanvas("C3DT","C3DT",800,600);
        track3D->Draw();
        if (saveHistograms){
            fStreamNameS = fStreamName + "track3D.root";
            const char *nameC3DT = fStreamNameS.c_str(); //WARNING http://www.cplusplus.com/forum/general/100714/
            C3DT->SaveAs(nameC3DT);
            delete C3DT;
        }
  }
    
 cout<<"\n\nEnd of program reached, encountered "<<globalNumberErrors<<" non-terminal errors.\nThere were "<<globalNumberOfPUevents<<" rejected pile-up events.\nThere were "<<globalNumberOfLongEvents<<" rejected long events.\n";
 delete preTBF;
 if (saveHistograms){
     if (drawAuxillaryChannels){ for (unsigned int divNumb=0;divNumb<6;divNumb++){ delete timeBucketFormAuxillary[divNumb];}}
     if (drawTBFeventRemanent){ delete TBFeventRemanent;}
     if (realignTracks){  delete TBFrealigned;}
     if (drawDerivatives){ delete whatsTheDerivative;}
     if (drawLengthVpeakHeight){  delete lengthVSheight;}
     if (drawLengthVtrackCharge){  delete lengthVScharge;}
     if (drawLengthVbackgrd){  delete lengthVbackground;}
     if (drawLengthVtac){  delete lengthVtac;}
     if (drawChargeVtac){  delete chargeVtac;}
     if (drawPeakHeightVtac){  delete peakHeightVtac;}
     if (drawChargeVpeakHeight){  delete chargeVpeakHeight;}
//      if (addSignalCharacOnPlots){  delete mLi,m,nLi,n,oLi,o;}
     if (draw3DTrack){delete track3D;}
     gApplication->Terminate();
 }
 return 0;
  
}

// /afs/crc.nd.edu/user/m/mrenaud1/Public/stockageGraphAnal/    /home/mrenaud1/Documents/Li8Analysis/requestedGraphs/
