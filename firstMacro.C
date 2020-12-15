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
* NOTE: because I'm an uncivilised Barbarian, this is supposed to be invoked as a root macro and run/compiled as such, not a regular compiled program. 
* I.e.: either 'root -b firstMacro.C' or 'root -b -L firstMacro.C++'
* 
* 
* NOTE HOWTO:at the beginning of the firstMacro(), a few self-explanatory boolians
* can be activated/deactivated depending on what you want the program to do.
*/

#include <stdlib.h>
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TBranch.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TFile.h"
#include "TF1.h"
#include "TMarker.h"
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

// #include "convolutionTA/convolution.c"  //NOTE: needs to be commented when not compiling


using namespace std;

Int_t globalNumberOfEvents;
Int_t globalNumberErrors;
Int_t globalNumberOfPUevents;
Int_t globalNumberOfLongEvents;
Int_t globalNumberOfTACchannelIDproblems;
Int_t globalNumberOfMassiveEvents;

Int_t whichAuxChanIsTAC(TH1I* histo, bool quiet){ //is the TAC in the auxillary channel 1 or 2 ?
    //locate the max and it's position: if in the back and over 3000, then it's MCP not TAC
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
        if (quiet == false) {cout<<"\n\n================ WARNING ! ================\nRevise your TAC-or-MCP conditions !\n================ WARNING ! ================\n\n";}
        globalNumberErrors++;
        globalNumberOfTACchannelIDproblems++;
        return 3;
    }
    else {return 1;}
}


//estimate pile-up by counting the number of TAC signals
Int_t pileUp(TH1I* histo){ 
    TF1 *fa0 = new TF1("fa0", "pol0(0)",0,511);
    histo->Fit(fa0,"0Q"); //options -> 0: not drawing, Q: quiet
    Float_t barem = fa0->GetParameter(0)*1.33;
    delete fa0;
    Int_t nbrOfEvent = 0;
    for (int i=0; i<511;i++){
       if (histo->GetBinContent(i) <= barem && histo->GetBinContent(i+1) > barem){nbrOfEvent++;}
    }
    return nbrOfEvent;
}


//locate peaks' positions in data to prepare ID-ing the one(s) that are signal
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
    //Sorting by Decreasing order
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

   
//gives a rough estimate of the derivative at ((startingPoint+numberOfBins)-startingPoint)/2
Float_t findDerivative(TH1F* histo_input, Double_t startingPoint, Int_t numberOfBins){
    if (numberOfBins > 0 && numberOfBins < 511/*histo_input->GetMaximumStored()*/){
        Float_t deri;
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




//=================================================================================================================================

Int_t firstMacro(){  
    
//END---------WHAT DO Y0U WISH TO DO HERE, MORTAL ? -------------------
    
    bool runOnTheCRC = true;               //deactivates mandatory user input & adapt file paths NOTE not optimal for the latter part.
    bool multiRunAnalysis = true;           //select between two lists of files (here 'list-preAnal.txt' and 'list-analysis.txt') for convenience in switching between debugging/actual analysis run NOTE: files can be commented in the list by putting an asterisk (*) in front
    
    bool fillIndividualChannels = false;    //if you want to look at a single event, channel by channel
    bool realignTracks = true;              //NOTE: BEST TRACK LENGTH FINDING ALGORYTHM IS LOCKED IN HERE; draws sum of all channels, summed for every event
    bool acceptPeaksBelowBeam = false;      // do you want to keep events that peak much further than the primary beam ? see 'pointOfPeakRejection'. Will exclude most events that are just noise and some disintegration events.
    bool useFusionGates = false;            //attempt to only keep events fulfilling certain conditions [hard-coded atm] TODO change that
    bool usePadsGates = true;              //only keep events with less than "maxNbrOfPadsFired" pads fired (NOTE not including 6 auxillary channels !)
    bool substractBckgrd = true;
    bool beQuiet = true;                   //reduced terminal output
    bool addSignalCharacOnPlots = true;     //plot expected values for fusion/contaminants on some of the data plots.
    bool trackBizarTACevents = false;       //Create a separate file listing events whose TAC channel behaviour are unexpected (i.e., generate a globalNumberOfTACchannelIDproblems error)
    bool deconvoluteThis = false;           //deconvolute data before any plotting WARNING currently not working
    
    bool drawWhichPadIsIt = false;
    bool drawIndividualChannels = false;
    bool drawSummedChannels = false;        //sums all channels of a single event into a track and draw.
    bool drawAuxillaryChannels = false;
    bool drawTBFeventRemanent = false;       //drawTimeBucketFormatEventRemanent = draw a summed profile of all tracks for all events
    bool drawDerivatives = false;           //NOTE : Hard-coded for only one event at a time
    bool drawLengthVpeakHeight = true;
    bool drawLengthVtrackCharge = true;
    bool drawLengthVbackgrd = true;
    bool drawLengthVtac = true;
    bool drawChargeVtac = true;
    bool drawPeakHeightVtac = true;
    bool drawChargeVpeakHeight = true;
    bool draw3DTrack = false;               //NOTE : Hard-coded for only one event at a time
    bool drawDoubleProjectTrack = false; 
    bool drawNbrOfPadsVcharge = true; 
    bool drawNbrOfPadsVpeakHeight = true;
    bool drawNbrOfPadsVtac = true;
    
    bool saveHistograms = true;            //NOTE WARNING: currently required to properly free up the memory used. Only falsefy when debugging
    
//END---THE CONTRACT HAS BEEN SEALED, MORTAL. I SHALL NOW PROCEED...-----
    
    
    
  const Int_t derivativeCalcWindowSize = 5;
  const Int_t estimatedWidth = 9.;              //estimated width of the track's beginning rising edge at high TimeBucket
  const Int_t pointOfPeakRejection = 75;
  const Int_t EVENT = 0;    //run_0090/!\ 11602;  // 14030->60; /!\ 14030 & 14041 & 14048 & 14053(noDipBelow); /!\ 14031 & 14032 & 14052 (1TAC2Peak) 14078 (1TAC3peaks); /!\ 4041 & 14056 (weirdMCP) /!\ 14030 & 14053 & 14059 (singleShort)
  char eRrOr;
  const Float_t meanBeamPeakHeight = 187.5;     //EDIT on 25 NOV 2020, prev.: 256.5; EDIT on 04 FEB 2020, prev.: 253.5;
  const Float_t meanBeamCharge = 7959;          //EDIT on 25 NOV 2020, prev.: 12095.5; EDIT on 04 FEB 2020, prev.: 11259.5;
  const Float_t driftV = 20.;                   //mm/us
  const Float_t timeBucketSize = 0.160;         //us
  const Float_t Li7ChargePercentage = 0.2537;   //charge dep by Li7 as percentage Li8 beam charge
  const Float_t Li7peakHeight = 0.9998;         //bragg peak height of Li7 as percentage Li8 bragg peak height
  const Float_t Li7avTrLenght = 57.5;           //Li7 bragg peak position in Time Buckets
  const int lengthOfConvolArray = 256;          //length of the (de)convolution array
  const Int_t maxNbrOfPadsFired = 8;            //user-defined limit to keep only events with straight tracks (no scattering, diffusion, decay, etc)


  
  //--ERRORS AND WARNINGS--------------------------------------------------
  if (fillIndividualChannels == false && drawIndividualChannels == true){cout<<"\n~~~~~~~~~~~~~~~~~~~~~WARNING !~~~~~~~~~~~~~~~~~~~\n\nYou need to read the channels to draw them.\nGo activate 'fillIndividualChannels' and try again.\n\n";return 0;}
  if (realignTracks == true && drawTBFeventRemanent == true && runOnTheCRC == false){
      cout<<"\nWARNING ! You asked for both the tracks to be re-aligned before drawing (realignTracks) AND for the tracks to be drawn without being re-aligned (drawTBFeventRemanent).\nDo you wish to proceed ? (y/n) : ";
      cin>>eRrOr;
      if (eRrOr == 'n') {return 0;}
      else if (eRrOr == 'y') {cout<<"\nProceeding..."<<endl;}
      else {cout<<"\n... you think you're being funny ?\n\n\n\tYou. Aren't.\n\n";return 0;}
  }
  if (useFusionGates == true && (drawLengthVpeakHeight == false || drawLengthVtrackCharge == false || realignTracks == false)){
      cout<<"\nWARNING\nFor now  please activate drawLengthVtrackCharge, drawLengthVpeakHeight & realignTracks to run useFusionGates.\n\n";
      return 0;
  }
  if (drawChargeVtac == true && drawLengthVtrackCharge == false){
      cout<<"\nWARNING ! 'drawChargeVtac == true && drawLengthVtrackCharge == false' is not a supported combo.\nPlease correct.\n\nThank you kindly.\n";
      return 0;
  }
  if (drawPeakHeightVtac == true && drawLengthVpeakHeight == false){
      cout<<"\nWARNING ! Asking for 'drawPeakHeightVtac == true && drawLengthVpeakHeight == false' is not a supported combo.\nPlease correct.\n\nThank you kindly.\n";
      return 0;
  }
  if ((draw3DTrack == true || drawDoubleProjectTrack == true) && realignTracks == false){
      cout<<"\nWell met, stranger. I hereby inform you that the options 'draw3DTrack' and 'drawDoubleProjectTrack' require option 'realignTracks' to be enabled, in order to be executed.\nThis message will now terminate your run.\n\nOur deepest apologies.\n";
      return 0;
  }
 //------------------------------------------------------------------------



  //--CALIBRATION ELEMENTS--------
  
  //Go fetch the mapping matrix
  const int numberOfPads = 2021;
  float mapp[numberOfPads][7];
  const char *fMapStreamChar;
  string fStreamName;
  int line = 0;
  if (runOnTheCRC == true){fStreamName = "/afs/crc.nd.edu/user/m/mrenaud1/Public/SEP-pATTPC-analysis/requiredFiles/TPC-SEP2019-fullMap.dat";}
  else{fStreamName = "requiredFiles/TPC-SEP2019-fullMap.dat";}
  fMapStreamChar = fStreamName.c_str(); //WARNING http://www.cplusplus.com/forum/general/100714/
  fstream mapStream(fMapStreamChar);
  if(mapStream.is_open()){//structure Cobo<< "\t" << Asad<< "\t" << Aget<< "\t" << Channel<< "\t" << PadNum << "\t" << xPos(mm) << "\t" << yPos(mm)
      while(mapStream>>mapp[line][0]>>mapp[line][1]>>mapp[line][2]>>mapp[line][3]>>mapp[line][4]>>mapp[line][5]>>mapp[line][6]){
          if (line == numberOfPads){
          cout<<"\n\nWARNING !\nThere are more lines in the 'mapStream' than space in the 'mapp' matrix.\nThis should be fixed.\n";
          return 0;
          }
          line++;
      }
  }
  else{
   cout<<"\nWarning ! Failed to open TPC-SEP2019-fullMap.dat, please start panicking at your earliest convenience.\n";
   return 0;
  }
  mapStream.close();
  
  const int numberOfEnergies = 16;
  float fusionCharac[numberOfEnergies][4];
  const char *fCharacStreamChar;
  if (runOnTheCRC == true){fStreamName = "/afs/crc.nd.edu/user/m/mrenaud1/Public/SEP-pATTPC-analysis/requiredFiles/signal-Li8Ar40to48Sc.dat";}
  else{fStreamName = "requiredFiles/signal-Li8Ar40to48Sc.dat";}
  fCharacStreamChar = fStreamName.c_str(); //WARNING http://www.cplusplus.com/forum/general/100714/
  fstream characStream(fCharacStreamChar);
  line = 0;
  if ((useFusionGates == true || addSignalCharacOnPlots == true) && characStream.is_open()){//energy[kev], charge[% of BeamE], peakHeight[fraction of beam peak], dist travelled [estim., mm]
      while (characStream>>fusionCharac[line][0]>>fusionCharac[line][1]>>fusionCharac[line][2]>>fusionCharac[line][3]){
          if (line == numberOfEnergies){
              cout<<"\n\nWARNING !\nThere are more lines in the 'characStream' than space in the 'fusionCharac' matrix.\nThis should be fixed.\n";
              return 0;
          }
          line++;
      }
  }
  else if (useFusionGates == true && characStream.is_open() == false){
      cout<<"\nWarning ! Failed to open signal-Li8Ar40To48Sc.dat, please start panicking at your earliest convenience.\n";
      return 0;
  }
  characStream.close();
  


  
  //Definition of objects
  const char *fListChar;
  if (runOnTheCRC == true){
	if (multiRunAnalysis == false) {fStreamName = "/afs/crc.nd.edu/user/m/mrenaud1/Public/SEP-pATTPC-analysis/requiredFiles/list-preAnal.txt";}
	else {fStreamName = "/afs/crc.nd.edu/user/m/mrenaud1/Public/SEP-pATTPC-analysis/requiredFiles/list-analysis.txt";}
  }
  else{
	if (multiRunAnalysis == false) {fStreamName = "requiredFiles/list-preAnal.txt";}
	else {fStreamName = "requiredFiles/list-analysis.txt";}
  }
  fListChar = fStreamName.c_str(); //WARNING http://www.cplusplus.com/forum/general/100714/
  ifstream List(fListChar);
  
  const char *fListTACproblems;
  if (runOnTheCRC == true){fStreamName = "/afs/crc.nd.edu/user/m/mrenaud1/Public/SEP-pATTPC-analysis/listOfTACproblemEvents.txt";}//create text file that will contain the run and event number for when 'globalNumberOfTACchannelIDproblems' is emitted
  else {fStreamName = "listOfTACproblemEvents.txt";}
  fListTACproblems = fStreamName.c_str(); //WARNING http://www.cplusplus.com/forum/general/100714/
  ofstream ListTAC(fListTACproblems);
  if (ListTAC.is_open() == false){
      cout<<"\nWarning ! Failed to create or open listOfTACproblemEvents.txt, please start grumbling at your earliest convenience.\n";
      return 0;
  }

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
  TH2F* track2Dfront[10];
  TH1F* trackProjection[10];
  TH2F* nbrOfPadsVtac;
  TH2F* nbrOfPadsVcharge;
  TH2F* nbrOfPadsVpeakHeight;
  TH1F* deconvolutedProfile;
  
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
  if (drawDoubleProjectTrack){
      for (int i=0;i<10;i++){
          track2Dfront[i] = new TH2F(Form("track2Dfront #%i",i),Form("Track #%i pad use;x [mm];y [mm]",i),244,-121,121,244,-121,121);
          track2Dfront[i]->Fill(0.,0.);
          trackProjection[i] = new TH1F(Form("Track #%i profile",i),Form("Track #%i profile;time bucket;Energy[arb.]",i),511,0,510);
      }
  }
  if (drawNbrOfPadsVcharge == true && substractBckgrd == true) nbrOfPadsVcharge = new TH2F("padsVcharge","#pads V. dep. charge;#fired pads;dep. charge [arb]",150,0,149,2500,0,25000);
  if (drawNbrOfPadsVcharge == true && substractBckgrd == false) nbrOfPadsVcharge = new TH2F("padsVcharge","#pads V. dep. charge;#fired pads;dep. charge [arb]",150,0,149,4500,0,90000);
  if (drawNbrOfPadsVpeakHeight == true && substractBckgrd == true) nbrOfPadsVpeakHeight = new TH2F("padsVpeakHeight","#pads V. peak Height;#fired pads;max peak height [arb]",150,0,149,700,0,700);
  if (drawNbrOfPadsVpeakHeight == true && substractBckgrd == false) nbrOfPadsVpeakHeight = new TH2F("padsVpeakHeight","#pads V. peak Height;#fired pads;max peak height [arb]",150,0,149,1100,0,1100);
  if (drawNbrOfPadsVtac) nbrOfPadsVtac = new TH2F("padsVtac","#pads V. TAC;#fired pads;tac [arb]",150,0,149,850,0,850);
  if (deconvoluteThis) deconvolutedProfile = new TH1F("deconvolutedProfile","Deconv. profile;time bucket;Energy[arb]",512,0,511);

  

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
  
  //   TH1F* timeBucketForm = new TH1F("TBF", "channel 1;time bucket;Energy[arb.]",512,0,511);
*/

  //--DEFINITION OF (DE)CONVOLUTION USTENSILS-------------------------------------
  float anArray1[lengthOfConvolArray],anArray3[lengthOfConvolArray], dataTest[lengthOfConvolArray], answerTest[2*lengthOfConvolArray];
  if (deconvoluteThis){
      for (int i = 0;i<lengthOfConvolArray;i++){
//           dataTest[i] = 2.0*TMath::Gaus(i,91,6);
          anArray1[i] = i;
          anArray3[i] = i;
//           anArray3[lengthOfConvolArray+i] = lengthOfConvolArray+i;
      }
  }
  
/*  const float *theArray1 = anArray1;
  const float *theArray3 = anArray3;
  const float *dataTest2 = dataTest;
  TGraph* grD = new TGraph(lengthOfConvolArray,anArray1,dataTest2);
  
  cout<<endl;
  convolution(dataTest,answerTest,1);
  
  const float *answerTest2 = answerTest;
  TGraph* grA = new TGraph(lengthOfConvolArray,anArray1,answerTest2);
  
  TCanvas* Ctestc = new TCanvas("Ctestc","Ctestc",600,600);
  grD->SetLineColor(kBlue);
  grD->Draw("AC");
  grA->SetLineColor(kGreen);
  grA->Draw("same");
  
  return 0;*/
  //--END DEFINITION OF (DE)CONVOLUTION USTENSILS---------------------------------
  
  
  
  //--DEFINITION OF VARIABLES FOR TREE /!\ variable have to match for root extraction, storing Int_t->Double_t won't fly
  Int_t num_hits;
  Int_t whatChan;
  Int_t data[2055][517];
  globalNumberOfEvents = 0;
  globalNumberErrors = 0;
  globalNumberOfPUevents = 0;
  globalNumberOfLongEvents = 0;
  globalNumberOfTACchannelIDproblems = 0;
  globalNumberOfMassiveEvents = 0;
  bool done = 0;
  Int_t auxChanDone = 0;
  bool checkCleared = 0;                //WARNING when 'no PileUp' condition gets relaxed, this check will be a problem
  Int_t doubleCheckPU = 0;              //stores calculated pile-up value for double checking event validity
  Double_t* daPeak;
  Int_t tacPeak;
  Double_t TACval = 0;
  Int_t daPeakElement = 0;
  Float_t storedDerivative[100], derivativeAverage[100];
  Int_t howManyInDerivArray = 0;
  Float_t maxDerivativePosition[2];
  Float_t temp[3];                      //variables for random use
  bool fusionTestValidated = 0;
  Float_t padsPosition[100][512][2];
  int checkTrack3D;
  bool did3DtrackGetMade = false;
  int projectionTrackLimit = 0;

  

  while (List >> WhichFile){
	
   FileName = WhichFile.c_str(); //TFile needs conversion to char
   if (FileName[0] != '*'){//NOTE filenames commented by adding a '*' in front won't be read
   cout<<"\nOpening next file, "<<FileName<<", processing...\n";
   
        
	//Reading & Initialising
	file = new TFile(FileName);
	if (file->IsOpen()){
    tree = (TTree*)(file->Get("T"));  
    tree->SetBranchAddress("br_len", &num_hits);
    tree->SetBranchAddress("myints", data);
    
    //NOTE IMPORTANT! structure: Cobo<< "\t" << Asad<< "\t" << Aget<< "\t" << Channel<< "\t" << PadNum
    
    //Recursively go through tree, getting "Fill" after "Fill"
    for (Int_t i=EVENT; i</*EVENT+10000*/tree->GetEntries(); i++){
        
        globalNumberOfEvents++;
        
        tree->GetEntry(i);
        if (beQuiet != 1) cout<<"\nEvent #"<<i<<", num_hits : "<<num_hits<<endl<<endl;
        if (num_hits-6 > 100){cout<<"\nWARNING, event "<<i<<" has too many entries.";globalNumberErrors++;globalNumberOfMassiveEvents++;}
        
        //resetting variables/histograms
        auxChanDone = 0;
        for (int i=0;i<512;i++){
            timeBucketFormAveraged->SetBinContent(i,0.);
            whatsTheDerivative->SetBinContent(i,0.);
            preTBF->SetBinContent(i,0.);
            if (projectionTrackLimit < 10){
                for (int j=0;j<100;j++){
                    padsPosition[j][i][0] = 0;padsPosition[j][i][1] = 0;
                }
            }
            if (i<100) {storedDerivative[i] = 0; derivativeAverage[i] = 0;}
            if (i<3) {temp[i] = 0;}
            if (i<2) {maxDerivativePosition[i] = 0;}
        }
        checkCleared = 0;
        howManyInDerivArray = 0; //how many non-zero entries in the array containing the derivative of the track
        doubleCheckPU = 0;
        daPeakElement = 0;
        fusionTestValidated = 0;
        TACval = 0;

        
        //loop through channels in event
        for (Int_t j=0; j<num_hits; j++) { //NOTE: num_hits = number of stored channels + 5 or 6 auxillary channels
            
            done = 0;
            
            for (Int_t k=0; k<517; k++){
                
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
                }// \end ID what the channel is
                if (done==1 && k>4 && data[j][0] != 2){//fill in time bucket format for non-auxillary channels
                    if (fillIndividualChannels == true && j<100){timeBucketForm[j]->SetBinContent(k-5,data[j][k]);}
                    timeBucketFormAveraged->AddBinContent(k-5,(float)(data[j][k])/max(num_hits-6,1)); //WARNING diff gains inner-outer
                }  
                if (done==1 && k>4 && data[j][0] == 2 && auxChanDone<6){ //fill in time bucket format for auxillary channels
                    timeBucketFormAuxillary[auxChanDone]->SetBinContent(k-5,data[j][k]);
                    if(k==516){auxChanDone++;} //move to next histogram
                }
                else if(auxChanDone >= 6){cout<<"\nError, which is sad.\nAnd 'sad' backwards is 'das', and das not good...\nLook at the Auxillaries.\n";if (trackBizarTACevents == true){ListTAC<<WhichFile<<"\t"<<i<<endl;} return 0;}
                if (k>5 && projectionTrackLimit <10 && (draw3DTrack == true || drawDoubleProjectTrack == true) /*&& did3DtrackGetMade == false*/ && num_hits-6 <100){
                    padsPosition[j][k-4][0] = mapp[checkTrack3D][5];
                    padsPosition[j][k-4][1] = mapp[checkTrack3D][6];
                }
            }
        }
        
        //section to see if event has pile-up & attempt clean-up/renorm
        num_hits = num_hits-(auxChanDone+1); //remove the auxillary channels from counting the number of fired pads
        temp[0] = whichAuxChanIsTAC(timeBucketFormAuxillary[1], beQuiet);
        if (temp[0] == 3){ //keeps track of events with problematic TAC or MCP in separate file
            if (trackBizarTACevents == true){ ListTAC<<WhichFile<<"\t"<<i<<endl;}
            temp[0] = 2;
        }
        Int_t pileUpInEvent = pileUp(timeBucketFormAuxillary[(int)temp[0]]);
        if (beQuiet != 1)cout<<"TAC is in aux. chan. "<<temp[0]<<endl;
        if (beQuiet != 1)cout<<"#ofEvents : "<<pileUpInEvent;
        if (pileUpInEvent == 1){ if (beQuiet != 1) cout<<"\t-> Accepted"<<endl;}
        else {if (beQuiet != 1){cout<<"\t-> Rejected"<<endl<<endl;} globalNumberOfPUevents++;}
        
        if (done == 1 && pileUpInEvent == 1){//if at end of event, and event was clean && no pile-up, perform last check; then add the average to the signal
            preTBF->Add(timeBucketFormAveraged);
            daPeak = locateTACsignal(timeBucketFormAveraged, pileUpInEvent+2); //stores position of peaks found in the averaged TimeBucket-stored output
            if (beQuiet != 1) cout<<"Watch : "<<daPeak[0]<<"\t"<<daPeak[1]<<"\t"<<daPeak[2]<<endl;
            for (int j=0;j<pileUpInEvent+2;j++){
                if (daPeak[j] >=1 && daPeak[j] <= 513){doubleCheckPU++;}
                else if (daPeak[j] >= 513){daPeakElement++;}//if you're here the real peak you want is not the first, so go fetch next entry, not the first "peak" -> daPeak[daPeakElement] instead of daPeak[0]
            }
            
            //Double-check pileUp
            if (pileUpInEvent == doubleCheckPU){
                if (acceptPeaksBelowBeam == true){checkCleared = 1;}
                else if (acceptPeaksBelowBeam == false && daPeak[daPeakElement] > pointOfPeakRejection){checkCleared =1;}
                else {if (beQuiet != 1) {cout<<"\n\tWARNING ! Peak beyond rejection point.\t-> Rejected"<<endl<<endl;} globalNumberErrors++;globalNumberOfLongEvents++;}
            }
            else if (doubleCheckPU == 0){if (beQuiet != 1) {cout<<"\n\tWARNING ! No signal !\t-> Rejected"<<endl<<endl;} globalNumberErrors++;}
            else{
                if (beQuiet != 1){cout<<"\n\tWARNING ! Multiple peaks for one TAC signal !\t-> Rejected"<<endl<<endl;}
                globalNumberErrors++;
                globalNumberOfPUevents++;
            }
            
            //Get TAC value
            tacPeak = timeBucketFormAuxillary[(int)temp[0]]->GetMaximumBin();
            TACval = timeBucketFormAuxillary[(int)temp[0]]->GetBinContent(tacPeak) - timeBucketFormAuxillary[(int)temp[0]]->GetBinContent(tacPeak-6); //get peak height, substract baseline
            if (beQuiet != 1) cout<<"TAC position : "<<tacPeak<<"\tTAC value : "<<TACval<<endl;
            
            //check number of fired Pads
            if (usePadsGates == true && (num_hits > maxNbrOfPadsFired)){
                checkCleared = 0;
                if (beQuiet != 1){cout<<"\n\tWARNING ! Too many fired pads : "<<num_hits<<" > "<<maxNbrOfPadsFired<<".\t-> Rejected"<<endl<<endl;}
            }
            
            if (checkCleared == 1){//actual clean & treatement
                temp[0]=0;
                //fit noise before track and substract it.
                if (beQuiet != 1) cout<<"location of peaks : "<<daPeak[daPeakElement]<<endl;
                TF1 *fa0 = new TF1("fa0", "pol0(0)",5,daPeak[daPeakElement]-20);
                timeBucketFormAveraged->Fit(fa0,"RQ"); //fit options: no draw, limited range, quiet
                if (drawLengthVbackgrd) temp[2] = fa0->GetParameter(0);
                if (beQuiet != 1) cout<<"fa0->GetParameter(0) = "<<fa0->GetParameter(0)<<endl;
                fa0->SetRange(0,511);
                if (substractBckgrd){preTBF->Add(fa0,-1,"");}
                preTBF->SetBinContent(0,0);
                
                //locate start of track and substract after track noise
                for (int bin=0;bin<(511.-daPeak[daPeakElement])/5.;bin++){//use the derivative of the track signal to help ID the end (sharp drop)
                    storedDerivative[bin] = findDerivative(preTBF,daPeak[daPeakElement]+bin*5,derivativeCalcWindowSize);
                    whatsTheDerivative->SetBinContent(daPeak[daPeakElement]+bin*5,storedDerivative[bin]);
                    howManyInDerivArray++;
                }
                fa0->SetRange(daPeak[daPeakElement]+10,505);
                preTBF->Fit(fa0,"RQ");
                if (drawLengthVbackgrd == true && substractBckgrd == true) temp[2] += fa0->GetParameter(0);
                else if (drawLengthVbackgrd == true && substractBckgrd == false) temp[2] = fa0->GetParameter(0);
                if (beQuiet != 1) cout<<"fa0->GetParameter(0), v2 = "<<fa0->GetParameter(0)<<endl;
                for (int j=daPeak[daPeakElement];j<511;j++){
                    if (preTBF->GetBinContent(j) <= max(fa0->GetParameter(0)*1.5,7.)){
                        for (int u=0; u<howManyInDerivArray-3; u++){
                            derivativeAverage[u]=AveragE(storedDerivative,3,3+u); //TODO optimise this section ! too many loops !
                        }
                        break;
                    }
                }
                
                for (int p = howManyInDerivArray-3;p>0;p--){// <- start from end of array
                    if (derivativeAverage[p] < -1.01){ // <- gives start of track
                        if (beQuiet != 1) cout<<"tentative : "<<p<<"\t"<<daPeak[daPeakElement]+(p+1)*5<<endl; //p +1(for stepping) 
                        fa0->SetRange(daPeak[daPeakElement]+(p+1)*5+10,505); // <- select plateau at the end of the time bucket form
                        preTBF->Fit(fa0,"RQ");
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
                            if (beQuiet != 1) cout<<"\nmaxDerivativePosition[0] = "<<maxDerivativePosition[0]<<", maxDerivativePosition[1] = "<<maxDerivativePosition[1]<<endl; //respectively position of start and end of track
                            
                            
                            if (drawLengthVbackgrd) lengthVbackground->Fill(maxDerivativePosition[0]+estimatedWidth-maxDerivativePosition[1],temp[2]); //temp[2] holds the background/baseline
                            
                            if (drawLengthVpeakHeight){
                                if (useFusionGates == false){
                                    lengthVSheight->Fill((maxDerivativePosition[0]+estimatedWidth-maxDerivativePosition[1]),preTBF->GetMaximum());
                                    if (drawPeakHeightVtac){peakHeightVtac->Fill(preTBF->GetMaximum(),TACval);}
                                }
                                else {//NOTE: Only keep events compatible with theoretical peak height. 
                                    for (line=15;line>=0;line--){
                                        if ((maxDerivativePosition[0]+estimatedWidth-maxDerivativePosition[1]) <= fusionCharac[line][3]/(driftV*timeBucketSize)){
                                            fusionTestValidated = 1;
                                            temp[0] = preTBF->Integral(maxDerivativePosition[1],maxDerivativePosition[0]+estimatedWidth,"");
                                            break;
                                        }
                                        globalNumberOfLongEvents++;
                                    }
                                    //WARNING hard-coded most conditions/gates. Current condition : <=92.5% of calculated charge
                                    if (fusionTestValidated == true && temp[0]<92.545455*(maxDerivativePosition[0]+estimatedWidth-maxDerivativePosition[1])-755.81819){
                                        lengthVSheight->Fill(maxDerivativePosition[0]+estimatedWidth-maxDerivativePosition[1],preTBF->GetMaximum());
                                    }
                                    if (drawPeakHeightVtac){peakHeightVtac->Fill(preTBF->GetMaximum(),TACval);}
                                }
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
                                did3DtrackGetMade = true;
                                draw3DTrack = false;
                            }
                            else if (draw3DTrack == true && num_hits >= 100){
                                if (!beQuiet) cout<<"\nWARNING ! Event incompatible with 'draw3DTrack' option. Dropping option.\n";
                                draw3DTrack = false; //TODO change to allow for next track to be selected
                            }
                            
                            if (drawDoubleProjectTrack==true && projectionTrackLimit <10){
                                for (int i=0;i<num_hits;i++){
                                    for (int j=maxDerivativePosition[1];j<maxDerivativePosition[0]+estimatedWidth;j++){
                                        track2Dfront[projectionTrackLimit]->Fill(padsPosition[i][j][0],padsPosition[i][j][1]);
                                    }
                                }
                                trackProjection[projectionTrackLimit]->Add(preTBF);
                                projectionTrackLimit++;
                            }
                            
                            if (drawNbrOfPadsVcharge){ nbrOfPadsVcharge->Fill(num_hits,preTBF->Integral/*AndError*/(maxDerivativePosition[1],maxDerivativePosition[0]+estimatedWidth,""));}
                            if (drawNbrOfPadsVtac){ nbrOfPadsVtac->Fill(num_hits,TACval);}
                            if (drawNbrOfPadsVpeakHeight){ nbrOfPadsVpeakHeight->Fill(num_hits,preTBF->GetMaximum());}
                            
                            if (deconvoluteThis){
                                for (int i=maxDerivativePosition[1];i<maxDerivativePosition[0]+estimatedWidth-maxDerivativePosition[1];i++){
                                    
                                }
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
                        if (realignTracks == false && drawNbrOfPadsVcharge == true){nbrOfPadsVcharge->Fill(num_hits,preTBF->Integral/*AndError*/(daPeak[daPeakElement]-10,daPeak[daPeakElement]+(p+1)*5,""));}
                        if (realignTracks == false && drawNbrOfPadsVpeakHeight == true){nbrOfPadsVpeakHeight->Fill(num_hits,preTBF->GetMaximum());}
                        if (realignTracks == false && drawNbrOfPadsVtac == true){nbrOfPadsVtac->Fill(num_hits,TACval);}

                        
                        p = 1; //break out of the "for (int p = howManyInDerivArray-3;p>0;p--)" loop
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
   else {cout << "\n\nWARNING ! FAILED TO OPEN FILE "<<WhichFile<<", ABORTING.\nYou should go find out what happened.\n\n";break;return 0;}
   
   
   file->Close();
   delete file;
   
   
   }//ENDFileName not '*'-ed
  }//ENDwhile loop
  
  List.close();
  ListTAC.close();
  
  
  
  
// -----------------------------DEF OF CANVAS AND PLOT---------------------------------------

  //Defining style of graph
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(1111111);//now showing integral under-and-overflow ! //old: 1001110
  gStyle->SetOptFit(1);
  gStyle->ToggleEventStatus();
  
  //save path definition
  string fStreamNameS;
  if (runOnTheCRC == false){fStreamName = "/home/maxime/Documents/Assistanat-KUL/Halo_Reactions-thesis/ExpNDSEP/Li8analysis/requestedGraphs/autoSaved-padTest-";}
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
        CauxChan->Divide(3,3);
        CauxChan->ToggleEventStatus();
        for (unsigned int divNumb=1;divNumb<7;divNumb++){
            CauxChan->cd(divNumb);
            timeBucketFormAuxillary[divNumb-1]->Draw();
        }
	if (saveHistograms){
	    fStreamNameS = fStreamName + "timeBucketFormAuxillary.root";
	    const char *nameCauxChan = fStreamNameS.c_str(); //WARNING http://www.cplusplus.com/forum/general/100714
	    CauxChan->SaveAs(nameCauxChan);
	    delete CauxChan;
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
  
  if (did3DtrackGetMade){
      TCanvas* C3DT = new TCanvas("C3DT","C3DT",800,600);
        track3D->Draw();
        if (saveHistograms){
            fStreamNameS = fStreamName + "track3D.root";
            const char *nameC3DT = fStreamNameS.c_str(); //WARNING http://www.cplusplus.com/forum/general/100714/
            C3DT->SaveAs(nameC3DT);
            delete C3DT;
        }
  }
  
  if (drawDoubleProjectTrack){
      TCanvas* CDPT = new TCanvas("CDPT","CDPT",1000,1000);
        CDPT->Divide(4,5); 
        CDPT->ToggleEventStatus();
        for (unsigned int divNumb=1;divNumb<11;divNumb++){
            CDPT->cd(divNumb);
            track2Dfront[divNumb-1]->Draw("colz");
            CDPT->cd(divNumb+10);
            trackProjection[divNumb-1]->Draw();
        }
        if (saveHistograms){
            fStreamNameS = fStreamName + "10Tracks.root";
            const char *nameCDPT = fStreamNameS.c_str(); //WARNING http://www.cplusplus.com/forum/general/100714/
            CDPT->SaveAs(nameCDPT);
            delete CDPT;
        }
  }
  
  if (drawNbrOfPadsVcharge){
      TCanvas* CNbrVQ = new TCanvas("CNbrVQ","CNbrVQ",600,600);
        CNbrVQ->SetLogz();
        nbrOfPadsVcharge->Draw("colz");
        if (saveHistograms){
            fStreamNameS = fStreamName + "nbrOfPadsVcharge.root";
            const char *nameCNbrVQ = fStreamNameS.c_str(); //WARNING http://www.cplusplus.com/forum/general/100714/
            CNbrVQ->SaveAs(nameCNbrVQ);
            delete CNbrVQ;
        }
  }
  
  if (drawNbrOfPadsVtac){
      TCanvas* CNbrVTAC = new TCanvas("CNbrVTAC","CNbrVTAC",600,600);
        CNbrVTAC->SetLogz();
        nbrOfPadsVtac->Draw("colz");
        if (saveHistograms){
            fStreamNameS = fStreamName + "nbrOfPadsVtac.root";
            const char *nameCNbrVTAC = fStreamNameS.c_str(); //WARNING http://www.cplusplus.com/forum/general/100714/
            CNbrVTAC->SaveAs(nameCNbrVTAC);
            delete CNbrVTAC;
        }
  }
  
  if (drawNbrOfPadsVpeakHeight){
      TCanvas* CNbrVPH = new TCanvas("CNbrVPH","CNbrVPH",600,600);
        CNbrVPH->SetLogz();
        nbrOfPadsVpeakHeight->Draw("colz");
        if (saveHistograms){
            fStreamNameS = fStreamName + "nbrOfPadsVpeakHeight.root";
            const char *nameCNbrVPH = fStreamNameS.c_str(); //WARNING http://www.cplusplus.com/forum/general/100714/
            CNbrVPH->SaveAs(nameCNbrVPH);
            delete CNbrVPH;
        }
  }
    
 cout<<"\n\nEnd of program reached, handled "<<globalNumberOfEvents<<" events.\nEncountered "<<globalNumberErrors<<" non-terminal errors ("<<100.*globalNumberErrors/globalNumberOfEvents<<"% of total).\nThere were "<<globalNumberOfPUevents<<" rejected pile-up events ("<<100.*globalNumberOfPUevents/globalNumberOfEvents<<"% of total).\nThere were "<<globalNumberOfLongEvents<<" rejected long events ("<<100.*globalNumberOfLongEvents/globalNumberOfEvents<<"% of total).\nThere were "<<globalNumberOfTACchannelIDproblems<<" events with unexpected TAC channel behaviour ("<<100.*globalNumberOfTACchannelIDproblems/globalNumberOfEvents<<"% of total).\nThere were "<<globalNumberOfMassiveEvents<<" events with over a 100 fired channels ("<<100.*globalNumberOfMassiveEvents/globalNumberOfEvents<<"% of total).\n";
 
 
 if (saveHistograms){
     //TODO do you need the 'if's here ? shouldn't every object get deleted regardless ?
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
     if (did3DtrackGetMade){delete track3D;}
     if (drawDoubleProjectTrack){for (unsigned int divNumb=0;divNumb<10;divNumb++){delete track2Dfront[divNumb]; delete trackProjection[divNumb];}}
     if (drawNbrOfPadsVcharge){delete nbrOfPadsVcharge;}
     if (drawNbrOfPadsVpeakHeight){delete nbrOfPadsVpeakHeight;}
     if (drawNbrOfPadsVtac){delete nbrOfPadsVtac;}
     delete preTBF;
     gApplication->Terminate();
 }
 return 0;
  
}

// /afs/crc.nd.edu/user/m/mrenaud1/Public/stockageGraphAnal/  
