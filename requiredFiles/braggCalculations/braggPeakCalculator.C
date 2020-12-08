#include <stdlib.h>
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TBranch.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TF1.h"
#include "TF2.h"
#include "TMath.h"
#include <stdio.h>
#include <random>
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

using namespace std;


int braggPeakCalculator(){
//SELECT YOUR OPTIONS----------------------------------------
	Int_t gas = 4; // 0: pure Ar ; 1: P10 ; 2: 0.9*Ar+0.1*H2 ; 3: like 2, but looking at 9Be daughter; 4: looking at He+Ar->44Ca contamination; 5: looking at Li7; 6: looking at He6 BROKEN
	Double_t NbrNuclFinal, NbrNuclBeam;    //number of nucleons
	bool drawRatio = 0;                    //also draw the ratio between the signal and the beam at that depth.
	bool drawRatioFromMax = 1;             //as above, but only divide by maximum from 8Li B-peak
	bool theWholeRatio = 1;                //0: not all, stop after 48Sc bragg peak ; 1: all, go the full length
	
	if (drawRatio == 1 && drawRatioFromMax == 1){
		 cout<<"\n\nYou're asking for both 'drawRatio' and 'drawRatioFromMax'. They're mutually exclusive.\nPlease pick one.\n\n";return 0;
	}
	
	 //convert
	 double conver;                                // stopp pwr from MeV/(mg/cm2) -> MeV/mm
	 double converEl = 1./0.03;                    //stopp power in energy to stopp pwr in #e-
	 Double_t chargeElectron = 1.602176634*TMath::Power(10,-19);
    double beamCharge = 1.54867e-09;
	 

//------------------------------------------------------------
    
	 //DEFINITION OF HISTOGRAMS
	 
    //for the bragg peaks in (dE/dx)(E), inverted scale
    string ListTwo;  //post-fusion
	 if (gas==0) {ListTwo = "sc48Ar235torr.dat"; conver = 0.051498; NbrNuclFinal = 48; NbrNuclBeam = 8;}
	 else if (gas==1) {ListTwo = "sc48P10-235torr.dat"; conver = 0.048398; NbrNuclFinal = 48; NbrNuclBeam = 8;}
	 else if (gas==2) {ListTwo = "sc48ArH235torr2.dat"; conver = 0.050198; NbrNuclFinal = 48; NbrNuclBeam = 8;}
	 else if (gas==3) {ListTwo = "be9ArH235torr2.dat"; conver = 0.050198; NbrNuclFinal = 9; NbrNuclBeam = 8;}
	 else if (gas==4) {ListTwo = "ca44ArH235torr2.dat"; conver = 0.050198; NbrNuclFinal = 44; NbrNuclBeam = 4;}
	 else if (gas==5) {ListTwo = "sc47ArH235torr2.dat"; conver = 0.050198; NbrNuclFinal = 47; NbrNuclBeam = 7;}
	 else if (gas==6) {ListTwo = "he6ArH235torr2.dat"; conver = 0.050198; NbrNuclFinal = 6; NbrNuclBeam = 6;}
	 else {cout<<"\n\n\nUnexpected value for 'gas', please read comment next to the variable.\nEnding now.";return 0;}

	 ifstream List2(ListTwo);
	 double stoppPwrCharacteristics[90][3], stoppPwrCharacteristics2[90][2];
	 int io=0;
    TH1F* sc48BraggPeak = new TH1F("sc48 B-peak","stopping power;E[keV];dE/dx[keV/um]",2600,0,26000);
	 if (List2.is_open()){
		 while (List2>>stoppPwrCharacteristics[io][0]>>stoppPwrCharacteristics[io][1]>>stoppPwrCharacteristics[io][2]){io++;}  //stoppPwrCharacteristics[x][0] = energy in keV; stoppPwrCharacteristics[x][1] = electric stopping power dE/dx; stoppPwrCharacteristics[x][2] = nuclear stopping power dE/dx;
	 }
	 for (int i=0;i<90;i++){
		 sc48BraggPeak->SetBinContent((int)(stoppPwrCharacteristics[i][0]/10.),conver*(stoppPwrCharacteristics[i][1]+stoppPwrCharacteristics[i][2]));
	 }
	 

    string ListOne;  //pre-fusion
    if (gas==0) {ListOne = "li8Ar235torr.dat";}
	 else if (gas==1) {ListOne = "li8P10-235torr.dat";}
	 else if (gas==2 || gas==3) {ListOne = "li8ArH235torr2.dat";}
	 else if (gas==4){ListOne = "he4ArH235torr2.dat";}
	 else if (gas==5){ListOne = "li7ArH235torr2.dat";}
	 else if (gas==6){ListOne = "he6ArH235torr2.dat";}
	 ifstream List(ListOne);
	 TH1F* li8BraggPeak = new TH1F("Li8 B-peak","stopping power;E[keV];dE/dx[keV/um]",2600,0,26000);
	 io=0;
	 double trash;
    if (List.is_open()){
		 while (List>>trash>>stoppPwrCharacteristics2[io][0]>>stoppPwrCharacteristics2[io][1]){io++;}
	 }
	 for (int i=0;i<90;i++){
		 li8BraggPeak->SetBinContent((int)(stoppPwrCharacteristics[i][0]/10.),conver*(stoppPwrCharacteristics2[i][0]+stoppPwrCharacteristics2[i][1]));
	 }
	 
	 //WARNING: those function are drawn at the bottom of the program to avoid empty canvas during run
	
//--------------------------------------------------------
	 
	 //interactive part: creation of the signal
	 
	 TCanvas* c1=new TCanvas("c1","c1",800,800);
    c1->SetGridy();
	 c1->ToggleEventStatus();
	 
	 cout<<"\n\n\nWELCOME, STRANGER.\n";
	 if (gas==0){cout<<"The gas target is 40 Argon at 235 torr.\n";}
	 else if (gas ==1){cout<<"The gas target is P10 at 235 torr.\n";}
	 else if (gas ==2 || gas ==3 || gas ==4 || gas ==5){cout<<"The gas target is Ar 90% and H2 10% at 235 torr.\n";}
	 if (drawRatio==1){cout<<"The purple line (datRatio) is the ratio between the signal and the beam at that depth.\n\n";}
	 if (drawRatioFromMax==1){cout<<"The purple line (datRatio) is the ratio between the signal and the maximum from beam Bragg peak.\n\n";}
	 if (gas != 5){cout<<"The beam is 8 Lithium, at 25 MeV.\n";}
	 else if (gas == 5){cout<<"The beam is 7 Lithium, at 12 MeV.\n";}
	 if (gas != 3 && gas != 4){cout<<"(Fusion reaction unlikely to be seen under 6 MeV => <1mb)";}
	 cout<<"\nPlease select the energy at which ";
	 if (gas==3) {cout<<"8Li+1H -> 9Be (in keV): \n";}
	 else if (gas==4){cout<<"4He+40Ar -> 44Ca (in keV): \n";}
	 else if (gas==5){cout<<"7Li+40Ar -> 47Sc (in keV): \n";}
	 else if (gas==6){cout<<"8Li replaced by 6He, please enter 11 :\n";}
	 else {cout<<"8Li+40Ar -> 48Sc (in keV): \n";}
	 double whatEnergy; //user-defined energy at which fusion happens, in lab-frame
	 char doYou;
	 bool cont =0;
	 cin>>whatEnergy;
	 for (int i=1;i<89;i++){//exclude 10 kev to avoid prob with stoppPwrCharacteristics[i-1], exclude 25 MeV bcs out of scale
		 if (whatEnergy==stoppPwrCharacteristics[i][0]){whatEnergy=i;cont=1; break;}
		 else if (whatEnergy<=stoppPwrCharacteristics[i][0]){whatEnergy=i-1;break;}
	 }
	 if(cont!=1){
		 cout<<"\n\nFor calculation purposes, your input has been changed to : "<<stoppPwrCharacteristics[(int) whatEnergy][0]<<" keV, are you fine with that (y/n)?\n";
		 cin>>doYou;
		 if (doYou == 'n'){cout<<"\nQuitting program.\n";return 0;}
		 else{cout<<"\nPerfect.\n";}
	 }
	 
	 TH1F* depthSignal = new TH1F("Signal","Stopping power;dist[mm];#e-/mm",500,0,500);
	 TH1F* depthLithium = new TH1F("Lith","Ref Stopp Pwr;dist [mm]; #e-/mm",500,0,500);
	 TH1F* datRatio = new TH1F("Ratio","Signal ratio;dist [mm]; Sc/Li",500,0,500);
	 double x=0;
	 double stepX=0;
	 double elecYield, elecYieldRef, memo, generatedCharge;
	 double xref=0;
	 double stepXRef=0;
	 cont=0;
    generatedCharge = 0;
	 int fromWhatEnergy;
	 if (gas != 4 && gas != 5 && gas != 6){fromWhatEnergy = 88;} //don't start from end of array, i.e. 89, because the distance travelled in each step (stepX) is calculated using the difference between 'energy at array point j' - 'energy at array point j-1', so starting from the end you try to access entries that don't exist.
	 else if (gas == 4){fromWhatEnergy = 80;} //why you start so late for He I don't remember... ANSWER: because otherwise He just punches through, but they will during the run.
	 else if (gas == 5){fromWhatEnergy = 81;} //Li7 energy at about 12MeV
	 else if (gas == 6){fromWhatEnergy = 83;} //He6 energy at about 14MeV
	 for (int j = fromWhatEnergy;j>1;j--){
		 if (j>=whatEnergy){
			 x += stepX/2.;
			 stepX = (stoppPwrCharacteristics[j][0]-stoppPwrCharacteristics[j-1][0])/(conver*1000*(stoppPwrCharacteristics2[j][0]+stoppPwrCharacteristics2[j][1]));
			 elecYield = converEl*(conver*1000.*(stoppPwrCharacteristics2[j][0]+stoppPwrCharacteristics2[j][1]));
          generatedCharge += chargeElectron*converEl*conver*1000.*stoppPwrCharacteristics[j][0]*stepX;
			 x += stepX/2.;
			 depthSignal->SetBinContent((int) x, elecYield);
		 }
		 else if (io!=0){
			 x += stepX/2.;
			 if (cont!=1){
				 memo = x;
				 cont = 1;
				 for (int i=1; i<89;i++){//finding fusion energy's ID in array, then continue loop from that energy (accounting conserv of momentum)
					 if (NbrNuclBeam*stoppPwrCharacteristics[(int)whatEnergy][0]/NbrNuclFinal == stoppPwrCharacteristics[i][0]){io=i; break;}
					 else if (NbrNuclBeam*stoppPwrCharacteristics[(int)whatEnergy][0]/NbrNuclFinal <= stoppPwrCharacteristics[i][0]){io=i-1;break;}
				 }
				 
			 } //continue with reduced recoil energy
			 stepX = (stoppPwrCharacteristics[io][0]-stoppPwrCharacteristics[io-1][0])/(conver*1000.*(stoppPwrCharacteristics[io][1]+stoppPwrCharacteristics[io][2]));
			 elecYield = converEl*(conver*1000.*(stoppPwrCharacteristics[io][1]+stoppPwrCharacteristics[io][2]));
          generatedCharge += chargeElectron*converEl*conver*1000.*stoppPwrCharacteristics[j][0]*stepX;
			 x += stepX/2.;
			 depthSignal->SetBinContent((int) x, elecYield);
			 io--;
		 }
		 xref += stepXRef/2.;
		 stepXRef = (stoppPwrCharacteristics[j][0]-stoppPwrCharacteristics[j-1][0])/(conver*1000*(stoppPwrCharacteristics2[j][0]+stoppPwrCharacteristics2[j][1]));
		 elecYieldRef = converEl*(conver*1000.*(stoppPwrCharacteristics2[j][0]+stoppPwrCharacteristics2[j][1]));
		 xref += stepXRef/2.;
		 depthLithium->SetBinContent((int) xref, elecYieldRef);
	 }
	 
	 //print generated charge---
	 cout<<"\nGenerated electric charge (unamplified) : "<<generatedCharge<<" C.\n\tThat is "<<100.*generatedCharge/beamCharge<<" % of the beam charge\n";
    
    //print distance traveled--
	 cout<<"\nDistance traveled : "<<x<<" mm.\n";
    
    
	 //fill datRatio -----------
	 bool remanence = 0;
	 int binMax;
	 double sigNal, sigLith, maxLith;
	 binMax = depthLithium->GetMaximumBin();
	 maxLith = depthLithium->GetBinContent(binMax);
	 for (int depth = 35;depth<500;depth++){
		 if (theWholeRatio == 1){remanence = 1;}
		 if (depthSignal->GetBinContent(depth) != 0){sigNal = depthSignal->GetBinContent(depth);remanence = 1;}
		 if (drawRatio == 1){//draws the ratio of sigNal/sigLith, but does not take into account the different position of bragg peaks
			 if (depthLithium->GetBinContent(depth) != 0){sigLith = depthLithium->GetBinContent(depth);}
			 if (sigNal/sigLith !=0 && isnan(sigNal/sigLith) == false && remanence == 1){
				 datRatio->SetBinContent(depth,1000*sigNal/sigLith);
				 remanence = 0;
			}
		 }
		 if (sigNal/maxLith !=0 && isnan(sigNal/maxLith) == false){
			 datRatio->SetBinContent(depth,1000*sigNal/maxLith);
		 }
	 }
	 
	 //get and print max ratio -------
	 cout<<"\n\nMaximum ratio reached by datRatio : "<<datRatio->GetBinContent(datRatio->GetMaximumBin())/1000.<<endl<<endl;
// 	 return 0;
    
	 //set drawing options -----------
	 depthSignal->SetMarkerStyle(kFullCircle);
	 depthSignal->SetMarkerColor(3);
	 depthSignal->Draw("P");
	 depthLithium->SetMarkerStyle(kOpenSquare);
	 depthLithium->SetMarkerColor(2);
	 depthLithium->Draw("same P");
	 if (drawRatio == 1 || drawRatioFromMax == 1){
		 datRatio->SetMarkerStyle(kStar);
		 datRatio->SetMarkerColor(9);
		 datRatio->Draw("same P");
	 }
	 
	 
	 auto theLegend27 = new TLegend(0.1,0.8,0.55,0.9);
	 theLegend27->SetTextFont(342);
	 theLegend27->SetHeader(Form("Fusion happens at %f keV, or %f mm",stoppPwrCharacteristics[(int) whatEnergy][0],memo),"C");
	 theLegend27->AddEntry(depthSignal,"Signal with fusion","P");
	 theLegend27->AddEntry(depthLithium,"8Li bragg peak","P");
	 theLegend27->AddEntry(datRatio,"Signal/B-peak * 1000 (so 1000 = 1)","P");
	 theLegend27->Draw();
	 
	 
// 	 -----------------------------------------------------
	 //Draw the dE/dx functions for both
	 TCanvas* c=new TCanvas("c","c",800,800);
    c->SetGridy();
	 c->ToggleEventStatus();
	 sc48BraggPeak->SetMarkerStyle(kFullCircle);
	 sc48BraggPeak->SetMarkerColor(3);
	 sc48BraggPeak->Draw("P");
	 
	 li8BraggPeak->SetMarkerStyle(kFullSquare);
	 li8BraggPeak->SetMarkerColor(2);
	 li8BraggPeak->Draw("same P");
	 
// 	 -------------------------------------------------------
	 
	 //legend making
	 char* ss;
	 if (gas==3){ss="9Be";}
	 else if (gas==4){ss="44Ca";}
	 else if (gas==5){ss="47Sc";}
	 else {ss="48Sc";}
	 auto legend = new TLegend(0.1,0.8,0.4,0.9);
	 legend->AddEntry(sc48BraggPeak,Form("%s stopping power",ss),"P");
    legend->AddEntry(li8BraggPeak,"8Li stopping power","P");
	 legend->Draw();
	 
	 
// 	 --------------------------------------------------------

	 
	 return 0;
}
