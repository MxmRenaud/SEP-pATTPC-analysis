/*==========================================================
 * 
 * FileName:     map.C
 * 
 * Description:  Take both .csv files of the pATTPC pad mapping, and creat a single .dat FileName
 * 
 * Created:      1/10(OCT)/2019
 * Author:       Maxime Renaud
 * 
 * =========================================================
 */
#include <iostream>
#include <fstream>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <cstdio>
#include <cstring>
using namespace std;

// // fstream& GotoLine(std::fstream& file, unsigned int num){
// //     file.seekg(std::ios::beg);
// //     for(int i=0; i < num - 1; ++i){
// //         file.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
// //     }
// //     return file;
// // }

int main(){
	
	ifstream channelMap("FlatLookupND2019.dat"); //tried with the .csv originals, got VERY freaky output. Attempt with caution.
	ifstream spaceMap("PadPositionsND2019.dat");
	
	fstream outputMap("TPC-SEP2019-fullMap.dat", ios::out);
	
	int Cobo, Asad, Aget, Channel, PadNum, counting = 0;
	double posX, posY;

	//start reading; store PadNum(ber) from first stream, as this is the line we need to read in the second
	if(channelMap.is_open()){
		while(channelMap>>Cobo>>Asad>>Aget>>Channel>>PadNum){
			if(PadNum < 2015 && spaceMap.is_open()){ //avoid trying to read more than file size
				while(spaceMap>>posX>>posY){
					if (counting == PadNum){//stop and reset
						counting =0;
						spaceMap.clear(); //those two lines to make sure you start reading from the beginning again
						spaceMap.seekg(0, ios::beg);
						break;
					}
					counting++;
				}
			}
			else if (spaceMap.is_open() == false){cout <<"\nWe're closed. Come back later.\n"; return 0;}
			else {cout << "\nEr is een probleempje met PadNum, vriend.\nJammer, heh ?\n";return 0;}
			outputMap << Cobo<< "\t" << Asad<< "\t" << Aget<< "\t" << Channel<< "\t" << PadNum<< "\t" << posX<< "\t" << posY <<endl;
		}
	}
	else if (channelMap.is_open() ==  false){cout<<"\nYou have an itty-bitty file problem, me friend\n";}
	return 0;
	
}