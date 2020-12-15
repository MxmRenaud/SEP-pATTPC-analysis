#include <iostream>
#include <fstream>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <cstdio>
#include <cstring>
using namespace std;

int makeList(){
    int startRun = 66;
    int endRun = 110;
    fstream outputList("list-analysis.txt", ios::out);
    
    for (int n=startRun;n<endRun+1;n++){
        if (n<10) {outputList <<"/home/mrenaud1/Public/converted/run_000"<<n<<".root"<<endl;}
        else if (10 <= n && n <100) {outputList <<"/home/mrenaud1/Public/converted/run_00"<<n<<".root"<<endl;}
        else if (100 <= n && n <1000) {outputList <<"/home/mrenaud1/Public/converted/run_0"<<n<<".root"<<endl;}
        
//         if (n<10) {outputList <<"/afs/crc/group/nsl/activetarget/data/TPC_Sepyember_2019/rootConverted/run_000"<<n<<".root"<<endl;}
//         else if (10 <= n && n <100) {outputList <<"/afs/crc/group/nsl/activetarget/data/TPC_Sepyember_2019/rootConverted/run_00"<<n<<".root"<<endl;}
//         else if (100 <= n && n <1000) {outputList <<"/afs/crc/group/nsl/activetarget/data/TPC_Sepyember_2019/rootConverted/run_0"<<n<<".root"<<endl;}
    }
    return 0;
}