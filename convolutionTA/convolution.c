/*
* ================================================================
*
*	Filename:    convolution.c
*	
*	Description: convolution and de-convolution program. WARNING: compilation mandatory !
*	             input isign = 1/-1 for convolution/deconvolution
*	             arguments should be flaot data[LENGTH], isign, float ans[2*LENGTH], so code accordingly
*
*	Modified:    17 FEB 2019
*               
*	Author:      Tan AHN, NSH, UND
*	Modifier:    Maxime RENAUD, NSH, KUL-UND
*
* ================================================================
* NOTE if run independently, one needs to create binary <foo>.o for all five of the "include files", and comment those "include".
*/ 



#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "convlv.c"
#include "four1.c"
#include "nrutil.c"
#include "realft.c"
#include "twofft.c"

#define LENGTH 256

// int main(){                                                  //COMMENT ME WHEN RUNNING WITH 'firstMacro.C'
int convolution(float data[], float ans[], int isign = 1){      //UN-COMMENT ME WHEN RUNNING WITH 'firstMacro.C'    
  
  unsigned int i;
  unsigned long m, n;
  float norm, sigma, data_sum, conv_sum;
//   int isign = -1;                                             //COMMENT ME WHEN RUNNING WITH 'firstMacro.C'
//   float data[LENGTH], ans[2*LENGTH];                         //COMMENT ME WHEN RUNNING WITH 'firstMacro.C'
  float respns[LENGTH];

  FILE *data_out, *respns_out, *ans_out;

  void convlv(float data[], unsigned long n, float respns[], unsigned long m, int isign, float ans[]);

  n = LENGTH;
  m = 79;

  data_sum = 0.0;
  for(i=0; i<LENGTH; i++){
//     data[i] = 2.0*exp(-pow((i-91.0)/6.0,2)/2.0);       //60/5 80/5 5/6?! 41/5 41/5 41/6 40/6  //COMMENT ME WHEN RUNNING WITH 'firstMacro.C'
    data_sum += data[i];                            //90/5-256 91/6-256 100/5-256
  }

  sigma = 2.6539;//1us of shaping time, 160ns TB -> 1us = FWHM = 2.355*sigma -> sigma = 6.25 Time buckets/2.355 = 2.6539   //old value: 10.0;

  for(i=0; i<(m-1)/2; i++){
    respns[i] = 1.0/(sigma*sqrt(2.0*M_PI))*exp(-pow(i/sigma,2));
  }
  for(i=m-1; i>=(m-1)/2+1; i--){
    respns[i] = 1.0/(sigma*sqrt(2.0*M_PI))*exp(-pow(-(i-m)/sigma,2));
  }
  
//   cout<<"\nNotice me Senpai !\n";

  convlv(data-1, n, respns-1, m, isign, ans-1);
  
//   cout<<"\nKIAAA ! Senpai noticed me !\n";

  conv_sum = 0.0;
  for(i=0; i<LENGTH; i++){
    conv_sum += ans[i];
  }

  norm = data_sum/conv_sum;

  data_out = fopen("exp1.dat", "w");
  respns_out = fopen("exp2.dat", "w");
  ans_out = fopen("exp3.dat", "w");

  for(i=0; i<LENGTH; i++){
    fprintf(data_out, "%d %f\n", i, data[i]);
    fprintf(respns_out, "%d %f\n", i, respns[i]);
    ans[i] = norm*ans[i];
    fprintf(ans_out, "%d %f\n", i, ans[i]);
  }

  fclose(data_out);
  fclose(respns_out);
  fclose(ans_out);

  return 0;
}
