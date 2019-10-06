/*
 Fred_moments_v3.c
 Sinead Farell, CPOM, 17/09/05

 What does code do?
 - Reads in first input data file containing 48 parameters: GLAS Tx_pulse((RXwf(k,j),k=1,48) 
 - Reads in second input data file containing 200 parameters: GLAS Rx_pulse((RXwf(k,j),k=1,200) 
 - Reverses GLAS Rx_pulse, so that it is in the correct direction.
 - Code calculates noise threshold for each Tx and Rx echo using the criteria:  mean_noise + 3Stdev_noise
 - The search criteria for which data is above nTH is: If Tx[j]-nTH>0 and Tx[j-1]-nTH>0 and Tx[j-2]-nTH>0 
   and Tx[j+1]-nTH>0 and Tx[j+2]-nTH>0, then I use Tx[j]-nTH.
 - Subtracts the noise level from each echo
 - Calculates waveform moments: 
   Calcultes the 1st,2nd,3rd,4th Moments of the GLAS Tx_pulse and Rx_pulse    
   (i.e. the Mean, Standard Deviation, Skewness, Kurtosis )
   Equations used to calculate waveform statistics are taken from:
   GLAS ATBD V4.0 "Derivation of range and range distributions from laser pulse waveform analysis for surface elevation, roughness, slope and vegetation
   heights", A.C. Brenner et al., NASA GSFC, Sept. 2003.
 
 Modifications:
 v1 10/10/05  (No longer a deconv component in this code). Also prints out Tx-nTHnorm and Rx-nTHnorm arrays as output files.
    11/10/05  Sorted out segmentation fault. problem was calculating the noise threshold (specifically "*x_mean_noise"). changed all variable in section to doubles.
 v2 10/12/05 Modified how we select the data above Threshold for Tx array (changed the loop min/max so it doesn't run outside the loop)
 v3 20/12/05 Modified the values for tmp* indices of the Tx_new and Rx_new arrays
   
*/ 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <unistd.h>

#define LEN1 48 /*the number of values in a line of Tx data  "j" */
#define LEN2 200 /*the number of values in a line of RX data "j" */
#define tmp LEN1-2 
#define tmp1 LEN1-1 
#define tmp2 LEN2-2 
#define tmp3 LEN2-1

//
//  Start main. 
//
int main(int argc, char *argv[])
{

  FILE *file1, *file2;
  char *filename1, *filename2;
  int count=0, count1=0, count2=0, count3=0, count4=0, count5=0;
  char dum[10000],dum1[1000000];
  int i,j,k;
  int *Tx,*Rx;
  double Tx_new[48];
  double Tx_norm[48];
  int Rx_rev[200];
  double Rx_new[200];
  double Rx_norm[200];
	  
/*check number of files entered*/
  if(argc < 3) 
    { printf("no input files given at command line\n");
      exit(0);    }

/*allocate memory to file names*/
  if((filename1 = (char *)malloc(255*sizeof(char))) == NULL)
    { printf("Memory allocation error A\n");
      exit(1);    }

 
  if((filename2 = (char *)malloc(255*sizeof(char))) == NULL)
    { printf("Memory allocation error B\n");
      exit(1);    }


/*read in data from two input files*/
  strcpy(filename1,argv[1]);
  if((file1=fopen(filename1,"r")) == NULL){ 
    printf("unable to read file1: %s\n ",filename1);
   }else{
    printf("Reading file : %s\n",filename1);
   }
   
  strcpy(filename2,argv[2]);
  if((file2=fopen(filename2,"r")) == NULL){ 
    printf("unable to read file2: %s\n ",filename2);
   }else{
    printf("Reading file : %s\n",filename2);
   }
 

/* Print out time at start of run to screen */ 
  time_t ttime;
  ttime=time(&ttime);
  fprintf(stderr,"\n             start time: %s\n",asctime(localtime(&ttime)));


/* Open the output files */
    FILE *fouta=fopen("TXmoments.dat","w");
    FILE *foutb=fopen("RXmoments.dat","w");
    FILE *foutc=fopen("TX-nTHnorm_array.dat","w");
    FILE *foutd=fopen("RX-nTHnorm_array.dat","w");


/* Find out the lenght of the TXpulse.dat and RXpulse.dat lines 
   i.e. find the EOF
   And check that they are the same length
   Call this lenght "count"					*/	

/* get each line in the file. It counts the no. of lines until either 10,000
   characters in the line are read or it reaches the end of the line */
  while (fgets(dum,10000,file1)) count++;
  rewind(file1);

  while (fgets(dum1,1000000,file2)) count1++;
  rewind(file2);
  
  printf("\n number of lines in %s = %d, number of lines in %s = %d\n",filename1,count,filename2,count1); 

/* Write a "infile_lines.flag" file if the number of lines in each file is not the same */  
  if(count!=count1) {
    FILE *foute=fopen("infile_lines.flag","w");
    fprintf(foute,"number of lines in %s = %d, number of lines in %s = %d\n",filename1,count,filename2,count1); 
    fclose(foute);
  }


fprintf(stderr,"-------------------------------------------------------------------------------------------------------\n");
fprintf(stderr,"Main Loop Begins here\n");


/* looping through each line of the input files */
    
    for(i=1; i<=count; i++){

//       printf("\n********************** \n"); 
//       printf("  Echo number=%d  \n",i); 
//       printf("********************** \n"); 


      	/*allocate dynamic memory to the input arrays*/
      	if((Tx = (int *)malloc(LEN1*sizeof(int))) == NULL)
      	  { printf("unable to allocate memory to Tx array\n");
	    exit(1); }

      	if((Rx = (int *)malloc(LEN2*sizeof(int))) == NULL)
      	  { printf("unable to allocate memory to Rx array\n");
	    exit(1); }

	//-----------------------------------------------------------------------------------------------------//
	/* working on the Tx array */
// 	printf("\n  Working on Tx array \n"); 
// 	printf(" --------------------- \n"); 
		
	/*scan in the Tx data array from input file1*/
      	for(j=0; j<LEN1; j++){ 
	  fscanf(file1,"%d ",Tx+j);
	}

	/* Calculating the Noise Threshold of Tx array*/
    	int l=0;
	double Tx_nvals=0.0;
    	double Tx_sumA=0.0;
	double Tx_sumB=0.0;
	double Tx_mean_noise=0.0,Tx_sd_noise=0.0,Tx_nTH=0.0;
	double Tx_maxA=0.0;
	
	/* I am using elements 0-9 for this calculation */
	for ( l=0; l<=9; l++) {
    	  Tx_sumA+=Tx[l];
    	  Tx_nvals++;
    	  Tx_mean_noise=(Tx_sumA/Tx_nvals);
	 }

      	for ( l=0; l<=9; l++) { 	
    	  Tx_sumB+=((Tx[l]-Tx_mean_noise)*(Tx[l]-Tx_mean_noise));
 	  Tx_sd_noise=(sqrt(Tx_sumB/Tx_nvals));
  	  Tx_nTH=(Tx_mean_noise+(3*Tx_sd_noise));
	}

// 	fprintf(stderr,"%8.5f elements used for noise calculation \n",Tx_nvals); 
//     	fprintf(stderr,"%8.5f sum of element values \n",Tx_sumA);
//     	fprintf(stderr,"%8.5f mean noise \n",Tx_mean_noise);
//     	fprintf(stderr,"%8.5f sumB \n",Tx_sumB);
// 	fprintf(stderr,"%8.5f standard deviation of noise \n",Tx_sd_noise);
// 	fprintf(stderr,"%8.5f noise threshold for Tx array \n",Tx_nTH);
 
 	/* Select data above the threshold */
//  	for(j=0; j<LEN1; j++) {
// 	/* search for data in bins where element value of bin-1 and of bin+1 is above 0 */
// 	  if ((Tx[j-2]-Tx_nTH)>0 && (Tx[j-1]-Tx_nTH)>0 && (Tx[j]-Tx_nTH)>0 && (Tx[j+1]-Tx_nTH)>0 && (Tx[j+2]-Tx_nTH)>0){
//     	    Tx_new[j]=(Tx[j]-Tx_nTH);
//     	  }else{    			
// 	    Tx_new[j]=0.0;
// 	  }
//  	}
 	for(j=2; j<(LEN1-2); j++) {
	/* search for data in bins where element value of bin-1 and of bin+1 is above 0 */
	  if ((Tx[j-2]-Tx_nTH)>0 && (Tx[j-1]-Tx_nTH)>0 && (Tx[j]-Tx_nTH)>0 && (Tx[j+1]-Tx_nTH)>0 && (Tx[j+2]-Tx_nTH)>0){
    	    Tx_new[j]=(Tx[j]-Tx_nTH);
    	  }else{    			
	    Tx_new[j]=0.0;
	  }
 	}
	Tx_new[0]=0.0;
	Tx_new[1]=0.0;
//	int tmp = (LEN1-2);
	Tx_new[tmp]=0.0;
//	int tmp1 = (LEN1-1);
	Tx_new[tmp1]=0.0;

	/* Find the maximum value in the Tx_new array */
     	for (j=0;j<LEN1;j++) {
          if (Tx_new[j]>Tx_maxA){
      	    Tx_maxA=Tx_new[j];
          }	
     	}
//     	fprintf(stderr,"Max value of Tx_new=%4.4f \n",Tx_maxA);


	/* Now we perform statistical calculations on data above the threshold (i.e. on array Tx_new) */
	/* calculating ATBD_defined MEAN */
	int Tx_nvalsB=0;
	double Tx_sumC=0.0,Tx_sumD=0.0,Tx_sumE=0.0,Tx_sumF=0.0,Tx_sumG=0.0;
	double Tx_mean=0.0;
	double Tx_prodA=0.0,Tx_prodB=0.0,Tx_prodC=0.0,Tx_prodD=0.0;
	double Tx_imm_sq=0.0,Tx_imm_cube=0.0,Tx_imm_four=0.0;
	double Tx_sd=0.0,Tx_sd_sq=0.0,Tx_sd_cube=0.0,Tx_sd_four=0.0;
	double Tx_skew=0.0,Tx_kurt=0.0;

	for(j=0; j<LEN1; j++) {
	  if (Tx_new[j]>0){
    	    Tx_sumC+=Tx_new[j];
    	    Tx_prodA=(j*Tx_new[j]);
	    Tx_sumD+=Tx_prodA;
	    Tx_nvalsB++;
    	    Tx_mean=(Tx_sumD/Tx_sumC);
	  }
	}

	/* calculating ATBD_defined STANDARD DEVIATION */
 	for(j=0; j<LEN1; j++) {
 	  if (Tx_new[j]>0){
    	    Tx_imm_sq=((j-Tx_mean)*(j-Tx_mean));
	    Tx_prodB=(Tx_imm_sq*Tx_new[j]);
	    Tx_sumE+=Tx_prodB;
	    Tx_sd_sq=(Tx_sumE/Tx_sumC);
	    Tx_sd=(sqrt(Tx_sumE/Tx_sumC));
	  }
	}

	/* calculating ATBD_defined SKEWNESS */
	for(j=0; j<LEN1; j++) {
	  if (Tx_new[j]>0){
     	    Tx_imm_cube=((j-Tx_mean)*(j-Tx_mean)*(j-Tx_mean));
	    Tx_prodC=(Tx_imm_cube*Tx_new[j]);
	    Tx_sumF+=Tx_prodC;
	    Tx_sd_cube=(Tx_sd_sq*Tx_sd);
	    Tx_skew=((1/Tx_sd_cube)*(Tx_sumF/Tx_sumC));
	  }
	}

	/* calculating ATBD_defined KURTOSIS */
	for(j=0; j<LEN1; j++) {
	  if (Tx_new[j]>0){
     	    Tx_imm_four=((j-Tx_mean)*(j-Tx_mean)*(j-Tx_mean)*(j-Tx_mean));
	    Tx_prodD=(Tx_imm_four*Tx_new[j]);
	    Tx_sumG+=Tx_prodD;
	    Tx_sd_four=(Tx_sd_cube*Tx_sd);
	    Tx_kurt=(((1/Tx_sd_four)*(Tx_sumG/Tx_sumC))-3);
	  }
	}

// 	/* Printing stats to screen */
// 	fprintf(stderr,"\n");
//     	fprintf(stderr,"%d values in Tx_Distribution above nTH \n",Tx_nvalsB); 
// 	fprintf(stderr,"\n");
//     	fprintf(stderr,"%8.5f sum of element values above threshold \n",Tx_sumC);
//  	fprintf(stderr,"%8.5f product of element value times element number for the last element in this sub-array \n",Tx_prodA);
//     	fprintf(stderr,"%8.5f sum of the products \n",Tx_sumD);
// 	fprintf(stderr,"%8.5f	 ATBD_defined mean \n",Tx_mean);
// 	fprintf(stderr,"\n");
//   	fprintf(stderr,"%8.5f square of element number minus ATBD_def_mean for the last element in this sub-array \n",Tx_imm_sq);
//     	fprintf(stderr,"%8.5f product of element value times i-mean_sq for the last element in this sub-array \n",Tx_prodB);
// 	fprintf(stderr,"%8.5f sum of the products \n",Tx_sumE);
// 	fprintf(stderr,"%8.5f square of the ATBD_defined StDev \n",Tx_sd_sq);
// 	fprintf(stderr,"%8.5f 	ATBD_defined StDev \n",Tx_sd);
// 	fprintf(stderr,"\n");
//   	fprintf(stderr,"%12.5f cube of element number minus ATBD_def_mean for the last element in this sub-array \n",Tx_imm_cube);
//     	fprintf(stderr,"%12.5f product of element value times i-mean_cube for the last element in this sub-array \n",Tx_prodC);
// 	fprintf(stderr,"%12.5f sum of the products \n",Tx_sumF);
// 	fprintf(stderr,"%12.5f cube of the ATBD_defined StDev \n",Tx_sd_cube);
// 	fprintf(stderr,"%12.5f 	ATBD_defined Skewness \n",Tx_skew);
// 	fprintf(stderr,"\n");
//    	fprintf(stderr,"%16.5f element number minus ATBD_def_mean to the four,for the last element in this sub-array \n",Tx_imm_four);
//     	fprintf(stderr,"%16.5f product of element value times i-mean_four for the last element in this sub-array \n",Tx_prodD);
// 	fprintf(stderr,"%16.5f sum of the products \n",Tx_sumG);
// 	fprintf(stderr,"%16.5f ATBD_defined StDev to the four \n",Tx_sd_four);
// 	fprintf(stderr,"%16.5f 	ATBD_defined Kurtosis \n",Tx_kurt);

	/*Normalise the new array */
     	double Tx_maxval=0.0;
	double Tx_sumH=0.0;
	int Tx_maxval_BinPos=0;
       	for (j=0;j<LEN1;j++) Tx_norm[j]=0.0;

     	for ( j=0; j<LEN1; j++) {
     	  if (Tx_maxA==0.0){
	    Tx_norm[j]=((double) Tx_new[j]);
	  }else{
     	    Tx_norm[j]=(((double) Tx_new[j])/((double) Tx_maxA));
	    Tx_sumH+=Tx_norm[j];
	  }
       	}
	
     	for ( j=0; j<LEN1; j++) {
          if (Tx_norm[j]>Tx_maxval){
      	    Tx_maxval=Tx_norm[j];
	    Tx_maxval_BinPos=j;
          }	
     	}
//     	fprintf(stderr,"Max value of Tx_norm=%4.4f, Bin Pos of Max Val=%d, Sum of Tx_norm=%4.4f \n",Tx_maxval,Tx_maxval_BinPos,Tx_sumH);

//  	/* Writing output to screen */
//   	for(j=0; j<LEN1; j++) {
//     	  fprintf(stderr,"Tx[%2i]:%4d,  noise threshold: %8.5f,   Tx_new[%2i]: %8.5f,   Tx_norm[%2i]: %8.5f \n",j,Tx[j],Tx_nTH,j,Tx_new[j],j,Tx_norm[j]);
//   	}

	/* writing output file A */
// 	/* output format: line no., Tx_norm_maxval, Tx_norm_maxval_BinPos, Tx_mean_noise, Tx_sd_noise, Tx_nTH,
// 	   Tx_nvals_abv_nTH, Tx_mean, Tx_sd, Tx_skew, Tx_kurt, Tx_norm[array] */
// 	fprintf(fouta,"%d %4.6f %d %4.6f %4.6f %4.6f %d %4.6f %4.6f %4.6f %4.6f",i,Tx_maxval,Tx_maxval_BinPos,Tx_mean_noise,Tx_sd_noise,Tx_nTH,Tx_nvalsB,Tx_mean,Tx_sd,Tx_skew,Tx_kurt);
// 	for(j=0; j<LEN1; j++) {
// 	  fprintf(fouta," %4.6f",Tx_norm[j]);
// 	}
// 	/* place a carriage return after printing each line of the output file */    
//     	fprintf(fouta,"\n");
	/* output format: line no., Tx_norm_maxval, Tx_norm_maxval_BinPos, Tx_mean_noise, Tx_sd_noise, Tx_nTH,
	   Tx_nvals_abv_nTH, Tx_mean, Tx_sd, Tx_skew, Tx_kurt */
	fprintf(fouta,"%d %4.6f %d %4.6f %4.6f %4.6f %d %4.6f %4.6f %4.6f %4.6f",i,Tx_maxval,Tx_maxval_BinPos,Tx_mean_noise,Tx_sd_noise,Tx_nTH,Tx_nvalsB,Tx_mean,Tx_sd,Tx_skew,Tx_kurt);
	/* place a carriage return after printing each line of the output file */    
    	fprintf(fouta,"\n");

	/* writing output file C */
	/* output format: Tx_norm[array] */
	for(j=0; j<LEN1; j++) {
	  fprintf(foutc," %4.6f",Tx_norm[j]);
	}
	/* place a carriage return after printing each line of the output file */    
    	fprintf(foutc,"\n");

	//-----------------------------------------------------------------------------------------------------//
	/* working on the Rx array */
// 	printf("\n  Working on Rx array \n"); 
// 	printf(" --------------------- \n"); 
		
	/*scan in the Rx data array from input file1*/
      	for(k=0; k<LEN2; k++){ 
	  fscanf(file2,"%d ",Rx+k);
	}

	/* Reverse the Rx array so that it is in the correct orientation (i.e increasing range runs with increasing bin # to the RH) */
	for(k=0; k<LEN2; k++){ 
    	  Rx_rev[k]=Rx[199-k];
	}
	//  	/* Writing output to screen */
	//   	for(k=0; k<LEN2; k++) {
	//     	  fprintf(stderr,"Rx[%2i]:%4d,  Rx_rev[%2i]: %4d \n",k,Rx[k],k,Rx_rev[k]);
	//   	}

	/* Calculating the Noise Threshold of Rx_rev array*/
    	int m=0;
	double Rx_nvals=0.0;
    	double Rx_sumA=0.0,Rx_sumB=0.0;
	double Rx_mean_noise=0.0,Rx_sd_noise=0.0,Rx_nTH=0.0;
	double Rx_maxA=0.0;
	
	/* I am using elements 0-49 for this calculation */
	for ( m=0; m<=49; m++) {
    	  Rx_sumA+=Rx_rev[m];
    	  Rx_nvals++;
    	  Rx_mean_noise=(Rx_sumA/Rx_nvals);
     	}

      	for ( m=0; m<=49; m++) { 	
    	  Rx_sumB+=((Rx_rev[m]-Rx_mean_noise)*(Rx_rev[m]-Rx_mean_noise));
	  Rx_sd_noise=(sqrt(Rx_sumB/Rx_nvals));
 	  Rx_nTH=(Rx_mean_noise+(3*Rx_sd_noise));
	}
	
// 	fprintf(stderr,"%8.5f elements used for noise calculation \n",Rx_nvals); 
//     	fprintf(stderr,"%8.5f sum of element values \n",Rx_sumA);
//     	fprintf(stderr,"%8.5f mean noise \n",Rx_mean_noise);
//     	fprintf(stderr,"%8.5f sumB \n",Rx_sumB);
// 	fprintf(stderr,"%8.5f standard deviation of noise \n",Rx_sd_noise);
// 	fprintf(stderr,"%8.5f noise threshold for Rx array \n",Rx_nTH);	
	
	
	/* Select data above the threshold */
//  	for(k=0; k<LEN2; k++) {
// 	/* search for data in bins where element value of bin-1 and of bin+1 is above 0 */
// 	  if ((Rx_rev[k-2]-Rx_nTH)>0 && (Rx_rev[k-1]-Rx_nTH)>0 && (Rx_rev[k]-Rx_nTH)>0 && (Rx_rev[k+1]-Rx_nTH)>0 && (Rx_rev[k+2]-Rx_nTH)>0){
//     	    Rx_new[k]=(Rx_rev[k]-Rx_nTH);
//     	  }else{    			
// 	    Rx_new[k]=0.0;
// 	  }
//  	}
 	for(k=2; k<(LEN2-2); k++) {
	/* search for data in bins where element value of bin-1 and of bin+1 is above 0 */
	  if ((Rx_rev[k-2]-Rx_nTH)>0 && (Rx_rev[k-1]-Rx_nTH)>0 && (Rx_rev[k]-Rx_nTH)>0 && (Rx_rev[k+1]-Rx_nTH)>0 && (Rx_rev[k+2]-Rx_nTH)>0){
    	    Rx_new[k]=(Rx_rev[k]-Rx_nTH);
    	  }else{    			
	    Rx_new[k]=0.0;
	  }
 	}
	Rx_new[0]=0.0;
	Rx_new[1]=0.0;
//	int tmp2 = (LEN2-2);
	Rx_new[tmp2]=0.0;
//	int tmp3 = (LEN2-1);
	Rx_new[tmp3]=0.0;

	/* Find the maximum value in the Rx_new array */
     	for (k=0;k<LEN2;k++) {
          if (Rx_new[k]>Rx_maxA){
      	    Rx_maxA=Rx_new[k];
          }	
     	}
//     	fprintf(stderr,"Max value of Rx_new=%4.4f \n",Rx_maxA);


	/* Now we perform statistical calculations on data above the threshold (i.e. on array Rx_new) */
	/* calculating ATBD_defined MEAN */
	int Rx_nvalsB=0;
	double Rx_sumC=0.0,Rx_sumD=0.0,Rx_sumE=0.0,Rx_sumF=0.0,Rx_sumG=0.0;
	double Rx_mean=0.0;
	double Rx_prodA=0.0,Rx_prodB=0.0,Rx_prodC=0.0,Rx_prodD=0.0;
	double Rx_imm_sq=0.0,Rx_imm_cube=0.0,Rx_imm_four=0.0;
	double Rx_sd=0.0,Rx_sd_sq=0.0,Rx_sd_cube=0.0,Rx_sd_four=0.0;
	double Rx_skew=0.0,Rx_kurt=0.0;

	for(k=0; k<LEN2; k++) {
	  if (Rx_new[k]>0){
    	    Rx_sumC+=Rx_new[k];
    	    Rx_prodA=(k*Rx_new[k]);
	    Rx_sumD+=Rx_prodA;
	    Rx_nvalsB++;
    	    Rx_mean=(Rx_sumD/Rx_sumC);
	  }
	}

	/* calculating ATBD_defined STANDARD DEVIATION */
 	for(k=0; k<LEN2; k++) {
 	  if (Rx_new[k]>0){
    	    Rx_imm_sq=((k-Rx_mean)*(k-Rx_mean));
	    Rx_prodB=(Rx_imm_sq*Rx_new[k]);
	    Rx_sumE+=Rx_prodB;
	    Rx_sd_sq=(Rx_sumE/Rx_sumC);
	    Rx_sd=(sqrt(Rx_sumE/Rx_sumC));
	  }
	}

	/* calculating ATBD_defined SKEWNESS */
	for(k=0; k<LEN2; k++) {
	  if (Rx_new[k]>0){
     	    Rx_imm_cube=((k-Rx_mean)*(k-Rx_mean)*(k-Rx_mean));
	    Rx_prodC=(Rx_imm_cube*Rx_new[k]);
	    Rx_sumF+=Rx_prodC;
	    Rx_sd_cube=(Rx_sd_sq*Rx_sd);
	    Rx_skew=((1/Rx_sd_cube)*(Rx_sumF/Rx_sumC));
	  }
	}

	/* calculating ATBD_defined KURTOSIS */
	for(k=0; k<LEN2; k++) {
	  if (Rx_new[k]>0){
     	    Rx_imm_four=((k-Rx_mean)*(k-Rx_mean)*(k-Rx_mean)*(k-Rx_mean));
	    Rx_prodD=(Rx_imm_four*Rx_new[k]);
	    Rx_sumG+=Rx_prodD;
	    Rx_sd_four=(Rx_sd_cube*Rx_sd);
	    Rx_kurt=(((1/Rx_sd_four)*(Rx_sumG/Rx_sumC))-3);
	  }
	}

// 	/* Printing stats to screen */
// 	fprintf(stderr,"\n");
//     	fprintf(stderr,"%d values in Rx_Distribution above nTH \n",Rx_nvalsB); 
// 	fprintf(stderr,"\n");
//     	fprintf(stderr,"%8.5f sum of element values above threshold \n",Rx_sumC);
//  	fprintf(stderr,"%8.5f product of element value times element number for the last element in this sub-array \n",Rx_prodA);
//     	fprintf(stderr,"%8.5f sum of the products \n",Rx_sumD);
// 	fprintf(stderr,"%8.5f	 ATBD_defined mean \n",Rx_mean);
// 	fprintf(stderr,"\n");
//   	fprintf(stderr,"%8.5f square of element number minus ATBD_def_mean for the last element in this sub-array \n",Rx_imm_sq);
//     	fprintf(stderr,"%8.5f product of element value times i-mean_sq for the last element in this sub-array \n",Rx_prodB);
// 	fprintf(stderr,"%8.5f sum of the products \n",Rx_sumE);
// 	fprintf(stderr,"%8.5f square of the ATBD_defined StDev \n",Rx_sd_sq);
// 	fprintf(stderr,"%8.5f 	ATBD_defined StDev \n",Rx_sd);
// 	fprintf(stderr,"\n");
//   	fprintf(stderr,"%12.5f cube of element number minus ATBD_def_mean for the last element in this sub-array \n",Rx_imm_cube);
//     	fprintf(stderr,"%12.5f product of element value times i-mean_cube for the last element in this sub-array \n",Rx_prodC);
// 	fprintf(stderr,"%12.5f sum of the products \n",Rx_sumF);
// 	fprintf(stderr,"%12.5f cube of the ATBD_defined StDev \n",Rx_sd_cube);
// 	fprintf(stderr,"%12.5f 	ATBD_defined Skewness \n",Rx_skew);
// 	fprintf(stderr,"\n");
//    	fprintf(stderr,"%16.5f element number minus ATBD_def_mean to the four,for the last element in this sub-array \n",Rx_imm_four);
//     	fprintf(stderr,"%16.5f product of element value times i-mean_four for the last element in this sub-array \n",Rx_prodD);
// 	fprintf(stderr,"%16.5f sum of the products \n",Rx_sumG);
// 	fprintf(stderr,"%16.5f ATBD_defined StDev to the four \n",Rx_sd_four);
// 	fprintf(stderr,"%16.5f 	ATBD_defined Kurtosis \n",Rx_kurt);

	/*Normalise the new array */
     	double Rx_maxval=0.0;
	double Rx_sumH=0.0;
	int Rx_maxval_BinPos=0;
       	for (k=0;k<LEN2;k++) Rx_norm[k]=0.0;

     	for ( k=0; k<LEN2; k++) {
     	  if (Rx_maxA==0.0){
	    Rx_norm[k]=((double) Rx_new[k]);
	  }else{
     	    Rx_norm[k]=(((double) Rx_new[k])/((double) Rx_maxA));
	    Rx_sumH+=Rx_norm[k];
	  }
       	}
	
     	for ( k=0; k<LEN2; k++) {
          if (Rx_norm[k]>Rx_maxval){
      	    Rx_maxval=Rx_norm[k];
	    Rx_maxval_BinPos=k;
          }	
     	}
//     	fprintf(stderr,"Max value of Rx_norm=%4.4f, Bin Pos of Max Val=%d, Sum of Rx_norm=%4.4f \n",Rx_maxval,Rx_maxval_BinPos,Rx_sumH);
 
//  	/* Writing output to screen */
//   	for(k=0; k<LEN2; k++) {
//     	  fprintf(stderr,"Rx[%2i]:%4d, Rx_rev[%2i]:%4d,  noise threshold: %8.5f,   Rx_new[%2i]: %8.5f,   Rx_norm[%2i]: %8.5f \n",k,Rx[k],k,Rx_rev[k],Rx_nTH,k,Rx_new[k],k,Rx_norm[k]);
//   	}

	/* writing output file B */
// 	/* output format: line no., Rx_norm_maxval, Rx_norm_maxval_BinPos, Rx_mean_noise, Rx_sd_noise, Rx_nTH,
// 	   Rx_nvals_abv_nTH, Rx_mean, Rx_sd, Rx_skew, Rx_kurt, Rx_norm[array] */
// 	fprintf(foutb,"%d %4.6f %d %4.6f %4.6f %4.6f %d %4.6f %4.6f %4.6f %4.6f",i,Rx_maxval,Rx_maxval_BinPos,Rx_mean_noise,Rx_sd_noise,Rx_nTH,Rx_nvalsB,Rx_mean,Rx_sd,Rx_skew,Rx_kurt);
// 	for(k=0; k<LEN2; k++) {
// 	  fprintf(foutb," %4.6f",Rx_norm[k]);
// 	}
// 	/* place a carriage return after printing each line of the output file */    
//     	fprintf(foutb,"\n");
	/* output format: line no., Rx_norm_maxval, Rx_norm_maxval_BinPos, Rx_mean_noise, Rx_sd_noise, Rx_nTH,
	   Rx_nvals_abv_nTH, Rx_mean, Rx_sd, Rx_skew, Rx_kurt */
	fprintf(foutb,"%d %4.6f %d %4.6f %4.6f %4.6f %d %4.6f %4.6f %4.6f %4.6f",i,Rx_maxval,Rx_maxval_BinPos,Rx_mean_noise,Rx_sd_noise,Rx_nTH,Rx_nvalsB,Rx_mean,Rx_sd,Rx_skew,Rx_kurt);
	/* place a carriage return after printing each line of the output file */    
    	fprintf(foutb,"\n");

	/* writing output file D */
	/* output format: line no.,  Rx_norm[array] */
	for(k=0; k<LEN2; k++) {
	  fprintf(foutd," %4.6f",Rx_norm[k]);
	}
	/* place a carriage return after printing each line of the output file */    
    	fprintf(foutd,"\n");


	//-----------------------------------------------------------------------------------------------------//
	/*free the Tx and Rx arrays so the next lines can be read*/
      	free(Tx);
       	free(Rx);
   }
     
fprintf(stderr,"-------------------------------------------------------------------------------------------------------\n");
fprintf(stderr,"Main Loop Ends here\n");

/* closing the output files having written to them*/
  fclose(fouta);
  fclose(foutb);
  fclose(foutc);
  fclose(foutd);

/* reopening the output files to count the number of lines in them - to read only! */
    fouta=fopen("TXmoments.dat","r");
    foutb=fopen("RXmoments.dat","r");
    foutc=fopen("TX-nTHnorm_array.dat","r");
    foutd=fopen("RX-nTHnorm_array.dat","r");


/* count each line in the out files. It counts the no. of lines until either 1,000,000
   characters in the line are read or it reaches the end of the line */
  while (fgets(dum1,1000000,fouta)) count2++;

  while (fgets(dum1,1000000,foutb)) count3++;

  while (fgets(dum1,1000000,foutc)) count4++;

  while (fgets(dum1,1000000,foutd)) count5++;

  printf("# lines in TXmoments.dat: %d, # lines in RXmoments.dat: %d, # lines in TX-nTHnorm_array.dat: %d, # lines in RX-nTHnorm_array.dat: %d\n",count2,count3,count4,count5); 

/* Write a "outfile_lines.flag" file if the number of lines in each file is not the same */  
  if(count2!=count3 || count2!=count4 || count2!=count5) {
    FILE *foutg=fopen("outfile_lines.flag","w");
    fprintf(foutg,"lines in TXmoments.dat: %d, lines in RXmoments.dat: %d, lines in TX-nTHnorm_array.dat: %d, lines in RX-nTHnorm_array.dat: %d\n",count2,count3,count4,count5); 
    fclose(foutg);
  }

/* closing the output files (again!) */
  fclose(fouta);
  fclose(foutb);
  fclose(foutc);
  fclose(foutd);


/*free the Tx and Rx arrays so the next lines can be read*/
  free(filename1);
  free(filename2);
  
/* Print out time at end of run to screen */
  ttime=time(&ttime);
  fprintf(stderr,"                 end time: %s\n",asctime(localtime(&ttime)));
   
}
