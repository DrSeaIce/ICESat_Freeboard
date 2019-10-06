/*
 Fred_xcorrel_v3.c
 Sinead Farell, CPOM, 10/10/05

 What does code do?
 - Reads in first input data file containing 48 parameters: GLAS TX-nTHnorm.dat (48 element array) 
 - Reads in second input data file containing 200 parameters: GLAS RX-nTHnorm.dat (200 element array)
 - The arrays have had their noise removed using Fred_v1.c (criteria:  mean_noise + 3Stdev_noise)
 - The X-correlation between the Tx pulse and the Rx pulse is calculated.
 - The maximum value of the Xcorrel is found and it's corresponding lag (will show how the Tx and Rx peaks are offset).
 
 
 Modifications:
 v1 10/10/05  No longer a deconv component in this code.
     	      Also prints out the following files: xcorrel_TX_RX_array.dat, xcorrel_TX_RX_rar_array.dat,
	       xcorrel_TX_RX_stats.dat, TX_padded_array.dat, RX_padded_array.dat. 
  	      Now padding out arrays to size 402 before calculating the X-correlation. This means the peaks of Tx and Rx can be 
	      at 200 and there will still be enough room for the rest of the Rx array without having an array over-run!!! 
 v2 21/12/05  Sorting out the array declarations and array sizes
 v3 04/01/06  Finding bug in v2 since it crashes when run with NAG error 2 (due to incorrect set up of NAG command line)
              Also changing the amount of data in the xcorrel_xy and xcorrel_xy_rar arrays that is printed to output files 
              Also no longer outputting the Tx and Rx padded arrays 
 
*/ 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <unistd.h>

#define LEN1 48 /*the number of values in a line of Tx data  "j" */
#define LEN2 200 /*the number of values in a line of RX data "j" */
#define N 402  /*size of Txp, Rxp, Xcorrel_xy, Xcorrel_xy_rar arrays */
#define NL 401  /*size of OUTxy array */

/* declare the fortran library NAG routines here */
/* G13BCF:
   calculates the cross correlation between two time series */
extern void 	g13bcf_();		

//
//  Start main. 
//
int main(int argc, char *argv[])
{
  FILE *file1, *file2;
  char *filename1, *filename2;
  int count=0, count1=0;
  char dum[10000],dum1[100000];
  int i,j,k;
  double *Tx,*Rx;
  double *Txp,*Rxp;
  double *OUTxy;

  int n_new=0,nl_new=0,n2=0,n2m1=0;
  /* Set and print local scalars */
  n_new = N;
  nl_new = NL;
  n2 = (N/2);
  n2m1 = (n2-1);
  fprintf(stderr,"\n N=%d,  n_new=%d,  n2=%d,  n2m1=%d,	NL=%d,  nl_new=%d \n",N,n_new,n2,n2m1,NL,nl_new);

  double Xcorrel_xy[402];
  double Xcorrel_xy_rar[402];

  /*check number of files entered*/
  if(argc < 3){
    printf("no input file given at command line\n");
    exit(0);
  }
      
  /*allocate memory to file names*/
  if((filename1 = (char *)malloc(255*sizeof(char))) == NULL){
    printf("Memory allocation error A\n");
    exit(1);
  }
 
  if((filename2 = (char *)malloc(255*sizeof(char))) == NULL){
    printf("Memory allocation error B\n");
    exit(1);
  }

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
  //  FILE *fouta=fopen("TX_padded_array.dat","w");
  //  FILE *foutb=fopen("RX_padded_array.dat","w");
  FILE *foutc=fopen("xcorrel_TX_RX_array.dat","w");
  FILE *foutd=fopen("xcorrel_TX_RX_rar_array.dat","w");
  FILE *foute=fopen("xcorrel_TX_RX_stats.dat","w");
  
  /* Find out the length of the TXpulse.dat and RXpulse.dat lines 
   i.e. find the EOF
   And check that they are the same length
   Call this lenght "count"					*/	
  /* get each line in the file. It counts the no. of lines until either 10,000
  characters in the line are read or it reaches the end of the line */
  while (fgets(dum,10000,file1)) count++;
  rewind(file1);
  while (fgets(dum1,100000,file2)) count1++;
  rewind(file2);
  printf("\n number of lines in %s = %d, number of lines in %s = %d\n",filename1,count,filename2,count1); 


  /* Write a "infile_lines.flag" file if the number of lines in each file is not the same */  
  if(count!=count1) {
    FILE *foutz=fopen("infile_lines.flag","w");
    fprintf(foutz,"number of lines in %s = %d, number of lines in %s = %d\n",filename1,count,filename2,count1); 
    fclose(foutz);
  }


//  fprintf(stderr,"-------------------------------------------------------------------------------------------------------\n");
//  fprintf(stderr,"Main Loop Begins here\n");


  /* looping through each line of the input files */    
  for(i=1; i<=count; i++){
    //    printf("\n********************** \n"); 
    //    printf("  Echo number=%d  \n",i); 
    //    printf("********************** \n"); 
    
    /*allocate dynamic memory to the arrays*/
    if((Tx = (double *)malloc(LEN1*sizeof(double))) == NULL){
      printf("unable to allocate memory to Tx array\n");
      exit(1);
    }

    if((Rx = (double *)malloc(LEN2*sizeof(double))) == NULL){
      printf("unable to allocate memory to Rx array\n");
      exit(1);
    }

    if((Txp = (double *)malloc(N*sizeof(double))) == NULL){
      printf("unable to allocate memory to Txp array\n");
      exit(1);
    }
    //  for (j=0;j<N;j++) Txp[j]=0.0;
    
    if((Rxp = (double *)malloc(N*sizeof(double))) == NULL){
      printf("unable to allocate memory to Rxp array\n");
      exit(1);
    }
    //  for (j=0;j<N;j++) Rxp[j]=0.0;

    if((OUTxy = (double *)malloc(NL*sizeof(double))) == NULL){
      printf("unable to allocate memory to OUTxy array\n");
      exit(1);
    }
    for (j=0;j<NL;j++) OUTxy[j]=0.0;
    
    //-----------------------------------------------------------------------------------------------------//
    /* working on the Tx array */
    //    printf("\n  Working on Tx array \n"); 
    //    printf(" --------------------- \n"); 
	
    /*scan in the Tx data array from input file1*/
    for(j=0; j<LEN1; j++){ 
      fscanf(file1," %le",Tx+j);
    }

    /*Find the max value of the Tx array and it's bin position */
    double Tx_maxval=0.0;
    int Tx_maxval_BinPos=0;
    double tx_sum=0.0, txp_sum=0.0;
    for ( j=0; j<LEN1; j++) {
      if (Tx[j]>Tx_maxval){
	Tx_maxval=Tx[j];
	Tx_maxval_BinPos=j;
      }	
    }

    //    /*Writing output to screen*/
    //    for(j=0; j<LEN1; j++) {
    //      fprintf(stderr,"Tx[%3d]: %10.6f \n",j,Tx[j]);
    //    }
    //    fprintf(stderr,"Max value of Tx=%10.6f, Bin Pos of Max Val=%d \n",Tx_maxval,Tx_maxval_BinPos);
    //    printf("\n"); 

    /*padding out Tx array to 402-element array and placing peak at bin #200 */
    for(j=0; j<N; j++) {
      if ( j>=(200-Tx_maxval_BinPos) && j<=(200+(47-Tx_maxval_BinPos)) ) {
	Txp[j] = Tx[ j-(200-Tx_maxval_BinPos) ];
      }else{
	Txp[j]=0.0;
      }
      //    fprintf(stderr,"Txp[%3d]: %10.6f \n",j,Txp[j]);
    }

    /* checksum: check the sum of both the Tx and Txp arrays are the same */
    for(j=0; j<LEN1; j++) {
      tx_sum+=Tx[j];
    }
    for(j=0; j<N; j++) {
      txp_sum+=Txp[j];
    }
    //    fprintf(stderr,"Tx_sum: %10.6f	Txp_sum: %10.6f \n",tx_sum,txp_sum);
    if ( tx_sum != txp_sum ) {
      fprintf(stderr,"!!!!!!! Array sums Differ !!!!!!! \n");
      fprintf(stderr,"Tx_sum: %10.6f	Txp_sum: %10.6f \n",tx_sum,txp_sum);
      break;
    }

    //-----------------------------------------------------------------------------------------------------//
    /* working on the Rx array */
    //    printf("\n  Working on Rx array \n"); 
    //    printf(" --------------------- \n"); 

    /*scan in the Rx data array from input file1*/
    for(j=0; j<LEN2; j++){ 
      fscanf(file2," %le",Rx+j);
    }

    /*Find the max value of the Rx array and it's bin position */
    double Rx_maxval=0.0;
    int Rx_maxval_BinPos=0;
    double rx_sum=0.0, rxp_sum=0.0;
	
    for ( j=0; j<LEN2; j++) {
      if (Rx[j]>Rx_maxval){
	Rx_maxval=Rx[j];
	Rx_maxval_BinPos=j;
      }	
    }

    //    /* Writing output to screen */
    //    for(j=0; j<LEN2; j++) {
    //      fprintf(stderr,"Rx[%3d]: %10.6f \n",j,Rx[j]);
    //    }
    //    fprintf(stderr,"Max value of Rx=%10.6f, Bin Pos of Max Val=%d \n",Rx_maxval,Rx_maxval_BinPos);
    //    printf("\n"); 

    /*padding out Rx array to 402-element array and placing peak at bin #200 */
    for(j=0; j<N; j++) {
      if ( j>=(200-Rx_maxval_BinPos) && j<=(200+(199-Rx_maxval_BinPos)) ) {
	Rxp[j] = Rx[j-(200-Rx_maxval_BinPos)];
      }else{
	Rxp[j]=0.0;
      }
      //    fprintf(stderr,"Rxp[%3d]: %10.6f \n",j,Rxp[j]);
    }
	
    /* checksum: check the sum of both the Rx and Rxp arrays are the same */
    for(j=0; j<LEN2; j++) {
      rx_sum+=Rx[j];
    }
    for(j=0; j<N; j++) {
      rxp_sum+=Rxp[j];
    }
    //    fprintf(stderr,"Rx_sum: %10.6f	Rxp_sum: %10.6f \n",rx_sum,rxp_sum);
    if ( rx_sum != rxp_sum ) {
      fprintf(stderr,"!!!!!!! Array sums Differ !!!!!!! \n");
      fprintf(stderr,"Rx_sum: %10.6f	Rxp_sum: %10.6f \n",rx_sum,rxp_sum);
      break;
    }

    //-----------------------------------------------------------------------------------------------------//
    /* calculating the x-correl between Txp and Rxp */
    double SDxy=0.0, OUT0xy=0.0, statxy=0.0;
    int IFAIL=0;

    /* writing out the data that is being read into NAG routine (for debugging) */
//     for(j=0; j<NL; j++) {
//       fprintf(stderr,"OUTxy[%3d]: %10.6f \n",j,OUTxy[j]);
//     }
//     fprintf(stderr,"\n");
//     
    //    fprintf(stderr,"n_new:%3d  nl_new:%3d  SDxy:%10.6f  OUT0xy:%10.6f statxy:%10.6f IFAIL:%d \n",n_new,nl_new,SDxy,OUT0xy,statxy,IFAIL);
//     
//     for(j=0; j<N; j++) {
//       fprintf(stderr,"Txp[%3d]: %10.6f \n",j,Txp[j]);
//     }
//     fprintf(stderr,"\n");
//     
//     for(j=0; j<N; j++) {
//       fprintf(stderr,"Rxp[%3d]: %10.6f \n",j,Rxp[j]);
//     }
//     fprintf(stderr,"\n");
//     
//     fprintf(stderr,"Rx_maxval:%10.6f \n",Rx_maxval);

    /* First check if Rx array contains only zeroes. 
       If it does then set X_correl array to zero and don't calculate x-correlation.
       Otherwise calculate x-correlation */

    if (Rx_maxval==0.0) {
      //printf("\n Rx array contains only zeroes => Not calculating the Cross-Correlation \n");
      OUT0xy=0.0;  
      SDxy=0.0;
      statxy=0.0;
      for (j=0;j<NL;j++) OUTxy[j]=0.0;
    }else{
      //printf("\n Calculating Cross-Correlation between Tx and Rx \n");
      //      fprintf(stderr," Error Indicator BEFORE NAG Routine is called: IFAIL=%d \n",IFAIL);
      /* using NAG routine G13BCF */
      g13bcf_(Txp, Rxp, &n_new, &nl_new, &SDxy, &OUT0xy, OUTxy, &statxy, &IFAIL);
      //      fprintf(stderr," Error Indicator AFTER NAG Routine is called: IFAIL=%d \n",IFAIL);
    }

    //    /* printing out x-correl results */
    //    printf("\n\n Results of Cross Correlation between Tx and Rx \n");
    //    fprintf(stderr," Standard Deviation Ratio: %7.4f	Test statistic: %12.6f\n",SDxy,statxy);
    //    fprintf(stderr,"\n Cross correlation at lag 0: %11.8f	\n",OUT0xy);
    //    for(j=0; j<NL; j++) {    
    //    fprintf(stderr,"OUTxy[%3d] %11.8f \n",j,OUTxy[j]);
    //    } 

    /* Generating Xcorrel array */
    //    printf("\n Making Cross-Correlation array \n");
    Xcorrel_xy[0]=OUT0xy;
    for(j=1; j<N; j++) {
      Xcorrel_xy[j] = OUTxy[j-1];
    }
    //    for(j=0; j<N; j++) {    
    //      fprintf(stderr,"Xcorrel_xy[%3d] Cross Correlation at lag %3d	%11.8f	\n",j,j,Xcorrel_xy[j]);
    //    } 

    //    printf("\n Finding the maximum value of the Cross-Correlation array \n");
    double Xcorrel_maxval=0.0;
    int Xcorrel_maxval_BinPos=0;
    for( j=0; j<N; j++) {
      if(Xcorrel_xy[j]>Xcorrel_maxval){
	Xcorrel_maxval=Xcorrel_xy[j];
	Xcorrel_maxval_BinPos=j;
      }	
    }
    //fprintf(stderr,"echo:%d  Xcorrel_maxval: %9.6f 	Xcorrel_maxval_BinPos: %3d	\n",i,Xcorrel_maxval,Xcorrel_maxval_BinPos);

    //    printf("\n Re-arranging Cross-Correlation array \n");
    for(j=0; j<n2; j++) {
      Xcorrel_xy_rar[j] = OUTxy[(n2m1-1)-j];
      Xcorrel_xy_rar[n2m1]=OUT0xy;
    }
    for(j=n2; j<N; j++) {
      Xcorrel_xy_rar[j] = OUTxy[j-n2];
    }
    //    for(j=0; j<N; j++) {
    //      fprintf(stderr,"Xcorrel_xy_rar[%3d] Cross Correlation at lag %3d	%11.8f	\n",j,j,Xcorrel_xy_rar[j]);
    //    } 
  
    //-----------------------------------------------------------------------------------------------------//
    /* Now we perform statistical calculations on data in Xcorrel_xy_rar array for data values above zero */

    /* calculating ATBD_defined MEAN */
    int k=0, Xcorrel_nvalsB=0;
    double Xcorrel_sumC=0.0,Xcorrel_sumD=0.0,Xcorrel_sumE=0.0,Xcorrel_sumF=0.0,Xcorrel_sumG=0.0;
    double Xcorrel_mean=0.0;
    double Xcorrel_prodA=0.0,Xcorrel_prodB=0.0,Xcorrel_prodC=0.0,Xcorrel_prodD=0.0;
    double Xcorrel_imm_sq=0.0,Xcorrel_imm_cube=0.0,Xcorrel_imm_four=0.0;
    double Xcorrel_sd=0.0,Xcorrel_sd_sq=0.0,Xcorrel_sd_cube=0.0,Xcorrel_sd_four=0.0;
    double Xcorrel_skew=-999.0,Xcorrel_kurt=-999.0;
    //    fprintf(stderr,"Xcorrel_skew: %9.6f 	Xcorrel_kurt: %9.6f	\n",Xcorrel_skew,Xcorrel_kurt);

    for(k=0; k<N; k++) {
      if (Xcorrel_xy_rar[k]>=0){
	Xcorrel_sumC+=Xcorrel_xy_rar[k];
	Xcorrel_prodA=(k*Xcorrel_xy_rar[k]);
	Xcorrel_sumD+=Xcorrel_prodA;
	Xcorrel_nvalsB++;
	Xcorrel_mean=(Xcorrel_sumD/Xcorrel_sumC);
	//	Xcorrel_mean=(Xcorrel_sumC/Xcorrel_nvalsB);
      }
    }

    if (Xcorrel_maxval==0) {
      Xcorrel_nvalsB=0;
      Xcorrel_mean=0.0;
    }
    //    fprintf(stderr,"Xcorrel_nvalsB: %d	\n",Xcorrel_nvalsB);
    //    fprintf(stderr,"Xcorrel_mean: %9.6f \n",Xcorrel_mean);

    /* calculating ATBD_defined STANDARD DEVIATION */
    for(k=0; k<N; k++) {
      if (Xcorrel_xy_rar[k]>=0){
	Xcorrel_imm_sq=((k-Xcorrel_mean)*(k-Xcorrel_mean));
	Xcorrel_prodB=(Xcorrel_imm_sq*Xcorrel_xy_rar[k]);
	Xcorrel_sumE+=Xcorrel_prodB;
	Xcorrel_sd_sq=(Xcorrel_sumE/Xcorrel_sumC);
	Xcorrel_sd=(sqrt(Xcorrel_sumE/Xcorrel_sumC));
      }
    }
    if (Xcorrel_maxval==0) {
      Xcorrel_sd=0.0;
    }
    //    fprintf(stderr,"Xcorrel_sd: %9.6f \n",Xcorrel_sd);

    /* calculating ATBD_defined SKEWNESS */
    for(k=0; k<N; k++) {
      if (Xcorrel_xy_rar[k]>=0){
	Xcorrel_imm_cube=((k-Xcorrel_mean)*(k-Xcorrel_mean)*(k-Xcorrel_mean));
	Xcorrel_prodC=(Xcorrel_imm_cube*Xcorrel_xy_rar[k]);
	Xcorrel_sumF+=Xcorrel_prodC;
	Xcorrel_sd_cube=(Xcorrel_sd_sq*Xcorrel_sd);
	Xcorrel_skew=((1/Xcorrel_sd_cube)*(Xcorrel_sumF/Xcorrel_sumC));
      }
    }
    if (Xcorrel_maxval==0) {
      Xcorrel_skew=0.0;
    }
    //    fprintf(stderr,"Xcorrel_skew: %9.6f \n",Xcorrel_skew);
	
    /* calculating ATBD_defined KURTOSIS */
    for(k=0; k<N; k++) {
      if (Xcorrel_xy_rar[k]>=0){
	Xcorrel_imm_four=((k-Xcorrel_mean)*(k-Xcorrel_mean)*(k-Xcorrel_mean)*(k-Xcorrel_mean));
	Xcorrel_prodD=(Xcorrel_imm_four*Xcorrel_xy_rar[k]);
	Xcorrel_sumG+=Xcorrel_prodD;
	Xcorrel_sd_four=(Xcorrel_sd_cube*Xcorrel_sd);
	Xcorrel_kurt=(((1/Xcorrel_sd_four)*(Xcorrel_sumG/Xcorrel_sumC))-3);
      }
    }
    if (Xcorrel_maxval==0) {
      Xcorrel_kurt=0.0;
    }
    //    fprintf(stderr,"Xcorrel_kurt: %9.6f \n",Xcorrel_kurt);

    /* Printing stats to screen */
    // 	fprintf(stderr,"\n");
    //  fprintf(stderr,"%d values in Xcorrel Distribution >= zero \n",Xcorrel_nvalsB); 
    // 	fprintf(stderr,"\n");
    //  fprintf(stderr,"%8.5f sum of element values above threshold \n",Xcorrel_sumC);
    //  fprintf(stderr,"%8.5f product of element value times element number for the last element in this sub-array \n",Xcorrel_prodA);
    //  fprintf(stderr,"%8.5f sum of the products \n",Xcorrel_sumD);
    // 	fprintf(stderr,"%8.5f	 ATBD_defined mean \n",Xcorrel_mean);
    // 	fprintf(stderr,"\n");
    //  fprintf(stderr,"%8.5f square of element number minus ATBD_def_mean for the last element in this sub-array \n",Xcorrel_imm_sq);
    //  fprintf(stderr,"%8.5f product of element value times i-mean_sq for the last element in this sub-array \n",Xcorrel_prodB);
    // 	fprintf(stderr,"%8.5f sum of the products \n",Xcorrel_sumE);
    // 	fprintf(stderr,"%8.5f square of the ATBD_defined StDev \n",Xcorrel_sd_sq);
    // 	fprintf(stderr,"%8.5f 	ATBD_defined StDev \n",Xcorrel_sd);
    // 	fprintf(stderr,"\n");
    //  fprintf(stderr,"%12.5f cube of element number minus ATBD_def_mean for the last element in this sub-array \n",Xcorrel_imm_cube);
    //  fprintf(stderr,"%12.5f product of element value times i-mean_cube for the last element in this sub-array \n",Xcorrel_prodC);
    // 	fprintf(stderr,"%12.5f sum of the products \n",Xcorrel_sumF);
    // 	fprintf(stderr,"%12.5f cube of the ATBD_defined StDev \n",Xcorrel_sd_cube);
    // 	fprintf(stderr,"%12.5f 	ATBD_defined Skewness \n",Xcorrel_skew);
    // 	fprintf(stderr,"\n");
    //  fprintf(stderr,"%16.5f element number minus ATBD_def_mean to the four,for the last element in this sub-array \n",Xcorrel_imm_four);
    //  fprintf(stderr,"%16.5f product of element value times i-mean_four for the last element in this sub-array \n",Xcorrel_prodD);
    // 	fprintf(stderr,"%16.5f sum of the products \n",Xcorrel_sumG);
    // 	fprintf(stderr,"%16.5f ATBD_defined StDev to the four \n",Xcorrel_sd_four);
    // 	fprintf(stderr,"%16.5f 	ATBD_defined Kurtosis \n",Xcorrel_kurt);

    //-----------------------------------------------------------------------------------------------------//
    /* writing output files */
	
    //    /* writing output file A */
    //    /* output format: Txp[array]  */
    //    for(j=0; j<N; j++) {
    //      fprintf(fouta," %9.6f",Txp[j]);
    //    }
    //    /* place a carriage return after printing each line of the output file */    
    //    fprintf(fouta,"\n");

    //    /* writing output file B */
    //    /* output format: Rxp[array]  */
    //    for(j=0; j<N; j++) {
    //      fprintf(foutb," %9.6f",Rxp[j]);
    //    }
    //    /* place a carriage return after printing each line of the output file */    
    //    fprintf(foutb,"\n");

    /* writing output file C */
    /* output format: Xcorrel_xy[array]  */
    /* only printing out the first 200 values in the Xcorrel_xy array */
    for(j=0; j<200; j++) fprintf(foutc," %9.6f",Xcorrel_xy[j]);
    /* place a carriage return after printing each line of the output file */    
    fprintf(foutc,"\n");
   
    /* writing output file D  */  
    /* output format: Xcorrel_xy_rar[array]  */
    /* only printing out the central 200 values in the Xcorrel_xy_rar array */
    for(j=100; j<299; j++) fprintf(foutd," %9.6f",Xcorrel_xy_rar[j]);
    /* place a carriage return after printing each line of the output file */    
    fprintf(foutd,"\n");

    /* writing output file E */
    /* output format:line no.,standard dev ratio,test statistic,x_correl@lag0,x_correl_maxval,x_correl_maxval_binpos,
       number of Xcorrel values >=0,Xcorrel_mean,Xcorrel_sd,Xcorrel_skew,Xcorrel_kurt  */
    /* write the stats line */
    fprintf(foute,"%7d %9.6f %14.6f %9.6f %9.6f %3d %3d %9.1f %9.6f %9.6f %9.6f",i,SDxy,statxy,OUT0xy,Xcorrel_maxval,Xcorrel_maxval_BinPos,Xcorrel_nvalsB,Xcorrel_mean,Xcorrel_sd,Xcorrel_skew,Xcorrel_kurt);
    /* place a carriage return after printing each line of the output file */    
    fprintf(foute,"\n");

    //-----------------------------------------------------------------------------------------------------//
    /*free the Tx and Rx arrays so the next lines can be read*/
    free(Tx);
    free(Rx);
    free(Txp);
    free(Rxp);
    free(OUTxy);

  }
     
  //fprintf(stderr,"-------------------------------------------------------------------------------------------------------\n");
  //fprintf(stderr,"Main Loop Ends here\n");

  /* closing the output files having written to them*/
  //  fclose(fouta);
  //  fclose(foutb);
  fclose(foutc);
  fclose(foutd);
  fclose(foute);

  /*free the memory allocated to the arrays containing the filenames*/
  free(filename1);
  free(filename2);
  
  /* Print out time at end of run to screen */
  ttime=time(&ttime);
  fprintf(stderr,"                 end time: %s\n",asctime(localtime(&ttime)));

}
