/* This driver file is written by Parisa Babaheidarian. The function computes PCA decomposition of sinogram measurement in CT simulations. 
Five bases of PCA are used which can be changed. "n" is the number of variables and is equal to the number of basis functions used. 
The line decomposition function uses L_BFGS-B algorithm. This source code also includes the MATLAB wrapper. */
#include<mex.h>
#include<matrix.h>
#include<math.h>
#include "lbfgsb.h"
#include "JacFPCA5.h"
/* Computational routine*/

void gLineDecompositionPCA(double bin1,double bin2,double bin3,double bin4,double bin5,
double bin6,double bin7,double bin8,double bin9,double bin10, double x[5],double * R)
{
    integer i__1;
    /* Local variables */
    static double f, *g;
    static double gg=0;
    g=&gg;
   // static double x[5];
    //Necessary constants
    static double B [5][100];
    static double Weight[10][100];
    FILE *myFile;
   // myFile= fopen("Bjuly.txt", "r"); // 10 by 100 Wbrute
   myFile= fopen("FPCA.txt", "r"); // 5 by 100 normalized PEC
    FILE *myFile2;
    myFile2 = fopen("WeightJuly.txt", "r");
    
    //Assigning values to Basis and Weight matrices
    for (int i=0; i<100; i++){
     for (int j=0; j<10; j++){
        fscanf(myFile2, "%lf", &Weight[j][i]); 
     }
    }
     for (int i=0; i<100; i++){
     for (int j=0; j<5; j++){
        fscanf(myFile, "%lf", &B[j][i]); 
     }
    }
    printf("We is: %lf ",Weight[0][0]);
    fclose(myFile);
    fclose(myFile2);
    static integer i__;
    static double l[5];
    static integer m, n;
    static double u[5], t1, t2, wa[43251];
    static integer nbd[5], iwa[3072];
/*     static char task[60]; */
    static integer taskValue;
    static integer *task=&taskValue; /* must initialize !! */
/*      http://stackoverflow.com/a/11278093/269192 */
    static double factr;
/*     static char csave[60]; */
    static integer csaveValue;
    static integer *csave=&csaveValue;
    static double dsave[29];
    static integer isave[44];
    static logical lsave[4];
    static double pgtol;
    static integer iprint;
   
/*     We wish to have output at every iteration. */
    iprint = 1; 
/*     iprint = 101; */
/*     We specify the tolerances in the stopping criteria. */
    factr = 1e7;
    pgtol = 1e-5;
/*     We specify the dimension n of the sample problem and the number */
/*        m of limited memory corrections stored.  (n and m should not */
/*        exceed the limits nmax and mmax respectively.) */
    n = 5;
    m = 10;
/*     We now provide nbd which defines the bounds on the variables: */
/*                    l   specifies the lower bounds, */
/*                    u   specifies the upper bounds. */
/*     First set bounds on the odd-numbered variables. */
    i__1 = n;
    for (i__ = 1; i__ <= i__1; i__ += 1) {
        nbd[i__ - 1] = 0;
       // l[i__ - 1] = 0.;
    }
    /*     We now define the starting point. */
    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
        x[i__ - 1] = 0.;
        /* L14: */
    }
//    x0=x[0];
  //  x1=x[1];
   // x2=x[2];
   // x3=x[3];
   // x4=x[4];
    printf("     Solving decomposition problem (Speck basis).\n");
    printf("      (f = 0.0 at the optimal solution.)\n");

    /*     We start the iteration by initializing task. */

    *task = (integer)START;
/*     s_copy(task, "START", (ftnlen)60, (ftnlen)5); */
    /*        ------- the beginning of the loop ---------- */
L111:
    /*     This is the call to the L-BFGS-B code. */
setulb(&n, &m, x, l, u, nbd, &f, g, &factr, &pgtol, wa, iwa, task, &
            iprint, csave, lsave, isave, dsave);
          //  x0=x[0];
          //  x1=x[1];
          //  x2=x[2];
          //  x3=x[3];
           // x4=x[4];
/*     if (s_cmp(task, "FG", (ftnlen)2, (ftnlen)2) == 0) { */
    if ( IS_FG(*task) ) {
        /*        the minimization routine has returned to request the */
        /*        function f and gradient g values at the current x. */
        /* 
        Compute function value f for the problem. */
        Essentials myfuns(x,B,Weight); 
        f=myfuns.FunMain(bin1,bin2,bin3,bin4,bin5,bin6,bin7,bin8,bin9,bin10);
        if (f>4)
        {printf("Fvalue is: %lf \n",f);}
        i__1 = n;
       
        /*        Compute gradient g for the problem. */
        
        g=myfuns.Grad(bin1,bin2,bin3,bin4,bin5,bin6,bin7,bin8,bin9,bin10);
       /* printf("gvalue is: %lf \n",g[0]);*/
        R=x;
        /*          go back to the minimization routine. */
        goto L111;
    }

/*     if (s_cmp(task, "NEW_X", (ftnlen)5, (ftnlen)5) == 0) { */
    if ( *task==NEW_X ) {
        goto L111;
    }
    
}

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])   { 

double * bin1; double *bin2; double * bin3; double * bin4; double * bin5; double * bin6;
double * bin7; double * bin8; double * bin9; double * bin10;
bin1=mxGetPr(prhs[0]); bin2=mxGetPr(prhs[1]);
bin3=mxGetPr(prhs[2]);bin4=mxGetPr(prhs[3]);
bin5=mxGetPr(prhs[4]); bin6=mxGetPr(prhs[5]);
bin7=mxGetPr(prhs[6]);bin8=mxGetPr(prhs[7]);
bin9=mxGetPr(prhs[8]); bin10=mxGetPr(prhs[9]);
double * x0=mxGetPr(prhs[10]); 
//double * x1=mxGetPr(prhs[11]);
//double * x2=mxGetPr(prhs[12]); 
//double * x3=mxGetPr(prhs[13]);
//double * x4=mxGetPr(prhs[14]);
double * R;

plhs[0]=mxCreateDoubleMatrix(5,1,mxREAL);
R=mxGetPr(plhs[0]);
gLineDecompositionPCA(*bin1,*bin2,*bin3,*bin4,*bin5,*bin6,*bin7,*bin8,*bin9,*bin10,x0,R);


}

/* Testing compilation with g++ with main function */
/*
int main (){
double bin[10]={1.,2.,1.,1.,1.,1.,1.,1.,1.,1.};
double * x;
x=LineDecomposition(bin[0],bin[1],bin[2],bin[3],bin[4],bin[5],bin[6],bin[7],bin[8],bin[9]);
printf("X[9] is: %a \n",x[9]);
return 0;
}*/
