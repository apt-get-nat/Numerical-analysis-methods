/* File: RK.cpp
 * Authors: Nathaniel Mathews, Michael Potter
 * Description: A collection of Runge-Kutta methods
 */
 
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <math.h>
#include <float.h>

#define STEPS 0
#define ERROR 1

//Returns [# steps required, error]
double* RK4(double (*funcp)(double, double), double yInit, double tInit, double trueVal){
	double* result = (double*) calloc(2,sizeof(double));
	result[STEPS]=1;
	result[ERROR]=DBL_MAX;
	double solution;
	double curTime;
	double lowBound = 1;
	double highBound = DBL_MAX;
	while(true){
		highBound == DBL_MAX ? result[STEPS]=2.0*result[STEPS] : result[STEPS] = ceil((lowBound + highBound)/2.0);
		solution = yInit;
		curTime = tInit;
		double stepSize = 1.0/result[STEPS];
		for(int i=0; i<(int)result[STEPS]; i++){
			double k1 = stepSize * funcp(curTime, solution);
			double k2 = stepSize * funcp(curTime+0.5*stepSize,solution+0.5*k1);
			double k3 = stepSize * funcp(curTime+0.5*stepSize,solution+0.5*k2);
			double k4 = stepSize * funcp(curTime+stepSize,solution+k3);
			solution=solution+(1.0/6.0)*k1+(1.0/3.0)*k2+(1.0/3.0)*k3+(1.0/6.0)*k4;
			curTime+=stepSize;
		}
		result[ERROR]=absolute(solution-trueVal);
		if(result[ERROR]<1e-6 && result[ERROR]>9e-7){
			break;
		}else{

			if(result[ERROR]<1e-6){ //Too accurate
				highBound = result[STEPS];
			}
			else{ //Not accurate Enough
				lowBound = result[STEPS];
			}
			if(highBound!=DBL_MAX){
				if(lowBound+1>=highBound && result[ERROR]<1e-6){ //Can't hit sweet spot
					break;
				}
			}
			
		}
	}
	return result;
}

double AdamsMoulton(double (*f1)(double, double), double h, double* y, double t0, double tEnd) {
	double ytemp;
	for(double t = t0; t <= tEnd; t+= h) {
		ytemp = y[0] + (h/720.0) * (
			  1901.0*f1(t    ,y[0])
			- 2774.0*f1(t-  h,y[1])
			+ 2616.0*f1(t-2*h,y[2])
			- 1274.0*f1(t-3*h,y[3])
			+  251.0*f1(t-4*h,y[4])
		);
		ytemp = y[0] + h/720.0*(
			  251.0*f1(t  +h, ytemp)
			+ 646.0*f1(t    , y[0])
			- 264.0*f1(t  -h, y[1])
			+ 106.0*f1(t-2*h, y[2])
			-  19.0*f1(t-3*h, y[3])
		);
		y[4]=y[3];y[3]=y[2];y[2]=y[1];y[1]=y[0];y[0]=ytemp;
	}
	return y[0];
}

// Where Y = dPsi/dx (ie, set up to solve second-order ODEs).
// Uses Dormand-Prince method.
long double RK45( long double (*dPsi)(long double, long double), long double (*dY)(long double, long double), 
			 long double xInit, long double psiInit, long double dPsiInit, long double epsilon, long double xEnd, 
			 long double h) {
    long double x = xInit;
    long double psi5 = psiInit;
    long double psi4 = psiInit;
    long double dPsi5 = dPsiInit;
    long double dPsi4 = dPsiInit;
    
    long double psi5temp = psi5;
    long double dPsi5temp = dPsi5;
    long double error = 0;
    long double q = 0;
    int exitFlag = 0;
    int count = 0;
    int redoCount = 0;
    
    long double k1psi = dPsi(x,dPsi5);
    long double k2psi = 0;
    long double k3psi = 0;
    long double k4psi = 0;
    long double k5psi = 0;
    long double k6psi = 0;
    long double k7psi = 0;

	long double k1y = dY(x,psi5);
    long double k2y = 0;
    long double k3y = 0;
    long double k4y = 0;
    long double k5y = 0;
    long double k6y = 0;
    long double k7y = 0;

    // long double c1 = 0;
    long double c2 = 1.0/5.0;
    long double c3 = 3.0/10.0;
    long double c4 = 4.0/5.0;
    long double c5 = 8.0/9.0;
    // long double c6 = 1;
    // long double c7 = 1;


    // long double a0 = 0;
    long double a11 = (1.0/5.0);
    long double a21 = (3.0/40.0);
    long double a22 = (9.0/40.0);
    long double a31 = (44.0/45.0);
    long double a32 = (-56.0/15.0);
    long double a33 = (32.0/9.0);
    long double a41 = (19372.0/6561.0);
    long double a42 = (-25360.0/2187.0);
    long double a43 = (64448.0/6561.0);
    long double a44 = (-212.0/729.0);
    long double a51 = (9017.0/3168.0);
    long double a52 = (-355.0/33.0);
    long double a53 = (46732.0/5247.0);
    long double a54 = (49.0/176.0);
    long double a55 = (-5103.0/18656.0);
    long double a61 = (35.0/384.0);
    // long double a62 = 0;
    long double a63 = (500.0/1113.0);
    long double a64 = (125.0/192.0);
    long double a65 = (-2187.0/6784.0);
    long double a66 = (11.0/84.0);

    long double b14 = (5179.0/57600.0);
    // long double b24 = 0;
    long double b34 = (7571.0/16695.0);
    long double b44 = (393.0/640);
    long double b54 = (-92097.0/339200.0);
    long double b64 = (187.0/2100.0);
    long double b74 = (1.0/40.0);

    long double b15 = (35.0/384.0);
    // long double b25 = 0;
    long double b35 = (500.0/1113.0);
    long double b45 = (125.0/192.0);
    long double b55 = (-2187.0/6784.0); 
    long double b65 = (11.0/84.0);
    // long double b76 = 0; 


    //Method
    while(x<=xEnd){
        //If willing to overstep the end then must be safe to step only to the end
        if(x+h>xEnd){
            h=xEnd-x;
            exitFlag = 1;
        }
        k2y   = dY(  x+c2*h,  psi5+h*a11*k1psi);
        k2psi = dPsi(x+c2*h, dPsi5+h*a11*k1y);
        
        k3y   = dY(  x+c3*h,  psi5+h*(a21*k1psi+a22*k2psi));
        k3psi = dPsi(x+c3*h, dPsi5+h*(a21*k1y+a22*k2y));
         
        k4y   = dY(  x+c4*h, psi5+h*(a31*k1psi+a32*k2psi+a33*k3psi));
        k4psi = dPsi(x+c4*h, dPsi5+h*(a31*k1y+a32*k2y+a33*k3y));
        
        k5y   = dY(  x+c5*h,  psi5+h*(a41*k1psi+a42*k2psi+a43*k3psi+a44*k4psi));
        k5psi = dPsi(x+c5*h, dPsi5+h*(a41*k1y+a42*k2y+a43*k3y+a44*k4y));
        
        k6y   =   dY(x+h,  psi5+h*(a51*k1psi+a52*k2psi+a53*k3psi+a54*k4psi+a55*k5psi));
        k6psi = dPsi(x+h, dPsi5+h*(a51*k1y+a52*k2y+a53*k3y+a54*k4y+a55*k5y));
        
        k7y   =   dY(x+h,  psi5+h*(a61*k1psi+a63*k3psi+a64*k4psi+a65*k5psi+a66*k6psi));
        k7psi = dPsi(x+h, dPsi5+h*(a61*k1y+a63*k3y+a64*k4y+a65*k5y+a66*k6y));
        
        dPsi4     = dPsi5 + h*(b14*k1y + b34*k3y + b44*k4y + b54*k5y + b64*k6y + b74*k7y);
        dPsi5temp = dPsi5 + h*(b15*k1y + b35*k3y + b45*k4y + b55*k5y + b65*k6y);
        
        psi4     = psi5 + h*(b14*k1psi + b34*k3psi + b44*k4psi + b54*k5psi + b64*k6psi + b74*k7psi);
        psi5temp = psi5 + h*(b15*k1psi + b35*k3psi + b45*k4psi + b55*k5psi + b65*k6psi);
        
        
        
        count++;
        //Get out if you're done
        if(exitFlag){
            psi5=psi5temp;
            x=x+h;
            break;
        }
        //Approximate the accumulated error to change step size
        error = max(absolute(psi5temp-psi4),absolute(dPsi5temp-dPsi4));
        error == 0 ? q = 2.0 : q = pow(epsilon*h/error,0.25);
        if (q < epsilon) q = epsilon;
        if(q>=1.0){ //Error is small enough
        	printf("#h=%.20Lf good: q=%Lf\n", h, q);
            psi5 = psi5temp;
            dPsi5 = dPsi5temp;
            printf("%.20Lf %.20Lf\n", x, psi5);
            x=x+h;
            q/pow(2.0,0.25) > 1.0 ? h = h*q/pow(2.0,0.25) : h = h;
            k1psi = k7psi;
            k1y = k7y;
        }else{ //Error too big, try again
        	printf("#h=%.20Lf bad: q=%Lf\n", h, q);
            h=h*q;
            redoCount++;
        }
    }
    printf("#Loop Count: %d\n",count);
    printf("#Wasted Loops: %d\n",redoCount);
    return psi5;
}

