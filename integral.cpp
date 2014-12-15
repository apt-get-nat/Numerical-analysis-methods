/* File: romberg.cpp
 * Author: Nathaniel Mathews
 * Description: An implementation of Simpson's Method
 */
 
#include <stdlib.h>
#include <math.h>
#include <stdio.h>

double SimpsonsMethod(double (*fH)(double), double a, double b, int n){
    double range = b-a;
    double result = fH(a) + fH(b);
    for(int i=1; i<n;i+=2) {
		result += 4*fH(range / (double)n * i + a);
	}
	for(int i=2; i<n-1;i+=2) {
		result += 2*fH(range / (double)n * i + a);
	}
	result = result * (double)range/(3.0*(double)n);
	return result;
}

