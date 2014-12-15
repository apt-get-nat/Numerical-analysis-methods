/* File: romberg.cpp
 * Author: Nathaniel Mathews
 * Description: A function to compute Romberg Integrals
 */
 
#include <stdlib.h>
#include <math.h>
#include <stdio.h>

double h(double a, double b, int n) {
    return (b-a)/pow(2.0,n);
}

double R(double (*f)(double), int n, int m, double a, double b){
    double result = 0;
    if(n==0 && m==0){
        return (b-a)/2.0 * (f(a)+f(b));
    } else if(m==0){
        for(int i = 1; i<=pow(2,n-1); i++){
           result+=f(a+((double)(2*i-1))*h(a,b,n)); 
        }
        result *= h(a,b,n);
        result += 0.5*R(n-1,0,a,b);
        return result;
    } else {
        return R(n, m-1, a, b) + 1/(pow(4, m)-1) * (R(n, m-1, a, b) - R(n-1, m-1,a,b));
    }
}