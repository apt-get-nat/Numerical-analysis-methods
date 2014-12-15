/* File: chebyshev.cpp
 * Author: Nathaniel Mathews
 * Description: Originally built as a module in a Relativity simulation.
 *              Does chebyshev decomposition and expansion.
 */
 
#include <math.h>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <iostream>

using namespace std;

double pi = 3.14159265;

vector < double > collocation_1D ( int n ) {
// Arguments: The integer representing the number of terms to
//            be used in the approximation.
// Returns: A vector of length n with the collocation points in -1,1.
// Preconditions: None.
// Postconditions: None.
	vector < double > col(n);
	for ( int k = 0; k < n; k++) {
		col.at(k) = cos( pi/n * (k+0.5) );
	}
	return col;
}

vector < double > cheb_coefficients ( vector < double > col_points, vector < double > col_val ) {
// Arguments: A vector the values at collocation points
// Returns: A vector of chebyshev coefficients necessary to model such points
// Preconditions: None.
// Postconditions: None.
	double alpha;
	int n_term = col_points.size();
	vector < double > coef(n_term) ;
	for ( int i = 0; i < n_term; i++ ) {
		if ( i == 0 ) { alpha = 1.0 / (double)n_term; }
		else { alpha = 2.0 / (double)n_term; }
		coef.at( i ) = 0;
		for ( int j = 0; j < n_term; j++ ) {
			coef.at( i ) += col_val.at(j) * cos ( i * acos ( col_points.at(j) ) );
		}
		coef.at(i) = alpha * coef.at(i);
	}
	return coef;
}

double cheb_reconstruct ( vector < double > coefficients, double xval ) {
// Arguments: A vector of coefficients of size n_term and the place we wish to approximate a value of
// Returns: A double approximating the value at the given point.
// Preconditions: None.
// Postconditions: None.
	double sum = 0;
	int n_term = coefficients.size();
	for ( int i = 0; i < n_term; i++ ) {
		sum += coefficients.at(i) * cos ( i * acos ( xval ) );
	}
	return sum;
}