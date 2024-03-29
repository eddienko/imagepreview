
/*
 * The following code is public domain.
 * Algorithm by Torben Mogensen, implementation by N. Devillard.
 * This code in public domain.
 */

#include <math.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include "fitsio.h"

float torben(float m[], int n)
{
    int         i, less, greater, equal;
    float  min, max, guess, maxltguess, mingtguess;
    
    min = max = m[0] ;
    for (i=1 ; i<n ; i++) {
        if (m[i]<min) min=m[i];
        if (m[i]>max) max=m[i];
    }
    
    while (1) {
        guess = (min+max)/2;
        less = 0; greater = 0; equal = 0;
        maxltguess = min ;
        mingtguess = max ;
        for (i=0; i<n; i++) {
            if (m[i]<guess) {
                less++;
                if (m[i]>maxltguess) maxltguess = m[i] ;
            } else if (m[i]>guess) {
                greater++;
                if (m[i]<mingtguess) mingtguess = m[i] ;
            } else equal++;
        }
        if (less <= (n+1)/2 && greater <= (n+1)/2) break ;
        else if (less>greater) max = maxltguess ;
        else min = mingtguess;
    }
    if (less >= (n+1)/2) return maxltguess;
    else if (less+equal >= (n+1)/2) return guess;
    else return mingtguess;
}

float mad(float m[], int n)
{
	int i;
	float median, mad;
	float *m2;
	
	m2 = (float *) calloc(n, TFLOAT);
	
	
	median = torben(m, n);
    printf("%f %f\n", median, m[200]);
	
	for (i=0; i<n; i++) m2[i]=fabs(m[i]-median);
	
	mad = torben(m2, n)*1.4826;
	
	return mad;
}

float *zscale(float m[], int n)
{
	static float retbuf[2];
	int i;
	float median, mad;
	float *m2;
	
	m2 = (float *) calloc(n, TFLOAT);
	median = torben(m, n);
	for (i=0; i<n; i++) m2[i]=fabs(m[i]-median);
	mad = torben(m2, n)*1.4826;
	
	retbuf[0]=median;
	retbuf[1]=mad;
	return retbuf;
}
