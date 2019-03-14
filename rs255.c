
/************************************************************************
  This implementation of the famous Reed-Solomon codes provides an 
  encoder and a decoder function for each. They both operate on arrays of
  255 char. Each char representing an 8-bit value in the GF(256).

  The encoder, rs255encode, takes the (255-k) input "bytes" stored in the 
  first (255-k) entries and calculates the k parity symbols, which are 
  then stored in the remaining k entries. The value of k is provided as
  input to the function.

  The decoder, rs255decode, has two input arrays and one output array.
  The first input array is the received codeword with errors. The 
  second array indicates the erasure positions. Any non-zero element of
  this array marks the corresponding element of the received codeword as
  unreliable. If no erasure information is available simply set this input
  to all zeros. The output array returns the corrected codeword.

  For more information about the return values see below.

  Note that a smaller codeword can easily be constructed by setting the
  "unused" input bytes to zeros when calling the encoder. Before calling
  the decoder the same entries are set to zeros.

  The decoder corrects all error patterns of weight k/2 and all erasure
  patterns of weight k. Weight is equivalent to the number of byte
  positions in error. Any combination of errors and erasures is corrected
  as long as (2*error + erasure) <= k.

  If the number of errors exceeds the limit, then the decoder detects the
  situation and reports an ucorrectable error pattern (with a very high
  probability).

  Author:

  Jens Jørgen Nielsen
  e-mail:jjn2970@gmail.com
  http:\\www.corix.dk
  
  Terms of use:
  
  To use this code in a project or product please ask the author for 
  a license. For non-commercial use it is currently (2015) provided for 
  free (2015). This however may change, so please ask.

  To distribute the source code please do NOT change anything in this file. 
  Otherwise feel free to pass it on.
*************************************************************************/

#include <stdio.h>
#include <stdlib.h>

#include "gf256.h"
#include "rs255.h"

static unsigned char debugFlag;

unsigned char polyEval(unsigned char a, unsigned char b[], long n)
{
    long            i;
    unsigned char   x;
	
	x = 0;
    for (i=(n-1); i>=0; i--) {
		x = gf256mul(a,x) ^ b[i];
	}
	return(x);
}


static void polyPrint(unsigned char b[], long n)
{
    long            i,k;
	
	k = 0;
    for (i=0; i<n; i++) {
		if (k==0) {
			printf("\n %3ld ",i);
			k = 20;
		}
		printf("%3hhx",b[i]);
		k = k-1;
	}
}

/* Find the roots in sigma. If the degree of sigma and the number of roots
   are equal each of the roots correspond to an error position. The value of
   each error is calculated. The return value is the number of roots found in
   sigma if it is equal to the degree of sigma, otherwise -1 is returned
*/
long rs255errors(unsigned char sigma[],    //error locator poly
                 unsigned char omega[],    //error value poly
                 unsigned char position[], //error positions
                 unsigned char value[],    //error value
                 long maxDegree)           //max sigma degree
{
    long i, sigmaDegree, deltaSigmaDegree, errorCount;
    unsigned char   x,y,y1,y2,y3;
    unsigned char   deltaSigma[256];
	

    errorCount = 0;
    sigmaDegree = maxDegree;
    while ((sigma[sigmaDegree]==0) && (sigmaDegree>0)) sigmaDegree = sigmaDegree - 1;
//if (maxDegree>31)    printf("\n Degree is : %3ld",sigmaDegree);

    i=0;
    while ((2*i+1)<=sigmaDegree) {
        deltaSigma[i] = sigma[2*i+1];
        i = i+1;
    }
    deltaSigmaDegree = i-1;

//if (maxDegree>31)    polyPrint(sigma, sigmaDegree+1);
//    polyPrint(deltaSigma, deltaSigmaDegree+1);
//    polyPrint(omega, maxDegree);

    for (i=0; i<=255; i++) {
        x = (unsigned char) i;
        y = polyEval(x, sigma, (sigmaDegree+1));
        if (y==0) {
            position[errorCount] = 255-gf256loga(x);
            y1 = polyEval(gf256mul(x,x),deltaSigma,deltaSigmaDegree+1);
            y2 = polyEval(x,omega,maxDegree);
            y3 = gf256x111(x);
            value[errorCount] = gf256mul(y3,gf256mul(y2,gf256inv(y1)));
//			printf("\n %3ld is a root. Pos : %3ld Val : %3ld", x, position[errorCount], value[errorCount]);
            errorCount = errorCount + 1;
		}
	}
    if (errorCount != sigmaDegree) {
//        printf("\nMore than %3ld errors. Sigma degree: %3ld, number of roots: %3ld",maxDegree, sigmaDegree, errorCount);
        errorCount = -1;
    }
//	printf("\nError count %3ld ", errorCount);
    return(errorCount);

}

void makeRS255generator(unsigned char b[], long n)
{
    long            i,j;
	unsigned char   x, a;

	a = rs255const_a;
	x = rs255const_a111;
    for (i=0; i<=n; i++) {
		b[i] = 0;
	}
	b[0] = 1;
    for (i=0; i<n; i++) {
		x = gf256mul(x,a);
		for (j=i; j>=0; j--) {
			b[j+1] = gf256mul(x,b[j+1]) ^ b[j];
		}
		b[0] = gf256mul(x,b[0]);
	}
}

void rs255syndroms(unsigned char b[], unsigned char s[], long n)
{
    long            i;
	unsigned char   x,a,flag;

	a = rs255const_a;
	x = rs255const_a111;
	flag = 1;
    for (i=0; i<n; i++) {
		x = gf256mul(x,a);
		s[i] = polyEval(x, b, 255);
	}
}

/*Find r(x) and t(x) so that s(x)*x^n + t(x)*b(x) = r(x), deg(t(x))<=k+((n-k)/2)-1*/
long rs255euclid(unsigned char b[], unsigned char r[], unsigned char t[], long n, long k)
{
    long            i,j, x;
    long            maxRemainderDegree;
	unsigned char   q,flag;
	unsigned char   r0[256],r1[256],t0[256],t1[256];
	long            r0deg, r1deg, t0deg, t1deg, tDegree;

    for (i=0; i<256; i++) {
		r0[i] = 0;
		r1[i] = 0;
		t0[i] = 0;
		t1[i] = 0;
	}

    r0deg = n-1;
    for (i=0; i<n; i++) {
		r0[i] = b[i];
	}

    maxRemainderDegree = k + ((n+1-k)/2) - 1;
    r1deg = n;
	r1[n] = 1;
	x = r1deg;

//	printf("\nmaxRemainderDegree %3ld ", maxRemainderDegree);

	t0[0] = 1;
    t0deg = 0;
	t1deg = 0;

    flag = 1;
	do {
		switch(flag) {
            case 0 : { 
                while (r1[r1deg]==0) {r1deg = r1deg - 1;}   /*find r1 poly's degree*/
                if (r1deg<0) r1deg = 0;
                if (r1deg >= r0deg) {
                    q = gf256mul(r1[r1deg], gf256inv(r0[r0deg]));
					for (j=r0deg; j>=0; j--) {
                        r1[j+(r1deg-r0deg)] = r1[j+(r1deg-r0deg)] ^ gf256mul(q, r0[j]);
                    }
					for (j=(n-1); j>=0; j--) {
                        t1[j+(r1deg-r0deg)] = t1[j+(r1deg-r0deg)] ^ gf256mul(q, t0[j]);
                    }
                }
                else {
                    flag = 1;
                    x = r1deg;
                }
                break;
			}
            case 1 : { 
                while (r0[r0deg]==0) {r0deg = r0deg - 1;}   /*find r0 poly's degree*/
                if (r0deg<0) r0deg = 0;
                if (r0deg >= r1deg) {
                    q = gf256mul(r0[r0deg], gf256inv(r1[r1deg]));
					for (j=r1deg; j>=0; j--) {
                        r0[j+(r0deg-r1deg)] = r0[j+(r0deg-r1deg)] ^ gf256mul(q, r1[j]);
                    }
					for (j=(n-1); j>=0; j--) {
                        t0[j+(r0deg-r1deg)] = t0[j+(r0deg-r1deg)] ^ gf256mul(q, t1[j]);
                    }
                }
                else {
                    flag = 0;
                    x = r0deg;
                }
                break;
			}
		}
/*
        if (abs(r0deg-r1deg)>3)
            if ((r0deg>16) && (r1deg>16)) {
                printf("\nr0deg, r1deg : %3ld %3ld",r0deg,r1deg);
                debugFlag = 1;
            }
        if (k==100) {
            printf("\n flag : %3ld ",flag);
            printf("\n degr : %3ld ",x);
            printf("\n q    : %3ld ",q);
            printf("\n r0deg: %3ld ",r0deg);
            printf("\n r1deg: %3ld ",r1deg);
            printf("\n r0: ");
            polyPrint(r0, n+1);
            printf("\n r1: ");
            polyPrint(r1, n+1);
            printf("\n t0: ");
            polyPrint(t0, n+1);
            printf("\n t1: ");
            polyPrint(t1, n+1);
            printf("\nEnter 0-3 : ");
            scanf( "%ld", &i );
        }

*/
	} while (x > maxRemainderDegree);

	switch(flag) {
        case 0 : { 
			for (j=(n-1); j>=0; j--) {
                r[j] = r0[j];
                t[j] = t0[j];
            }
            break;
		}
        case 1 : { 
			for (j=(n-1); j>=0; j--) {
                r[j] = r1[j];
                t[j] = t1[j];
            }
            break;
		}
    }
    tDegree = (n-1);
    while ((t[tDegree]==0) && (tDegree>0)) tDegree = tDegree - 1;
    return(tDegree);

}

void rs255encode(unsigned char b[], long n)
{
    long            i,j;
	unsigned char   a,g[256], p[256];

	makeRS255generator(g, n);
	for (i=(254-n); i>=0; i--) {
		p[i+n] = b[i];
		b[i+n] = b[i];
	}
	for (i=(n-1); i>=0; i--) {
		p[i] = 0;
	}
	for (i=254; i>=n; i--) {
		a = p[i];
	    for (j=0; j<=n; j++) {
            p[i-j] = p[i-j] ^ gf256mul(a,g[n-j]);
		}
	}
	for (i=(n-1); i>=0; i--) {
		b[i] = p[i];
	}
}

long rs255erasure(unsigned char u[], unsigned char s[])
{
    long            i,j, count;
    unsigned char   x;
    unsigned char   ep[256];

    count=0;
    for (i=0; i<255; i++) {
        if (u[i]!=0) {
            ep[count] = (unsigned char) i;
            count++;
        }
        s[i]=0;
    }
    s[0]=1;
    for (i=0; i<count; i++) {
		x = gf256expa((unsigned char)(255-ep[i]));
		for (j=i; j>=0; j--) {
			s[j+1] = gf256mul(x,s[j+1]) ^ s[j];
		}
		s[0] = gf256mul(x,s[0]);
    }
    return(count);
}

long polyMultiply(unsigned char p1[], unsigned char p2[], unsigned char p3[], long p1deg, long p2deg)
{
    long            i,j;

    for (i=0; i<255; i++) {
        p3[i]=0;
    }
    for (i=0; i<=p1deg; i++) {
        for (j=0; j<=p2deg; j++) {
            p3[i+j] = p3[i+j] ^ gf256mul(p1[i],p2[j]);
        }
    }
    return(p1deg+p2deg);
}

/*
  Title: rs255decode

  Description: 
    Decodes an RS255 codeword, 'r', that may or may not have errors and/or erasures. Erasure positions
    are indicated by non-zero elements in the 'u' array. All errors and erasures are corrected if
    (2*errorCount + erasureCount) <= syndromeCount, the number of parity symbols in the codeword.
  
  Input:
    r[]           255bytes RS255 codeword + errors + erasures
    u[]           255bytes Erasure positions indicated by non-zero elements
    syndromeCount long     Number of parity symbols in the codeword

  Output:
    rs255decode   long     Number of errors/erasures corrected, -1 if decoding error
    c[]           255bytes Resulting codeword from the decoding process

  Local variables:
    erasureCount  long     Number of erasure symbols (as indicated by u[])
	errorCount    long     Number of error symbols as detected by the RS decoding algorithm
	errorEraseCount 
	              long     The sum of erasure and error symbols
	errorSigma[]  255bytes The error sigma polynomial. Output from the Euclids algorithm
	erasureSigma[] 
	              255bytes The erasure sigma polynomial calculated from u[]
	erasureSyndrome[] 
	              255bytes The product of the syndrome and the erasure sigma polynomials
	syndrome[]    255bytes The syndrome values of r[]
	sigma[]       255bytes The product of the erasureSigma and the errorSigma polynomials.
	omega[]       255bytes The error locator polynomial
	ePos[]        255bytes The error positions (indices)
	eVal[]        255bytes The error values

  Author: Jens Joergen Nielsen

  Rev. Date       Init. Comment
  1.00 2001-08-19 JJN   First release


*/
long rs255decode(unsigned char r[], unsigned char u[], unsigned char c[], long syndromeCount)
{
    long            i,erasureCount, errorCount, errorEraseCount;
    unsigned char   errorSigma[256], erasureSigma[256], erasureSyndrome[256];
	unsigned char   syndrome[256], omega[256], sigma[256], ePos[256], eVal[256];

    erasureCount = rs255erasure(u, erasureSigma);
    rs255syndroms(r,syndrome,syndromeCount);
    polyMultiply(erasureSigma, syndrome, erasureSyndrome, erasureCount, syndromeCount);
    errorCount = rs255euclid(erasureSyndrome, omega, errorSigma, syndromeCount, erasureCount);
    polyMultiply(erasureSigma, errorSigma, sigma, erasureCount, errorCount);
//	polyPrint(erasureSigma, erasureCount+1);
//	polyPrint(errorSigma, errorCount+1);
//	polyPrint(sigma, errorCount+erasureCount+1);
    errorEraseCount = rs255errors(sigma, omega, ePos, eVal, (errorCount+erasureCount));

	for (i=0; i<255; i++) c[i] = r[i];
	for (i=0; i<errorEraseCount; i++) {
		c[ePos[i]] = r[ePos[i]] ^ eVal[i];
	}
    return(errorEraseCount);
}

