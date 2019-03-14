//
//  test.c
//  
//
//  Created by Varshant Dhar on 3/1/19.
//

#include <math.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "rs255.h"
#include "gf256.h"

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
 **/

/*Compare the two arrays and return the number of non equal elements*/
long rs255compare(unsigned char x[], unsigned char y[], long n)
{
    long            i;
    unsigned char   a;
    
    a = 0;
    for (i=0; i<n; i++) {
        if (x[i] != y[i]) {
            a++;
           // printf("\nerror in position %3ld %3hhu %3hhu ",i,x[i],y[i]);
        }
    }
    return(a);
}

/*Return an array of 255 elements, each element different from 0 with a probability of
 (px+py). If different from 0 all elements are equally likely. A fraction of the non-
 zero elements (corresponding to the probability py) are termed erasures and a "1" in
 the y-array points to the "erased" symbol.
 */
long rs255channel(unsigned char x[], unsigned char y[], float px, float py)
{
    long            i,k;
    unsigned char   a;
    int             n,m;
    
    k = 0;
    a = 0;
    for (i=0; i<255; i++) {
        x[i]=0;
        y[i]=0;
        n = rand();
        if (n < (32768*px)) {
            m = rand();
            x[i] = m / 128;
            if (x[i]>0) k = k + 1;
        }
        if (((32768*px) <= n) && (n < (32768*(px+py)))) {
            m = rand();
            x[i] = m / 128;
            if (x[i]>0) {
                y[i]=1;
            }
        }
    }
    return(k);
}


void append(char* s, char c)
{
    int len = strlen(s);
    s[len] = c;
    s[len+1] = '\0';
}

void randomErrorErasureTest()
{
    long            i,j,k;
    float           errorProbability, erasureProbability;
    long            paritySymbolCount, correctableErrorCount;
    unsigned char   x[256],y[256],z[256],e[256],u[256];
    long            errorsInserted, erasuresInserted;
    long            errorIndex, erasureIndex;
    long            errorsCorrected, errorsFound;
    long            errorEvent;
    long            correctableEventsCorrected;
    long            correctableEventsNotCorrected;
    long            notCorrectableEventsCorrected;
    long            notCorrectableEventsNotCorrected;
    unsigned char   c[256];
    float           alpha;
    long            testCount;
    FILE            *fp;
    
    char str[256] = "zdiscovery";
    testCount = 0;
    fp = fopen("test.csv", "r+");
    memcpy(x, str, strlen(str));
    
    for (int i = strlen("zdiscovery"); i < 201; i++) {
        errorsFound = 0;
        paritySymbolCount = 256 - i;
        
        erasureProbability = (float) 1.0;
        errorProbability   = (float) 0.0;
        
        correctableEventsCorrected = 0;
        correctableEventsNotCorrected = 0;
        notCorrectableEventsCorrected = 0;
        notCorrectableEventsNotCorrected = 0;
        
        for (int k = 0; k < 20; k++) {
            rs255encode(x, paritySymbolCount);
            
            errorsInserted = rs255channel(e, u, errorProbability, erasureProbability);
            
            erasuresInserted=0;
            for (j=0; j<255; j++) {
                y[j] = e[j] ^ x[j];
                if (u[j] !=0)
                    erasuresInserted++;
            }
            
            errorEvent = errorsInserted+erasuresInserted;
            
            errorsCorrected = rs255decode(y, u, z, paritySymbolCount);
            
            errorsFound += rs255compare(x,z,255);
            testCount++;
        
            if (errorsFound != 0) {
        
                if ((2*errorsInserted+erasuresInserted) <= paritySymbolCount) {
                    correctableEventsNotCorrected++;
                }
                else
                    notCorrectableEventsNotCorrected++;
            }
            else {

                if ((2*errorsInserted+erasuresInserted) <= paritySymbolCount)
                    correctableEventsCorrected++;
                else
                    notCorrectableEventsCorrected++;
            }
        }
        strcat(str, "C");
        memcpy(x, str, strlen(str));
        
        long RS_err = correctableEventsNotCorrected + notCorrectableEventsCorrected;
        fprintf(fp, "%ld,%ld,%ld\n", RS_err, testCount, strlen((char *) x) - 1);
        testCount = 0;
        RS_err = 0;
    }
    fclose(fp);
}


int main()
{
    randomErrorErasureTest();
    return 0;
    
}

