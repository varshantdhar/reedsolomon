
/************************************************************************
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

#include "gf256.h"
#include "rs255.h"

static unsigned char gf256invTable[256];
static unsigned char gf256mulTable[256][256];

/*Calculate x=a*b using table look-up*/
unsigned char gf256mul(unsigned char a, unsigned char b)
{
    return(gf256mulTable[a][b]);
}


/*Calculate x=a*b using polynomial multiplication mod q*/
unsigned char gf256mulOld(unsigned char a, unsigned char b)
{
    long            i;
    unsigned char   x, y, q, msb;
	
	x = 0;
	y = b;
	q = gf256generator;
    for (i=7; i>=0; i--) {
		msb = x & 0x80;
		x = x + x;
		if (msb == 0x80) {
			x = x ^ q;
		}
		msb = y & 0x80;
		if (msb == 0x80) {
			x = x ^ a;
		}
		y = y + y;
	}
	return(x);
}

/*Calculate the inverse using table look-up*/
unsigned char gf256inv(unsigned char a)
{
	return(gf256invTable[a]);
}

/*Make 256 by 256 multiplication table*/
void initGF256mulTable()
{
    long   i,j;
	
    for (i=0; i<=255; i++) {
        for (j=0; j<=255; j++) {
            gf256mulTable[i][j] = gf256mulOld((unsigned char) i,(unsigned char) j);
        }
    }
}

/*Calculate the inverse using 1/x = x^254, 254=1111.1110b*/
void initGF256invTable()
{
    long            i,j;
    unsigned char   x,a;
	
    for (j=1; j<=255; j++) {
        a = (unsigned char) j;
        x = a;
        for (i=6; i>=1; i--) {  
    		x = gf256mul(x,x);
	    	x = gf256mul(x,a);
        }
        gf256invTable[a] = gf256mul(x,x);
    }
    gf256invTable[0] = 0;
}

void initGF256tables()
{
    initGF256mulTable();
    initGF256invTable();
}

/*Calculate log(x) with base rs255const_a*/
unsigned char gf256loga(unsigned char x)
{
    long            i;
    unsigned char   a, b, y;
	
	a = rs255const_a;
    b = 1;
    y = 255;
    i = 0;
    for (i=0; i<=255; i++) {
        if (b == x) y = (unsigned char) i;
		b = gf256mul(b,a);
	}
	return(y);
}

/*Calculate exp(x) with base rs255const_a,  exp(x)=a^x)*/
unsigned char gf256expa(unsigned char x)
{
    long            i;
    unsigned char   a, b;
	
	a = rs255const_a;
    b = 1;
    for (i=0; i<x; i++) {
		b = gf256mul(b,a);
	}
	return(b);
}

/*Calculate x^111*/
unsigned char gf256x111(unsigned char x)
{
    long            i;
    unsigned char   a,y;
	
	y = 1;
    a = 111;
    for (i=7; i>=0; i--) {
		y = gf256mul(y,y);
        if ((a & 0x80) == 0x80) y = gf256mul(x,y);
        a = a + a;
	}
	return(y);
}

