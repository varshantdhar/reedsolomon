#define  rs255const_a     0xac; /*p(x)=x^7+x^5+x^3+x^2+1 a primitive element in GF(256)*/
#define  rs255const_a111  0x0f; /*a^111*/

void rs255encode(unsigned char b[], long n);
long rs255decode(unsigned char r[], unsigned char u[], unsigned char c[], long nParitySymbols);
