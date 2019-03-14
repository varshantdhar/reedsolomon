#define  gf256generator   0x87; /*GF(256) field generator polynomium g(x) = x^8 + x^7 + x^2 + x + 1*/

unsigned char gf256mul(unsigned char a, unsigned char b);
unsigned char gf256inv(unsigned char a);
unsigned char gf256loga(unsigned char x);
unsigned char gf256expa(unsigned char x);
unsigned char gf256x111(unsigned char x);
void initGF256tables();
