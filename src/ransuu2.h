/* 
   A C-program for MT19937, with initialization improved 2002/1/26.
   Coded by Takuji Nishimura and Makoto Matsumoto.

   Before using, initialize the state by using init_genrand(seed)  
   or init_by_array(init_key, key_length).

   Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
   All rights reserved.                          

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions
   are met:

     1. Redistributions of source code must retain the above copyright
        notice, this list of conditions and the following disclaimer.

     2. Redistributions in binary form must reproduce the above copyright
        notice, this list of conditions and the following disclaimer in the
        documentation and/or other materials provided with the distribution.

     3. The names of its contributors may not be used to endorse or promote 
        products derived from this software without specific prior written 
        permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


   Any feedback is very welcome.
   http://www.math.keio.ac.jp/matumoto/emt.html
   email: matumoto@math.keio.ac.jp
*/

// #include <stdio.h>

/* Period parameters */  
#define MTN 624
#define MTM 397
#define MATRIX_A 0x9908b0dfUL   /* constant vector a */
#define UPPER_MASK 0x80000000UL /* most significant w-r bits */
#define LOWER_MASK 0x7fffffffUL /* least significant r bits */

//	#define PI 3.141592654  /* personal option */

static unsigned long mt[MTN]; /* the array for the state vector  */
static int mti=MTN+1; /* mti==MTN+1 means mt[MTN] is not initialized */

/* initializes mt[MTN] with a seed */
void init_genrand(unsigned long s)
{
    mt[0]= s & 0xffffffffUL;
    for (mti=1; mti<MTN; mti++) {
        mt[mti] = 
	    (1812433253UL * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti); 
        /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
        /* In the previous versions, MSBs of the seed affect   */
        /* only MSBs of the array mt[].                        */
        /* 2002/01/09 modified by Makoto Matsumoto             */
        mt[mti] &= 0xffffffffUL;
        /* for >32 bit machines */
    }
}

/* initialize by an array with array-length */
/* init_key is the array for initializing keys */
/* key_length is its length */
void init_by_array(unsigned long init_key[], long key_length)
{
    int i, j, k;
    init_genrand(19650218UL);
    i=1; j=0;
    k = (MTN>key_length ? MTN : key_length);
    for (; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1664525UL))
          + init_key[j] + j; /* non linear */
        mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++; j++;
        if (i>=MTN) { mt[0] = mt[MTN-1]; i=1; }
        if (j>=key_length) j=0;
    }
    for (k=MTN-1; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1566083941UL))
          - i; /* non linear */
        mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++;
        if (i>=MTN) { mt[0] = mt[MTN-1]; i=1; }
    }

    mt[0] = 0x80000000UL; /* MSB is 1; assuring non-zero initial array */ 
}


/* generates a random number on [0,0xffffffff]-interval */
unsigned long genrand_int32(void)
{
    unsigned long y;
    static unsigned long mag01[2]={0x0UL, MATRIX_A};
    /* mag01[x] = x * MATRIX_A  for x=0,1 */

    if (mti >= MTN) { /* generate MTN words at one time */
        int kk;

        if (mti == MTN+1)   /* if init_genrand() has not been called, */
            init_genrand(5489UL); /* a default initial seed is used */	/* SEED !!!!!!!! */

        for (kk=0;kk<MTN-MTM;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+MTM] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        for (;kk<MTN-1;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+(MTM - MTN)] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        y = (mt[MTN-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
        mt[MTN-1] = mt[MTM-1] ^ (y >> 1) ^ mag01[y & 0x1UL];

        mti = 0;
    }
  
    y = mt[mti++];

    /* Tempering */
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return y;
}

/* generates a random number on [0,0x7fffffff]-interval */
long genrand_int31(void)
{
    return (long)(genrand_int32()>>1);
}

/* generates a random number on [0,1]-real-interval */
double genrand_real1(void)
{
    return genrand_int32()*(1.0/4294967295.0); 
    /* divided by 2^32-1 */ 
}

/* generates a random number on [0,1)-real-interval */
double genrand_real2(void)
{
    return genrand_int32()*(1.0/4294967296.0); 
    /* divided by 2^32 */
}

/* generates a random number on (0,1)-real-interval */
double genrand_real3(void)
{
    return (((double)genrand_int32()) + 0.5)*(1.0/4294967296.0); 
    /* divided by 2^32 */
}

/* generates a random number on [0,1) with 53-bit resolution*/
double genrand_res53(void) 
{ 
    unsigned long a=genrand_int32()>>5, b=genrand_int32()>>6; 
    return(a*67108864.0+b)*(1.0/9007199254740992.0); 
} 

	//  for Poisson    (personal)

	// ---------- Poisson random number -------------------
	//				using MT19937ar

int poidev(float xm)
{
	float gammln(float xx);
	static float sq,alxm,g,oldm=(-1.0);
	float em,t,y;

	if (xm < 12.0) {
	    if (xm != oldm) {
		oldm=xm;
		g=exp(-xm);
	    }
	    em = -1;
	    t=1.0;
	    do {
		++em;
		t *= genrand_real1();
	    } while (t > g);
	} else {
	    if (xm != oldm) {
		oldm=xm;
		sq=sqrt(2.0*xm);
		alxm=log(xm);
		g=xm*alxm-gammln(xm+1.0);
	    }
	    do {
		do {
		    y=tan(PI*genrand_real1());
		    em=sq*y+xm;
	        } while (em < 0.0);
		em=floor(em);
		t=0.9*(1.0+y*y)*exp(em*alxm-gammln(em+1.0)-g);
	    } while (genrand_real1() > t);
	}
	return em;
}

float gammln(float xx)
{
	double x,y,tmp,ser;
	static double cof[6]={76.18009172947146,-86.50532032941677,
	    24.01409824083091,-1.231739572450155,
	    0.1208650973866179e-2,-0.5395239384953e-5};
	int j;

	y=x=xx;
	tmp=x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser=1.000000000190015;
	for (j=0;j<=5;j++) ser += cof[j]/++y;
	return -tmp+log(2.5066282746310005*ser/x);
}

//
//
//


int multirand(int k, double p0[])
{
	int i,x;
	double xx, x0, x1, xk;
	
	x0=0.0;xk=1.0;x=0;
	xx=genrand_real1();
	if(xx==0.0)
	{
		x=0;
	}else
	if(xx<0.5)
	{
		for(i=0;i<k;i++)
		{
			x1=x0+p0[i];
			if((xx>x0)&&(xx<=x1))
			{
				x=i;
				break;
			}else
			{
				x0=x1;
			};
		};
	}else
	{
		for(i=k-1;i>0;i--)
		{
			x1=xk-p0[i];
			if((xx>x1)&&(xx<=xk))
			{
				x=i;
				break;
			}else
			{
				xk=x1;
			};
		};
	};
	
	return x;
}

int multirand2(int k, double p0[])
{
	int i,x=0;
	double xx, x1;
	
	xx=genrand_real1();
	if(xx<0.5)
	{
		x=0;
		x1=0.0;i=0;
		while(xx>x1)
		{	
			x=i;	
			x1 += p0[i];
			i++;
		};
	}else
	{
		x1=1.0;i=k-1;
		while(xx<=x1)
		{
			x=i;
			x1 -= p0[i];
			i--;
		};
	};
	
	return x;
}
