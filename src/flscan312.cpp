#include <Rcpp.h>
using namespace Rcpp;

// flscan312.cpp
// Modified by Suzuryo, 2011/9/22. l.1302-1316
// Added by Suzuryo, 2011/9/28
// Modified by Otani, 2019

/* #pragma inline */

#include <stdlib.h>
//#include <stdio.h>
#include <math.h>
#include <time.h>
#include <memory.h>
#include <string.h>
//#pragma hdrstop

//#include "com.h"
//#include "ransuu2.h"
//#include "ranlib.h"

#define ErrMemory (1)
#define ErrFile0 (2)
#define ErrFile1 (3)
#define ErrFile2 (4)
#define ErrFile3 (5)
#define ErrFile4 (6)
#define ErrFile9 (7)
#define ErrNearest (8)
#define ErrSimLim (9)
#define ErrFile1ID (11)
#define ErrFile2ID (12)
#define ErrFile4ID (14)
#define ErrFile10 (20)
#define ErrFile11 (21)
//***Added by suzuryo, 2011.9.28
#define ErrFile12 (22)
//***end
//***Added by suzuryo, 2011.10.6
#define ErrFile13 (23)
//***end
//***Added by Otani, 2019
#define ErrFile14 (24)
//***end
#define ErrData (30)
#define Err99 (99)

#define KLIM (64)       /* limit of K */
#define SIMLIM (1000000) /* limit of Monte Carlo simulation */
#define PAI (3.141592653589793238)

char	*LAMBDAFILE = NULL;
char	*CLSTGEO = NULL;
//***Added by suzuryo, 2011/9/28.
char	*NODEFILE = NULL;
//***end
//***Added by suzuryo, 2011/10/6.
char	*EDGEFILE = NULL;
//***end
//***Added by Otani, 2019
char  *RFILE = NULL;
//***end

int     SCANMETHOD = 0; /* 0=Flexible, 1=Circular */

int     DEFSEED = 5489; /* seed of random number */
int     SIMCOUNT = 999; /* number of Monte Carlo simulation */
int     SIM;
typedef short areaidx;  /* type of area-index */

int		RANTYPE = 0;		/* 0=Multinomial, 1=Poisson */
int		EXTYPE = 0;		/* 0=Internal, 1=External */

int		MODEL = 0;		/* 0=Poisson, 1=Binomial */
int		STATTYPE = 0;		/* 0=LLR, 1=LLR (restricted PVALTYPE<RALPHA) */

float	RALPHA = 0.2;

int		CARTESIAN = 0;	/* 0:latitude&longitude, 1:Euclid distance */
double  R_EARTH = 6370;

int     K = 15;         /* maximum length of connected areas */
int		K2;
int     N;              /* number of all areas */

areaidx	*w;             /* w[] */
areaidx	z[KLIM];        /* z[] */
int     z_length;       /* length of z[] */

int		**a;            /* connection matrix a[N][N] */

FILE    *fp0, *fp1, *fp2, *fp3, *fp4, *fp9, *fp10, *fp11;
//***Added by suzuryo, 2011.9.28
FILE	*fp12;
//***end
//***Added by suzuryo, 2011.10.6
FILE	*fp13;
//***end
//***Added by otani, 2019
FILE  *fp14;
//***end

char    *filename1, *filename2, *filename3, *filename4, *filename9;

struct  TArea {         /* area information */
int     index;  /* index number (0..N-1) */
char    *id;    /* ID of the area */
double  popul;  /* no. of population */
int     cases;  /* no. of cases */
double  m, l;   /* coordinates */
double  dist;   /* distance from the center area */
};
TArea   *area;          /* area[N] */
TArea   **area_sorted;   /* sorted area[] */

int     **cases;        /* no. of cases for Monte Carlo simulations */
double  *popul;         /* copy of area[].popul */
int     *detectedarea;  /* flags to detect exclusive clusters */
double  *pp;			/* pp[i] is a probability for multinomial, i=0..N-1 */
int	*rtmp;

double  *maxstat;		/* test statistic */

int     nZ[SIMLIM];
double mZ;         /* sum of population and cases in a selected cluster */
int     *nG;
double mG;         /* nG, mG(cases) */
int		nGmax;
int     simworksize;

struct TLikelyCluster { /* information for the most likely cluster */
areaidx	z[KLIM];
  int     z_length;
  int     nZ;
  double mZ;
  double	lambda;
};
TLikelyCluster lkc;

int     phase;          /* calculation phase */

/* for Poisson model */

int		lors = 1;

double  **minmZ;        /* mimmZ[0..SIMCOUNT][0..snG[s]] */
int     *minmZ_zlength;
areaidx **minmZ_z;

int     MLC_zlength;
areaidx *MLC_z;

char	comment[256];

int     *calen;
areaidx **ca;
int     cazlen;
areaidx *caz;
short   *masksw;

/* for LLR with restriction */
double	**pv0;			/* pv0[0..N-1][0..SIMCOUNT] */

double	*Lpoi0;

/* for Binomial model */
double	*Lbin0;

time_t  tm1, tm2;

/* missing */
int		misarea;

/*---------------------------------------------------------------------------------*/
/* Sort maxstat[] in descending sequence */
int     sort_func3(const void *a, const void *b) {
  double  aa = *(double *)a, bb = *(double *)b;
  return((aa > bb) ? -1 : (aa < bb) ? 1 : 0);
}

/*---------------------------------------------------------------------------------*/
/* Sort areas by distance in ascending sequence */
int     sort_func0(const void *a, const void *b) {
  double  aa = (*(TArea **)a)->dist, bb = (*(TArea **)b)->dist;
  return((aa > bb) ? 1 : (aa < bb) ? -1 : 0);
}
/*---------------------------------------------------------------------------------*/
/* Sort z[] by area index in ascending sequence */
int sort_func1(const void *a, const void *b) {
  return(*(areaidx *)a - *(areaidx *)b);
}
/*---------------------------------------------------------------------------------*/
double  distance(double ma, double la, double mb, double lb) {
  double  xa, ya, za, xb, yb, zb, d, tab;
  
  if (CARTESIAN == 0) {
    if (ma == mb && la == lb)
      return(0.0);
    xa = cos(PAI / 180 * la) * cos(PAI / 180 * ma);
    ya = cos(PAI / 180 * la) * sin(PAI / 180 * ma);
    za = sin(PAI / 180 * la);
    xb = cos(PAI / 180 * lb) * cos(PAI / 180 * mb);
    yb = cos(PAI / 180 * lb) * sin(PAI / 180 * mb);
    zb = sin(PAI / 180 * lb);
    tab = xa * xb + ya * yb + za * zb;
    d = (tab == 1) ? 0 : R_EARTH * acos(tab);
  } else {
    d = sqrt((ma - mb) * (ma - mb) + (la - lb) * (la - lb));
  };
  return(d);
}
/*---------------------------------------------------------------------------------*/
double  MaxDistance(areaidx *z, int zlen, areaidx *z1, areaidx *z2) {
  double  maxdist = 0;
  double  m1, l1, d;
  int     i, j;
  
  if (zlen == 1) {
    *z1 = *z2 = z[0];
    return(maxdist);
  };
  for (i = 0; i < zlen; ++i) {
    m1 = area[z[i]].m;
    l1 = area[z[i]].l;
    for (j = 0; j < zlen; ++j) {
      d = distance(m1, l1, area[z[j]].m, area[z[j]].l);
      if (d > maxdist) {
        maxdist = d;
        *z1 = z[i];
        *z2 = z[j];
      };
    };
  };
  return(maxdist);
}
/*---------------------------------------------------------------------------------*/
void	ScanNearestNeighbours(int center) {
  int     i;
  double  m1, l1;
  
  m1 = area[center].m;
  l1 = area[center].l;
  for (i = 0; i < N; ++i)
    area[i].dist = distance(m1, l1, area[i].m, area[i].l);
  for (i = 0; i < N; ++i)
    area_sorted[i] = area+i;
  qsort((void *)area_sorted, N, sizeof(*area_sorted), sort_func0);
  for (i = 0; i < K; ++i)
    w[i] = area_sorted[i]->index;
  if (w[0] != center)
    Rcpp::stop("ERROR! Code:", ErrNearest);
}
/*---------------------------------------------------------------------------------*/
/* some functions for Restriction */
/*---------------------------------------------------------------------------------*/
/* Pr{X<z} foe N(0,1) */
double p_nor(double z) {
  int i;
  double z2, prev, p, t;
  
  if (z < -12.0) {
    return 0.0;
  } else if (z > 12.0) {
    return 1.0;
  } else {
    z2 = z * z;
    t = p = z * exp(-0.5 * z2) / sqrt(2 * PAI);
    for (i = 3; i < 200; i += 2) {
      prev = p;
      t *= z2 / i;
      p += t;
      if (p == prev) return 0.5 + p;
    }
    return (z > 0 ? 1.0 : 0.0);
  };
}
/*---------------------------------------------------------------------------------*/
/* n! */
double kaijo(int n) {
  int		s;
  double	tt;
  
  if (n == 0) {
    return 1.0;
  } else {
    tt = 0;
    for (s = 1; s <= n; ++s)
      tt += (double)log((double)s);
    return(exp(tt));
  };
}
/*---------------------------------------------------------------------------------*/
/* Peizer & Platt No.2 for Poisson */
double pplattP(int xx, double lambda) {
  double	x2, d1, d2, u2;
  
  
  x2 = (double)xx + 0.5;
  d1 = x2 - lambda + (1.0 / 6.0);
  d2 = d1 + 0.02 / (double)(xx + 1);
  u2 = (d2 / (fabs(x2 - lambda))) * sqrt(2.0 * x2 * log(x2 / lambda) + 2.0 * (lambda - x2));
  
  return(p_nor(u2));
}
/*---------------------------------------------------------------------------------*/
/* calculation of Pr{X>x}+0.5Pr{X=x} on Poisson Distribution */
double Ppfm(int ax, double ex) {
  double tmp0, tmp1;
  
  tmp0 = 1 - pplattP(ax, ex);
  
  if (ax > 0)
    tmp1 = pplattP(ax, ex) - pplattP(ax - 1, ex);
  else
    tmp1 = pplattP(0, ex);
  
  tmp0 = tmp0 + 0.5 * tmp1;
  
  return tmp0;
}
/*---------------------------------------------------------------------------------*/
/* Peizer & Platt No.2 for Binomial */
double pplattB(int xx, int bn, double bp) {
  double	x2, d1, d2, u2, xt0, xt1, xt2;
  double	bq = 1.0 - bp;
  
  x2 = (double)xx + 0.5;
  d1 = x2 - bn * bp  + (bq - bp) / (6.0);
  d2 = d1 + 0.02 * ((bq / (double)(xx + 1)) - (bp / (double)(bn - xx)) + ((bq - bp) / (2.0 * (bn + 1))));
  xt0 = log((double)(x2)) - log((double)(bn * bp));
  xt1 = log((double)(bn - x2)) - log((double)(bn * bq));
  xt2 = 2.0 / (1.0 + (1 / (6.0 * bn)));
  u2 = (d2 / (fabs(x2 - bn * bp))) * sqrt(xt2 * (x2 * xt0 + (bn - x2) * xt1));
  
  return(p_nor(u2));
}
/*---------------------------------------------------------------------------------*/
/* calculation of Pr{X>x}+0.5Pr{X=x} on Binomial Distribution */
double Pbfm(int ax, int bn, double bp) {
  double	tmp0, tmp1;
  
  if (ax > bn - 1)
    tmp0 = 0.0;
  else
    tmp0 = 1 - pplattB(ax, bn, bp);
  
  if (ax > bn - 1) {
    tmp1 = 1.0 - pplattB(ax - 1, bn, bp);
  } else if (ax > 0) {
    tmp1 = pplattB(ax, bn, bp) - pplattB(ax - 1, bn, bp);
  } else {
    tmp1 = pplattB(0, bn, bp);
  };
  
  tmp0 = tmp0 + 0.5 * tmp1;
  
  return tmp0;
}

/*---------------------------------------------------------------------------------*/
/* for original LLR / Poisson small nG */
/*---------------------------------------------------------------------------------*/
/* calculation of lambda and Monte Carlo simulation */
void	CalcLambda0s() {
  int     s;
  int     nZ1, nZ0;
  double  c1, c2, c3, lambda, nGf, mGf, nZf, mZf;
  
  nZ0 = -1;
  for (s = 0; s <= SIM; ++s) {
    nGf = (double)nG[s];
    mGf = (double)mG;
    c3 = log(nGf / mGf) * nGf;
    maxstat[s] = 0;
    for (nZ1 = 1; nZ1 <= nG[s]; ++nZ1) {
      if (minmZ[s][nZ1] == mG)
        continue;
      mZf = (double)minmZ[s][nZ1];
      nZf = (double)nZ1;
      c1 = nZf / mZf;
      c2 = (nGf - nZf) / (mGf - mZf);
      if (c1 > c2) {
        c1 = log(c1) * nZf;
        c2 = (c2 == 0) ? 0 : log(c2) * (nGf - nZf);
        lambda = c1 + c2 - c3;
        if (lambda > maxstat[s]) {
          maxstat[s] = lambda;
          if (s == 0)
            nZ0 = nZ1;
        };
      };
    };
  };
  
  if (nZ0 == (-1)) {
    lkc.z_length = 0;
    return;
  };
  
  lkc.z_length = minmZ_zlength[nZ0];
  for (s = 0; s < lkc.z_length; ++s)
    lkc.z[s] = minmZ_z[nZ0][s];
  lkc.lambda = maxstat[0];
  lkc.nZ = nZ0;
  lkc.mZ = minmZ[0][nZ0];
  
  qsort((void *)(&lkc.z[0]), lkc.z_length, sizeof(lkc.z[0]), sort_func1);
}
/*---------------------------------------------------------------------------------*/
void FlexibleScan0s(int zlen) {
  short int i, j;
  int     s;
  int     cazlensav;
  areaidx r;
  int     *cases_r;
  double  mZbak;
  
  if (detectedarea[z[zlen-1]] != 0)
    return;
  
  /* check lambda */
  for (s = 0; s <= SIM; ++s) {
    if (mZ < minmZ[s][nZ[s]]) {
      minmZ[s][nZ[s]] = mZ;
      if (s == 0) {
        for (i = 0; i < zlen; ++i)
          minmZ_z[nZ[s]][i] = z[i];
        minmZ_zlength[nZ[s]] = zlen;
      };
    };
  };
  
  if (zlen == 1) {
    for (i = 0; i < N; ++i)
      masksw[i] = -2;
    for (i = 0; i < K2; ++i)
      masksw[w[i]] = 0;
    masksw[z[0]] = -1;
    for (cazlen = 0; cazlen < calen[z[0]]; ++cazlen)
      caz[cazlen] = ca[z[0]][cazlen];
  };
  
  /* if zlen==K2 then end */
  if (zlen == K2)
    return;
  
  /* add each connected area caz[] to z[] */
  for (i = 0; i < cazlen; ++i) {
    /* if already included or finished then skip. */
    if (masksw[caz[i]] != 0) /* -1:included, 1,2,3..:finished */
    continue;
    /* add an area to z[] */
    z[zlen] = r = caz[i];
    masksw[r] = -1;
    /* update connected areas to this new cluster z[] */
    cazlensav = cazlen;
    
    for (j = 0; j < calen[r]; ++j) {
      if (masksw[ca[r][j]] == 0)
        caz[cazlen++] = ca[r][j];
    };
    
    mZbak = mZ;
    mZ += popul[r];
    cases_r = cases[r];
    for (s = 0; s <= SIM; ++s)
      nZ[s] += cases_r[s];
    /* reentrant */
    FlexibleScan0s(zlen + 1);
    
    mZ = mZbak;
    for (s = 0; s <= SIM; ++s)
      nZ[s] -= cases_r[s];
    
    /* mark finished */
    masksw[caz[i]] = zlen + 1;
    cazlen = cazlensav;
  };
  /* reset finished marks */
  for (i = 0; i < cazlen; ++i)
    if (masksw[caz[i]] == zlen + 1)
      masksw[caz[i]] = 0;
}
/*---------------------------------------------------------------------------------*/
void	CircularScan0s(int zlen) {
  int     i, r;
  int     s;
  int     *cases_r;
  double  mZbak;
  
  if (detectedarea[w[zlen-1]] != 0)
    return;
  
  for (s = 0; s <= SIM; ++s) {
    if (mZ < minmZ[s][nZ[s]]) {
      minmZ[s][nZ[s]] = mZ;
      if (s == 0) {
        for (i = 0; i < zlen; ++i)
          minmZ_z[nZ[s]][i] = w[i];
        minmZ_zlength[nZ[s]] = zlen;
      };
    };
  };
  
  if (zlen == K2)
    return;
  
  r = w[zlen];
  mZbak = mZ;
  mZ += popul[r];
  cases_r = cases[r];
  for (s = 0; s <= SIM; ++s)
    nZ[s] += cases_r[s];
  
  CircularScan0s(zlen + 1);
  
  mZ = mZbak;
  for (s = 0; s <= SIM; ++s)
    nZ[s] -= cases_r[s];
}
/*---------------------------------------------------------------------------------*/
/* for LLR with Restriction / Poisson model small nG*/
/*---------------------------------------------------------------------------------*/
void FlexibleScan1s(int zlen, int ss) {
  short int i, j;
  int		s;
  int     cazlensav;
  areaidx r;
  double  mZbak;
  
  s = ss;
  
  if (detectedarea[z[zlen-1]] != 0)
    return;
  
  /* check lambda */
  
  if (mZ < minmZ[s][nZ[s]]) {
    minmZ[s][nZ[s]] = mZ;
    if (s == 0) {
      for (i = 0; i < zlen; ++i)
        minmZ_z[nZ[s]][i] = z[i];
      minmZ_zlength[nZ[s]] = zlen;
    };
  };
  
  if (zlen == 1) {
    for (i = 0; i < N; ++i)
      masksw[i] = -2;
    for (i = 0; i < K2; ++i)
      masksw[w[i]] = 0;
    masksw[z[0]] = -1;
    for (cazlen = 0; cazlen < calen[z[0]]; ++cazlen)
      caz[cazlen] = ca[z[0]][cazlen];
  };
  
  /* if zlen==K2 then end */
  if (zlen == K2)
    return;
  
  /* add each connected area caz[] to z[] */
  for (i = 0; i < cazlen; ++i) {
    /* if already included or finished then skip. */
    if (masksw[caz[i]] != 0) /* -1:included, 1,2,3..:finished */
    continue;
    /* add an area to z[] */
    z[zlen] = r = caz[i];
    masksw[r] = -1;
    /* update connected areas to this new cluster z[] */
    cazlensav = cazlen;
    
    for (j = 0; j < calen[r]; ++j) {
      if (masksw[ca[r][j]] == 0)
        caz[cazlen++] = ca[r][j];
    };
    
    mZbak = mZ;
    mZ += popul[r];
    nZ[s] += cases[r][s];
    
    if (pv0[r][s] < RALPHA)
      FlexibleScan1s(zlen + 1, s);
    
    mZ = mZbak;
    nZ[s] -= cases[r][s];
    
    /* mark finished */
    masksw[caz[i]] = zlen + 1;
    cazlen = cazlensav;
  };
  /* reset finished marks */
  for (i = 0; i < cazlen; ++i)
    if (masksw[caz[i]] == zlen + 1)
      masksw[caz[i]] = 0;
}
/*---------------------------------------------------------------------------------*/
void	CircularScan1s(int zlen, int ss) {
  int     i, r;
  int     s;
  double  mZbak;
  
  s = ss;
  
  if (detectedarea[w[zlen-1]] != 0)
    return;
  
  if (mZ < minmZ[s][nZ[s]]) {
    minmZ[s][nZ[s]] = mZ;
    if (s == 0) {
      for (i = 0; i < zlen; ++i)
        minmZ_z[nZ[s]][i] = w[i];
      minmZ_zlength[nZ[s]] = zlen;
    };
  };
  
  if (zlen == K2)
    return;
  
  r = w[zlen];
  mZbak = mZ;
  mZ += popul[r];
  nZ[s] += cases[r][s];
  
  if (pv0[r][s] < RALPHA)
    CircularScan1s(zlen + 1, s);
  
  mZ = mZbak;
  nZ[s] -= cases[r][s];
}


/*---------------------------------------------------------------------------------*/
/* for original LLR / Poisson large nG */
/*---------------------------------------------------------------------------------*/
double calcstatP0(int bnZ, double bmZ, int bnG, double bmG) {
  double	c1, c2, c3, c4, c5;
  double	plmb;
  
  c1 = bnZ / bmZ;
  c2 = (double)(bnG - bnZ);
  c3 = c2 / (bmG - bmZ);
  c4 = log(c1) * bnZ;
  c5 = (c2 == 0) ? 0 : log(c3) * c2;
  plmb = c4 + c5;
  
  return(plmb);
}
/*---------------------------------------------------------------------------------*/
void FlexibleScan0l(int zlen) {
  short int i, j;
  int     s;
  int     cazlensav;
  areaidx r;
  int     *cases_r;
  double  mZbak;
  double	st0;
  
  
  if (detectedarea[z[zlen-1]] != 0)
    return;
  
  for (s = 0; s <= SIM; ++s) {
    if ((nZ[s] / mZ) > ((nG[s] - nZ[s]) / (mG - mZ)))
      st0 = calcstatP0(nZ[s], mZ, nG[s], mG) - Lpoi0[s];
    else
      st0 = 0;
    if (st0 > maxstat[s]) {
      maxstat[s] = st0;
      if (s == 0) {
        for (i = 0; i < zlen; ++i)
          MLC_z[i] = z[i];
        MLC_zlength = zlen;
      };
    };
  };
  
  if (zlen == 1) {
    for (i = 0; i < N; ++i)
      masksw[i] = -2;
    for (i = 0; i < K2; ++i)
      masksw[w[i]] = 0;
    masksw[z[0]] = -1;
    for (cazlen = 0; cazlen < calen[z[0]]; ++cazlen)
      caz[cazlen] = ca[z[0]][cazlen];
  };
  
  if (zlen == K2)
    return;
  
  /* add each connected area caz[] to z[] */
  for (i = 0; i < cazlen; ++i) {
    /* if already included or finished then skip. */
    if (masksw[caz[i]] != 0) /* -1:included, 1,2,3..:finished */
    continue;
    /* add an area to z[] */
    z[zlen] = r = caz[i];
    masksw[r] = -1;
    /* update connected areas to this new cluster z[] */
    cazlensav = cazlen;
    for (j = 0; j < calen[r]; ++j) {
      if (masksw[ca[r][j]] == 0)
        caz[cazlen++] = ca[r][j];
    };
    
    mZbak = mZ;
    mZ += popul[r];
    cases_r = cases[r];
    for (s = 0; s <= SIM; ++s)
      nZ[s] += cases_r[s];
    /* reentrant */
    FlexibleScan0l(zlen + 1);
    
    mZ = mZbak;
    for (s = 0; s <= SIM; ++s)
      nZ[s] -= cases_r[s];
    
    /* mark finished */
    masksw[caz[i]] = zlen + 1;
    cazlen = cazlensav;
  };
  /* reset finished marks */
  for (i = 0; i < cazlen; ++i)
    if (masksw[caz[i]] == zlen + 1)
      masksw[caz[i]] = 0;
}

/*---------------------------------------------------------------------------------*/
void	CircularScan0l(int zlen) {
  int     i, r;
  int     s;
  int		*cases_r;
  double  mZbak;
  double	st0;
  
  
  if (detectedarea[w[zlen-1]] != 0)
    return;
  
  for (s = 0; s <= SIM; ++s) {
    if ((nZ[s] / mZ) > ((nG[s] - nZ[s]) / (mG - mZ)))
      st0 = calcstatP0(nZ[s], mZ, nG[s], mG) - Lpoi0[s];
    else
      st0 = 0;
    if (st0 > maxstat[s]) {
      maxstat[s] = st0;
      if (s == 0) {
        for (i = 0; i < zlen; ++i)
          MLC_z[i] = w[i];
        MLC_zlength = zlen;
      };
    };
  };
  
  if (zlen == K2)
    return;
  
  r = w[zlen];
  mZbak = mZ;
  mZ += popul[r];
  cases_r = cases[r];
  for (s = 0; s <= SIM; ++s)
    nZ[s] += cases_r[s];
  
  CircularScan0l(zlen + 1);
  
  mZ = mZbak;
  for (s = 0; s <= SIM; ++s)
    nZ[s] -= cases_r[s];
}

/*---------------------------------------------------------------------------------*/
/* for LLR with Restriction / Poisson model large nG*/
/*---------------------------------------------------------------------------------*/
void FlexibleScan1l(int zlen, int ss) {
  short int i, j;
  int     s;
  int     cazlensav;
  areaidx r;
  double  mZbak;
  double	st0;
  
  s = ss;
  
  if (detectedarea[z[zlen-1]] != 0)
    return;
  
  if ((nZ[s] / mZ) > ((nG[s] - nZ[s]) / (mG - mZ)))
    st0 = calcstatP0(nZ[s], mZ, nG[s], mG) - Lpoi0[s];
  else
    st0 = 0;
  if (st0 > maxstat[s]) {
    maxstat[s] = st0;
    if (s == 0) {
      for (i = 0; i < zlen; ++i)
        MLC_z[i] = z[i];
      MLC_zlength = zlen;
    };
  };
  
  if (zlen == 1) {
    for (i = 0; i < N; ++i)
      masksw[i] = -2;
    for (i = 0; i < K2; ++i)
      masksw[w[i]] = 0;
    masksw[z[0]] = -1;
    for (cazlen = 0; cazlen < calen[z[0]]; ++cazlen)
      caz[cazlen] = ca[z[0]][cazlen];
  };
  
  if (zlen == K2)
    return;
  
  /* add each connected area caz[] to z[] */
  for (i = 0; i < cazlen; ++i) {
    /* if already included or finished then skip. */
    if (masksw[caz[i]] != 0) /* -1:included, 1,2,3..:finished */
    continue;
    /* add an area to z[] */
    z[zlen] = r = caz[i];
    masksw[r] = -1;
    /* update connected areas to this new cluster z[] */
    cazlensav = cazlen;
    for (j = 0; j < calen[r]; ++j) {
      if (masksw[ca[r][j]] == 0)
        caz[cazlen++] = ca[r][j];
    };
    
    mZbak = mZ;
    mZ += popul[r];
    nZ[s] += cases[r][s];
    
    /* reentrant */
    if (pv0[r][s] < RALPHA)
      FlexibleScan1l(zlen + 1, s);
    
    mZ = mZbak;
    nZ[s] -= cases[r][s];
    
    /* mark finished */
    masksw[caz[i]] = zlen + 1;
    cazlen = cazlensav;
  };
  /* reset finished marks */
  for (i = 0; i < cazlen; ++i)
    if (masksw[caz[i]] == zlen + 1)
      masksw[caz[i]] = 0;
}

void	CircularScan1l(int zlen, int ss) {
  int     i, r;
  int     s;
  double  mZbak;
  double	st0;
  
  s = ss;
  
  if (detectedarea[w[zlen-1]] != 0)
    return;
  
  if ((nZ[s] / mZ) > ((nG[s] - nZ[s]) / (mG - mZ)))
    st0 = calcstatP0(nZ[s], mZ, nG[s], mG) - Lpoi0[s];
  else
    st0 = 0;
  if (st0 > maxstat[s]) {
    maxstat[s] = st0;
    if (s == 0) {
      for (i = 0; i < zlen; ++i)
        MLC_z[i] = w[i];
      MLC_zlength = zlen;
    };
  };
  
  if (zlen == K2)
    return;
  
  r = w[zlen];
  mZbak = mZ;
  mZ += popul[r];
  nZ[s] += cases[r][s];
  
  if (pv0[r][s] < RALPHA)
    CircularScan1l(zlen + 1, s);
  
  mZ = mZbak;
  nZ[s] -= cases[r][s];
}






/*---------------------------------------------------------------------------------*/
/* for original LLR / Binomial  */
/*---------------------------------------------------------------------------------*/
double calcstatB0(int bnZ, double bmZ, int bnG, double bmG) {
  double	c1, c2, c3, c4;
  double	blmb;
  
  c1 = bnZ / bmZ;
  c2 = (double)(bnG - bnZ);
  c3 = bmG - bmZ;
  
  if (c1 >= 1.0)
    c4 = 0.0;
  else
    c4 = bnZ * log(c1) + (bmZ - bnZ) * log(1.0 - c1);
  
  blmb = c4 + c2 * log(c2 / c3) + (c3 - c2) * log((c3 - c2) / c3);
  
  return(blmb);
}
/*---------------------------------------------------------------------------------*/
void FlexibleScanB0(int zlen) {
  short int i, j;
  int     s;
  int     cazlensav;
  areaidx r;
  int     *cases_r;
  double  mZbak;
  double	st0;
  
  
  if (detectedarea[z[zlen-1]] != 0)
    return;
  
  for (s = 0; s <= SIM; ++s) {
    if ((nZ[s] / mZ) > ((nG[s] - nZ[s]) / (mG - mZ)))
      st0 = calcstatB0(nZ[s], mZ, nG[s], mG) - Lbin0[s];
    else
      st0 = 0;
    if (st0 > maxstat[s]) {
      maxstat[s] = st0;
      if (s == 0) {
        for (i = 0; i < zlen; ++i)
          MLC_z[i] = z[i];
        MLC_zlength = zlen;
      };
    };
  };
  
  if (zlen == 1) {
    for (i = 0; i < N; ++i)
      masksw[i] = -2;
    for (i = 0; i < K2; ++i)
      masksw[w[i]] = 0;
    masksw[z[0]] = -1;
    for (cazlen = 0; cazlen < calen[z[0]]; ++cazlen)
      caz[cazlen] = ca[z[0]][cazlen];
  };
  
  if (zlen == K2)
    return;
  
  /* add each connected area caz[] to z[] */
  for (i = 0; i < cazlen; ++i) {
    /* if already included or finished then skip. */
    if (masksw[caz[i]] != 0) /* -1:included, 1,2,3..:finished */
    continue;
    /* add an area to z[] */
    z[zlen] = r = caz[i];
    masksw[r] = -1;
    /* update connected areas to this new cluster z[] */
    cazlensav = cazlen;
    for (j = 0; j < calen[r]; ++j) {
      if (masksw[ca[r][j]] == 0)
        caz[cazlen++] = ca[r][j];
    };
    
    mZbak = mZ;
    mZ += popul[r];
    cases_r = cases[r];
    for (s = 0; s <= SIM; ++s)
      nZ[s] += cases_r[s];
    /* reentrant */
    FlexibleScanB0(zlen + 1);
    
    mZ = mZbak;
    for (s = 0; s <= SIM; ++s)
      nZ[s] -= cases_r[s];
    
    /* mark finished */
    masksw[caz[i]] = zlen + 1;
    cazlen = cazlensav;
  };
  /* reset finished marks */
  for (i = 0; i < cazlen; ++i)
    if (masksw[caz[i]] == zlen + 1)
      masksw[caz[i]] = 0;
}

/*---------------------------------------------------------------------------------*/
void	CircularScanB0(int zlen) {
  int     i, r;
  int     s;
  int		*cases_r;
  double  mZbak;
  double	st0;
  
  if (detectedarea[w[zlen-1]] != 0)
    return;
  
  for (s = 0; s <= SIM; ++s) {
    if ((nZ[s] / mZ) > ((nG[s] - nZ[s]) / (mG - mZ)))
      st0 = calcstatB0(nZ[s], mZ, nG[s], mG) - Lbin0[s];
    else
      st0 = 0;
    if (st0 > maxstat[s]) {
      maxstat[s] = st0;
      if (s == 0) {
        for (i = 0; i < zlen; ++i)
          MLC_z[i] = w[i];
        MLC_zlength = zlen;
      };
    };
  };
  
  if (zlen == K2)
    return;
  
  r = w[zlen];
  mZbak = mZ;
  mZ += popul[r];
  cases_r = cases[r];
  for (s = 0; s <= SIM; ++s)
    nZ[s] += cases_r[s];
  
  CircularScanB0(zlen + 1);
  
  mZ = mZbak;
  for (s = 0; s <= SIM; ++s)
    nZ[s] -= cases_r[s];
  
}

/*---------------------------------------------------------------------------------*/
/* for LLR with restriction / Binomial  */
/*---------------------------------------------------------------------------------*/
void FlexibleScanB1(int zlen, int ss) {
  short int i, j;
  int     s;
  int     cazlensav;
  areaidx r;
  double  mZbak;
  double	st0;
  
  s = ss;
  
  if (detectedarea[z[zlen-1]] != 0)
    return;
  
  if ((nZ[s] / mZ) > ((nG[s] - nZ[s]) / (mG - mZ)))
    st0 = calcstatB0(nZ[s], mZ, nG[s], mG) - Lbin0[s];
  else
    st0 = 0;
  if (st0 > maxstat[s]) {
    maxstat[s] = st0;
    if (s == 0) {
      for (i = 0; i < zlen; ++i)
        MLC_z[i] = z[i];
      MLC_zlength = zlen;
    };
  };
  
  if (zlen == 1) {
    for (i = 0; i < N; ++i)
      masksw[i] = -2;
    for (i = 0; i < K2; ++i)
      masksw[w[i]] = 0;
    masksw[z[0]] = -1;
    for (cazlen = 0; cazlen < calen[z[0]]; ++cazlen)
      caz[cazlen] = ca[z[0]][cazlen];
  };
  
  if (zlen == K2)
    return;
  
  /* add each connected area caz[] to z[] */
  for (i = 0; i < cazlen; ++i) {
    /* if already included or finished then skip. */
    if (masksw[caz[i]] != 0) /* -1:included, 1,2,3..:finished */
    continue;
    /* add an area to z[] */
    z[zlen] = r = caz[i];
    masksw[r] = -1;
    /* update connected areas to this new cluster z[] */
    cazlensav = cazlen;
    for (j = 0; j < calen[r]; ++j) {
      if (masksw[ca[r][j]] == 0)
        caz[cazlen++] = ca[r][j];
    };
    
    mZbak = mZ;
    mZ += popul[r];
    nZ[s] += cases[r][s];
    
    /* reentrant */
    if (pv0[r][s] < RALPHA)
      FlexibleScanB1(zlen + 1, s);
    
    mZ = mZbak;
    nZ[s] -= cases[r][s];
    
    /* mark finished */
    masksw[caz[i]] = zlen + 1;
    cazlen = cazlensav;
  };
  /* reset finished marks */
  for (i = 0; i < cazlen; ++i)
    if (masksw[caz[i]] == zlen + 1)
      masksw[caz[i]] = 0;
}

void	CircularScanB1(int zlen, int ss) {
  int     i, r;
  int     s;
  double  mZbak;
  double	st0;
  
  s = ss;
  
  if (detectedarea[w[zlen-1]] != 0)
    return;
  
  if ((nZ[s] / mZ) > ((nG[s] - nZ[s]) / (mG - mZ)))
    st0 = calcstatB0(nZ[s], mZ, nG[s], mG) - Lbin0[s];
  else
    st0 = 0;
  if (st0 > maxstat[s]) {
    maxstat[s] = st0;
    if (s == 0) {
      for (i = 0; i < zlen; ++i)
        MLC_z[i] = w[i];
      MLC_zlength = zlen;
    };
  };
  
  if (zlen == K2)
    return;
  
  r = w[zlen];
  mZbak = mZ;
  mZ += popul[r];
  nZ[s] += cases[r][s];
  
  if (pv0[r][s] < RALPHA)
    CircularScanB1(zlen + 1, s);
  
  mZ = mZbak;
  nZ[s] -= cases[r][s];
}




/*---------------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------------*/
/*  Main Calculation */
/*---------------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------------*/
void	FlexScan() {
  int     i, j, center, s, ind_k;
  //TLikelyCluster *lc;
  double  maxdist;
  areaidx maxdistz1, maxdistz2;
  double	l1, l2, m1, m2;
  
  int		rnk;
  
  int     SECONDARY = 64;	/* number of secondary clusters */
int		NOMORE = 0;		/* 1:nomore secondary cluster */

if (MODEL == 0 && lors == 0)
  for (i = 0; i <= nGmax; ++i)
    minmZ_zlength[i] = 0;

for (phase = 0; phase <= SECONDARY; ++phase) {
  if (phase == 0) {
    SIM = SIMCOUNT;
    Rprintf("*** SCANNING MOST LIKELY CLUSTER with %d MONTE CARLO REPLICATIONS ***\n", SIM);
  } else {
    SIM = 0;
    Rprintf("*** SCANNING SECONDARY CLUSTERS (%d) ***\n", phase);
  };
  
  if (MODEL == 0 && lors == 0) {
    for (s = 0; s <= SIM; ++s)
      for (i = 0; i <= nGmax; ++i)
        minmZ[s][i] = mG;
  };
  
  for (s = 0; s <= SIM; ++s)
    maxstat[s] = -1.0;
  
  for (center = 0; center < N; ++center) {
    if (detectedarea[center] == -1)
      continue;
    
    if (phase == 0) {
      if (center == 0)
        Rprintf("Scanning areas around %s", area[center].id);
      else
        Rprintf(", %s", area[center].id);
    };
    
    ScanNearestNeighbours(center); /* and set the list to w[] */

ind_k = 0;
for (K2 = 0; K2 < N; ++K2) {
  if (detectedarea[w[K2]] != -1)
    ++ind_k;
  
  if (ind_k > K)
    break;
};

z[0] = center;

if (MODEL == 0 && lors == 0) {	/* Poisson Model small nG */
if (STATTYPE == 1) {	/* LLR with restriction */
for (s = 0; s <= SIM; ++s) {
  if (pv0[center][s] < RALPHA) {
    mZ = popul[center];
    nZ[s] = cases[center][s];
    
    if (SCANMETHOD == 1)
      CircularScan1s(1, s);
    else
      FlexibleScan1s(1, s);
  };
};
} else {	/* LLR original */
mZ = popul[center];
  for (s = 0; s <= SIM; ++s)
    nZ[s] = cases[center][s];
  
  if (SCANMETHOD == 1)
    CircularScan0s(1);
  else
    FlexibleScan0s(1);
};
} else if (MODEL == 0 && lors == 1) { /* Poisson Model large nG */
if (STATTYPE == 1) {	/* LLR with restriction */
for (s = 0; s <= SIM; ++s) {
  if (pv0[center][s] < RALPHA) {
    mZ = popul[center];
    nZ[s] = cases[center][s];
    if (SCANMETHOD == 1)
      CircularScan1l(1, s);
    else
      FlexibleScan1l(1, s);
  };
};
} else {	/* LLR original */
mZ = popul[center];
  for (s = 0; s <= SIM; ++s)
    nZ[s] = cases[center][s];
  
  if (SCANMETHOD == 1)
    CircularScan0l(1);
  else
    FlexibleScan0l(1);
};
} else if (MODEL == 1) {	/* Binomial Model */
if (STATTYPE == 1) {	/* LLR with restriction */
for (s = 0; s <= SIM; ++s) {
  if (pv0[center][s] < RALPHA) {
    mZ = popul[center];
    nZ[s] = cases[center][s];
    
    if (SCANMETHOD == 1)
      CircularScanB1(1, s);
    else
      FlexibleScanB1(1, s);
  };
};
} else {	/* LLR original */
mZ = popul[center];
  for (s = 0; s <= SIM; ++s)
    nZ[s] = cases[center][s];
  
  if (SCANMETHOD == 1)
    CircularScanB0(1);
  else
    FlexibleScanB0(1);
};
};

  };
  
  if (phase == 0)
    Rprintf("\n");
  
  if (MODEL == 0 && lors == 0) {
    CalcLambda0s();
  };
  
  if (MODEL == 1 || lors == 1) {
    lkc.z_length = MLC_zlength;
    lkc.nZ = 0;
    lkc.mZ = 0;
    for (j = 0; j < lkc.z_length; ++j) {
      lkc.z[j] = MLC_z[j];
      lkc.nZ += cases[lkc.z[j]][0];
      lkc.mZ += popul[lkc.z[j]];
    };
    lkc.lambda = maxstat[0];
    
    qsort((void *)(&lkc.z[0]), lkc.z_length, sizeof(lkc.z[0]), sort_func1);
  };
  
  
  for (j = 0; j < lkc.z_length; ++j)
    detectedarea[lkc.z[j]] = 1;
  
  if (phase == 0 && LAMBDAFILE != NULL) {
    if ((fp13 = fopen(NODEFILE, "wt")) == NULL)
      Rcpp::stop("ERROR! Code:", ErrFile13);
    if ((fp10 = fopen(LAMBDAFILE, "wt")) == NULL)
      Rcpp::stop("ERROR! Code:", ErrFile10);
    if ((fp14 = fopen(RFILE, "wt")) == NULL)
      Rcpp::stop("ERROR! Code:", ErrFile14);
    fprintf(fp10, "%d\n", SIMCOUNT + 1);
    for (i = 0; i <= SIMCOUNT; ++i)
      fprintf(fp10, "%.22g\n", maxstat[i]);
    l1 = l2 = area[0].l;
    m1 = m2 = area[0].m;
    for (i = 1; i < N; ++i) {
      if (area[i].l < l1) l1 = area[i].l;
      if (area[i].m < m1) m1 = area[i].m;
      if (area[i].l > l2) l2 = area[i].l;
      if (area[i].m > m2) m2 = area[i].m;
    };
    fprintf(fp10, "%.22g %.22g %.22g %.22g\n", l1, m1, l2, m2);
    for (i = 0; i < N; ++i) {
      for (j = i + 1; j < N; ++j)
        if (a[i][j]) {
          //***Buggy block starat
          //fprintf(fp10,"%d %.22lg %.22lg %.22lg %.22lg\n",
          //	(detectedarea[i] && detectedarea[j]) ? 1 : 0,
          //	area[i].l, area[i].m,
          //	area[j].l, area[j].m
          //);
          //***end
          
          //***Modified by suzuryo, 2011/9/22.
          fprintf(fp10,"%d %.22g %.22g %.22g %.22g\n",
                  (detectedarea[i]==1 && detectedarea[j]==1) ? 1 : 0,
                  area[i].l, area[i].m,
                  area[j].l, area[j].m
          );
          //***end
          
          //***Added by suzuryo, 2011/10/6.
          fprintf(fp12,"%d %d\n",area[i].index,area[j].index);
          //***end
        };
    };
    fprintf(fp10, "-\n");
    for (i = 0; i < N; ++i)
      fprintf(fp10, "%d %s %d %.22g %.22g\n", i, area[i].id, detectedarea[i], area[i].l, area[i].m);
    fclose(fp10);
  };
  
  /* output the results */
  if (phase == 0) {
    fprintf(fp9, "-------------------------------------------------------\n");
    fprintf(fp9, "        FleXScan ver3.1.2 -- purely spatial      \n");
    fprintf(fp9, "-------------------------------------------------------\n");
    fprintf(fp9, "\n");
    
    if (comment[0] != '\0')
      fprintf(fp9, "Comment: %s\n\n", comment);
    
    if (SCANMETHOD == 1)
      fprintf(fp9, "Scanning method: Circular spatial scan.\n\n");
    else
      fprintf(fp9, "Scanning method: Flexible spatial scan.\n\n");
    
    fprintf(fp9, "<STATISTICS>\n");
    fprintf(fp9, " %s.",
            (STATTYPE == 1) ? "Log likelihood ratio with restriction" : "Original log likelihood ratio");
    if (STATTYPE == 1)
      fprintf(fp9, "\n (Pr{X>x}+0.5*Pr{X=x}<%f)\n\n", RALPHA);
    else
      fprintf(fp9, "\n\n");
    
    fprintf(fp9, "Program run on: %s", ctime(&tm1));
    fprintf(fp9, "\n");
    fprintf(fp9, "Purely Spatial analysis\n");
    fprintf(fp9, "scanning for clusters with\n");
    fprintf(fp9, "high rates using the %s model.\n", (MODEL == 1) ? "Binomial" : "Poisson");
    fprintf(fp9, "\n");
    fprintf(fp9, "--------------------------------------------------------\n");
    fprintf(fp9, "\n");
    fprintf(fp9, "SUMMARY OF DATA\n");
    fprintf(fp9, "\n");
    fprintf(fp9, "Limit length of cluster: %d\n", K);
    fprintf(fp9, "Number of census areas.: %d\n", N - misarea);
    if (misarea != 0) {
      fprintf(fp9, "  (Number of excluded area: %d : ", misarea);
      for (i = 0; i < N; ++i) {
        if (detectedarea[i] == -1)
          fprintf(fp9, "%s ", area[i].id);
      };
      fprintf(fp9, "\n");
    };
    fprintf(fp9, "Total cases ...........: %d\n", nG[0]);
    if (MODEL == 1)
      fprintf(fp9, "Total population ......: %d\n", (int)mG);
    if (MODEL == 0) {
      if (EXTYPE == 1)
        fprintf(fp9,  " (expected number is adjusted by Total cases as *Total expects* = *Total cases*)\n");
    }
      fprintf(fp9, "\n");
      fprintf(fp9, "--------------------------------------------------------\n");
      fprintf(fp9, "\n");
      fprintf(fp9, "MOST LIKELY CLUSTER\n");
      fprintf(fp9, "\n");
      
      qsort((void *)(&maxstat[1]), SIMCOUNT, sizeof(maxstat[1]), sort_func3);
  } else if (phase == 1) {
    fprintf(fp9, "--------------------------------------------------------\n");
    fprintf(fp9, "\n");
    fprintf(fp9, "SECONDARY CLUSTERS\n");
    fprintf(fp9, "\n");
  };
  
  if (lkc.z_length == 0) {
    if (phase == 0) {
      Rprintf("*** There is no cluster ***\n");
      fprintf(fp9, "*** There is no cluster ***\n\n");
    } else {
      Rprintf("*** There are no more secondary clusters ***\n");
      fprintf(fp9, "*** There are no more secondary clusters ***\n\n");
    };
    
    NOMORE = 1;
    
  };
  
  if (NOMORE == 0) {
    fprintf(fp9, "%d.Census areas included .: ", phase + 1);
    fprintf(fp9, "%s", area[lkc.z[0]].id);
    fprintf(fp11, "%s\t%d\t%f\t%f\n", area[lkc.z[0]].id, phase + 1, area[lkc.z[0]].l, area[lkc.z[0]].m);
    //***added by suzuryo, 2011.9.28
    fprintf(fp13, "%d\t%d\t%f\t%f\n", lkc.z[0], phase + 1, area[lkc.z[0]].l, area[lkc.z[0]].m); // Write out numbers rather than IDs
    //***end
    for (j = 1; j < lkc.z_length; ++j) {
      fprintf(fp9, ", %s", area[lkc.z[j]].id);
      fprintf(fp11, "%s\t%d\t%f\t%f\n", area[lkc.z[j]].id, phase + 1, area[lkc.z[j]].l, area[lkc.z[j]].m);
      //***added by suzuryo, 2011.9.28
      fprintf(fp13, "%d\t%d\t%f\t%f\n", lkc.z[j], phase + 1, area[lkc.z[j]].l, area[lkc.z[j]].m); // Write out numbers rather than IDs
      //***end
    };
    fprintf(fp9, "\n");
    
    maxdist = MaxDistance(&(lkc.z[0]), lkc.z_length, &maxdistz1, &maxdistz2);
    
    if (CARTESIAN == 0)
      fprintf(fp9, "  Maximum distance.......: %g km (areas: %s to %s)\n",
              maxdist, area[maxdistz1].id, area[maxdistz2].id
      );
    else if (CARTESIAN == 1)
      fprintf(fp9, "  Maximum distance.......: %g    (areas: %s to %s)\n",
              maxdist, area[maxdistz1].id, area[maxdistz2].id
      );
    
    fprintf(fp9, "  Number of cases .......: %d\n", lkc.nZ);
    if (MODEL == 0) {
      fprintf(fp9, "  (Expected number of cases: %g)\n", (double)lkc.mZ);
      fprintf(fp9, "  Overall relative risk .: %g\n", ((double)lkc.nZ / (double)lkc.mZ) / ((double)nG[0] / (double)mG));
    } else if (MODEL == 1) {
      fprintf(fp9, "  Population ............: %d\n", (int)lkc.mZ);
    }
    fprintf(fp9, "  Statistic value .......: %g\n", lkc.lambda);
    
    for (j = 1; j <= SIMCOUNT; ++j)
      if (lkc.lambda > maxstat[j])
        break;
      
    fprintf(fp9, "  Monte Carlo rank ......: %d/%d\n", (rnk = j), SIMCOUNT + 1);
    fprintf(fp9, "  P-value ...............: %g\n", (double)j / (double)(SIMCOUNT + 1));
    //***added by suzuryo,2011.9.28
    fprintf(fp13, "***%g\n", (double)j / (double)(SIMCOUNT + 1));
    //***end
    if (rnk == SIMCOUNT + 1) {
      Rprintf("*** There are no more secondary clusters ***\n");
      fprintf(fp9, "\n*** There are no more secondary clusters ***\n");
      NOMORE = 1;
    };
    
    fprintf(fp9, "\n");
    
    //***added by Otani, 2019
    if (phase == 0) {
      if (MODEL == 0)
        fprintf(fp14, "max_dist from to n_case expected RR stats rank pval");
      else
        fprintf(fp14, "max_dist from to n_case population stats rank pval");
      
      for (i = 0; i < N; ++i)
        fprintf(fp14, " %s", area[i].id);
      
      fprintf(fp14, "\n");
    }
    
    if (MODEL == 0) {
      fprintf(fp14, "%g %s %s %d %g %g %g %d %g",
              maxdist, area[maxdistz1].id, area[maxdistz2].id,
              lkc.nZ, (double)lkc.mZ, ((double)lkc.nZ / (double)lkc.mZ) / ((double)nG[0] / (double)mG),
              lkc.lambda, rnk, (double)j / (double)(SIMCOUNT + 1));
    } else {
      fprintf(fp14, "%g %s %s %d %d %g %d %g",
              maxdist, area[maxdistz1].id, area[maxdistz2].id,
              lkc.nZ, (int)lkc.mZ,
              lkc.lambda, rnk, (double)j / (double)(SIMCOUNT + 1));
    }
    
    int *idx = (int*)calloc(N, sizeof(int));
    for (j = 0; j < lkc.z_length; ++j) {
      idx[lkc.z[j]] = 1;
    }
    for (j = 0; j < N; ++j) {
      fprintf(fp14, " %d", idx[j]);
    }
    fprintf(fp14, "\n");
    free(idx);
    
    //***end
  }
  
  if (NOMORE == 1) {
    fprintf(fp9, "--------------------------------------------------------\n");
    fprintf(fp9, "The statistic value required for an observed\ncluster to be significant at level\n");
    fprintf(fp9, "... 0.01: %g\n", maxstat[(int)(0.01*(SIMCOUNT+1))]);
    fprintf(fp9, "... 0.05: %g\n", maxstat[(int)(0.05*(SIMCOUNT+1))]);
    fprintf(fp9, "\n");
    fprintf(fp9, "--------------------------------------------------------\n");
    fprintf(fp9, "\n");
    fprintf(fp9, "PARAMETER SETTINGS\n");
    fprintf(fp9, "\n");
    fprintf(fp9, "Input Files\n");
    fprintf(fp9, "-----------\n");
    fprintf(fp9, "  Case File:        %s\n", filename1);
    fprintf(fp9, "  Coordinates File: %s\n", filename3);
    fprintf(fp9, "  Matrix File:      %s\n", filename4);
    fprintf(fp9, "\n");
    fprintf(fp9, "Analysis\n");
    fprintf(fp9, "--------\n");
    fprintf(fp9, "  Type of Analysis     : Purely Spatial\n");
    fprintf(fp9, "  Probability Model    : %s\n", (MODEL == 1) ? "Binomial" : "Poisson");
    
    fprintf(fp9, "  <STATISTICS>\n");
    fprintf(fp9, "   %s.",
            (STATTYPE == 1) ? "Log likelihood ratio with restriction" : "Original log likelihood ratio");
    if (STATTYPE == 1)
      fprintf(fp9, " (%s<%f)\n\n", "Pr{X>x}+0.5*Pr{X=x}", RALPHA);
    else
      fprintf(fp9, "\n\n");
    
    fprintf(fp9, "  Scan for Area with   : High Rates\n");
    fprintf(fp9, "  Coordinates          : %s\n", (CARTESIAN == 0) ? "Latitude/Longitude" : "Cartesian");
    if (CARTESIAN == 0)
      fprintf(fp9, "  Radius of Earth      : %g km\n", R_EARTH);
    fprintf(fp9, "\n");
    fprintf(fp9, "  Number of Replications : %d\n", SIMCOUNT);
    fprintf(fp9, "  Type of Random number  : %s\n", (RANTYPE == 0) ? "Multinomial" : ((MODEL == 1) ? "Binomial" : "Poisson"));
    fprintf(fp9, "  Seed of Random number  : %d\n", DEFSEED);
    fprintf(fp9, "\n");
    fprintf(fp9, "Scanning Method\n");
    fprintf(fp9, "---------------\n");
    fprintf(fp9, "  Scanning Method :  %s scan\n", (SCANMETHOD == 1) ? "Circular" : "Flexible");
    fprintf(fp9, "  Maximum Spatial Cluster Size : %d areas\n", K);
    fprintf(fp9, "\n");
    fprintf(fp9, "Output\n");
    fprintf(fp9, "------\n");
    fprintf(fp9, "  Results File : %s\n", filename9);
    fprintf(fp9, "  Criteria for Reporting Secondary Clusters : No Geographical Overlap\n");
    fprintf(fp9, "\n");
    fprintf(fp9, "--------------------------------------------------------\n");
    fprintf(fp9, "\n");
    
    //				fprintf(fp9,"Working Memory Size : %d bytes\n",simworksize);
    time(&tm2);
    fprintf(fp9, "Program completed  : %s", ctime(&tm2));
    fprintf(fp9, "Total Running Time : %d seconds\n", (int)(tm2 - tm1));
    return;
  };
};
}

int     LoadData() {
  int     i, j, c;
  char    buf[256], ids[128];
  char	buf2[256];
  int		k;
  double	mG_temp;
  float	ppbin;
  
  /* N=number of areas */
  if ((fp1 = fopen(filename1, "rt")) == NULL) {
    Rprintf("ErrFile %s.", filename1);
    Rcpp::stop("ERROR! Code:", ErrFile1);
  };
  
  
  for (N = 0; fgets(buf, 256, fp1) != NULL; ++N) {
    if (sscanf(buf, "%s", buf2) != 1)
      N--;
  }
    
    if ((area = (TArea *)calloc(N, sizeof(TArea))) == NULL) {
      Rprintf("ErrMemory - area");
      Rcpp::stop("ERROR! Code:", ErrMemory);
    };
    
    if ((area_sorted = (TArea **)calloc(N, sizeof(TArea *))) == NULL) {
      Rprintf("ErrMemory - area_sorted");
      Rcpp::stop("ERROR! Code:", ErrMemory);
    };
    
    if ((detectedarea = (int *)calloc(N, sizeof(int))) == NULL) {
      Rprintf("ErrMemory - detectedarea");
      Rcpp::stop("ERROR! Code:", ErrMemory);
    };
    for (i = 0; i < N; ++i)
      detectedarea[i] = 0;
    
    /* Read coordinates */
    if ((fp3 = fopen(filename3, "rt")) == NULL) {
      Rprintf("ErrFile %s.", filename3);
      Rcpp::stop("ERROR! Code:", ErrFile3);
    };
    for (i = 0; i < N; ++i) {
      if (fgets(buf, 256, fp3) == NULL) {
        Rprintf("Read Error - %s", filename3);
        Rcpp::stop("ERROR! Code:", ErrFile3);
      };
      
      if (sscanf(buf, "%s", buf2) != 1) {
        i--;
        continue;
      };
      
      area[i].index = i;
      
      if (sscanf(buf, "%s %lf %lf", ids, &(area[i].l), &(area[i].m)) != 3) {
        Rprintf("Read Error -- %s", filename3);
        Rcpp::stop("ERROR! Code:", ErrFile3);
      };
      area[i].id = strdup(ids);
    };
    fclose(fp3);
    
    if (K > N)
      K = N;
    
    /* Read cases */
    rewind(fp1);
    for (i = 0; i < N; ++i) {
      if (fgets(buf, 256, fp1) == NULL)
        Rcpp::stop("ERROR! Code:", ErrFile1);
      
      if (sscanf(buf, "%s", buf2) != 1) {
        i--;
        continue;
      };
      
      area[i].index = i;
      
      if (sscanf(buf, "%s %d %lf", ids, &(area[i].cases), &(area[i].popul)) != 3) {
        Rprintf("Read Error -- case file");
        Rcpp::stop("ERROR! Code:", ErrFile1);
      };
      if (strcmp(ids, area[i].id) != 0) {
        Rprintf("ErrFile1ID");
        Rcpp::stop("ERROR! Code:", ErrFile1ID);
      };
    };
    fclose(fp1);
    
    /* Allocate and read connection matrix a[0..N-1][0..N-1] */
    if ((fp4 = fopen(filename4, "rt")) == NULL) {
      Rprintf("ErrMemory -- fp4");
      Rcpp::stop("ERROR! Code:", ErrFile4);
    };
    if ((a = (int **)calloc(N, sizeof(int *))) == NULL) {
      Rprintf("ErrMemory -- a");
      Rcpp::stop("ERROR! Code:", ErrMemory);
    };
    for (i = 0; i < N; ++i) {
      if ((a[i] = (int *)calloc(N, sizeof(int))) == NULL) {
        Rprintf("ErrMemory -- a[i]");
        Rcpp::stop("ERROR! Code:", ErrMemory);
      };
      if (fscanf(fp4, "%s", ids) != 1) {
        Rprintf("ErrFile4");
        Rcpp::stop("ERROR! Code:", ErrFile4);
      };
      if (strcmp(ids, area[i].id) != 0) {
        Rprintf("ErrFile4ID");
        Rcpp::stop("ERROR! Code:", ErrFile4ID);
      };
      for (j = 0; j < N; ++j) {
        if (fscanf(fp4, "%d", &c) != 1)
          Rcpp::stop("ERROR! Code:", ErrFile4);
        a[i][j] = c;
      };
    };
    
    /* check symmetry */
    for (i = 0; i < N; ++i) {
      for (j = 0; j < N; ++j) {
        if (a[i][j] != a[j][i])
          Rcpp::stop("ERROR! Code:", ErrFile4);
      }
    }
        fclose(fp4);
        
        if ((calen = (int *)calloc(N, sizeof(int))) == NULL) {
          Rprintf("ErrMemory -- calen");
          Rcpp::stop("ERROR! Code:", ErrMemory);
        };
        if ((ca = (areaidx **)calloc(N, sizeof(areaidx *))) == NULL) {
          Rprintf("ErrMemory -- ca");
          Rcpp::stop("ERROR! Code:", ErrMemory);
        };
        for (i = 0; i < N; ++i) {
          if ((ca[i] = (areaidx *)calloc(N, sizeof(areaidx))) == NULL) {
            Rprintf("ErrMemory -- ca[i]");
            Rcpp::stop("ERROR! Code:", ErrMemory);
          }
        }
        if ((caz = (areaidx *)calloc(N * K, sizeof(areaidx))) == NULL) {
          Rprintf("ErrMemory -- caz");
          Rcpp::stop("ERROR! Code:", ErrMemory);
        };
        if ((masksw = (short *)calloc(N, sizeof(short))) == NULL) {
          Rprintf("ErrMemory -- masksw");
          Rcpp::stop("ERROR! Code:", ErrMemory);
        };
        for (i = 0; i < N; ++i) {
          calen[i] = 0;
          for (j = 0; j < N; ++j)
            if (a[i][j] && i != j)
              ca[i][calen[i]++] = j;
        };
        
        /* Allocate and initialize data for Monte Carlo simulation */
        
        /* popul[0..N-1] */
        if ((popul = (double *)calloc(N, sizeof(double))) == NULL) {
          Rprintf("ErrMemory - popul");
          Rcpp::stop("ERROR! Code:", ErrMemory);
        };
        /* cases[0..N-1][0..SIMCOUNT] */
        if ((cases = (int **)calloc(N, sizeof(int *))) == NULL) {
          Rprintf("ErrMemory - cases");
          Rcpp::stop("ERROR! Code:", ErrMemory);
        };
        for (i = 0; i < N; ++i) {
          if ((cases[i] = (int *)calloc(SIMCOUNT + 1, sizeof(int))) == NULL) {
            Rprintf("ErrMemory -- cases");
            Rcpp::stop("ERROR! Code:", ErrMemory);
          }
        }
          /* pp[0..N-1] */
          if ((pp = (double *)calloc(N, sizeof(double))) == NULL) {
            Rprintf("ErrMemory - pp");
            Rcpp::stop("ERROR! Code:", ErrMemory);
          };
          if ((rtmp = (int *)calloc(N, sizeof(int))) == NULL) {
            Rprintf("ErrMemory - rtemp");
            Rcpp::stop("ERROR! Code:", ErrMemory);
          };
          
          /* pv0[i][s] */
          if ((pv0 = (double **)calloc(N, sizeof(double *))) == NULL) {
            Rprintf("ErrMemory - pv0");
            Rcpp::stop("ERROR! Code:", ErrMemory);
          };
          for (i = 0; i < N; ++i)
            if ((pv0[i] = (double *)calloc(SIMCOUNT + 1, sizeof(double))) == NULL) {
              Rprintf("ErrMemory -- pv0");
              Rcpp::stop("ERROR! Code:", ErrMemory);
            };
            
            /* nG[0..SIMCOUNT] */
            if ((nG = (int *)calloc(SIMCOUNT + 1, sizeof(int))) == NULL) {
              Rprintf("ErrMemory - nG");
              Rcpp::stop("ERROR! Code:", ErrMemory);
            };
            
            /* Calculate total population and cases */
            mG_temp = nG[0] = 0;
            misarea = 0;
            for (i = 0; i < N; ++i) {
              nG[0] += (cases[i][0] = area[i].cases);
              if (area[i].popul <= 0.0) {
                Rprintf("- %s is excluded.\n", area[i].id);
                detectedarea[i] = -1;
                ++misarea;
                //				Rcpp::stop("ERROR! Code:", ErrData);
              } else
                mG_temp += (area[i].popul);
            };
            
            if (MODEL == 0) {
              if (abs((double)nG[0] - mG_temp) < 0.001)
                EXTYPE = 0;
            } else {
                EXTYPE = 1;
            }
              
              if (MODEL == 1) { /* Binomial model */
            for (i = 0; i < N; ++i)
              popul[i] = area[i].popul;
                mG = mG_temp;
              } else { /* Poisson model */
            for (i = 0; i < N; ++i)
              popul[i] = (nG[0] * area[i].popul) / ((double)mG_temp);
                mG = (double)nG[0];
              };
              
              /* Generate random numbers */
              
              Rprintf("\nsimulation data generating...");
              
              //ignlgi();
              //advnst(DEFSEED);
              
              nGmax = nG[0];
              
              if (MODEL == 1) { /* Binomial model */
                if (RANTYPE == 0) { /* Multinomial */
                  for (i = 0; i < N; ++i) {
                    if (popul[i] > 0.0)
                      pp[i] = (double)popul[i] / (double)mG;
                    else
                      pp[i] = 0.0;
                  };
                  
                  for (j = 1; j <= SIMCOUNT; ++j) {
                    k = 0;
                    //genmul(nG[0], pp, N, rtmp);
                    rmultinom(nG[0], pp, N, rtmp);
                    for (i = 0; i < N; ++i) {
                      if (rtmp[i] > popul[i]) {
                        ++k;
                        cases[i][j] = popul[i];
                      } else
                        cases[i][j] = rtmp[i];
                    }
                    nG[j] = nG[0] - k;
                  };
                } else {	/* Binomial */
                ppbin = nG[0] / (double)mG;
                  
                  for (j = 1; j <= SIMCOUNT; ++j) {
                    nG[j] = 0;
                    for (i = 0; i < N; ++i) {
                      //cases[i][j] = ignbin((long)popul[i], ppbin);
                      cases[i][j] = R::rbinom((double)popul[i], (double)ppbin);
                      nG[j] += cases[i][j];
                    };
                    if (nG[j] > nGmax)
                      nGmax = nG[j];
                  };
                }
              } else {	/* Poisson model */
              if (RANTYPE == 0) { /* Multinomial */
              for (i = 0; i < N; ++i) {
                if (popul[i] > 0.0)
                  pp[i] = (double)popul[i] / (double)mG;
                else
                  pp[i] = 0.0;
              };
                
                for (j = 1; j <= SIMCOUNT; ++j) {
                  //genmul(nG[0], pp, N, rtmp);
                  rmultinom(nG[0], pp, N, rtmp);
                  for (i = 0; i < N; ++i)
                    cases[i][j] = rtmp[i];
                  nG[j] = nG[0];
                };
              } else if (RANTYPE == 1) { /* Poisson */
              for (j = 1; j <= SIMCOUNT; ++j) {
                nG[j] = 0;
                for (i = 0; i < N; ++i) {
                  if (popul[i] > 0.0)
                    cases[i][j] = R::rpois(popul[i]); //cases[i][j] = ignpoi((float)popul[i]);
                  else
                    cases[i][j] = 0.0;
                  nG[j] += cases[i][j];
                };
                if (nG[j] > nGmax)
                  nGmax = nG[j];
              };
              };
              };
              Rprintf("done\n");
              
              /* maxstat[0..SIMCOUNT] */
              if ((maxstat = (double *)calloc(SIMCOUNT + 1, sizeof(double))) == NULL) {
                Rprintf("ErrMemory - maxstat");
                Rcpp::stop("ERROR! Code:", ErrMemory);
              };
              
              if (MODEL == 0) {
                if ((SIMCOUNT * nGmax) < 100000000)
                  lors = 0;
                else
                  lors = 1;
                
                if (lors == 0) {
                  /* minmZ_zlength[nGmax] */
                  if ((minmZ_zlength = (int *)calloc(nGmax + 1, sizeof(int))) == NULL) {
                    Rprintf("ErrMemory - minmZ_zlength");
                    Rcpp::stop("ERROR! Code:", ErrMemory);
                  };
                  
                  /* minmZ_z[0..nGmax][K-1] */
                  if ((minmZ_z = (areaidx **)calloc(nGmax + 1, sizeof(areaidx *))) == NULL) {
                    Rprintf("ErrMemory -- minmZ_z");
                    Rcpp::stop("ERROR! Code:", ErrMemory);
                  };
                  for (i = 0; i <= nGmax; ++i)
                    if ((minmZ_z[i] = (areaidx *)calloc(K, sizeof(areaidx))) == NULL) {
                      Rprintf("ErrMemory -- minmZ_z[i]");
                      Rcpp::stop("ERROR! Code:", ErrMemory);
                    };
                    
                    /* *minmZ[0..SIMCOUNT] */
                    if ((minmZ = (double **)calloc(SIMCOUNT + 1, sizeof(double *))) == NULL) {
                      Rprintf("ErrMemory -- minmZ");
                      Rcpp::stop("ERROR! Code:", ErrMemory);
                    };
                    /* minmZ[0..SIMCOUNT][0..snG] */
                    for (i = 0; i <= SIMCOUNT; ++i) {
                      if ((minmZ[i] = (double *)calloc(nGmax + 1, sizeof(double))) == NULL) {
                        Rprintf("ErrMemory -- minmZ[i]");
                        Rcpp::stop("ERROR! Code:", ErrMemory);
                      };
                    };
                } else {
                  if ((MLC_z = (areaidx *)calloc(K, sizeof(areaidx))) == NULL) {
                    Rprintf("ErrMemory - MLC_z");
                    Rcpp::stop("ERROR! Code:", ErrMemory);
                  };
                  
                  if ((Lpoi0 = (double *)calloc(SIMCOUNT + 1, sizeof(double))) == NULL) {
                    Rprintf("ErrMemory - Lbin0");
                    Rcpp::stop("ERROR! Code:", ErrMemory);
                  };
                  
                  for (j = 0; j <= SIMCOUNT; ++j)
                    Lpoi0[j] = nG[j] * log(nG[j] / (double)mG);
                };
              } else if (MODEL == 1) {
                /* MLC_z[K-1] */
                if ((MLC_z = (areaidx *)calloc(K, sizeof(areaidx))) == NULL) {
                  Rprintf("ErrMemory - MLC_z");
                  Rcpp::stop("ERROR! Code:", ErrMemory);
                };
                
                /* for Binomila Model */
                if ((Lbin0 = (double *)calloc(SIMCOUNT + 1, sizeof(double))) == NULL) {
                  Rprintf("ErrMemory - Lbin0");
                  Rcpp::stop("ERROR! Code:", ErrMemory);
                };
                
                for (j = 0; j <= SIMCOUNT; ++j)
                  Lbin0[j] = nG[j] * log(nG[j] / (double)mG) + (mG - nG[j]) * log((mG - nG[j]) / (double)(mG));
              };
              
              return(0);
}

/* --------------------------------------------------------------*/
char    *ffgets(char *buf, int len, FILE *fp) {
  char    *status = fgets(buf, len, fp);
  int     i;
  
  if (buf == NULL)
    return(NULL);
  i = strlen(buf);
  if (buf[i-1] == '\n')
    buf[i-1] = '\0';
  return(status);
}

/* -----------------------------------------------------------------*/
int     flexcore(int argc, char *argv[]) {
  char    buf[256];
  int     i, s;
  int		SIM2;
  double	totp;
  
  tm1 = time(NULL);
  
  // parameter file
  if (argc == 2 && argv[1][0] == '/') {
    if ((fp0 = fopen(&argv[1][1], "rt")) == NULL) {
      Rprintf("Cannot open the parameter file.");
      return(ErrFile0);
    }
    
    while (ffgets(buf, 256, fp0) != NULL) {
      if (strncasecmp(buf, "CASE=", 5) == 0)
        filename1 = strdup(&buf[5]);
      else if (strncasecmp(buf, "POPULATION=", 11) == 0)
        filename2 = strdup(&buf[11]);
      else if (strncasecmp(buf, "COORDINATES=", 12) == 0)
        filename3 = strdup(&buf[12]);
      else if (strncasecmp(buf, "MATRIX=", 7) == 0)
        filename4 = strdup(&buf[7]);
      else if (strncasecmp(buf, "RESULTS=", 8) == 0) {
        filename9 = strdup(&buf[8]);
        strcat(&buf[8], ".clt");
        CLSTGEO = strdup(&buf[8]);
      } else if (strncasecmp(buf, "LAMBDAFILE=", 11) == 0)
        LAMBDAFILE = strdup(&buf[11]);
      else if (strncasecmp(buf, "SIMCOUNT=", 9) == 0) {
        SIMCOUNT = atoi(&buf[9]);
        if (SIMCOUNT >= SIMLIM) {
          Rprintf("SIMCOUNT is too large!");
          Rcpp::stop("ERROR! Code:", ErrMemory);
        }
      } else if (strncasecmp(buf, "RANSEED=", 8) == 0)
        DEFSEED = atoi(&buf[8]);
      else if (strncasecmp(buf, "CLUSTERSIZE=", 12) == 0)
        K = atoi(&buf[12]);
      else if (strncasecmp(buf, "CARTESIAN=", 10) == 0)
        CARTESIAN = (atoi(&buf[10]) == 0) ? 0 : 1;
      else if (strncasecmp(buf, "RADIUS=", 7) == 0)
        R_EARTH = atof(&buf[7]);
      else if (strncasecmp(buf, "COMMENT=", 8) == 0)
        strcpy(comment, &buf[8]);
      else if (strncasecmp(buf, "MODEL=", 6) == 0)
        MODEL = (strcasecmp(&buf[6], "BINOMIAL") == 0) ? 1 : 0;
      else if (strncasecmp(buf, "RANTYPE=", 8) == 0)
        RANTYPE = (strcasecmp(&buf[8], "MULTINOMIAL") == 0) ? 0 : 1;
      else if (strncasecmp(buf, "SCANMETHOD=", 11) == 0)
        SCANMETHOD = (strcasecmp(&buf[11], "FLEXIBLE") == 0) ? 0 : 1;
      else if (strncasecmp(buf, "STATTYPE=", 9) == 0) {
        STATTYPE = atoi(&buf[9]);
        if (STATTYPE != 1)
          STATTYPE = 0;
      } else if (strncasecmp(buf, "RALPHA=", 7) == 0)
        RALPHA = atof(&buf[7]);
      //added by suzuryo, 2011/9/28.
      else if(strncasecmp(buf,"NODEFILE=",9)==0){
        NODEFILE=strdup(&buf[9]); // fprintf(stderr,"NODEFILE: %s\n",NODEFILE);
      }
      //end
      //added by suzuryo, 2011/10/6.
      else if(strncasecmp(buf,"EDGEFILE=",9)==0){
        EDGEFILE=strdup(&buf[9]); // fprintf(stderr,"EDGEFILE: %s\n",EDGEFILE);
      }
      //end
      //added by Otani, 2019
      else if (strncasecmp(buf, "RFILE=",6)==0){
        RFILE=strdup(&buf[6]); // fprintf(stderr,"RFILE: %s\n",RFILE);
      }
    };
    fclose(fp0);
  } else {
    Rprintf("---------------------------\n");
    Rprintf(" FleXScan ver.3.1.2\n");
    Rprintf("---------------------------\n");
    Rprintf("FLSCMD /settingfile\n");
    return(0);
  };
  
  Rprintf("<STATISTICAL MODEL>\n");
  Rprintf(" %s.\n", (MODEL == 1) ? "Binomial" : "Poisson");
  Rprintf("<SCANING METHOD>\n");
  Rprintf(" %s spatial scan by length.\n", (SCANMETHOD == 1) ? "circular" : "flexible");
  Rprintf("<STATISTICS>\n");
  Rprintf(" %s.",
          (STATTYPE == 1) ? "Log likelihood with restriction" : "Original log likelihood ratio");
  if (STATTYPE == 1)
    Rprintf(" (%s<%f)\n", "Pr{X>x}+0.5*Pr{X=x}", RALPHA);
  else
    Rprintf("\n");
  Rprintf("<INPUT FILES>\n");
  Rprintf(" Case file [%s].\n", filename1);
  Rprintf(" Coodinates file [%s].\n", filename3);
  Rprintf(" Matrix file [%s].\n", filename4);
  Rprintf("<OUTPUT FILE>\n");
  Rprintf(" Results file [%s].\n", filename9);
  Rprintf("<SETTINGS>\n");
  Rprintf(" Maximum area length = %d.\n", K);
  Rprintf(" Number of simulation = %d.\n", SIMCOUNT);
  Rprintf(" Seed of random number = %d.\n", DEFSEED);
  Rprintf(" Random number = %s.\n", (RANTYPE == 0) ? "Multinomial" : ((MODEL == 0) ? "Poisson" : "Binomial"));
  Rprintf(" Coordinates = %s.\n", (CARTESIAN == 0) ? "Latitude/Longitude" : "Cartesian");
  if (CARTESIAN == 0)
    Rprintf(" Radius of Earth = %g km.\n", R_EARTH);
  
  Rprintf("\nInitializing...\n");
  
  if (LoadData() != 0)
    return(-1);
  if ((w = (areaidx *)calloc(N, sizeof(areaidx))) == NULL)
    return(-1);
  
  Rprintf("\n--  CALCULATING  --\n");
  //Rprintf("<ctrl-C> to terminate\n");
  if ((fp9 = fopen(filename9, "wt")) == NULL)
    return(ErrFile9);
  if ((fp11 = fopen(CLSTGEO, "wt")) == NULL)
    return(ErrFile11);
  if ((fp12 = fopen(EDGEFILE, "wt")) == NULL)
    return(ErrFile12);
  
  SIM2 = SIMCOUNT;
  
  if (MODEL == 0) {
    if (STATTYPE == 1) {
      for (i = 0; i < N; ++i) {
        if (detectedarea[i] != -1)
          for (s = 0; s <= SIM2; ++s) {
            pv0[i][s] = Ppfm(cases[i][s], popul[i]);
          };
      };
    };
  } else if (MODEL == 1) {
    if (STATTYPE == 1) {
      totp = nG[0] / double(mG);
      for (i = 0; i < N; ++i) {
        if (detectedarea[i] != -1)
          for (s = 0; s <= SIM2; ++s) {
            pv0[i][s] = Pbfm(cases[i][s], (int)popul[i], totp);
          };
      };
    };
  };
  
  FlexScan();
  
  if (lors == 0)
    free(minmZ);
  
  Rprintf("-----\n");
  Rprintf("Completed in %d seconds.\n", tm2 - tm1);
  
  fclose(fp9);
  fclose(fp11);
  fclose(fp12);
  fclose(fp13);
  fclose(fp14);
  
  return(0);
}

// [[Rcpp::export]]
NumericVector runFleXScan(String settingfile) {
  int argc = 2;
  char *argv[2];
  
  argv[0] = (char*)malloc(256 * sizeof(char));
  strcpy(argv[0], "flscan312");
  argv[1] = (char*)malloc(256 * sizeof(char));
  settingfile.push_front("/");
  strcpy(argv[1], settingfile.get_cstring());
  
  int retval = flexcore(argc, argv);
  
  free(argv[0]);
  free(argv[1]);
  
  return retval;
}

