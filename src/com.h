#ifndef __COM_H_
#define __COM_H_

void advnst(long k);
void getsd(long *iseed1,long *iseed2);
long ignlgi(void);
void initgn(long isdtyp);
void inrgcm(void);
void setall(long iseed1,long iseed2);
void setant(long qvalue);
void setsd(long iseed1,long iseed2);

#endif // __COM_H_
