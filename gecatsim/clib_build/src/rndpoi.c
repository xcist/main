// Copyright 2020, General Electric Company. All rights reserved. See https://github.com/xcist/code/blob/master/LICENSE

#ifdef WIN32
#define DLLEXPORT __declspec(dllexport)
#else
#define DLLEXPORT
//#define LIMIT_LOAD
#endif

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>


void gsrgs(int getset,int *qvalue);
void gssst(int getset,int *qset);
void gscgn(int getset,int *g);
void inrgcm(void);

int Xm1,Xm2,Xa1,Xa2,Xcg1[32],Xcg2[32],Xa1w,Xa2w,Xig1[32],Xig2[32],Xlg1[32],
    Xlg2[32],Xa1vw,Xa2vw;
int Xqanti[32];

/*
 * int multMod(int a, int s, int m)
 * Returns (a*s)%m
*/
int multMod(int a, int s, int m)
{
    static const int h=32768;
    static int a0,a1,k,p,q,qh,rh;

    if(a <= 0 || a >= m || s <= 0 || s >= m)
    {
        fprintf(stderr, "Error: multMod requires: 0 < a < m; 0 < s < m\n");
        exit(1);
    }
    if (a < h)
    {
        a0 = a;
        p = 0;
    }
    else
    {
        a1 = a/h;
        a0 = a-h*a1;
        qh = m/h;
        rh = m-h*qh;
        if (a1 >= h)
        {
            a1 -= h;
            k = s/qh;
            p = h*(s-k*qh)-k*rh;
            while (p < 0)
            {
                p += m;
            }
        }
        else
        {
            p = 0;
        }
        if (a1 != 0)
        {
            q = m/a1;
            k = s/q;
            p -= (k*(m-a1*q));
            if(p > 0)
            {
                p -= m;
            }
            p += (a1*(s-k*q));
            while (p < 0)
            {
                p += m;
            }
            k = p/qh;
            p = h*(p-k*qh)-k*rh;
            while (p < 0)
            {
                p += m;
            }
        }
    }
    if (a0 != 0)
    {
        q = m/a0;
        k = s/q;
        p -= (k*(m-a0*q));
        if(p > 0)
        {
            p -= m;
        }
        p += (a0*(s-k*q));
        while (p < 0)
        {
            p += m;
        }
    }
    return p;
}

void gsrgs(int bSet,int *qvalue)
/*
**********************************************************************
     void gsrgs(int bSet,int *qvalue)
               Get/Set Random Generators Set
     Gets or sets whether random generators set (initialized).
     Initially (data statement) state is not set
     If bSet is 1 state is set to qvalue
     If bSet is 0 state returned in qvalue
**********************************************************************
*/
{
static int qinit = 0;

    if(bSet == 0)
    {
        *qvalue = qinit;
    }
    else
    {
        qinit = *qvalue;
    }
}

void gssst(int bSet,int *qset)
/*
**********************************************************************
     void gssst(int bSet,int *qset)
          Get or Set whether Seed is Set
     Initialize to Seed not Set
     If bSet is 1 sets state to Seed Set
     If bSet is 0 returns T in qset if Seed Set
     Else returns F in qset
**********************************************************************
*/
{
static int qstate = 0;
    if(bSet)
    {
        qstate = 1;
    }
    else
    {
        *qset = qstate;
    }
}

void gscgn(int bSet,int *g)
/*
**********************************************************************
     void gscgn(int bSet,int *g)
                         Get/Set GeNerator
     Gets or returns in G the number of the current generator
                              Arguments
     bSet --> 0 Get
              1 Set
     g <-- Number of the current random number generator (1..32)
**********************************************************************
*/
{
    static const int numg=32;
    static int curntg = 1;
    if(bSet == 0) *g = curntg;
    else
    {
        if(*g < 0 || *g > numg)
        {
            fprintf(stderr, "Error: Generator number out of range in GSCGN\n");
            exit(1);
        }
        curntg = *g;
    }
}

void inrgcm(void)
/*
**********************************************************************
     void inrgcm(void)
          INitialize Random number Generator CoMmon
                              Function
     Initializes common area  for random number  generator.  This saves
     the  nuisance  of  a  BLOCK DATA  routine  and the  difficulty  of
     assuring that the routine is loaded with the other routines.
**********************************************************************
*/
{
    static const int numg=32;
    static int T1;
    static int i;
/*
     V=20;                            W=30;
     A1W = MOD(A1**(2**W),M1)         A2W = MOD(A2**(2**W),M2)
     A1VW = MOD(A1**(2**(V+W)),M1)    A2VW = MOD(A2**(2**(V+W)),M2)
   If V or W is changed A1W, A2W, A1VW, and A2VW need to be recomputed.
    An efficient way to precompute a**(2*j) MOD m is to start with
    a and square it j times modulo m using the function MLTMOD.
*/
    Xm1 = 2147483563;
    Xm2 = 2147483399;
    Xa1 = 40014;
    Xa2 = 40692;
    Xa1w = 1033780774;
    Xa2w = 1494757890;
    Xa1vw = 2082007225;
    Xa2vw = 784306273;
    for(i=0; i<numg; i++)
    {
        Xqanti[i] = 0;
    }
    T1 = 1;
    gsrgs(1,&T1);
}

/*
**********************************************************************
     void initgn(int isdtyp)
          INIT-ialize current G-e-N-erator
     Reinitializes the state of the current generator
                              Arguments
     isdtyp -> The state to which the generator is to be set
          isdtyp = -1  => sets the seeds to their initial value
          isdtyp =  0  => sets the seeds to the first value of
                          the current block
          isdtyp =  1  => sets the seeds to the first value of
                          the next block
**********************************************************************
*/
void initGn(int isdtyp)
{
    static int g;
    static int qrgnin;

    gsrgs(0L,&qrgnin);
    if(!qrgnin)
    {
        fprintf(stderr,"%s\n",
            "INITGN called before random number generator  initialized -- abort!");
        exit(1);
    }
    gscgn(0,&g);
    switch (isdtyp)
    {
        case -1:
            Xlg1[g-1] = Xig1[g-1];
            Xlg2[g-1] = Xig2[g-1];
            break;
        case 0:
            break;
        case 1:
            Xlg1[g-1] = multMod(Xa1w,Xlg1[g-1],Xm1);
            Xlg2[g-1] = multMod(Xa2w,Xlg2[g-1],Xm2);
            break;
        default:
            fprintf(stderr,"%s\n","isdtyp not in range in INITGN");
            exit(1);
    }
    Xcg1[g-1] = Xlg1[g-1];
    Xcg2[g-1] = Xlg2[g-1];
}


void setall(int iseed1,int iseed2)
/*
**********************************************************************
     void setall(int iseed1,int iseed2)
               SET ALL random number generators
     Sets the initial seed of generator 1 to ISEED1 and ISEED2. The
     initial seeds of the other generators are set accordingly, and
     all generators states are set to these seeds.
     This is a transcription from Pascal to Fortran of routine
     Set_Initial_Seed from the paper
     L'Ecuyer, P. and Cote, S. "Implementing a Random Number Package
     with Splitting Facilities." ACM Transactions on Mathematical
     Software, 17:98-111 (1991)
                              Arguments
     iseed1 -> First of two integer seeds
     iseed2 -> Second of two integer seeds
**********************************************************************
*/
{
    static const int numg=32;
    static int T1;
    static int g,ocgn;
    static int qrgnin;
    T1 = 1;
/*
     TELL IGNLGI, THE ACTUAL NUMBER GENERATOR, THAT THIS ROUTINE
      HAS BEEN CALLED.
*/
    gssst(1,&T1);
    gscgn(0L,&ocgn);
/*
     Initialize Common Block if Necessary
*/
    gsrgs(0L,&qrgnin);
    if(!qrgnin) inrgcm();
    Xig1[0] = iseed1;
    Xig2[0] = iseed2;
    initGn(-1L);
    for(g=2; g<=numg; g++) {
        Xig1[g-1] = multMod(Xa1vw,Xig1[g-2],Xm1);
        Xig2[g-1] = multMod(Xa2vw,Xig2[g-2],Xm2);
        gscgn(1,&g);
        initGn(-1);
    }
    gscgn(1,&ocgn);
}

int ignlgi(void)
/*
**********************************************************************
     int ignlgi(void)
               GeNerate LarGe Integer
     Returns a random integer following a uniform distribution over
     (1, 2147483562) using the current generator.
     This is a transcription from Pascal to Fortran of routine
     Random from the paper
     L'Ecuyer, P. and Cote, S. "Implementing a Random Number Package
     with Splitting Facilities." ACM Transactions on Mathematical
     Software, 17:98-111 (1991)
**********************************************************************
*/
{
    static const int numg=32;
    static int ignlgi,curntg,k,s1,s2,z;
    static int qqssd,qrgnin;
/*
     IF THE RANDOM NUMBER PACKAGE HAS NOT BEEN INITIALIZED YET, DO SO.
     IT CAN BE INITIALIZED IN ONE OF TWO WAYS : 1) THE FIRST CALL TO
     THIS ROUTINE  2) A CALL TO SETALL.
*/
    gsrgs(0L,&qrgnin);
    if(!qrgnin)
    {
        inrgcm();
    }
    gssst(0,&qqssd);
    if(!qqssd)
    {
        //setall(1234567890,123456789);
/*
		To make seeds different for each runs, use time as the seeds. Mingye
*/
		time_t seed1=time(0); //current time
		time_t seed2=(int)(seed1/10);
		setall(seed1,seed2);
    }
/*
     Get Current Generator
*/
    gscgn(0L,&curntg);
    s1 = Xcg1[curntg-1];
    s2 = Xcg2[curntg-1];
    k = s1/53668;
    s1 = Xa1*(s1-k*53668)-k*12211;
    if(s1 < 0)
    {
        s1 += Xm1;
    }
    k = s2/52774;
    s2 = Xa2*(s2-k*52774)-k*3791;
    if(s2 < 0)
    {
        s2 += Xm2;
    }
    Xcg1[curntg-1] = s1;
    Xcg2[curntg-1] = s2;
    z = s1-s2;
    if(z < 1) z += (Xm1-1);
    if(Xqanti[curntg-1])
    {
        z = Xm1-z;
    }
    ignlgi = z;
    return ignlgi;
}

float fsign( float num, float sign )
/* Transfers sign of argument sign to argument num */
{
if ( ( sign>0.0f && num<0.0f ) || ( sign<0.0f && num>0.0f ) )
    return -num;
else return num;
}

float ranf(void)
/*
**********************************************************************
     float ranf(void)
                RANDom number generator as a Function
     Returns a random floating point number from a uniform distribution
     over 0 - 1 (endpoints of this interval are not returned) using the
     current generator
     This is a transcription from Pascal to Fortran of routine
     Uniform_01 from the paper
     L'Ecuyer, P. and Cote, S. "Implementing a Random Number Package
     with Splitting Facilities." ACM Transactions on Mathematical
     Software, 17:98-111 (1991)
**********************************************************************
*/
{
static float ranf;
/*
     4.656613057E-10 is 1/M1  M1 is set in a data statement in IGNLGI
      and is currently 2147483563. If M1 changes, change this also.
*/
    ranf = ignlgi()*4.656613057E-10;
    return ranf;
}


/*
 * expDistr()
 * Exponential distribution
*/
float expDistr()
{
    static float q[8] =
    {
        0.6931472,0.9333737,0.9888778,0.9984959,0.9998293,0.9999833,0.9999986,1.0
    };
    int i=0;
    float a,u,ustar,umin;

    a=-q[0];
    u=ranf();
    while (u <= 1.0)
    {
        a += q[0];
        u += u;
    }
    u -= 1.0;
    if (u <= q[0])
    {
        return a+u;
    }
    umin=ranf();
    do
    {
        ustar=ranf();
        if (ustar < umin)
        {
            umin = ustar;
        }
        i++;
    } while (u > q[i]);
    return a+umin*q[0];
}

/*
 * normDistr()
 * Normal distribution
*/
float normDistr()
{
    static float a[32] = {
        0.0,3.917609E-2,7.841241E-2,0.11777,0.1573107,0.1970991,0.2372021,0.2776904,
        0.3186394,0.36013,0.4022501,0.4450965,0.4887764,0.5334097,0.5791322,
        0.626099,0.6744898,0.7245144,0.7764218,0.8305109,0.8871466,0.9467818,
        1.00999,1.077516,1.150349,1.229859,1.318011,1.417797,1.534121,1.67594,
        1.862732,2.153875
    };
    static float d[31] = {
        0.0,0.0,0.0,0.0,0.0,0.2636843,0.2425085,0.2255674,0.2116342,0.1999243,
        0.1899108,0.1812252,0.1736014,0.1668419,0.1607967,0.1553497,0.1504094,
        0.1459026,0.14177,0.1379632,0.1344418,0.1311722,0.128126,0.1252791,
        0.1226109,0.1201036,0.1177417,0.1155119,0.1134023,0.1114027,0.1095039
    };
    static float t[31] = {
        7.673828E-4,2.30687E-3,3.860618E-3,5.438454E-3,7.0507E-3,8.708396E-3,
        1.042357E-2,1.220953E-2,1.408125E-2,1.605579E-2,1.81529E-2,2.039573E-2,
        2.281177E-2,2.543407E-2,2.830296E-2,3.146822E-2,3.499233E-2,3.895483E-2,
        4.345878E-2,4.864035E-2,5.468334E-2,6.184222E-2,7.047983E-2,8.113195E-2,
        9.462444E-2,0.1123001,0.136498,0.1716886,0.2276241,0.330498,0.5847031
    };
    static float h[31] = {
        3.920617E-2,3.932705E-2,3.951E-2,3.975703E-2,4.007093E-2,4.045533E-2,
        4.091481E-2,4.145507E-2,4.208311E-2,4.280748E-2,4.363863E-2,4.458932E-2,
        4.567523E-2,4.691571E-2,4.833487E-2,4.996298E-2,5.183859E-2,5.401138E-2,
        5.654656E-2,5.95313E-2,6.308489E-2,6.737503E-2,7.264544E-2,7.926471E-2,
        8.781922E-2,9.930398E-2,0.11556,0.1404344,0.1836142,0.2790016,0.7010474
    };
    static int i;
    static float snorm,u,s,ustar,aa,w,y,tt;

    u = ranf();
    s = 0.0;
    if(u > 0.5) s = 1.0;
    u += (u-s);
    u = 32.0*u;
    i = (int) (u);
    if(i == 32)
    {
        i = 31;
    }

    if(i == 0)
    {
        i = 6;
        aa = *(a+31);
        u += u;
        while (u < 1.0)
        {
           aa += d[i-1]; 
           i++;
           u += u;
        }
        u -= 1.0;
        while (1)
        {
            w = u*d[i-1];
            tt = (0.5*w+aa)*w;
            while (1)
            {
                ustar = ranf();
                if (ustar > tt)
                {
                    return (s == 1.0 ? -(aa+w) : aa+w);
                }
                u = ranf();
                if (ustar < u)
                {
                    break;
                }
                tt = u;
            }
            u = ranf();
        }
    }
    else
    {
        ustar = u-(float)i;
        aa = a[i-1];
        while (ustar <= t[i-1])
        {
            u = ranf();
            w = u*(a[i]-aa);
            tt = (0.5*w+aa)*w;
            while (1)
            {
                if (ustar > tt)
                {
                    return (s == 1.0 ? -(aa+w) : aa+w);
                }
                u = ranf();
                if (ustar < u)
                {
                    break;
                }
                tt = u;
                ustar = ranf();
            }
            ustar = ranf();
        }
    }
    w = (ustar-t[i-1])*h[i-1];
    return (s == 1.0 ? -(aa+w) : aa+w);
}

/*
 * int poissonDistr(float mu)
 * Poisson distribution
 * mu --> The mean of the Poisson distribution from which
            a random deviate is to be generated.
*/
int poissonDistr(float mu)
{
    static float a0 = -0.5;
    static float a1 = 0.3333333;
    static float a2 = -0.2500068;
    static float a3 = 0.2000118;
    static float a4 = -0.1661269;
    static float a5 = 0.1421878;
    static float a6 = -0.1384794;
    static float a7 = 0.125006;
    static float muold = 0.0;
    static float muprev = 0.0;
    static float fact[10] = {
        1.0,1.0,2.0,6.0,24.0,120.0,720.0,5040.0,40320.0,362880.0
    };
    static int ignpoi,j,k,kflag,l,m;
    static float b1,b2,c,c0,c1,c2,c3,d,del,difmuk,e,fk,fx,fy,g,omega,p,p0,px,py,q,s,
        t,u,v,x,xx,pp[35];
    int expSamplFlag=1;

    if (mu >= 10.0)
    {
        if (mu != muprev)
        {
            muprev = mu;
            s = sqrt(mu);
            d = 6.0*mu*mu;
            l = (int)(mu-1.1484);
        }
        g = mu+s*normDistr();
        if (g >= 0.0)
        {
            ignpoi = (int)(g);
            if (ignpoi >= l)
            {
                return ignpoi;
            }
            fk = (float)ignpoi;
            difmuk = mu-fk;
            u = ranf();
            if (d*u >= difmuk*difmuk*difmuk)
            {
                return ignpoi;
            }
        }
        if (mu != muold)
        {
            muold = mu;
            omega = 0.3989423/s;
            b1 = 4.166667E-2/mu;
            b2 = 0.3*b1*b1;
            c3 = 0.1428571*b1*b2;
            c2 = b2-15.0*c3;
            c1 = b1-6.0*b2+45.0*c3;
            c0 = 1.0-b1+3.0*b2-15.0*c3;
            c = 0.1069/mu;
        }
        if (g >= 0.0)
        {
            kflag = 0;
            expSamplFlag=0;
        }
        while (1)
        {
            if (expSamplFlag)
            {
                do
                {
                    e = expDistr();
                    u = ranf();
                    u += (u-1.0);
                    t = 1.8+fsign(e,u);
                } while (t <= -0.6744);
                ignpoi = (int) (mu+s*t);
                fk = (float)ignpoi;
                difmuk = mu-fk;
                kflag = 1;
            }
            if(ignpoi < 10)
            {
                px = -mu;
                py = pow((double)mu,(double)ignpoi)/ *(fact+ignpoi);
            }
            else
            {
                del = 8.333333E-2/fk;
                del -= (4.8*del*del*del);
                v = difmuk/fk;
                if (fabs(v) <= 0.25)
                {
                    px = fk*v*v*(((((((a7*v+a6)*v+a5)*v+a4)*v+a3)*v+a2)*v+a1)*v+a0)-del;
                }
                else
                {
                    px = fk*log(1.0+v)-difmuk-del;
                }
                py = 0.3989423/sqrt(fk);
            }
            x = (0.5-difmuk)/s;
            xx = x*x;
            fx = -0.5*xx;
            fy = omega*(((c3*xx+c2)*xx+c1)*xx+c0);
            if (kflag <= 0)
            {
                if (fy-u*fy <= py*exp(px-fx))
                {
                    return ignpoi;
                }
            }
            else
            {
                if (c*fabs(u) <= py*exp(px+e)-fy*exp(fx+e))
                {
                    return ignpoi;
                }
            }
            expSamplFlag=1;
        }
    }
    else
    {
        muprev = 0.0;
        if (mu != muold)
        {
            muold = mu;
            m = (int)(mu) > 1 ? (int)(mu) : 1;
            l = 0;
            p = exp(-mu);
            q = p0 = p;
        }
        while (1)
        {
            u = ranf();
            ignpoi = 0;
            if (u <= p0) 
            {
                return ignpoi;
            }
            if (l != 0)
            {
                j = 1;
                if (u > 0.458)
                {
                    j = l < m ? l : m;
                }
                for(k=j; k <= l; k++)
                {
                    if (u <= pp[k-1])
                    {
                        break;
                    }
                }
                if (k <= l)
                {
                    break;
                }
                if (l == 35)
                {
                    continue;
                }
            }
            l++;
            for (k=l; k <= 35; k++)
            {
                p = p*mu/(float)k;
                q += p;
                pp[k-1] = q;
                if (u <= q)
                {
                    break;
                }
            }
            if (k <= 35)
            {
                l=k;
                break;
            }
            l=35;
        }
    }
    ignpoi = k;
    return ignpoi;
}
//#endif
//long ignpoi(float mu);

DLLEXPORT void rndpoi(float *lambda, int len)
{
    int i;
    for (i=0; i < len; i++)
    {
        lambda[i]=(float)poissonDistr(lambda[i]);
//        lambda[i]=(float)ignpoi(lambda[i]);
    }
}

#if 0
int main(int argc, char **argv)
{
    int i, j, k;
//    int Xa1w = 1033780774L;
//    int iseed1=1234567890L;
//    int Xm1 = 2147483563L;
//    printf("res=%ld (%ld)\n", mltmod(Xa1w, iseed1, Xm1), (Xa1w*iseed1)%Xm1);
    for (j=1; j < 10000; j++)
    {
        printf("### mn %d\n", j);
        for (i=0; i < 10000; i++)
        {
//            printf("%d: %d\n", i, poissonDistr((j+i)%111+1));
//            printf("%d: %d\n", i, ignpoi((j+i)%111+1));
            printf("%d: %d\n", i, poissonDistr(j));
//            printf("%d: %d\n", i, ignpoi(j));
/*            for (k=0; k < 32; k++)
            {
                printf("k=%d: Xlg1=%ld, Xlg2=%ld, Xig1=%ld, Xig2=%ld\n",
                        k, Xlg1[k], Xlg2[k], Xig1[k], Xig2[k]);
            } */
        }
    }
    return 0;
}
#endif

