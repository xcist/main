#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define pi 3.14159265358979

typedef struct TestStruct {
    float    ScanR;
    float    DistD;
    int       YL;
    int       ZL;
    float    dectorYoffset;
    float    dectorZoffset;
    float    XOffSet;
    float    YOffSet;
    float    ZOffSet;
    float    phantomXOffSet;
    float    phantomYOffSet;
    float    phantomZOffSet;
    float    DecFanAng;
    float    DecHeight;
    float    DecWidth;
    float    dx;
    float    dy;
    float    dz;
    float    h;
    float    BetaS;
    float    BetaE;
    int       AngleNumber;
    int       N_2pi;
    float    Radius;
    int       RecSize;
    int       RecSizeZ;
//    int       FOILength;
//    int       FOIWidth;
//    int       FOIHeight;
    float    delta;
    float    HSCoef;
    float    k1;
    float    ***GF;
    float    ***RecIm;
} TestStruct;

extern void fbp(TestStruct *t) {

	float ScanR, DistD,DecL,DecHeight,DecWidth,ObjR,dectorYoffset,dectorZoffset;
	float dx,dy,dz,dYL,dZL,DeltaFai,YLC,ZLC,RadiusSquare,XOffSet,YOffSet,ZOffSet;
	float XNC, YNC, ZNC,startangle;
	float h,h1, BetaE, BetaS, delta, HSCoef, k1;
	int YL,ZL,PN,RecSize, RecSizeZ, N_2pi,N_pi,XN, YN, ZN;
//	int FOILength,FOIWidth,FOIHeight;

	ScanR = t->ScanR;        /*source object distance*/
    DistD = t->DistD;
	DecL = t->DecFanAng;     /*Project Angle*/
	YL = t->YL;              /*Detection Number */
    ZL = t->ZL;
    DecHeight = t->DecHeight;
    DecWidth = t->DecWidth;
    h1 = t->h;
	ObjR = t->Radius;      /*Radius of the field of view*/
	RecSize = t->RecSize;  /*reconstruction size x*/
	RecSizeZ = t->RecSizeZ;  /*reconstruction size x*/
	delta = t->delta;
	HSCoef = t->HSCoef;
	k1 = t->k1;

    BetaS = t->BetaS;
    N_2pi = t->N_2pi;
    PN = t->AngleNumber;
    dx = t->dx;   
	dy = t->dy;
    dz = t->dz;
    dectorYoffset = t->dectorYoffset;
//    dectorZoffset = t->dectorZoffset;
//    startangle = t->startangle;
    XOffSet = t->XOffSet;
    YOffSet = t->YOffSet;
    ZOffSet = t->ZOffSet;
    XN = RecSize;
    XNC = (XN-1)*0.5;
    YN = RecSize;
    YNC = (YN-1)*0.5;
    ZN = RecSizeZ;
    ZNC = (ZN-1)*0.5;
    h = h1*DecHeight;

	dYL= DecL/YL;       /*Each move Angle of projection*/
	dZL= DecHeight/(ZL);       /*Each move Angle of projection*/
	YLC    = (YL-1)*0.5;
	ZLC    = (ZL-1)*0.5 + dectorZoffset;

	RadiusSquare= ObjR*ObjR;
	DeltaFai = 2*pi/N_2pi;
    N_pi = N_2pi/2;

    dYL = DecWidth/YL;
    dZL = DecHeight/ZL;
    DeltaFai = 2*pi/N_2pi;

    float *w;
    w     = (float*)malloc(sizeof(float)*N_2pi);

     //////begin of the  main code
     float x,y,z,Dey,Dez,touying,UU,U1,V1,Beta0,Yr,Zr,View,weight,weight1,weight2,Gama,Gama_C,m1,m2;
     int ProjInd,xi,yi,zi,U,V,s0,s1,s2,d1,d2,L,Shift;

	 for(zi = 0; zi<ZN; zi++)  // create two slices for testing
	 {
         printf("   recon slice %d/%d...\n", zi, ZN);
		 ///compute the projection position for every grid on the image plane
         z = (zi-ZNC) * dz+ZOffSet;
         Beta0 = 2 * pi * z / h;
         s0 = ceil((Beta0-BetaS) / DeltaFai-0.5);
         s1 = s0-ceil(N_pi*HSCoef);       //s1 = ceil((Beta0 + N_Circle * pi - pi ) /DeltaFai);
         s2 = s0+ceil(N_pi*HSCoef)-1;     //s2 = ceil((Beta0 + N_Circle * pi + pi ) /DeltaFai);

         if ((s1<PN)||(s2>0))
         {
           if (s1 < 0)  {s1 = 0;}
           if (s2 > PN-1) {s2 = PN-1;}
          //////////////////////////////////////////
         ////Producing the weighting function
         for (int k=0;k<N_2pi;k++)
           { w[k] = 0;}
         L = s2-s1+1;
         Shift = N_pi - (s0-s1);

         if (L<2*delta)
         {
             for (int k=0;k<L;k++)
               w[k+Shift]= pow(cos((pi/2)*(2*k-L+1)/L),2);
         }
         else
         {
           for (int k=0;k<L;k++)
           {
             if (0 <= k && k<delta)
                 w[k+Shift]= pow(cos((pi/2)*(delta-k-0.5)/delta),2);
                
             else if(L-delta<=k && k < L)
                 w[k+Shift]= pow(cos((pi/2)*(k-(L-delta)+0.5)/delta),2);
             else
                 w[k+Shift] = 1;
          }
         }//if(L<2*delta)
         ///////////////////////////////////////////
         for (ProjInd = s1; ProjInd <= s2; ProjInd++ )
         {
             View = BetaS + ProjInd * DeltaFai;
             d1   = N_pi-(s0-ProjInd); //d1 = ProjInd;
             if (ProjInd < s0)
             {
                 d2 = d1+N_pi;
             }
             else //(ProjInd >= s0)
             {
                 d2 = d1-N_pi;
             }

             for(yi=0;yi<YN;yi++)
             {
                 y = -(yi-YNC)*dy-YOffSet;
                 #pragma omp parallel for private(xi,x, UU, Yr, Zr, U1, U,V1,V, Dey,Dez,touying,weight1,weight2,Gama,Gama_C,m1,m2,weight)
                 for(xi=0;xi<XN;xi++)
                 {
                    x  = -(xi-XNC)*dx-XOffSet;
                    UU = -x*cos(View)-y*sin(View);
				    Yr = -x*sin(View)+y*cos(View);
                    Zr = (z-h*(View+asin(Yr/ScanR))/(2.0*pi))*(DistD)/(sqrt(ScanR*ScanR-Yr*Yr)+UU);///03/05/23 Yu
                    U1 = Yr/dYL+YLC;
                    U  = ceil(U1);
                    V1 = Zr/dZL+ZLC;
                    V  = ceil(V1);
                    Dey = U-U1;
                    Dez = V-V1;
                    //Linear interploate
                    if ((U>0)&&(U<YL)&&(V>0)&&(V<ZL))
                    {
                          touying = Dey*Dez*t->GF[U-1][V-1][ProjInd]
                                    +Dey*(1-Dez)*t->GF[U-1][V][ProjInd]
                                   +(1-Dey)*Dez*t->GF[U][V-1][ProjInd]
                                   +(1-Dey)*(1-Dez)*t->GF[U][V][ProjInd];

                           weight1 = w[d1];
                           weight2 = w[d2];
                       
                           Gama   = fabs((z-h*View/(2.0*pi))/(sqrt(ScanR*ScanR-Yr*Yr)+UU));
                           if (ProjInd < s0)
                           {
                               Gama_C = fabs((z-h*(View+pi)/(2.0*pi))/(sqrt(ScanR*ScanR-Yr*Yr)-UU));
                           }
                           else
                           {
                               Gama_C = fabs((z-h*(View-pi)/(2.0*pi))/(sqrt(ScanR*ScanR-Yr*Yr)-UU));
                           }
                           m1 = pow(Gama,  k1);    //m1     = std::real(std::pow(Gama,k1*h));
                           m2 = pow(Gama_C, k1);  //m2     = std::real(std::pow(Gama_C,k1*h));
                           weight = (weight1*m2)/(weight2*m1+weight1*m2);
                           
                           t->RecIm[yi][xi][zi]=t->RecIm[yi][xi][zi]+weight*touying*DeltaFai;
                    }
                   }	//yi
			 }//xi
         }//ProjInd
         }//  if ((s1<PN)||(s2>0))
	 } //zi
     //////end of the main code
    free(w);
 }
