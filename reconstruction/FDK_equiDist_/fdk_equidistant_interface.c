#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define PI 3.14159265358979

typedef struct TestStruct {
    double       ScanR;
    double     DecWidth;
    double     DecHeigh;
    int       YL;
    int       ZL;
    double    YOffSet;
    double    ZOffSet;
    int       ProjScale;
    double       DistD;
    double       Radius;
    int       RecSize;
    int       centerX;
    int       centerY;
    int       centerZ;
    int       FOILength;
    int       FOIWidth;
    int       FOIHeigh;
    double    ***GF;
    double    ***RecIm;
} TestStruct;

extern void fbp(TestStruct *t) {

	double ScanR, DecWidth,DecHeigh,Radius;
	double DeltaY,DeltaZ,YCtr,ZCtr,DeltaR,DeltaL,RCtr,RadiusSquare,YOffSet,ZOffSet;
	int YL,ZL,ProjScale,RecSize, DistD;
	int centerX,centerY,centerZ,FOILength,FOIWidth,FOIHeigh,xl,yl,zl;

	ScanR = t->ScanR;
	DecWidth = t->DecWidth;
	DecHeigh = t->DecHeigh;
	YL = t->YL;
	ZL = t->ZL;
	ProjScale = t->ProjScale;
	DistD = t->DistD;
	Radius = t->Radius;
	RecSize = t->RecSize;

    centerX = t->centerX;
    centerY = t->centerY;
    centerZ = t->centerZ;
    FOILength = t->FOILength;
    FOIWidth = t->FOIWidth;
    FOIHeigh = t->FOIHeigh;
    YOffSet = t->YOffSet;
    ZOffSet = t->ZOffSet;

	DeltaY  = DecWidth/YL;
	DeltaZ  = DecHeigh/ZL ;
	YCtr    = (YL-1)*0.5 + YOffSet;
	ZCtr    = (ZL-1)*0.5 + ZOffSet;
	DeltaR  = 2*Radius/RecSize;
    RCtr    = (RecSize-1)*0.5;
	RadiusSquare= Radius*Radius;
	DeltaL = 2*PI/ProjScale;

    xl = (int)(centerX - FOILength*0.5);
    yl = (int)(centerY - FOIWidth*0.5);
    zl = (int)(centerZ - FOIHeigh*0.5);

	/*Compute the coordinates of scanning locus and its corresponding local coordinate*/
	double  *VectorS, *VectorE,*xCor,*yCor,*zCor;
    int  loop;
	double temp;

	VectorS  = (double*)malloc(sizeof(double)*2*ProjScale);
	VectorE = (double*)malloc(sizeof(double)*2*ProjScale);
	xCor     = (double*)malloc(sizeof(double)*RecSize);
	yCor     = (double*)malloc(sizeof(double)*RecSize);
	zCor     = (double*)malloc(sizeof(double)*RecSize);

	for(loop=0;loop<ProjScale;loop++)
	{
		temp = (loop+180)*DeltaL;
		VectorS[loop*2  ] = ScanR*cos(temp);
		VectorS[loop*2+1] = ScanR*sin(temp);
		VectorE[loop*2  ]= cos(temp);
		VectorE[loop*2+1]= sin(temp);
	}

	for(loop=0;loop<RecSize;loop++)
		xCor[loop] = (loop-RCtr)*DeltaR;
	for(loop=0;loop<RecSize;loop++)
		yCor[loop] = (loop-RCtr)*DeltaR;
	for(loop=0;loop<RecSize;loop++)
		zCor[loop] = (loop-RCtr)*DeltaR;

	/* argument about object:  ObjR RecMX RecMY RecMZ*/

	int i,j,k,ProjIndex,UU,UL,VL,VV;
	double xcor,ycor,zcor,theta,cost,sint,dis,tpdata,Dlocal;
	double m1,m2,UCor,VCor,alfa,beta, W2;

    for (i=0; i<FOILength; i++)
        {
        for (j=0;j<FOIWidth;j++)
            {
                for (k=0;k<FOIHeigh;k++)
//                for (k=48;k<50;k++)
                {
                    if( (xCor[i]*xCor[i] + yCor[j]*yCor[j]) < 2*RadiusSquare)
                    {
                        t->RecIm[i][j][k] = 0;
                        for(ProjIndex=0;ProjIndex<ProjScale;ProjIndex++)
                        {
                            m1= xCor[i]-VectorS[ProjIndex*2];
                            m2= yCor[j]-VectorS[ProjIndex*2+1];
                            Dlocal= sqrt(m1*m1+m2*m2+zCor[k]*zCor[k]);
                            dis   = fabs(xCor[i]*VectorE[ProjIndex*2]+yCor[j]*VectorE[ProjIndex*2+1]-ScanR);
                            UCor  = -(m1*VectorE[ProjIndex*2+1]-m2*VectorE[ProjIndex*2])*DistD/dis;
                            VCor  = zCor[k]*DistD/dis;
                            W2    = sqrt(DistD*DistD+UCor*UCor+VCor*VCor)/Dlocal;
                            UCor  = UCor/DeltaY+YCtr+0;
                            VCor  = VCor/DeltaZ+ZCtr+0;
                            UL    = (int)UCor;
                            VL    = (int)VCor;
                            UU    = UL+1;
                            VV    = VL+1;
                            alfa  = UU-UCor;
                            beta  = VV-VCor;
                            if((UL>0)&(UU<YL)&(VL>0)&(VV<ZL))
                            {
                             t->RecIm[i][j][k] = (t->GF[UU][VV][ProjIndex]*(1-alfa)*(1-beta)+
                                                  t->GF[UL][VV][ProjIndex]*alfa*(1-beta)+
                                                  t->GF[UU][VL][ProjIndex]*(1-alfa)*beta+
                                                  t->GF[UL][VL][ProjIndex]*alfa*beta)*W2+t->RecIm[i][j][k];
                            }
                        }//for(projindex=0;Projindex<ProjNum;Projindex++)
                        t->RecIm[i][j][k] = t->RecIm[i][j][k]/(2*ProjScale);
                    }//if(tpdata<0)
                }//for (k=97;k<98;k++)
            }//for (j=0;j<RecSize;j++)
        }//for (i=0; i<RecSize; i++)
        free(VectorS);
	    free(VectorE);
	    free(xCor);
	    free(yCor);
	    free(zCor);
 }