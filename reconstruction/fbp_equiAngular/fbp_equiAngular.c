#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define PI 3.14159265358979

typedef struct TestStruct {
    double       ScanR;
    double     DecFanAng;
    int       YL;
    int       AngleNumber;
    double       Radius;
    int       RecSizeX;
    int       RecSizeY;
    int       centerX;
    int       centerY;
    int       FOILength;
    int       FOIWidth;
    double    **GF;
    double    **RecIm;
} TestStruct;

extern void fbp(TestStruct *t) {

	double ScanR, DecL,Radius,DeltaY,DeltaL,YCtr,DeltaRX,DeltaRY,RCtrX,RCtrY,RadiusSquare;
	int YL,ProjNum,RecSizeX,RecSizeY,centerX,centerY,FOILength,FOIWidth,xl,xr,yl,yr;

	ScanR = t->ScanR;
	DecL = t->DecFanAng;
	YL = t->YL;
	ProjNum = t->AngleNumber;
	Radius = t->Radius;
	RecSizeX = t->RecSizeX;
	RecSizeY = t->RecSizeY;

    centerX = t->centerX;
    centerY = t->centerY;
    FOILength = t->FOILength;
    FOIWidth = t->FOIWidth;

	DeltaY  = DecL/YL;
	YCtr    = (YL-1)*0.5;
	DeltaRX  = 2*Radius/RecSizeX;
	DeltaRY  = 2*Radius/RecSizeY;
    RCtrX    = (RecSizeX-1)*0.5;
    RCtrY    = (RecSizeY-1)*0.5;
	RadiusSquare= Radius*Radius;
	DeltaL = 2*PI/ProjNum;

    xl = (int)(centerX - FOILength*0.5);
    xr = (int)(centerX + FOILength*0.5);
    yl = (int)(centerY - FOIWidth*0.5);
    yr = (int)(centerY + FOIWidth*0.5);

    /*Compute the coordinates of scanning locus and its corresponding local coordinate*/

	double  *VectorS, *VectorE,*xCor,*yCor,*zCor;
    int  loop;
	double temp;

	VectorS  = (double*)malloc(sizeof(double)*2*ProjNum);
	VectorE = (double*)malloc(sizeof(double)*2*ProjNum);
	xCor     = (double*)malloc(sizeof(double)*RecSizeX);
	yCor     = (double*)malloc(sizeof(double)*RecSizeY);


	for(loop=0;loop<ProjNum;loop++)
	{
		temp = (loop+0)*DeltaL;
		VectorS[loop*2  ] = ScanR*cos(temp);
		VectorS[loop*2+1] = ScanR*sin(temp);
		VectorE[loop*2  ]= cos(temp);
		VectorE[loop*2+1]= sin(temp);
	}

	for(loop=0;loop<RecSizeX;loop++)
		xCor[loop] = (loop-RCtrX+xl)*DeltaRX;
	for(loop=0;loop<RecSizeY;loop++)
		yCor[loop] = -(loop-RCtrY+yl)*DeltaRY;

    /* argument about object*/

	int i,j,ProjIndex,UU,UL;
	double xcor,ycor,theta,cost,sint,dis,tpdata,Dlocal;
	double DSX[2],UCor,alfa;


	for (i=0;i<FOILength;i++)
		{
			for (j=0;j<FOIWidth;j++)
			{
				tpdata = xcor*xcor+ycor*ycor-RadiusSquare;
				if(tpdata<0)
				{
                    for(ProjIndex=0;ProjIndex<ProjNum;ProjIndex++)
					{
	                    DSX[0]= xCor[i]-VectorS[ProjIndex*2];
						DSX[1]= yCor[j]-VectorS[ProjIndex*2+1];
						Dlocal = sqrt(DSX[0]*DSX[0]+DSX[1]*DSX[1]);
						dis   = fabs(xCor[i]*VectorE[ProjIndex*2]+yCor[j]*VectorE[ProjIndex*2+1]-ScanR);
						UCor  = atan((DSX[0]*VectorE[ProjIndex*2+1]-DSX[1]*VectorE[ProjIndex*2])/dis);
						UCor  = UCor/DeltaY+YCtr+0;
						UL    = (int)UCor;
						UU    = UL+1;
						alfa  = UU-UCor;
                        if((UL>0)&(UU<YL))
                        {
                         t->RecIm[i][j] = (t->GF[UU][ProjIndex]*(1-alfa)+t->GF[UL][ProjIndex]*alfa)/Dlocal+t->RecIm[i][j];
                        }
					}//for(projindex=0;Projindex<ProjNum;Projindex++)
					t->RecIm[i][j] = -t->RecIm[i][j]/1440;
				}//if(tpdata<0)
			}//for (y=0;y<RecSize;y++)
		}//for (x=0;x<RecSize;x++)
	free(VectorS);
	free(VectorE);
	free(xCor);
	free(yCor);
 }
