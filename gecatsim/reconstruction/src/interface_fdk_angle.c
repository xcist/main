#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define PI 3.14159265358979

typedef struct TestStruct {
    float    ScanR;
    float    DecFanAng;
    float    startangle;
    int      rotdir; //rotation direction
    float    DecHeight;
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
    int       AngleNumber;
    float    DistD;
    float    Radius;
    int       RecSize;
    float    sliceThickness;
    int       FOILength;
    int       FOIWidth;
    int       FOIHeight;
    float    ***GF;
    float    ***RecIm;
} TestStruct;

extern void fbp(TestStruct *t) {

	float ScanR, DecL,DecHeight,DistD,Radius,sliceThickness,dectorYoffset,dectorZoffset;
	float DeltaY,DeltaZ,DeltaL,YCtr,ZCtr,DeltaR,RCtr,RadiusSquare,XOffSet,YOffSet,ZOffSet;
	float centerX,centerY,centerZ,startangle,phantomXOffSet,phantomYOffSet,phantomZOffSet;
	int YL,ZL,ProjNum,RecSize, rotdir;
	int FOILength,FOIWidth,FOIHeight;

	ScanR = t->ScanR;        /*source object distance*/
	DecL = t->DecFanAng;     /*Project Angle*/
	YL = t->YL;              /*Detection Number */
	ProjNum = t->AngleNumber;/*Projection Number */
	Radius = t->Radius;      /*Radius of the field of view*/
	RecSize = t->RecSize;  /*reconstruction size x*/
    sliceThickness = t->sliceThickness;
    DistD = t->DistD;
    DecHeight = t->DecHeight;
    ZL = t->ZL;
    dectorYoffset = t->dectorYoffset;
    dectorZoffset = t->dectorZoffset;
    phantomXOffSet = t->phantomXOffSet;
    phantomYOffSet = t->phantomYOffSet;
    phantomZOffSet = t->phantomZOffSet;

    startangle = t->startangle;
    rotdir = t->rotdir;
    FOILength = t->FOILength;
    FOIWidth = t->FOIWidth;
    FOIHeight = t->FOIHeight;
    XOffSet = t->XOffSet;
    YOffSet = t->YOffSet;
    ZOffSet = t->ZOffSet;
    centerX = (FOILength-1)*0.5;
    centerY = (FOIWidth-1)*0.5;
    centerZ = (FOIHeight-1)*0.5;
	DeltaY  = DecL/YL;       /*Each move Angle of projection*/
	DeltaZ  = DecHeight/(ZL);       /*Each move Angle of projection*/
	YCtr    = (YL-1)*0.5+dectorYoffset;
	ZCtr    = (ZL-1)*0.5+dectorZoffset;
	DeltaR  = 2*Radius/RecSize;
    RCtr    = (RecSize-1)*0.5;
	RadiusSquare= Radius*Radius;
	DeltaL = 2*PI/ProjNum;

	float  *VectorS, *VectorE,*xCor,*yCor,*zCor;
    int  loop;
	float temp;

	VectorS  = (float*)malloc(sizeof(float)*2*ProjNum);
	VectorE  = (float*)malloc(sizeof(float)*2*ProjNum);
	xCor     = (float*)malloc(sizeof(float)*RecSize);
	yCor     = (float*)malloc(sizeof(float)*RecSize);
	zCor     = (float*)malloc(sizeof(float)*RecSize);

	for(loop=0;loop<ProjNum;loop++)
	{
		temp = (rotdir*loop+startangle)*DeltaL;
		VectorS[loop*2  ] = ScanR*cos(temp);
		VectorS[loop*2+1] = ScanR*sin(temp);
		VectorE[loop*2  ]= cos(temp);
		VectorE[loop*2+1]= sin(temp);
	}

	for(loop=0;loop<RecSize;loop++) {
		xCor[loop] = (loop-centerX)*DeltaR-XOffSet-phantomXOffSet;}
	for(loop=0;loop<RecSize;loop++) {
		yCor[loop] = (loop-centerY)*DeltaR-YOffSet-phantomYOffSet;}
    for(loop=0;loop<RecSize;loop++) {
        zCor[loop] = (loop-centerZ)*sliceThickness-ZOffSet-phantomZOffSet;}

	/* argument about object:  ObjR RecMX RecMY RecMZ*/

	int i,j,k,ProjIndex,UU,UL,VV,VL;//,IR;
	float dis,Dlocal, DSX[2],UCor,VCor,alfa,beta;

    //IR = (int)(RecSize/10);

#pragma omp parallel for collapse(3) private(i,j,k,ProjIndex, DSX, Dlocal, dis, UCor, VCor, UL, VL, UU, VV, alfa, beta)
	for (i=0;i<FOILength;i++)
		{
            //if ((i%IR)==0)
            //printf("Reconstruction process:  %d %% \n", 10*i/IR);
			for (j=0;j<FOIWidth;j++)
			{
                for (k=0;k<FOIHeight;k++)

                {
//				if( (xCor[i]*xCor[i] + yCor[j]*yCor[j])< 2*RadiusSquare)
//			    	{
					t->RecIm[i][j][k] = 0;
                        for(ProjIndex=0;ProjIndex<ProjNum;ProjIndex++)
                        {
                            DSX[0]= xCor[i]-VectorS[ProjIndex*2];
                            DSX[1]= yCor[j]-VectorS[ProjIndex*2+1];
                            Dlocal = sqrt(DSX[0]*DSX[0]+DSX[1]*DSX[1]+zCor[k]*zCor[k]);
                            dis   = fabs(xCor[i]*VectorE[ProjIndex*2]+yCor[j]*VectorE[ProjIndex*2+1]-ScanR);
                            UCor  = atan((DSX[0]*VectorE[ProjIndex*2+1]-DSX[1]*VectorE[ProjIndex*2])/dis);
                            VCor  = (zCor[k]*DistD)/dis;
                            VCor  = VCor/DeltaZ+ZCtr;
                            UCor  = UCor/DeltaY+YCtr;
                            UL    = (floor)(UCor);
                            VL    = (floor)(VCor);
                            UU    = UL+1;
                            VV    = VL+1;
                            alfa  = UU-UCor;
                            beta  = VV-VCor;
                            if((UL>=0)&(UU<YL)&(VL>=0)&(VV<ZL))
                            {
                             t->RecIm[i][j][k] = (t->GF[UU][VL][ProjIndex]*(1-alfa)*beta+
                                                  t->GF[UU][VV][ProjIndex]*(1-alfa)*(1-beta)+
                                                  t->GF[UL][VL][ProjIndex]*alfa*beta+
                                                  t->GF[UL][VV][ProjIndex]*alfa*(1-beta)
                                                  )/(Dlocal*Dlocal)+t->RecIm[i][j][k];
                            }
                        }//for(projindex=0;Projindex<ProjNum;Projindex++)
					t->RecIm[i][j][k] = -ScanR*PI*t->RecIm[i][j][k]/(ProjNum);

//				}
            }//for (j=0;j<RecSize;j++)
        }//for (i=0;i<RecSize;i++)
    }//for (k=0;i<RecSize;k++)
    free(VectorS);
	free(VectorE);
	free(xCor);
	free(yCor);
	free(zCor);
 }
