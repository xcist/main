// Copyright 2020, General Electric Company. All rights reserved. See https://github.com/xcist/code/blob/master/LICENSE

#include <stdio.h>

/* colnum, rownum, xspot, zspot, viewnum */
/* start = rownum*linesize + xspot*rownum*colnum + zspot*xcount*rownum*colnum*/

#ifdef WIN32
__declspec(dllexport)
#endif
void ExtractSino(char *filename, float *sino, 
		 int colcount, int rowcount, int xycount, int zcount, int viewcount,
		 int rownum, int znum) {
  FILE *fp;
  int i, j;
  fp = fopen(filename,"rb");
  for (i=0;i<viewcount;i++) {
    for (j=0;j<xycount;j++) {
      fseek(fp,sizeof(float)*(rownum*colcount+j*rowcount*colcount+znum*xycount*rowcount*colcount+i*colcount*rowcount*xycount*zcount),SEEK_SET);
      fread(sino+i*xycount*colcount+j*colcount,sizeof(float),colcount,fp);
    }
  }
  fclose(fp);
}

#ifdef WIN32
__declspec(dllexport)
#endif
void DistributeSino(char *filename, float *sino, int rowcount, int rownum, int colcount, int viewcount) {
  FILE *fp;
  char buffer[1000];
  int i;
  for (i=0;i<viewcount;i++) {
    sprintf(buffer,filename,i);
    fp = fopen(buffer,"r+");
    fseek(fp,(i*rowcount*colcount+rownum*colcount)*sizeof(float),SEEK_SET);
    fwrite(sino+i*colcount,sizeof(float),colcount,fp);
    fclose(fp);
  }
}


