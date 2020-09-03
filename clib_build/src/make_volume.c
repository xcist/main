// Copyright 2020, General Electric Company. All rights reserved. See https://github.com/xcist/code/blob/master/LICENSE

#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#define MIN(a,b) (((a) > (b)) ? (b) : (a))
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

double* datavec[8] = {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
double* *Xdata=NULL;
int* XdataL=NULL;
int datavecL[8];
int* idatavec[8] = {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
int idatavecL[8];
double* vector;

double solve_cubic(double *a)
     
{
     double a12,Q,Q3,R,D,theta,srD,S,T;

     a[0]=a[0]*(1.0/3.0);
     a12=a[0]*a[0];
     Q=(1.0/3.0)*a[1]-a12;
     Q3=Q*Q*Q;
     R=0.5*(a[0]*a[1]-a[2])-a12*a[0];
     D=Q3+R*R;
     if (D<0) {
       theta=acos(R/sqrt(-Q3)); 
       return 2*sqrt(-Q)*cos(theta*(1.0/3.0))-a[0];
     } else {
       srD=sqrt(D);
       S=cbrt(R+srD);
       T=cbrt(R-srD);
       return S+T-a[0];
     }

}

double magnitude(double r, double i)

{
  return sqrt(r*r+i*i);
}

void sqrtm(double in_r,double in_i,double *out_r,double *out_i)

{
  double ang,mag;
  ang=atan2(-in_i,-in_r)/2.0;
  mag=sqrt(magnitude(in_r,in_i));
  out_r[0]=mag*cos(ang);
  out_i[0]=mag*sin(ang);
}

void complex_multiply(double r1, double i1, double r2, double i2,double *real_out, double *imag_out)

{
  real_out[0]=r1*r2-i1*i2;
  imag_out[0]=r1*i2+r2*i1;
}

void solve_cubic_all(double *a,double *zr,double *zi)

{
  double a0,a1,a2,p,q,tmp,tmpr,tmpi,ur[3],ui[3],mag,ang,rmag2,pi=3.141592653589793;
  int i;
  //Uses Cardano's method
  
  a2=a[0];
  a1=a[1];
  a0=a[2];
  
  p=a1-a2*a2/3.0;
  q=a0+(2.0*a2*a2*a2-9.0*a2*a1)/27.0;
  
  tmp=q*q/4.0+p*p*p/27.0;
  //printf("P Q TMP: %1.12lf %1.12lf %1.12lf\n\r",p,q,tmp);
  if (q>0.0)
    if(tmp>=0)
      {ur[0]=cbrt(q/2.0+sqrt(tmp));ui[0]=0;mag=ur[0];ang=0;}
    else
      {
	tmpr=q/2.0;tmpi=sqrt(-tmp);
	mag=cbrt(sqrt(tmpr*tmpr+tmpi*tmpi));
	ang=atan2(tmpi,tmpr)/3.0;
	ur[0]=mag*cos(ang);
	ui[0]=mag*sin(ang);
      }
  else
    if(tmp>=0)
      {ur[0]=cbrt(q/2.0-sqrt(tmp));ui[0]=0;mag=ur[0];ang=0;}
    else
      {
	tmpr=q/2.0;tmpi=-sqrt(-tmp);
	mag=cbrt(magnitude(tmpr,tmpi));
	ang=atan2(tmpi,tmpr)/3.0;
	ur[0]=mag*cos(ang);
	ui[0]=mag*sin(ang);
      }
  ang+=2.0*pi/3.0;
  ur[1]=mag*cos(ang);
  ui[1]=mag*sin(ang);

  ang+=2.0*pi/3.0;
  ur[2]=mag*cos(ang);
  ui[2]=mag*sin(ang);

  //  for(i=0;i<3;i++) printf("%1.12lf %1.12lf %1.12lf\r\n",ur[i],ui[i],magnitude(ui[i],ur[i]));
  
  rmag2=1.0/(mag*mag);

  for(i=0;i<3;i++) {
    zr[i]=p/3.0*(ur[i]*rmag2)-ur[i]-a2/3;
    zi[i]=p/3.0*(-ui[i]*rmag2)-ui[i];
  }
  // for(i=0;i<3;i++) printf("%1.12lf %1.12lf %1.12lf\r\n",zr[i],zi[i],magnitude(zi[i],zr[i]));
}

int compare_doubles (const void *a, const void *b)

{
  const double *da = (const double *) a;
  const double *db = (const double *) b;
  
  return (*da > *db) - (*da < *db);
}

int comp (const void *a, const void *b)

{
  const int da = * (const int *) a;
  const int db = * (const int *) b;
  
  return (vector[da]>vector[db]) - (vector[da]<vector[db]);
}

int solve_quartic2(double *a,double *z)

{
// Faucette's method:
//   Faucette, W. M. "A Geometric Interpretation of the Solution of the General Quartic Polynomial." 
//   Amer. Math. Monthly 103, 51-57, 1996. 

  double yr[3],yi[3],a0,a1,a2,a3,p,q,r,cub[3],sr[3],si[3],zr[4],zi[4],errr,erri,err,errt; 
  int num_real,i;
  
  a3=a[0];
  a2=a[1];
  a1=a[2];
  a0=a[3];
  
  //  First we solve the resolvent cubic:
  p=a2-0.375*a3*a3;
  q=a1-0.5*a2*a3+0.125*a3*a3*a3;
  r=a0-0.25*a1*a3+0.0625*a2*a3*a3-(3.0/256.0)*a3*a3*a3*a3;
  cub[0]=-2.0*p;
  cub[1]=p*p-4.0*r;
  cub[2]=q*q;

  //  for(i=0;i<3;i++) printf("input: %1.12lf\n\r",cub[i]);

  solve_cubic_all(cub,yr,yi);

  //  for(i=0;i<3;i++) printf("Y%d :  %1.12lf +i %1.12lf\n\r",i,yr[i],yi[i]);

  for(i=0;i<3;i++) sqrtm(yr[i],yi[i],&sr[i],&si[i]);

  //for(i=0;i<3;i++) printf("S%d :  %1.12lf +i %1.12lf\n\r",i,sr[i],si[i]);
  
  zr[0]=(sr[0]+sr[1]-sr[2])*0.5;
  zr[1]=(sr[2]+sr[0]-sr[1])*0.5;
  zr[2]=(sr[2]+sr[1]-sr[0])*0.5;
  zr[3]=-zr[0]-zr[1]-zr[2];

  zi[0]=(si[0]+si[1]-si[2])*0.5;
  zi[1]=(si[2]+si[0]-si[1])*0.5;
  zi[2]=(si[2]+si[1]-si[0])*0.5;
  zi[3]=-zi[0]-zi[1]-zi[2];
  
  for(i=0;i<4;i++){
    complex_multiply(zr[i],zi[i],zr[i],zi[i],&errr,&erri);
    errr+=p;
    complex_multiply(errr,erri,zr[i],zi[i],&errr,&erri);
    errr+=q;
    complex_multiply(errr,erri,zr[i],zi[i],&errr,&erri);
    errr+=r;
    errr=magnitude(errr,erri);
    errt+=errr*errr;
  }
  errt=sqrt(errt);
  
  if (errt>1e-5){
    zr[0]+=sr[2];
    zr[1]-=sr[2];
    zr[2]-=sr[2];
    zr[3]+=sr[2];
    zi[0]+=si[2];
    zi[1]-=si[2];
    zi[2]-=si[2];
    zi[3]+=si[2];
  }
  for(i=0;i<4;i++) zr[i]-=a3/4;
  num_real=0;
  for(i=0;i<4;i++) if ((fabs(zi[i])/fabs(zr[i]))<1e-10) {z[num_real]=zr[i];num_real++;}
  qsort(z,num_real,sizeof(double),compare_doubles);
  return num_real;
}

int solve_quartic(double *a,double *z)

{
  // Equations found in: 
  //   Weisstein, Eric W. "Quartic Equation." From MathWorld--A Wolfram Web Resource. 
  //   http://mathworld.wolfram.com/QuarticEquation.html 
  
  double tmp,y1,tmp2,R,alpha,beta,D,E,Dt,Et,cub[3];
  int i;
  
  tmp=a[0];
  for(i=0;i<4;i++) a[i]=a[i+1]/tmp;
  
  tmp=0.25*a[0]*a[0]-a[1];
  cub[0]=-a[1];
  cub[1]=a[2]*a[0]-4.0*a[3];
  cub[2]=4.0*a[1]*a[3]-a[2]*a[2]-a[0]*a[0]*a[3];
  y1=solve_cubic(cub);

  
  //  Next we compute R, D, and E:
  tmp2=tmp+y1;
  //printf("y1  = %1.12lf\n\rtmp = %1.12lf     %1.12lf\n\r",y1,tmp,fabs(tmp2)/MAX(1.0,fabs(tmp)));
  if (fabs(tmp2)/MAX(1.0,fabs(tmp))>1e-4){
    if (tmp2>0){
      R=sqrt(tmp2);
      alpha=0.75*a[0]*a[0]-R*R-2.0*a[1];
      beta=0.25*(4.0*a[0]*a[1]-8.0*a[2]-a[0]*a[0]*a[0])/R;
      Dt=alpha+beta;
      Et=alpha-beta;
      //  printf("Dt = %1.12f\n\rEt = %1.12f\n\r",Dt,Et);

      if (Dt>0){
	D=sqrt(Dt);
	z[1]=-0.25*a[0]+0.5*R+0.5*D;
	z[0]=-0.25*a[0]+0.5*R-0.5*D;
	//printf("                          E = %1.12lf %1.12lf %1.12lf\n\r",D,R,z[0]);
	if (Et>0){
	  E=sqrt(Et);
	  z[3]=-0.25*a[0]-0.5*R+0.5*E;
	  z[2]=-0.25*a[0]-0.5*R-0.5*E;
	  if (z[2]<z[1]){
	    tmp=z[1];z[1]=z[2];z[2]=tmp;
	    if (z[1]<z[0]){tmp=z[0];z[0]=z[1];z[1]=tmp;}
	    if (z[3]<z[2]){tmp=z[2];z[2]=z[3];z[3]=tmp;if (z[2]<z[1]){tmp=z[1];z[1]=z[2];z[2]=tmp;}}
	  }
	  return 4;
	} else return 2;
      } else{
	if (Et>0){
	  E=sqrt(Et);
	  z[1]=-0.25*a[0]-0.5*R+0.5*E;
	  z[0]=-0.25*a[0]-0.5*R-0.5*E;
	  return 2;
	} else return 0;
      }
    } 
    else return 0;
  }
  else {
    //for(i=0;i<4;i++) printf("input: %1.12lf\n\r",a[i]);
    return solve_quartic2(a,z);
    //    return -1;
  }
}

double quadratic_form(double *vec1,double *matrix, double *vec2)
{
  int i;
  double x=0;
  for(i=0;i<3;i++){x+=(matrix[i]*vec2[0]+matrix[i+3]*vec2[1]+matrix[i+6]*vec2[2])*vec1[i];}
  return x;
}

int quartic_intersect(double *a0,double *alpha,double *tc,int obj)

{
  int i,out;
  double scale,displ,a0t[3],tmp[3],aa,b,c,d,e,f,C[5],*Ql,*Qr,sh;
  
  Ql=&datavec[4][obj*18];
  Qr=&datavec[4][obj*18+9];
  sh=datavec[3][obj];
  scale=1.0;
  displ=0;
  for(i=0;i<3;i++) displ+=a0[i]*a0[i];
  displ=sqrt(displ);                  // scaling/displacement done only f0r numerical reasons

  for(i=0;i<3;i++) a0t[i]=(a0[i]+displ*alpha[i])/scale;

  // The quartic formula is of the form    A t^4 + B t^3 + C t^2 + D t + E = 0
  aa=quadratic_form(alpha,Ql,alpha);
  b=2.0*quadratic_form(a0t,Ql,alpha);
  //c=1-sh*sh+quadratic_form(a0t,Ql,a0t);
  c=sh+quadratic_form(a0t,Ql,a0t);
  d=4.0*quadratic_form(alpha,Qr,alpha);
  e=8.0*quadratic_form(a0t,Qr,alpha);
  f=4.0*quadratic_form(a0t,Qr,a0t);

  C[0]=aa*aa;
  C[1]=2*aa*b;
  C[2]=b*b+2*aa*c-d;
  C[3]=2*b*c-e;
  C[4]=c*c-f;

  out=solve_quartic(C,tc);
  for(i=0;i<out;i++) tc[i]=tc[i]*scale+displ;
  return out;
}

int clip_all(double *a,double *alpha,double maxD,double *tc2,int out,double *st_list,double *en_list,double *den_list,int *pri_list,int *mat_list,int num_int,int i)

{
  double b[3],s1,s2,tcrit,tmin,tmax,*eta,*s,den;
  int j,cp,ddd;

  ddd=datavec[2][i];
  eta=&datavec[0][ddd*3];
  s=&datavec[1][ddd];
  cp=idatavec[0][i];
  den=datavec[6][i];
  tmin=-1e50;tmax=1e50;
  for(j=0;j<3;j++) b[j]=a[j]+maxD*alpha[j];
  
  for(j=0;j<cp;j++){
    s1=(eta[j*3]*a[0]+eta[j*3+1]*a[1]+eta[j*3+2]*a[2]-s[j]);
    s2=(eta[j*3]*b[0]+eta[j*3+1]*b[1]+eta[j*3+2]*b[2]-s[j]);
    if(s1*s2<0){
      tcrit=fabs(s1)/(fabs(s1)+fabs(s2))*maxD;
      if (s1<s2) 
	{if (tcrit<tmax) tmax=tcrit;}
      else 
	{if (tcrit>tmin) tmin=tcrit;}
    }
    else 
    if ((s1+s2)>0){
      tmin=0;tmax=0;
      break;
    }
  }
  
  if(out>2) {
    if(tmax<tc2[2])
      out=2;
    else 
      if(tmin>tc2[1]) {
	tc2[0]=tc2[2];
	tc2[1]=tc2[3];
      }
      else{
	tc2[3]=MIN(tmax,tc2[3]);
	st_list[num_int]=tc2[2];
	en_list[num_int]=tc2[3];
	den_list[num_int]=den;
	pri_list[num_int]=i;
	mat_list[num_int]=1;  // needs to be changed
	num_int=num_int+1;
      }
  }
  tc2[1]=MIN(tmax,tc2[1]);
  tc2[0]=MAX(tmin,tc2[0]);
  if (tc2[1]>tc2[0]){
    st_list[num_int]=tc2[0];
    en_list[num_int]=tc2[1];
    den_list[num_int]=den;
    pri_list[num_int]=i;
    mat_list[num_int]=1;  // needs to be changed
    num_int=num_int+1;
  }
  return num_int;
}

int quadratic_intersect(double *a0,double *alpha,int pars11,double *tc2,int obj)

{
  int out;
  double A,B,C,tmp,*Q,k;
    // The quadratic formula is of the form    A t^2 + B t + C = 0
  Q=&datavec[4][obj*18];
  k=datavec[3][obj];
  C=quadratic_form(a0,Q,a0)-k;
  B=2.0*quadratic_form(alpha,Q,a0);
  A=quadratic_form(alpha,Q,alpha);
  if ((B*B)>(4.0*A*C)) {                    // If sign of determinant is positive there is an intersection with the
                                            //  non-clipped quadratic
    if(A>=0.0){
      tmp=sqrt(B*B-4.0*A*C);
      tc2[0]=(-B-tmp)/(2.0*A);
      tc2[1]=(-B+tmp)/(2.0*A);
      out=2;
    } else {                       // For cones and hyperboloids the segment bounded by the roots may be the part
                                   //  outside the object instead of the part inside the object (i.e., when A<0)
      tmp=sqrt(B*B-4.0*A*C);
      tc2[0]=-1e50;
      tc2[1]=(-B+tmp)/(2.0*A);
      tc2[2]=(-B-tmp)/(2.0*A);
      tc2[3]=1e50;
      out=4;
    }
  }
  else 
    if((pars11==5)&&(C<0)){
      tc2[0]=-1e50;                // For hyperboloid of one sheet the entire ray may be inside the quadratic surface
      tc2[1]=1e50;
      out=2;
    }
    else out=0;
  return out;
}


void SetXData(int length, double *data,int whichVec,int numObj) {
  int i;
  if(whichVec==0){ 
    if(Xdata!=NULL){
      for(i=0;i<numObj;i++)
	if(Xdata[i]!=NULL) {free(Xdata[i]);Xdata[i]=NULL;}
      free(Xdata);
      Xdata=NULL;
    }
    if(XdataL!=NULL) {free(XdataL);XdataL=NULL;}
    Xdata=malloc(sizeof(double *)*numObj);
    XdataL=malloc(sizeof(int)*numObj);
    for(i=0;i<numObj;i++)
      Xdata[i]=NULL;
  }
  if (Xdata[whichVec]!=NULL) {free(Xdata[whichVec]);Xdata[whichVec]=NULL;}
  Xdata[whichVec] = malloc(sizeof(double)*length*3);
  memcpy(Xdata[whichVec],data,length*sizeof(double)*3);
  XdataL[whichVec] = length;
}

void SetDataVec(int length, double *data,int whichVec) {
  if (datavec[whichVec]!=NULL) {free(datavec[whichVec]);datavec[whichVec]=NULL;}
  datavec[whichVec] = malloc(sizeof(double)*length);
  memcpy(datavec[whichVec],data,length*sizeof(double));
  datavecL[whichVec] = length;
}

void PrintDataVec(int whichVec) {
  int i;
  if (datavec[whichVec]!=NULL)
    for(i=0;i<datavecL[whichVec];i++) printf("data: %1.12lf\r\n",datavec[whichVec][i]);
}

void SetIDataVec(int length, int *data,int whichVec) {
  if (idatavec[whichVec]!=NULL) {free(idatavec[whichVec]);idatavec[whichVec]=NULL;}
  idatavec[whichVec] = malloc(sizeof(int)*length);
  memcpy(idatavec[whichVec],data,length*sizeof(int));
  idatavecL[whichVec] = length;
}

void intersections_full_list(int *objlist,int lol,double *a,double *alpha,double maxD,double *t_ends,int *matls,double *dens,int *num_segs)

{
  int i,n,p11,j,pri_list[lol*2],mat_list[lol*2],num_int=0,out,*ibuff,sp,ep,lm,*ien,*ist,ibuffs,p,ip,nn,ii,seg=0;
  double tc2[4],a0[3],st_list[lol*2],en_list[lol*2],den_list[lol*2],lc,ld,*el,*sl;  //[lol*2] used to be [50]

  for(n=0;n<lol;n++){
    i=objlist[n];         // Loop over objects of interest
    p11=idatavec[1][i];
    for(j=0;j<3;j++) a0[j]=(a[j]-datavec[5][i*3+j]);              // vector pointing from obj. center to src
    for(j=0;j<4;j++) tc2[j]=0.0;
    if((p11!=3)&&(p11!=7))                  // If not a torus or vessel segment
      out=quadratic_intersect(a0,alpha,p11,tc2,i);
    else                                   // The torus ca$e requires solving a quartic
      out=quartic_intersect(a0,alpha,tc2,i);
    if(out)        // At least part of the ray is inside the quadratic/quartic surface
      num_int=clip_all(a,alpha,maxD,tc2,out,st_list,en_list,den_list,pri_list,mat_list,num_int,i);
  }                            // Loop over objects is done


  if(num_int){
    ien=malloc(sizeof(double)*num_int);
    ist=malloc(sizeof(double)*num_int);
    el=malloc(sizeof(double)*num_int);
    sl=malloc(sizeof(double)*(num_int+1));
    ibuff=malloc(sizeof(int)*num_int);
    for(i=0;i<num_int;i++) ibuff[i]=0;
    if (num_int>1){
      for(i=0;i<num_int;i++) {ien[i]=i;ist[i]=i;}
      vector=en_list;
      qsort(ien,num_int,sizeof(int),comp);
      vector=st_list;
      qsort(ist,num_int,sizeof(int),comp);
      for(i=0;i<num_int;i++) 
	{
	  el[i]=en_list[ien[i]];sl[i]=st_list[ist[i]];
	}
    }
    else {
      sl[0]=st_list[0];
      el[0]=en_list[0];
      ien[0]=0;
      ist[0]=0;
    }
    sl[num_int]=el[num_int-1]+1;
    sp=0;ep=0;lc=0;ld=0;lm=1;ibuffs=0;                 // sp points to the current interval start
                                                               // ep points to the current interval End
							       // lc : last change (transition from one object to another)
							       // ld : last density
							       // lm : last material
    while(ep<num_int){
      if (sl[sp]<el[ep]){
	p=ist[sp];
	// push into ibuff
	ip=0;
	while(ibuffs>ip){
	  if(pri_list[p]>pri_list[ibuff[ip]]) 
	    break;
	  ip++;
	}
	for(i=ibuffs;i>ip;i--) ibuff[i]=ibuff[i-1];
	ibuff[ip]=p;
	ibuffs++;
	if(ip==0) {
	  //	  mb[lm-1]=mb[lm-1]+ld*(st_list[p]-lc);
	  t_ends[seg]=st_list[p];
	  matls[seg]=lm;
	  dens[seg]=ld;
	  seg++;
	  lc=st_list[p];
	  ld=den_list[p];
	  lm=mat_list[p];
	}
	sp++;
      }
      else {
	p=ien[ep];
	// pull out of ibuff
	ip=0;
	while(ibuff[ip]!=p) ip++;
	for(i=ip;i<(ibuffs-1);i++) ibuff[i]=ibuff[i+1];
	ibuffs--;
	if(ip==0){
	  //	  mb[lm-1]=mb[lm-1]+ld*(en_list[p]-lc);
	  t_ends[seg]=en_list[p];
	  matls[seg]=lm;
	  dens[seg]=ld;
	  seg++;
	  lc=en_list[p];
	  if(ibuffs) {
	    ld=den_list[ibuff[0]];
	    lm=mat_list[ibuff[0]];
	  }
	  else
	    ld=0;  // lm shouldn't matter when ld=0 so we leave it alone
	}
	ep++;
      }
    }
    free(ien);
    free(ist);
    free(el);
    free(sl);
    free(ibuff);
  }
  num_segs[0]=seg;
}


// gcc -shared -fPIC -O3 -mfpmath=sse -msse -o make_volume.so make_volume.c


