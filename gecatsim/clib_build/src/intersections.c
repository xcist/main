// Copyright 2020, General Electric Company. All rights reserved. See https://github.com/xcist/code/blob/master/LICENSE

void intersections(int *objlist,int lol,double *a,double *alpha,double maxD,double *mb)

{
  int i,n,p11,j,pri_list[lol*2],mat_list[lol*2],num_int=0,out,*ibuff,sp,ep,lm,*ien,*ist,ibuffs,p,ip,nn,ii;
  double tc2[4],a0[3],st_list[lol*2],en_list[lol*2],den_list[lol*2],lc,ld,*el,*sl;  //[lol*2] used to be [50]

  //  printf("objlist: %d lol: %d a=[%lf %lf %lf] alpha=[%lf %lf %lf] maxD: %lf mb[0]: %lf\n",objlist[0],lol,a[0],a[1],a[2],alpha[0],alpha[1],alpha[2],maxD,mb[0]);
  for(n=0;n<lol;n++){
    i=objlist[n];         // Loop over objects of interest
    //    if((i==4)&&(prnt)) pnt=1; else pnt=0;
    p11=idatavec[1][i];
    for(j=0;j<3;j++) a0[j]=(a[j]-datavec[5][i*3+j]);              // vector pointing from obj. center to src
    for(j=0;j<4;j++) tc2[j]=0.0;
    if((p11!=3)&&(p11!=7))                  // If not a torus or vessel segment
      out=quadratic_intersect(a0,alpha,p11,tc2,i);
    else                                   // The torus ca$e requires solving a quartic
      out=quartic_intersect(a0,alpha,tc2,i);
    //if((pnt)&&(tc2[0]>100.0))
    //  printf("Intersection points for object 5, module 181:\r\n   %1.12lf   %1.12lf\r\n\r\n",tc2[0],tc2[1]);
    //if (out!=0) {for(ii=0;ii<out;ii++) printf("%1.12lf  ",tc2[ii]); printf("\r\n");}
    if(out)        // At least part of the ray is inside the quadratic/quartic surface
      num_int=clip_all(a,alpha,maxD,tc2,out,st_list,en_list,den_list,pri_list,mat_list,num_int,i);
  }                            // Loop over objects is done
  //printf("%d\r\n",num_int);

  if(num_int){
    ien=malloc(sizeof(double)*num_int);
    ist=malloc(sizeof(double)*num_int);
    el=malloc(sizeof(double)*num_int);
    sl=malloc(sizeof(double)*(num_int+1));
    ibuff=malloc(sizeof(int)*num_int);
    for(i=0;i<num_int;i++) ibuff[i]=0;
  //printf("a%d\r\n",num_int);
    if (num_int>1){
      for(i=0;i<num_int;i++) {ien[i]=i;ist[i]=i;}
      vector=en_list;//for(i=0;i<num_int;i++) printf("en_vector: %1.12f\r\n",vector[i]);
      qsort(ien,num_int,sizeof(int),comp);
      vector=st_list;//for(i=0;i<num_int;i++) printf("st_vector: %1.12f\r\n",vector[i]);
      qsort(ist,num_int,sizeof(int),comp);
      for(i=0;i<num_int;i++) 
	{
	  el[i]=en_list[ien[i]];sl[i]=st_list[ist[i]];
	  //printf("sl %1.12lf  el %1.12lf  ist: %d  ien: %d\r\n",sl[i],el[i],ist[i],ien[i]);
	  //if (i>0) if(sl[i]<sl[i-1]) printf("ERROR!! %1.12lf\r\n",a0[10000]);
	}
      //printf("b%d\r\n",num_int);
    }
    else {
      sl[0]=st_list[0];
      el[0]=en_list[0];
      ien[0]=0;
      ist[0]=0;
    }
    //printf("c%d\r\n",num_int);
    sl[num_int]=el[num_int-1]+1;
    sp=0;ep=0;lc=0;ld=0;lm=1;ibuffs=0;                 // sp points to the current interval start
                                                               // ep points to the current interval End
							       // lc : last change (transition from one object to another)
							       // ld : last density
							       // lm : last material
    //printf("d%d\r\n",num_int);
    while(ep<num_int){
      //printf("E  %1.12lf  %1.12lf  %1.12lf  %d  %d\r\n",sl[0],el[0],sl[1],sp,ep);
      if (sl[sp]<el[ep]){
	//printf("f%d\r\n",num_int);
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
	//for(i=0;i<ibuffs;i++) printf("Ibuff: %d\n\r",ibuff[i]);
	if(ip==0) {
	  mb[lm-1]=mb[lm-1]+ld*(st_list[p]-lc);
	  lc=st_list[p];
	  ld=den_list[p];
	  lm=mat_list[p];
	}
	sp++;
      }
      else {
	//printf("g%d\r\n",num_int);
	//printf("g%d%d\r\n",num_int,ibuff[0]);
	p=ien[ep];
	// pull out of ibuff
	ip=0;
  //printf("h%d %d\r\n",p,ibuffs);
  //for(i=0;i<ibuffs;i++) printf("ibuff/w: %d\n\r",ibuff[i]);
	while(ibuff[ip]!=p) ip++;
  //printf("i%d\r\n",num_int);
	for(i=ip;i<(ibuffs-1);i++) ibuff[i]=ibuff[i+1];
  //printf("j%d\r\n",num_int);
	ibuffs--;
	//for(i=0;i<ibuffs;i++) printf("ibuff/wo: %d\n\r",ibuff[i]);
  //printf("k%d\r\n",num_int);
	if(ip==0){
  //printf("l   %1.12lf %d %1.12lf %1.12lf %d %1.12lf\r\n",mb[0],lm,ld,en_list[0],p,lc);
	  mb[lm-1]=mb[lm-1]+ld*(en_list[p]-lc);
	  //	  printf("ld: %1.12lf diff: %1.12lf\r\n",ld,en_list[p]-lc);
  //printf("m%d\r\n",num_int);
	  lc=en_list[p];
  //printf("n%d\r\n",num_int);
	  if(ibuffs) {
  //printf("o%d\r\n",num_int);
	    ld=den_list[ibuff[0]];
  //printf("p%d\r\n",num_int);
	    lm=mat_list[ibuff[0]];
  //printf("q%d\r\n",num_int);
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
}


