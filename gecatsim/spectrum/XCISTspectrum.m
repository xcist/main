function [newE,newI] = XCISTspectrum(kVp, angle, dE)

  % ------------------------------------------------------------------------------
  % Created: B. De Man, August 27, 2020, 2020 - GE Research, Niskayuna, NY 12065, USA
  % Last modified: B. De Man, August 27, 2020, 2020 - GE Research, Niskayuna, NY 12065, USA
  % modified by Jiayong Zhang at June 2023: added the correct scaling factor for spectrum
  %
  % Copyright 2024, GE Precision HealthCare
  % Code Developed Using US Government Funding via National Cancer Institute Grant 1U01CA231860-01A1 (ITCR)
  % All rights reserved. See https://github.com/xcist/main/tree/master/license
  %
  % Purpose
  %    Function to generate normalized, unfiltered X-ray CT spectrum for tungsten anode for tube voltages between
  %    80 kVp and 140 kVp and for target angles between 5 and 9 degrees.
  %    Spectra were initialized based using SpekPY (physics-based tool) and iteratively updated based on empirical
  %    CT data and gradient-ascent minimization of the attenuation least squares error (post-log).
  %    The spectra are compactly parametrized using a 17 parameter vector. The average RMS
  %    attenuation errors were 0.026, 0.011, 0.007 and 0.006 at 80, 100, 120 and 140 kVp respectively.
  %    A detailed description of the spectrum estimation methodology is reported in :
  %        "Spectrum estimation for X-ray computed tomography"
  %        Paul FitzGerald, Stephen Araujo, Mingye Wu, Bruno De Man
  %        Submitted to Medical Physics
  %
  % Inputs
  %    kVp: X-ray tube voltage
  %    angle: X-ray tube target angle (or takeoff angle of the x-ray beam relative to the target surface)
  %    dE: desire spectrum sampling density (in keV)
  %
  % Outputs
  %    newE: array of all keV values at which the spectrum is defined
  %    newI: array of all normalized photon count values of the spectrum
  %       (note that I scales inversely proportional with dE to preserve total counts)
  %
  % Example usage
  %    [E,I]=XCISTspectrum(120);
  %    figure; plot(E,I); grid on;
  %
  % Limitations
  %    Only valid for Tungsten targets
  %    Input parameter space is limited to 80-140 kVp and 5-9 degrees
  %    Does not include any inherent tube filtration.
  % ------------------------------------------------------------------------------

  %-------------------------------------------------------------------
  % validate inputs
  %-------------------------------------------------------------------
  if kVp < 80 | kVp > 140
    error('Tube voltage has to be between 80 and 140 kVp !')
  end
  if nargin < 2
    angle=7.0;
  end
  if angle < 5 | angle > 9
    error('Angle has to be between 5 and 9 degrees !')
  end
  if nargin < 3
    dE=1.0;
  end
  if dE < 0.5 | dE > 20
    error('Sampling density angle has to be between 0.5 and 20 keV !')
  end

  %-------------------------------------------------------------------
  % spectrum parameters
  %-------------------------------------------------------------------
  allparams=zeros(19,3,4);

  allparams(:,:,1)=...
  [      0         0         0
         0         0         0
	 9.2280    9.0751    8.4031
	 62.1915   64.5380   69.3968
	 109.3572  110.8911  116.9230
	 120.2025  119.3448  121.0521
	 84.8994   83.3843   80.0781
	 47.5492   46.5216   43.1936
	 31.1946   31.4719   30.0455
	 20.2388   19.9671   19.4418
	 -20.2388  -19.9671  -19.4418
	 -60.7163  -59.9014  -58.3255
	 -121.4327 -119.8027 -116.6511
	 -202.3878 -199.6712 -194.4185
         0         0         0
	 38.4669   41.5628   42.7531
	 67.8298   73.3348   74.9963
	 18.4094   20.5120   19.9288
         0         0         0];

  allparams(:,:,2)=...
  [ 0         0         0
    1.4528    1.4060    2.1868
    16.0157   16.4726   21.1581
    42.3437   43.5510   49.8557
    63.3678   64.3889   68.6747
    73.3139   73.2596   73.9680
    65.9847   64.7651   63.1847
    50.6813   49.0172   47.0623
    34.3827   35.3012   35.0264
    29.7647   30.1220   29.5544
    19.1145   18.9863   18.1651
    7.3835    7.2983    6.7053
    -14.7670  -14.5966  -13.4106
    -44.3011  -43.7897  -40.2317
    0         0         0
    109.8538  117.4195  122.2282
    197.4694  210.8765  219.0544
    73.0435   77.4888   79.7090
    15.6381   16.5770   17.0522];

  allparams(:,:,3)=...
  [0         0         0
   0.8564    0.8994    1.4183
   9.8912   10.5994   14.0769
   27.7893   29.1783   34.4463
   44.0825   45.1210   49.2789
   54.6745   54.4745   55.8673
   53.4040   52.0372   51.2659
   44.9039   43.1173   41.7155
   28.7341   30.0178   30.4752
   26.6774   27.4320   27.5315
   21.3651   21.4246   21.0232
   15.5955   15.4494   14.7377
   6.8805    6.8382    6.0445
   -6.8805   -6.8382   -6.0445
   0         0         0
   142.3841  150.0002  156.9327
   256.7763  270.3153  282.1586
   97.3396  101.9380  105.2753
   22.1465   23.1369   23.7903];

  allparams(:,:,4)=...
  [0         0         0
   0.6769    0.7496    1.1367
   7.5893    8.3309   11.1569
   21.3809   22.7056   27.4418
   34.4374   35.3683   39.6081
   43.8696   43.5788   45.6035
   44.5953   43.1870   43.0696
   39.3997   37.5810   36.5738
   23.4315   24.8254   25.6625
   22.6010   23.5537   24.0224
   19.8653   20.1196   20.0579
   16.4030   16.3043   15.9226
   11.0143   10.8408   10.1404
   4.0105    4.0080    3.3360
   0         0         0
   157.5624  163.9847  172.8835
   284.6359  296.0568  311.1753
   109.3829  113.3104  117.3900
   25.3486   26.1938   26.9541];

  %-------------------------------------------------------------------
  % Interpolate in angle
  %-------------------------------------------------------------------
  angleparams=zeros(19,4);
  for i=1:19
    for j=1:4
      pp=polyfit([5,7,9],[allparams(i,:,j)],2);
      angleparams(i,j)=(polyval(pp,angle));
    end
  end
  % interpolate between angle and kVp
  group_ang = 8.73:-0.495:5.26;
  rescale = 0.01.*[ %0.01 is due to the spectrum unit used in testing was cm^2
  11681 11135   10812   11722   11672   10199   10152   10263
  22428 21704   21307   22676   22581   19954   19581   19226
  36252 35615   34373   36867   36385   32570   31924   31214
  51002 49395   48234   51116   50536   45081   44094   41995];
  % http://www.ece.northwestern.edu/local-apps/matlabhelp/techdoc/ref/interp2.html
  allkVp = 80:20:140;
  allangle = group_ang;
  valid_idx = [1 2 3 6 7 8]; % note the center two 8-row groups are not used because abnormal photons/mA
  thisrescale = interp2(allangle(valid_idx), allkVp, rescale(:, valid_idx), angle, kVp, 'spline');
  
  %-------------------------------------------------------------------
  % Interpolate in tube voltage
  %-------------------------------------------------------------------
  params=zeros(19,1);
  for i=1:19
    params(i)=max(spline([80,100,120,140],angleparams(i,:),kVp),0);
  end

  %-------------------------------------------------------------------
  % Define spline left and right of K-edge and characteristic peaks
  %-------------------------------------------------------------------
  % Define spline left of K-edge
  ee=[15,20,26,32,38,46,56,66]; % spline knot energies below K-edge
  ii=params(1:8);

  % Define spline right of K-edge
  ee2=[70.25,75,85,95,110,130]; % spline knot energies above K-edge
  index=find(ee2<(kVp-10));
  if ~isempty(index)
    ee2=[ee2(index),kVp];
  else
    ee2=[69.75,kVp];
  end
  ii2=[params(9:9+length(ee2)-2);0];

  % Define characteristic peaks
  eek=[57.982,59.318,67.244,69.100]; % 4 characteristic peaks (drop the fifth one at 69.517)
  iik=params(16:19);

  %-------------------------------------------------------------------
  % Generate new spectrum with sampling distance dE
  %-------------------------------------------------------------------
  % Interpolate from knots left of K-edge
  newE=dE:dE:kVp; %define the new energy coordinates
  indexl=find(newE<=69.5); %indices where we want to interpolate from left size of K-edge
  xx=newE(indexl);
  yy=spline([ee],[ii],xx)*dE/0.5; %re-scale to correct absolute count depending on keV sampling density
  yy(find(xx<ee(1)))=0; % set to zero left of first knot

  % Eliminate oscillations below 30 keV
  therefix=find(xx<30);
  f1=0.98^(dE/0.5);
  f2=0.8^(dE/0.5);
  for i=1:length(therefix)
    i2=therefix(end-i+1);
    if yy(i2)>f1*yy(i2+1)
      yy(i2)=f2*yy(i2+1);
    end
  end

  % Interpolate right side of K-edge
  indexr=find(newE>69.5); %indices where we want to interpolate from right size of K-edge
  xx2=newE(indexr); %re-scale to correct absolute count depending on keV sampling density
  yy2=spline(ee2,ii2,xx2)*dE/0.5;
  yy2(find(xx2>kVp))=0; % set to zero right of last knot

  % Stack left and right sides together
  newI=max([yy,yy2],0);
  

  % Add counts for characteristic peak #1
  index=find(newE>eek(1));indexh=index(1); indexl=indexh-1; %find coords that straddle the peak
  dxx=newE(indexh)-newE(indexl);
  newI(indexl)=newI(indexl)+iik(1)*(newE(indexh)-eek(1))/dxx; %spread counts across two nearest energies
  newI(indexh)=newI(indexh)+iik(1)*(eek(1)-newE(indexl))/dxx;

  % Add counts for characteristic peak #2
  index=find(newE>eek(2));indexh=index(1); indexl=indexh-1;
  dxx=newE(indexh)-newE(indexl);
  newI(indexl)=newI(indexl)+iik(2)*(newE(indexh)-eek(2))/dxx;
  newI(indexh)=newI(indexh)+iik(2)*(eek(2)-newE(indexl))/dxx;

  % Add counts for characteristic peak #3
  index=find(newE>eek(3));indexh=index(1); indexl=indexh-1;
  dxx=newE(indexh)-newE(indexl);
  newI(indexl)=newI(indexl)+iik(3)*(newE(indexh)-eek(3))/dxx;
  newI(indexh)=newI(indexh)+iik(3)*(eek(3)-newE(indexl))/dxx;

  % Add counts for characteristic peak #4
  index=find(newE>eek(4));indexh=index(1); indexl=indexh-1;
  dxx=newE(indexh)-newE(indexl);
  newI(indexl)=newI(indexl)+iik(4)*(newE(indexh)-eek(4))/dxx;
  newI(indexh)=newI(indexh)+iik(4)*(eek(4)-newE(indexl))/dxx;
  newI = newI.*thisrescale;

  %outvar = [newE, newI];
  outname = strcat('xcist_kVp', num2str(kVp), '_tar', num2str(angle), '_bin', num2str(dE), '.mat');
  save(outname,'newE','newI');
