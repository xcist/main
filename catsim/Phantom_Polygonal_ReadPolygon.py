#% -----------------------------------------------------------------------
#%   Program Name: Phantom_Polygonal_ReadPolygon.m                                     
#%   Author:  Samit Basu, Bruno De Man, Jed Pack (GE Global Research)
#%   Organization: GE Global Research and GE Healthcare, General Electric Company
#%   Version:  6.0.3
#%   Date:  Feb 3, 2015
#%   Class: GE Confidential. General Electric Proprietary Data (c) 2012 General Electric Company
#%
#% History: 
#%   2012-10-12 Paul FitzGerald (GE Global Research)
#%              Renamed - was load_polygon.m.
#%              Changed bin_directory to BinDirectory.
#% -----------------------------------------------------------------------
import os
import numpy as np

'''
log:
    verified by comparing std, sum of Vx, nV
'''

def Phantom_Polygonal_ReadPolygon(Verts, debug=False):

    BinDirectory = "/home/dp058565/work/catsim/base/bin/"
    ddir = BinDirectory
    filename = 'poly{}'.format(Verts)
    with open(os.path.join(ddir, filename),'rb') as fid:
        data_array = np.fromfile(fid, dtype=np.int32, count=4)
        #breakpoint()
        nv_sz1, nv_sz2, vx_sz1, vx_sz2 = data_array
        #breakpoint()
        #nV = np.fromfile(fid, dtype=np.int64, count=nv_sz1*nv_sz2, offset=4*4)
        #Vx = np.fromfile(fid, dtype=np.float64, count=vx_sz1*vx_sz2, offset=(4*4+8*nv_sz1*nv_sz1))
        tmp = np.fromfile(fid, dtype=np.float64)
        nV = tmp[:nv_sz1*nv_sz2].reshape(nv_sz2,nv_sz1).T
        #nV = nV.astype(int)
        Vx = tmp[nv_sz1*nv_sz2:].reshape(vx_sz2,vx_sz1).T
        #nV = tmp[nv_sz1-1:nv_sz2]
        #Vx = tmp[vx_sz1-1:vx_sz2]
        #nV = fread(fid,[nv_sz1 nv_sz2],'double';
        #Vx = fread(fid,[vx_sz1 vx_sz2],'double');
    
    if debug: breakpoint()

    return Vx,nV

if __name__=="__main__":
    test = Phantom_Polygonal_ReadPolygon(20, debug=True)
