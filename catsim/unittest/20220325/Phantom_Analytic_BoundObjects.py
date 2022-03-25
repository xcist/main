# Generated with SMOP  0.41
# Phantom_Analytic_BoundObjects.m

    # -----------------------------------------------------------------------
#   Program Name: Phantom_Analytic_BoundObjects.m                                     
#   Authors:  Samit Basu, Bruno De Man, Jed Pack (GE Global Research)
#   Organization: GE Global Research and GE Healthcare, General Electric Company
#   Version:  6.0.3
#   Date:  Feb 3, 2015
#   Class: GE Confidential. General Electric Proprietary Data (c) 2012 General Electric Company
    
    # Aim
#   This function defines a polyhedron that bounds (completely encompasses) objects in an analytic phantom,
#   and passes the boundary to C.
    
    # Inputs
#   fields in phantObject:  eta,s,k,Q,xo,cp,typ
    
    # Outputs
#   phantObject is modified to include:
#        vertices (x,y,z)
#        neighboring vertices (vertex indicies in CCW (counterclockwise)
#        order when looking from the outside of the polyhedron)
#        neighbor start vec
    
    # History:
#   2012-10-12 Paul FitzGerald (GE Global Research)
#              Renamed - was YYpoly_bound.
#              Changed call to set_bounding_info to call wrapper rather than calling library function.
#              Added "Verbose" output.
# -----------------------------------------------------------------------

'''
log:
    20220322: tested by comparing std, sum
'''

import numpy as np
from catsim.Phantom_Polygonal_ReadPolygon import Phantom_Polygonal_ReadPolygon
from catsim.CommonTools import *
from ctypes import *
from numpy.ctypeslib import ndpointer

#NOTE by jiayong: has been verified to be correct as of current version
def Phantom_Analytic_BoundObjects(phantObject=None,params=None, debug=False):
    
    print('Bounding objects in the ANALYTIC phantom.')
    speedup = 0
    
    Good = [4,4,4,4,5,6,7,9,9,10,12,12,14,14,17,17,17,20,20,20,21,22,23,24,25,26,28,28,32,32,32,32,33,34,36,36,38,38,40,40,41,44,44,44,49,49,49,49,49,50,52,52,55,55,55,58,58,58,60,60,62,62,66,66,66,66,68,68,71,71,71,75,75,75,75,77,77,79,79,86,86,86,86,86,86,86,92,92,92,92,92,92,96,96,96,96,100,100,100,100]
    numObjects = len(params)

    if debug: all_verts = []
    for i in range(numObjects):
        pars = params[i]
        typ = pars[10]
        if typ == 1:
            c = 2.5
            # prod calc production per column TODO: need to rewrite
            # For FB, only 1 and 2 are involved
            Volume = np.prod(pars[3:6])
            Verts = Good[min(100, int(np.ceil(Volume**(4/15)*c)))-1]
            Vx, nV = Phantom_Polygonal_ReadPolygon(Verts)
            tmp = np.vstack(([pars[3],0,0],[0,pars[4],0],[0,0,pars[5]]))
            Vx = phantObject[i]['R'].T@tmp@Vx + pars[:3, None]@np.ones((1,Verts))
        elif typ == 2:
            c = 2
            Volume = np.prod(pars[3:6])
            nodes=max(round(c*Volume ** (2 / 15)),3)
            ang = np.linspace(0, 2*np.pi ,nodes + 1)
            ang = ang[:-1]
            sc = 1/np.cos(np.pi / nodes)
            x = np.cos(ang)*sc
            y = np.sin(ang)*sc
            z = x*0
            Vx = np.stack((np.concatenate((x,x), axis=0),
                                np.concatenate((y,y), axis=0),
                                np.concatenate((z-1,z+1), axis=0)), axis=0)
            Vx_old = Vx.copy()
            tmp1 = (np.arange(1,nodes+1) - 2)%nodes + 1
            tmp2 = (np.arange(1,nodes+1))%nodes + 1
            tmp3 = np.arange(1,nodes+1) + nodes
            #nV = np.array([[tmp1.T,tmp2.T + nodes], [tmp2.T,tmp1.T + nodes], [tmp3.T,tmp3.T - nodes]])
            nV = np.stack((np.concatenate((tmp1, tmp2 + nodes), axis=0),
                           np.concatenate((tmp2, tmp1 + nodes), axis=0),
                           np.concatenate((tmp3, tmp3 - nodes), axis=0)), axis=0)
            Verts = Vx.shape[1]
            tmp_mat = np.array(([pars[3],0,0],[0,pars[4],0],[0,0,pars[5]]))
            Vx = np.matmul(phantObject[i]['R'].T, np.matmul(tmp_mat,Vx)) + pars[:3, None]*np.ones((1,Verts))
        else: printf("new types need attention:", typ)
        if debug: all_verts.append(Verts)
        #if i==9: breakpoint()
        #elif typ == 3:
        #    c=2
        #    s=pars(13)
        #    Volume=dot(dot(prod(pars(arange(4,6))),(1 + s) ** 2),s)
        #    nodes=max(round(dot(c,Volume ** (2 / 15))),3)
        #    ang=linspace(0,dot(2,pi),nodes + 1)
        #    ang=ang(arange(1,end() - 1))
        #    sc=dot(1 / cos(pi / nodes),(1 + s))
        #    x=dot(cos(ang),sc)
        #    y=dot(sin(ang),sc)
        #    z=dot(x,0)
        #    Vx=concat([[x,x],[y,y],[z - s,z + s]])

        #    tmp1=mod(concat([arange(1,nodes)]) - 2,nodes) + 1
        #    tmp2=mod(concat([arange(1,nodes)]),nodes) + 1
        #    tmp3=concat([arange(1,nodes)]) + nodes
        #    nV=concat([concat([[tmp1.T],[tmp2.T + nodes]]),concat([[tmp2.T],[tmp1.T + nodes]]),concat([[tmp3.T],[tmp3.T - nodes]])])
        #    Verts=size(Vx,2)
        #    Vx=dot(dot(phantObject[i].R.T,concat([[pars(4),0,0],[0,pars(5),0],[0,0,pars(6)]])),Vx) + dot(pars(arange(1,3)).T,ones(1,Verts))
        #elif typ == 4:
        #    c=2
        #    alim=sort(pars(arange(14,15)))
        #    rat=alim(2) / pars(6)
        #    Volume=dot(dot(dot(pars(4),pars(5)),rat ** 2),alim(2)) / 3
        #    nodes=max(round(dot(c,Volume ** (2 / 15))),3)
        #    ang=linspace(0,dot(2,pi),nodes + 1)
        #    ang=ang(arange(1,end() - 1))
        #    sc=1 / cos(pi / nodes)
        #    x=dot(cos(ang),sc)
        #    y=dot(sin(ang),sc)
        #    z=dot(x,0)
        #    scl=max(alim(1) / alim(2),0.05)
        #    Vx=concat([[dot(x,scl),x],[dot(y,scl),y],[z + alim(1),z + alim(2)]])
        #    tmp1=mod(concat([arange(1,nodes)]) - 2,nodes) + 1
        #    tmp2=mod(concat([arange(1,nodes)]),nodes) + 1
        #    tmp3=concat([arange(1,nodes)]) + nodes
        #    nV=concat([concat([[tmp1.T],[tmp2.T + nodes]]),concat([[tmp2.T],[tmp1.T + nodes]]),concat([[tmp3.T],[tmp3.T - nodes]])])
        #    Verts=size(Vx,2)
        #    Vx=dot(dot(phantObject[i].R.T,concat([[dot(pars(4),rat),0,0],[0,dot(pars(5),rat),0],[0,0,1]])),Vx) + dot(pars(arange(1,3)).T,ones(1,Verts))
        #if typ == 5:
        #    c=2
        #    alim=sort(pars(arange(14,15)))
        #    mag=sqrt(1 + alim ** 2 / pars(6) ** 2)
        #    Volume=dot(dot(dot(pars(4),pars(5)),prod(mag)),(alim(2) - alim(1)))
        #    nodes=max(round(dot(c,Volume ** (2 / 15))),3)
        #    ang=linspace(0,dot(2,pi),nodes + 1)
        #    ang=ang(arange(1,end() - 1))
        #    sc=1 / cos(pi / nodes)
        #    x=dot(cos(ang),sc)
        #    y=dot(sin(ang),sc)
        #    z=dot(x,0)
        #    Vx=concat([[dot(x,mag(1)),dot(x,mag(2))],[dot(y,mag(1)),dot(y,mag(2))],[z + alim(1),z + alim(2)]])
        #    tmp1=mod(concat([arange(1,nodes)]) - 2,nodes) + 1
        #    tmp2=mod(concat([arange(1,nodes)]),nodes) + 1
        #    tmp3=concat([arange(1,nodes)]) + nodes
        #    nV=concat([concat([[tmp1.T],[tmp2.T + nodes]]),concat([[tmp2.T],[tmp1.T + nodes]]),concat([[tmp3.T],[tmp3.T - nodes]])])
        #    Verts=size(Vx,2)
        #    Vx=dot(dot(phantObject[i].R.T,concat([[pars(4),0,0],[0,pars(5),0],[0,0,1]])),Vx) + dot(pars(arange(1,3)).T,ones(1,Verts))
        #if typ == 6:
        #    c=2
        #    alim=sort(pars(arange(14,15)))
        #    lim=max(abs(alim))
        #    if (alim(1) > - pars(6)) and (alim(1) < pars(6)):
        #        alim[1]=min(pars(6),alim(2))
        #    if (alim(2) > - pars(6)) and (alim(2) < pars(6)):
        #        alim[2]=max(- pars(6),alim(1))
        #    mag=sqrt(- 1 + lim ** 2 / pars(6) ** 2)
        #    Volume=dot(dot(dot(pars(4),pars(5)),prod(mag)),(alim(2) - alim(1)))
        #    nodes=max(round(dot(c,Volume ** (2 / 15))),3)
        #    ang=linspace(0,dot(2,pi),nodes + 1)
        #    ang=ang(arange(1,end() - 1))
        #    sc=1 / cos(pi / nodes)
        #    x=dot(cos(ang),sc)
        #    y=dot(sin(ang),sc)
        #    z=dot(x,0)
        #    Vx=concat([[dot(x,mag),dot(x,mag)],[dot(y,mag),dot(y,mag)],[z + alim(1),z + alim(2)]])
        #    tmp1=mod(concat([arange(1,nodes)]) - 2,nodes) + 1
        #    tmp2=mod(concat([arange(1,nodes)]),nodes) + 1
        #    tmp3=concat([arange(1,nodes)]) + nodes
        #    nV=concat([concat([[tmp1.T],[tmp2.T + nodes]]),concat([[tmp2.T],[tmp1.T + nodes]]),concat([[tmp3.T],[tmp3.T - nodes]])])
        #    Verts=size(Vx,2)
        #    Vx=dot(dot(phantObject[i].R.T,concat([[pars(4),0,0],[0,pars(5),0],[0,0,1]])),Vx) + dot(pars(arange(1,3)).T,ones(1,Verts))
        #if typ == 7:
        #    C=pars(13)
        #    A=pars(6)
        #    T=pars(15)
        #    theta=pars(14)
        #    delta=pars(5)
        #    k=pars(4)
        #    if C <= 0:
        #        error('Torus shape value must be > 0.')
        #    #    st_ang=(pi/2-theta)*T;
        #    st_ang=pi / 4 + T - theta / 2
        #    en_ang=st_ang + theta
        #    deltas=dot(delta,cos(dot(2,concat([st_ang,en_ang]))))
        #    if sign1(sin(st_ang)) != sign1(sin(en_ang)):
        #        delta_max=copy(delta)
        #    else:
        #        delta_max=max(deltas)
        #    alphas=dot(2 / C ** 2,(1 + concat([deltas,delta_max]) / 2))
        #    s1=sqrt(alphas - k + sqrt(alphas ** 2 - dot(dot(2,alphas),k)))
        #    s2=sqrt(alphas - k - sqrt(alphas ** 2 - dot(dot(2,alphas),k)))
        #    rmax=((s1(3) - s2(3)) / 2 + sqrt(alphas(3) / 2)) / A
        #    r=(- (s1(3) - s2(3)) / 2 + sqrt(alphas(arange(1,2)) / 2)) / A / rmax
        #    s=(s1(3) - s2(3)) / 2 / A
        #    c=2
        #    shp=(s1(3) - s2(3)) / 2 / sqrt(alphas(3) / 2)
        #    Volume=dot(dot(theta,C ** 3),shp ** 2)
        #    nodes=max(round(dot(c,Volume ** (2 / 15))),3)
        #    ang=linspace(st_ang,en_ang,nodes)
        #    sc=dot(1 / cos(theta / 2 / (nodes - 1)),(rmax))
        #    x=dot(cos(ang),sc)
        #    y=dot(sin(ang),sc)
        #    z=dot(x,0)
        #    if min(r) <= 0.1:
        #        Vx=concat([[x,0,x,0],[y,0,y,0],[z - s,- s,z + s,s]])
        #        nodes=nodes + 1
        #    else:
        #        Vx=concat([[concat([x,dot(r(2),x(end())),dot(r(1),x(1)),x,dot(r(2),x(end())),dot(r(1),x(1))])],[concat([y,dot(r(2),y(end())),dot(r(1),y(1)),y,dot(r(2),y(end())),dot(r(1),y(1))])],[concat([z - s,- s,- s,z + s,s,s])]])
        #        nodes=nodes + 2
        #    tmp1=mod(concat([arange(1,nodes)]) - 2,nodes) + 1
        #    tmp2=mod(concat([arange(1,nodes)]),nodes) + 1
        #    tmp3=concat([arange(1,nodes)]) + nodes
        #    nV=concat([concat([[tmp1.T],[tmp2.T + nodes]]),concat([[tmp2.T],[tmp1.T + nodes]]),concat([[tmp3.T],[tmp3.T - nodes]])])
        #    Verts=size(Vx,2)
        #    Vx=dot(phantObject[i].R.T,Vx) + dot(pars(arange(1,3)).T,ones(1,Verts))

        #speedup is set to 0 for now
        if speedup:
            for j in arange(phantObject[i]['cp']):
                #concat([i,j])
                Vx,nV=clip_polygon(Vx,nV,phantObject[i]['eta'][:,j],phantObject[i]['s'][j])
        phantObject[i]['Vx'] = Vx
    
    Vx_all = np.hstack([this['Vx'] for this in phantObject])
    nm1 = [this['Vx'].shape[1] for this in phantObject]
    
    cum = np.cumsum(nm1)
    #cum = np.insert(cum, 0, 0)
    cum = np.concatenate(([0], cum))
    if debug: breakpoint()
    #TODO: need to finish C_Phantom_Analytic_SetBoundary.py file
    #cum, Vx_all = feval('C_Phantom_Analytic_SetBoundary', numObjects, cum,Vx_all, cum[-1])
    clib = load_C_lib()
    func = clib.set_bounding_info
    func.argtypes = [c_int, ndpointer(c_int), ndpointer(c_double), c_int]
    func.restype = None
    cum = cum.astype(np.int32)
    #cum1 = np.copy(cum)
    #Vx_all1 = np.copy(Vx_all)
    #breakpoint()
    func(numObjects, cum, Vx_all, cum[-1])
    #breakpoint()

    return phantObject, cum, Vx_all
    
#NOTE: not tested
def clip_polygon(vX=None,nV=None,clipVec=None,clipS=None):

    # This function clips the bounding polyhedron with the clipping planes.
	# This helps to reduce the volume of the bounding polyhedron and thereby improves efficiency.
	# In some (rare) cases, this clipping is also necessary to prevent the source from being inside the polyhedron.
	# If the source is inside the bounding polyhedron, no horizon can be computed and the simulator will fail.
    
    numVert = vX.shape[1]
    newVerts = 0
    thresh = 0.01
    distOutside = clipVec.T*vX - clipS
    if max(distOutside) > thresh:
        outside = (distOutside > 0)
        for i in range(1,numVert+1):
            # This won't happen when we are at index 1 so we start with index 2.
            numNeighbors = sum(np.logical_not(np.isinf(nV[i])))
            for j in range(numNeighbors):
                neighbor = nV[i,j]
                if neighbor < i:
                    if outside[i] + outsize[neighbor] == 1:
                        newVerts += 1
                        wt = abs(distOutside(neighbor)) / (abs(distOutside(neighbor)) + abs(distOutside(i)))
                        vX = [vX, wt*vX[:,i] +(1 - wt)*vX[:, neighbor]]
                        if outside[i]:
                            ind = find(nV(neighbor,arange()) == i)
                            nV[neighbor,ind]=numVert + newVerts
                            nV[numVert + newVerts,arange()]=inf
                            nV[numVert + newVerts,1]=neighbor
                        else:
                            nV[i,j]=numVert + newVerts
                            nV[numVert + newVerts,arange()]=inf
                            nV[numVert + newVerts,1]=i
        # Next, we connect the dots on the new vertices  
		# Find a point on the face
        newFaceCenter=mean(vX(arange(),numVert + concat([arange(1,newVerts)])),2)
        vec1=vX(arange(),numVert + 1) - newFaceCenter
        vec1=vec1 / norm_cs(vec1)
        vec2=cross(clipVec,vec1)
        vec2=vec2 / norm_cs(vec2)
        coords=dot(concat([[vec1.T],[vec2]]),(vX(arange(),numVert + concat([arange(1,newVerts)])) - dot(newFaceCenter,ones(1,newVerts))))
        angles=atan2(coords(2,arange()),coords(1,arange()))
        angles,order=sort(angles,nargout=2)
        order=order + numVert
        nV[order,2]=order(concat([arange(2,end()),1]))
        nV[order,3]=order(concat([end(),arange(1,end() - 1)]))
        outside[numVert + concat([arange(1,newVerts)])]=0
        ind=find(outside == 0)
        newIndex[ind]=concat([arange(1,length(ind))])
        nV=nV(ind,arange())
        vX=vX(arange(),ind)
        ii=find(logical_not(isinf(nV)))
        nV[ii]=newIndex(nV(ii))
        if sum(nV == 0):
            sys.exit('Unexpected.')
        # show_poly_stereographic(vX,nV,newFaceCenter,cross(vec1,vec2)',[vec1 vec2'])
        # FIXME : Need to eliminate any very small faces to prevent precision errors
    
    
    return vX,nV
    
    
#def show_poly_stereographic(vX=None,nV=None,newfacecenter=None,newfacevector=None,orthogonalvecs=None,*args,**kwargs):
#    varargin = show_poly_stereographic.varargin
#    nargin = show_poly_stereographic.nargin
#
#    focus=newfacecenter + dot(newfacevector,5)
#    vX2=dot(concat([orthogonalvecs,newfacevector]).T,(vX - repmat(focus,concat([1,size(vX,2)]))))
#    x=- vX2(1,arange()) / vX2(3,arange())
#    y=- vX2(2,arange()) / vX2(3,arange())
#    figure(1)
#    hold('off')
#    for i in arange(1,size(vX,2)).reshape(-1):
#        text(x(i),y(i),sprintf('%d',i))
#        hold('on')
#        nb=nV(i,find(logical_not(isinf(nV(i,arange())))))
#        for j in arange(1,length(nb)).reshape(-1):
#            h=plot(x(concat([i,nb(j)])),y(concat([i,nb(j)])))
#    
#    
#    return
    
if __name__ == '__main__':
    from scipy import io as sio
    mat_in = sio.loadmat("Phantom_Analytic_BoundObject_phantObject_in.mat")['phantObject']
    python_in = []
    for i in range(mat_in.shape[1]):
        python_in.append({})
        python_in[-1]['eta'] = mat_in[0][i][0][0][0]
        python_in[-1]['s'] = mat_in[0][i][0][0][1] #1d array
        python_in[-1]['k'] = mat_in[0][i][0][0][2][0,0] #scalar
        python_in[-1]['Q'] = mat_in[0][i][0][0][3]
        python_in[-1]['xo'] = mat_in[0][i][0][0][4] #1d
        python_in[-1]['cp'] = mat_in[0][i][0][0][5][0,0] #scalar
        python_in[-1]['typ'] = mat_in[0][i][0][0][6][0,0] #scalar
        python_in[-1]['matIndex'] = mat_in[0][i][0][0][7][0,0] #scalar
        python_in[-1]['R'] = mat_in[0][i][0][0][8] #2d

    params_in = sio.loadmat("Phantom_Analytic_BoundObject_params_in.mat")['params']
    #params_in = sio.loadmat("params_in_BoundObjects.mat")['params']
    #params_in = params_in.astype(int)
    #breakpoint()
    phantObjec, cum, Vx_all = Phantom_Analytic_BoundObjects(phantObject=python_in, params=params_in, debug=False)

    mat_out = sio.loadmat("Phantom_Analytic_BoundObjects_out.mat")
    #phantObject_mat = mat_out['phantObject']
    cum_mat = mat_out['cum']
    Vx_all_mat = mat_out['Vx_all']

    #breakpoint()
    try:
        assert np.allclose(cum, cum_mat, atol=1.E-8), 'cum'
        assert np.allclose(Vx_all, Vx_all_mat, atol=1.E-8), 'Vx_all'
    except AssertionError as err:
        print(err)
        breakpoint()
