# Copyright 2024, GE Precision HealthCare. All rights reserved. See https://github.com/xcist/main/tree/master/license

import re
import sys
from ctypes import *
from numpy.ctypeslib import ndpointer
from gecatsim.pyfiles.GetMu import GetMu
from gecatsim.pyfiles.CommonTools import *
    
def Phantom_Analytic(cfg):
    
    print("Starting to read ANALYTIC phantom...")

    phobj, phobject, numObjects, T, cumCP, NumCP, materialIndex, X, K, Q, Eta, S, D = Phantom_Analytic_Get(cfg)

    cfg.phantom.numberOfMaterials = len(phobject['materialList'])
    set_materials(cfg, phobject['materialList'])
    set_volume(cfg, numObjects, T, cumCP, NumCP, materialIndex, X, K, Q, Eta, S, D)

    print('... done reading phantom.')
    return cfg


def set_materials(cfg, materialList):
    Evec = cfg.spec.Evec
    nMat = len(materialList)
    Mus = np.zeros([Evec.size, nMat], dtype=np.float64)
    for i in range(nMat):
        Mus[:, i] = GetMu(materialList[i], Evec)/10

    # analytic_projector.c: void set_material_info(int materialCount, int eBinCount, double *muTable)
    fun = cfg.clib.set_material_info
    fun.argtypes = [c_int, c_int, ndpointer(c_double)]
    fun.restype = None
    fun(nMat, Evec.size, Mus)


def set_volume(cfg, numObjects, T, cumCP, NumCP, materialIndex, X, K, Q, Eta, S, D):
    X = X.T.ravel()
    K = K.T.ravel()
    Eta = Eta.T.ravel()

    # analytic_projector.c: void set_phantom_info(int numObjs, int *objType, int *clipStInd, int *nPlanes, int *matInd, double *objCent, double *shp, double *Qmat, double *clipNormVec, double *clipDist, double *den, int totalNumPlanes)
    # in matlabe: calllib(CatSimLib, 'set_phantom_info', numObjects, T, cumCP, NumCP, materialIndex, X, K, Q, Eta, S, D, cumCP_end);
    fun = cfg.clib.set_phantom_info
    fun.argtypes = [c_int, ndpointer(c_int), ndpointer(c_int), ndpointer(c_int), ndpointer(c_int), ndpointer(c_double), ndpointer(c_double), ndpointer(c_double), ndpointer(c_double), ndpointer(c_double), ndpointer(c_double), c_int]
    fun.restype = None
    fun(numObjects, T, cumCP, NumCP, materialIndex, X, K, Q, Eta, S, D, cumCP[-1])
    
'''
NOTE: phobject here is object in matlab [simply the parsed phantom], phantobject here is phobject in matlab
'''
def Phantom_Analytic_Get(cfg):
    ###----------- phantom file
    cfg.phantom.filename = my_path.find('phantom', cfg.phantom.filename, '')

    PhantomFilename=cfg.phantom.filename
    Scale=cfg.phantom.scale
    print('Reading phantom file {}...'.format(PhantomFilename))

    PhantomPath = os.path.dirname(PhantomFilename)
    BaseName = os.path.basename(PhantomFilename)
    PhantomBaseName, PhantomExtension = BaseName.split('.')
    if 'pp' == PhantomExtension:
        ppmPhantomFilename = PhantomFilename.replace('.pp','.ppm')
        if cfg.force_phantom_conversion or os.path.getmtime(PhantomFilename) > os.path.getmtime(ppmPhantomFilename):
            if os.path.exist(PhantomFilename):
                pass
                # Phantom_Analytic_pp_to_ppm(PhantomFilename.split('.')[0])
            else:
                sys.exit('Phantom file {} not found'.format(PhantomFilename))
    else:
        if 'ppm' == PhantomExtension:
            ppmPhantomFilename = PhantomFilename
    
    print('Reading phantom file {}.'.format(ppmPhantomFilename))
   
    phobject = parse_analytical_ppm(ppmPhantomFilename)
  
    #FIXME: will fix this part
    #if isfield(cfg,'phantom_rotation_y') and any(cfg.phantom_rotation_y) != 0:
    #    cfg.phantom_rotation_y_analytic = copy(cfg.phantom_rotation_y)
    #
    #if isfield(cfg,'phantom_rotation_y_analytic') and cfg.phantom_rotation_y_analytic != 0:
    #    ang=cfg.phantom_rotation_y_analytic
    #    sina=sin(dot(- ang / 180,pi))
    #    cosa=cos(dot(- ang / 180,pi))
    #    transmat=concat([[cosa,0,sina],[0,1,0],[- sina,0,cosa]])
    #    for ii in arange(1,length(object.type)).reshape(-1):
    #        if any(object.euler_angs(ii,arange(2,3)) != 0) or all(object.euler_angs(ii,1) != concat([0,90,- 90])):
    #            error('CatSim: phantom rotation currently only works without original rotation')
    #        object.center[ii,arange()]=dot(object.center(ii,arange()),transmat)
    #        object.euler_angs[ii,arange()]=concat([90,ang,90])
    #        if logical_not(isempty(object.clip[ii])):
    #            for jj in arange(1,size(object.clip[ii],1)).reshape(-1):
    #                if all(abs(object.clip[ii](jj,arange(1,3))) == concat([1,0,0])):
    #                    object.clip[ii][jj,arange(1,3)]=dot(object.clip[ii](jj,arange(1,3)),transmat)
    #
    # phantom position offset, Mingye
    if hasattr(cfg.phantom,'centerOffset') and not np.allclose(cfg.phantom.centerOffset, 0):
        setattr(cfg.phantom, "centerOffset_analytic", cfg.phantom.centerOffset)
    
    if hasattr(cfg.phantom, 'centerOffset_analytic') and not np.allclose(cfg.phantom.centerOffset_analytic, 0):
        for ii in range(len(phobject['type'])):
            phobject['center'][ii] = np.array(phobject['center'][ii])
            phobject['center'][ii] += cfg.phantom.centerOffset_analytic
            if len(phobject['clip'][ii]) > 0:
                # Mingye Wu, Sept 15 2017
                # get the delta D of clip distance
                for jj in range(len(phobject['clip'][ii])):
                    d = get_clip_dD(phobject['clip'][ii][jj][:3], cfg.phantom.centerOffset_analytic)
                    phobject['clip'][ii][jj][3] += d

    if len(phobject['materialList']) < max(phobject['material']):
        sys.exit('Phantom cannot reference a material that is not in the list')
    objs = {'params': [[] for _ in range(len(phobject['type']))],
            'clip': [[] for _ in range(len(phobject['type']))]}
    if np.any(np.array(phobject['type'])==100):
        #FIXME: finish this part
        #Phantom_Polygonal(object,Scale)
        pass
    else:
        print('Converting boxes to cylinders with clipping.')
        for i in range(len(phobject['type'])):
            if phobject['type'][i] == 8:
                phobject['type'][i] = 2
                Bt = Rmat(phobject['euler_angs'][i]).T
                phobject['clip'][i] = np.vstack((
                    np.hstack((Bt[:,0].T, phobject['half_axes'][i][0] + np.array(phobject['center'][i])@Bt[:,0])),
                    np.hstack((-Bt[:,0].T, phobject['half_axes'][i][0] - np.array(phobject['center'][i])@Bt[:,0])),
                    np.hstack((Bt[:,1].T, phobject['half_axes'][i][1] + np.array(phobject['center'][i])@Bt[:,1])),
                    np.hstack((-Bt[:,1].T, phobject['half_axes'][i][1] - np.array(phobject['center'][i])@Bt[:,1]))))
                phobject['half_axes'][i][:2]= np.array(phobject['half_axes'][i][:2])*np.sqrt(2)

        print('Scaling and compacting format.')
        indsToScale= list(range(6)) + [13,14]

        for i in range(len(phobject['type'])):
            if not phobject['type'][i] in [3,7]:
                tmp = list(phobject['center'][i]) + phobject['half_axes'][i] + phobject['euler_angs'][i] + [phobject['density'][i]] + [phobject['type'][i]] + [0] + [phobject['shape'][i]] + phobject['axial_lims'][i] + [phobject['material'][i]]
                tmp = np.array(tmp)
                tmp[indsToScale] = tmp[indsToScale]*Scale
                objs['params'][i] = tmp
                tmp2 = np.array(phobject['clip'][i])
                if len(tmp2)!=0:
                    tmp2[:,3]= tmp2[:,3]*Scale
                objs['clip'][i] = tmp2
            else:
                objs['params'][i] = np.vstack([phobject['center'][i],phobject['half_axes'][i],phobject['euler_angs'][i],phobject['density'][i],phobject['type'][i],0,phobject['shape'][i],phobject['axial_lims'][i],phobject['material'][i]])
                objs['clip'][i] = phobject['clip'][i]
        # Set the objects, and pass to C.
        objs['params'] = np.array(objs['params'])
        phantObject, numObjects, T, cumCP, NumCP, materialIndex, X, K, Q, Eta, S, D = Phantom_Analytic_SetObjects(objs)
        objs['params'] = Phantom_Analytic_ConvertTori(phantObject,objs['params'],objs['clip'])
        phantObject = Phantom_Analytic_BoundObjects(phantObject,objs['params'])
        print('... done with phantom.')
    
    return phantObject, phobject, numObjects, T, cumCP, NumCP, materialIndex, X, K, Q, Eta, S, D
    
    
def Rmat(p=None):
    p=np.array(p)*3.141592653589793 / 180
    p1=p[0]
    p2=p[1]
    p3=p[2]
    R=np.vstack(([np.cos(p1)*np.cos(p3) - np.sin(p1)*np.cos(p2)*np.sin(p3), np.cos(p1)*np.sin(p3) + np.sin(p1)*np.cos(p2)*np.cos(p3), np.sin(p1)*np.sin(p2)],
            [-np.sin(p1)*np.cos(p3) - np.cos(p1)*np.cos(p2)*np.sin(p3), -np.sin(p1)*np.sin(p3) + np.cos(p1)*np.cos(p2)*np.cos(p3),np.cos(p1)*np.sin(p2)],
            [np.sin(p2)*np.sin(p3), -np.sin(p2)*np.cos(p3), np.cos(p2)]))
    return R
    
# Mingye Wu, Sept 15 2017
# get the delta D of clip
# c0: original clip vector
# c1: phantom off-centering coordinates
def get_clip_dD(c0=None,c1=None):
    
    c0 = np.array(c0)
    c1 = np.array(c1)
    r0 = np.linalg.norm(c0)
    r1 = np.linalg.norm(c1)
    d = np.linalg.norm(c0 - c1)
    cosA = (r0**2 + r1**2 - d**2)/(2*r0*r1)
    cosA = np.nan_to_num(cosA)
    #cosA[cosA==np.inf] = 0.
    delta = r1*cosA
    
    return delta
    
def Phantom_Analytic_BoundObjects(phantObject=None,params=None):
    
    print('Bounding objects in the ANALYTIC phantom.')
    speedup = 0
    
    Good = [4,4,4,4,5,6,7,9,9,10,12,12,14,14,17,17,17,20,20,20,21,22,23,24,25,26,28,28,32,32,32,32,33,34,36,36,38,38,40,40,41,44,44,44,49,49,49,49,49,50,52,52,55,55,55,58,58,58,60,60,62,62,66,66,66,66,68,68,71,71,71,75,75,75,75,77,77,79,79,86,86,86,86,86,86,86,92,92,92,92,92,92,96,96,96,96,100,100,100,100]
    numObjects = len(params)

    for i in range(numObjects):
        pars = params[i]
        typ = pars[10]
        if typ == 1:
            c = 2.5
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
            nV = np.stack((np.concatenate((tmp1, tmp2 + nodes), axis=0),
                           np.concatenate((tmp2, tmp1 + nodes), axis=0),
                           np.concatenate((tmp3, tmp3 - nodes), axis=0)), axis=0)
            Verts = Vx.shape[1]
            tmp_mat = np.array(([pars[3],0,0],[0,pars[4],0],[0,0,pars[5]]))
            Vx = np.matmul(phantObject[i]['R'].T, np.matmul(tmp_mat,Vx)) + pars[:3, None]*np.ones((1,Verts))
        else: print("new types need attention:", typ)
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

        if speedup:
            for j in range(phantObject[i]['cp']):
                Vx,nV=clip_polygon(Vx,nV,phantObject[i]['eta'][:,j],phantObject[i]['s'][j])
        phantObject[i]['Vx'] = Vx
    
    Vx_all = np.hstack([this['Vx'] for this in phantObject])
    nm1 = [this['Vx'].shape[1] for this in phantObject]
    
    cum = np.cumsum(nm1)
    cum = np.concatenate(([0], cum))

    clib = load_C_lib()
    func = clib.set_bounding_info
    func.argtypes = [c_int, ndpointer(c_int), ndpointer(c_double), c_int]
    func.restype = None
    cum = cum.astype(np.int32)
    Vx_all = Vx_all.T.ravel()

    func(numObjects, cum, Vx_all, cum[-1])

    return phantObject, cum, Vx_all
    
#FIXME: not tested
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
    
def Phantom_Analytic_ConvertTori(phantObject=None,params=None,clip=None):
    pi = 3.141592653589793

    print('Converting tori in the ANALYTIC phantom.')
    for i in range(params.shape[0]):
        pars = params[i]
        clp = clip[i]
        if pars[10] == 3:
            if clp.shape[0] == 2:
                ax = phantObject[i]['R'][2].T
                dot = clp[:,:3]*ax
                center = pars[:3]
                dot2 = clp[:,:3]*center - clp[:,3]
                if max(abs(dot) + abs(dot2)) < 5e-05:
                    if np.std(pars[4:6]) == 0:
                        clip_ax_vec[:,1] = np.cross(ax,clp[1,:3]).T
                        clip_ax_vec[:,2] = np.cross(ax,clp[2,:3]).T
                        clip_ax_vec = phantObject[i]['R']*clip_ax_vec
                        ang = np.arctan2(clip_ax_vec[2], clip_ax_vec[1])
                        if (ang(2) - ang(1))%(2*pi) > pi:
                            ang[2] += pi
                        else:
                            ang[1] += pi
                            ang = [[ang[2]],[ang[1]]]
                        ang = ang%(2*pi)
                        l = (ang(2) - ang(1))%(2*pi)
                        relative_cent_ang = ang(1) + l / 2 - pi / 4
                        C = pars(6)
                        A = 1 / C ** 2
                        shp = pars(13)
                        pars[4] = A*(1 - shp ** 2)
                        pars[5] = 0
                        pars[6] = A
                        pars[13] = C
                        pars[14] = l
                        pars[15] = relative_cent_ang
                        pars[11] = 7
                        params[i] = pars
    
    # 
	# 	Info on parameter storage:
	# 	--------------------------
	# 
	# 	vessel seg.    par #   original objects
	# 	------------------------------------------
	# 	center		1:3   - object center
	# 	k,d,A		4:6   - object half axes
	# 	Euler angles	7:9   - Euler angles
	# 	intensity	10    - intensity
	# 	7		11    - obj. type  (1=ell,2=cyl,3=tor,4=cone,5=hyp1,6=hyp2,7=vessel seg.)
	# 	trans		12    - transparency (vis)
	# 	C              13    - torus shape
	# 	l,T		14,15 - ax. lim.
	# 	min OS		16    - minimum detector oversampling
	# 
	# 	vessel definition
	# 	-----------------
	# 
	# 	start_point				(x,y,z)
	# 	start_direction			(theta,phi)
	# 	start_curvature_direction		(gamma)
	# 	start_radius				(R)
	# 
	# 	R_scaling_factor			(vector [Num_segs])  default:1
	# 	length				(vector [Num_segs])  default:pi/4
	# 	ellipticity				(vector [Num_segs])  default:0
	# 	curvature				(vector [Num_segs])  default:20
	# 	torsion				(vector [Num_segs])  default:0
	# 	tapering nature (rel_cent_ang)	(vector [Num_segs])  default:0
	# 
	# 
	# 	Program reads in and converts to:
	# 
	# 	vessel seg store
	# 	----------------
	# 
	# 	center point
	# 	clipping planes (array[2 by 4])
	# 	k
	# 	ellipticity (d)
	# 	A
	# 	curvature (C)
	# 	length (l)   ( 0 < l < pi )
	# 	tapering nature (T)
	# 	Ql,Qr,R
	# 
	# 	Also can convert circular torus segments to the above form if there are exactly two clipping planes and they
	# 	both contain the major axis of the torus (C=radius of curvature;A=1/C^2; s=1-shape^2; k=A*s;d=0;T=0;)
   
    return params

def parse_analytical_ppm(ppmPhantomFilename):
    #NOTE: consider change obj to a class
    obj = dict.fromkeys(['materialList', 'center', 'half_axes', 'euler_angs', 'density', 'type', 'material', 'axial_lims', 'shape', 'clip'])
    with open(ppmPhantomFilename, 'r') as f:
        all_lines = f.readlines()

    key_indice = []
    for line in all_lines:
        if "materialList" in line:
            material_list = line.strip().rstrip(";").split('=')[-1].lstrip('{').rstrip('}').split()
            material_list = [tmp.strip().lstrip('{').rstrip('}').strip("'") for tmp in material_list]
            obj['materialList'] = material_list
        elif len(line) < 4: continue
        else:
            key_str, value_str = line.rstrip(";").split('=')
            key_idx = int(''.join(filter(str.isdigit, key_str)))
            key_indice.append(key_idx)

    key_indice = set(key_indice)
    for key in obj.keys():
        if key != 'materialList': obj[key] = [[] for _ in range(len(key_indice))]

    for line in all_lines:
        if len(line) < 4: continue
        elif "materialList" not in line:
            key_str, value_str = line.rstrip(";").split('=')
            key_word = key_str.split('.')[-1].split('(')[0].split('{')[0]
            key_word = str(key_word)
            key_idx = int(''.join(filter(str.isdigit, key_str)))

            if '[' in value_str:
                # extract substring in []
                value_str = re.findall(r'\[([^]]*)\]', value_str)[0]
                split_valstr = value_str.rstrip(';').split(';')
                if len(split_valstr) == 1:
                    this_str = split_valstr[0].strip().rstrip(";").split(' ')
                    if len(this_str)==1 and this_str[0] == '':
                        value = []
                    else:
                        value = [eval(x) for x in this_str]
                else:
                    value = []
                    for this_splvalstr in split_valstr:
                        _tmp = this_splvalstr.split(' ')
                        #breakpoint()
                        this_value = [eval(x) for x in _tmp]
                        value.append(this_value)
            else:
                this_str = value_str.strip().rstrip(";")
                value = eval(this_str)
            if 'clip' in key_word and len(value)>0 and len(np.array(value).shape)==1: value=[value]

            obj[key_word][key_idx-1] = value

    return obj

def Phantom_Analytic_SetObjects(objs=None, debug=False):
    numObjects = objs['params'].shape[0]
    print('Setting {} objects in the ANALYTIC phantom.'.format(numObjects))
    phantObject = []
    for i in range(numObjects):
        tmp_dict = dict.fromkeys(('eta', 's', 'k', 'Q', 'xo', 'cp', 'matIndex', 'R'))
        phantObject.append(tmp_dict)

        pars = objs['params'][i]
        cp = len(objs['clip'][i])
        if cp > 0:
            phantObject[i]['eta'] = objs['clip'][i][:,:3].T
            phantObject[i]['s'] = objs['clip'][i][:,3].T
        else:
            #TODO: consider change to []
            #phantObject[i]['eta'] = []
            #phantObject[i]['s'] = []
            phantObject[i]['eta'] = None
            phantObject[i]['s'] = None
        p1=pars[6]*np.pi / 180
        p2=pars[7]*np.pi / 180
        p3=pars[8]*np.pi / 180
        R = np.array([[np.cos(p1)*np.cos(p3) - np.sin(p1)*np.cos(p2)*np.sin(p3), np.cos(p1)*np.sin(p3) + np.sin(p1)*np.cos(p2)*np.cos(p3), np.sin(p1)*np.sin(p2)],
            [-np.sin(p1)*np.cos(p3) - np.cos(p1)*np.cos(p2)*np.sin(p3), -np.sin(p1)*np.sin(p3) + np.cos(p1)*np.cos(p2)*np.cos(p3), np.cos(p1)*np.sin(p2)],
            [np.sin(p2)*np.sin(p3), - np.sin(p2)*np.cos(p3), np.cos(p2)]])

        if (pars[10] != 3) and (pars[10] != 7):
            if pars[10] == 1:
                phantObject[i]['k'] = 1.
                d=np.array(pars[3:6]) ** 2
            elif pars[10] == 2:
                phantObject[i]['k'] = 1.
                d = np.concatenate((pars[3:5], [np.inf])) ** 2

                if phantObject[i]['eta'] is None:
                    phantObject[i]['eta'] = np.vstack((R[2], -R[2])).T
                else:
                    phantObject[i]['eta'] = np.hstack((phantObject[i]['eta'], R[2][:,None], -R[2][:,None])) 
                if phantObject[i]['s'] is None:
                    phantObject[i]['s'] = np.vstack((pars[5] + pars[:3]@R[2].T, pars[5] - pars[:3]@R[2].T))
                else:
                    phantObject[i]['s'] = np.hstack((phantObject[i]['s'], pars[5] + pars[:3]@R[2].T, pars[5] - pars[:3]@R[2].T))
            else:
                print("Error! New materials, need to fix", pars[10])
                sys.exit(1)
            #elif pars(11) == 4:
            #    phantObject[i].k = copy(0)
            #    d=concat([pars(arange(4,5)) ** 2,- pars(6) ** 2])
            #    phantObject[i].eta = copy(concat([phantObject[i].eta,R(3,arange()).T,- R(3,arange()).T]))
            #    phantObject[i].s = copy(concat([phantObject[i].s,pars(15) + dot(pars(arange(1,3)),R(3,arange()).T),- pars(14) - dot(pars(arange(1,3)),R(3,arange()).T)]))
            #elif pars(11) == 5:
            #    phantObject[i].k = copy(1)
            #    d=concat([pars(arange(4,5)) ** 2,- pars(6) ** 2])
            #    phantObject[i].eta = copy(concat([phantObject[i].eta,R(3,arange()).T,- R(3,arange()).T]))
            #    phantObject[i].s = copy(concat([phantObject[i].s,pars(15) + dot(pars(arange(1,3)),R(3,arange()).T),- pars(14) - dot(pars(arange(1,3)),R(3,arange()).T)]))
            #elif pars(11) == 6:
            #    phantObject[i].k = copy(- 1)
            #    d=concat([pars(arange(4,5)) ** 2,- pars(6) ** 2])
            #    phantObject[i].eta = copy(concat([phantObject[i].eta,R(3,arange()).T,- R(3,arange()).T]))
            #    phantObject[i].s = copy(concat([phantObject[i].s,pars(15) + dot(pars(arange(1,3)),R(3,arange()).T),- pars(14) - dot(pars(arange(1,3)),R(3,arange()).T)]))
            S = np.zeros((3,3), dtype=np.float64)
            np.fill_diagonal(S, 1.0 / d) #in place ops
            phantObject[i]['Q'] = R.T@S@R
        elif pars[10] == 3:
            Sl=np.zeros((3,3))
            Sr=np.copy(Sl)
            parst = pars[3:6]
            np.fill_diagonal(S, 1.0 / parst**2)
            Sr[0,0]=1.0 / parst[0] ** 2
            Sr[1,1]=1.0 / parst[1] ** 2
            phantObject[i]['Q'] = np.vstack((R.T*Sl*R, R.T*Sr*R))
            phantObject[i]['k'] = 1. - pars[12] ** 2
        #elif pars[11] == 7:
        #    #    error('attempt to use incomplete code')
        #    C=pars(13)
        #    delta=pars(5)
        #    A=pars(6)
        #    k=pars(4)
        #    s=k / A
        #    shape=sqrt(1 - s)
        #    D_true=dot(1 / C ** 2,(1 + delta / 2))
        #    E_true=dot(1 / C ** 2,(1 - delta / 2))
        #    Sl=zeros(3)
        #    Sr=copy(Sl)
        #    Sl[concat([1,5,9])]=dot(concat([1,1,1]),A)
        #    Sr[concat([1,5])]=concat([D_true,E_true])
        #    phantObject[i].Q = copy(concat([dot(dot(R.T,Sl),R),dot(dot(R.T,Sr),R)]))
        #    phantObject[i].k = copy(s)
        else:
            sys.exit('Unknown object type.')

        phantObject[i]['xo'] = pars[:3].T
        if phantObject[i]['s'] is None:
            phantObject[i]['cp'] = 0
        else:
            phantObject[i]['cp'] = len(phantObject[i]['s'])
        phantObject[i]['typ'] = pars[10]
        phantObject[i]['matIndex'] = pars[15]
        phantObject[i]['R'] = R
    
    # Write some of the data in phantObject to the memory of C
    Q=[]
    qt=np.zeros((3,6), dtype=np.float64)
    Eta= np.hstack([x['eta'] for x in phantObject if x['eta'] is not None])
    S = []
    for x in phantObject:
        if x['s'] is not None:
            S += list(x['s'].ravel())
    NumCP = []
    for x in phantObject:
        if x['s'] is not None:
            NumCP.append(len(x['s']))
        else:
            NumCP.append(0)
    K = np.hstack([x['k'] for x in phantObject  if x['k'] is not None])
    X = np.vstack([x['xo'] for x in phantObject if x['xo'] is not None])
    C = np.hstack([x['cp'] for x in phantObject if x['cp'] is not None])
    T = np.hstack([int(x['typ']) for x in phantObject if x['typ'] is not None])
    D = np.hstack([x[9] for x in objs['params'] if x[9] is not None])
    materialIndex = np.hstack([x['matIndex'] for x in phantObject if x['matIndex'] is not None])
    for i in range(numObjects):
        if phantObject[i]['Q'].shape[1] == 3:
            Q.append(np.vstack((phantObject[i]['Q'], phantObject[i]['Q'])))
        else:
            Q.append(phantObject[i]['Q'])

    Q = np.vstack(Q).T
    
    materialCount = len(materialIndex)
    cumCP = np.cumsum(NumCP)
    cumCP= np.insert(cumCP, 0, 0)
    if debug:
        print("Q: sum={}, std={}\n".format(np.sum(Q), np.std(Q)))
        breakpoint()
    T = T.astype(np.int32)
    S = np.array(S, dtype=np.float64)
    X = X.T
    cumCP = cumCP.astype(np.int32)
    NumCP = np.array(NumCP).astype(np.int32)
    materialIndex = materialIndex.astype(np.int32)
    clib = load_C_lib()
    func = clib.set_phantom_info
    func.argtypes = [c_int, ndpointer(c_int), ndpointer(c_int), ndpointer(c_int), ndpointer(c_int), ndpointer(c_double), ndpointer(c_double), ndpointer(c_double), ndpointer(c_double), ndpointer(c_double), ndpointer(c_double), c_int]
    func.restype = None
    func(numObjects, T, cumCP, NumCP, materialIndex, X, K, Q, Eta, S, D, cumCP[-1])

    return phantObject, numObjects, T, cumCP, NumCP, materialIndex, X, K, Q, Eta, S, D

# NOTE:not verified yet
#def Phantom_Analytic_pp_to_ppm(PhantomFileBasename, Scale=1., debug=False):
#
#    print('Converting phantom file {}.pp to .ppm file...'.format(PhantomFileBasename))
#
#    ppPhantomFilename = PhantomFileBasename + '.pp'
#    tmpPhantomFilename = PhantomFileBasename + '.tmp'
#
#    #TODO: need to change this to python wrapper
#    feval('C_Phantom_Analytic_FORBILD_to_tmp',Scale,ppPhantomFilename,tmpPhantomFilename)
#
#    print('Reading and deleting {}.'.format(tmpPhantomFilename))
#
#    with open(tmpPhantomFilename, 'r') as f:
#        lns = f.readlines()
#
#    os.remove(tmpPhantomFilename)
#
#    ppmPhantomFilename = PhantomFileBasename + '.ppm'
#    print('Writing {}.'.format(ppmPhantomFilename))
#
#    f = open(ppmPhantomFilename, 'w')
#    ind = [i for i in range(len(lns)) if '35' in lns[i, 1]]
#    if len(ind) < 2:
#        sys.exit('material table not found')
#    else:
#        f.write('materialList = {')
#        for i in range(1,ind(2) - 2):
#            f.write('\'%s\' ', lns[i + 1])
#        f.write('};\n\n')
#        offset = ind[2]
#    
#    i=0
#    while (offset < size(lns,1) - 15):
#        i += 1
#        tmpo,tmpc,lns,offset,clip = read_object(lns, offset)
#        # TODO: first save string then flush all at once
#        f.write('object.center({%d},:) = [{%f} {%f} {%f}];\n'.format(i,tmpo.cent(1),tmpo.cent(2),tmpo.cent(3)));
#        f.write('object.half_axes({%d},:) = [{%f} {%f} {%f}];\n'.format(i,tmpo.hax(1),tmpo.hax(2),tmpo.hax(3)));
#        f.write('object.euler_angs({%d},:) = [{%f} {%f} {%f}];\n'.format(i,tmpo.ea(1),tmpo.ea(2),tmpo.ea(3)));
#        f.write('object.density({%d}) = {%f};\\n'.format(i,tmpo.density))
#        f.write('object.type({%d}) = {%d};\\n'.format(i,tmpo.type))
#        f.write('object.material({%d}) = {%d};\\n'.format(i,tmpo.material))
#        f.write('object.axial_lims({%d},:) = [0 0];\\n'.format(i))
#        f.write('object.shape(%d) = 0;\\n',i)
#        f.write('object.clip{%d} = ['.format(i))
#        for j in range(1,tmpc.shape[0]):
#            f.write('{%f} {%f} {%f} {%f};'.format(tmpc(j,1),tmpc(j,2),tmpc(j,3),tmpc(j,4)))
#        f.write('];\\n\\n')
#        clip[i]=tmpc
#        obj[i]=tmpo
#    
#    f.close()
#    print('... done writing {}'.format(ppmPhantomFilename))
    
#def read_object(lines=None,offset=None):
#
#    obj.type = lines[1 + offset]
#    for i in range(4):
#        A[i] = float(lines[i + 4 + offset])
#    
#    if np.linalg.norm(A - [0,0,0,1]) > 1e-05:
#        sys.exit('Unexpected transform')
#    
#    obj.cent = np.copy(A[:3,3].T)
#    B = A[0:3, 0:3]
#    obj.hax = np.sqrt(np.sum(B*B))
#    B /= np.ones(3,1)*obj.hax
#    B = B.T
#    obj.ea = euler_angs(B)*180/pi
#    obj.transform = copy(A)
#    tmp=sscanf(lines(15 + offset,arange()),'%s',1)
#    tmp = lines[15+offset]
#    #NOTE: we can step matlab code and compare python code
#    if tmp(1) == 'C':
#        i=1
#        while tmp(1) == 'C':
#            v = list(map(float, lines[16+offset, 21:]))
#            clip[i,0:3] = v
#            clip[i,4] = float(lines[16 + offset, 26:])
#            clip[i,:] = -clip(i)
#            offset += 3
#            tmp = lines[15+offset]
#            i += 1
#    else:
#        clip=[]
#    
#    if np.sum((obj.type[1:3]) == 'Cyl') == 3:
#        obj.type = 2
#        obj.hax[3]=obj.hax(3) / 2
#    else:
#        if sum(obj.type[:3] == 'Cub') == 3:
#            obj.type = 2
#            old=0
#            if old:
#                Bt=copy(B)
#            else:
#                Bt=B.T
#
#            clip = clip + [Bt[:,1].T, obj.hax(1)/2 + obj.cent*Bt[:,1]]
#            clip = clip + [-Bt[:,1].T, obj.hax(1)/2 - obj.cent*Bt[:,1]]
#            clip = clip + [Bt[:,2].T, obj.hax(2)/2 + obj.cent*Bt[:,2]]
#            clip = clip + [-Bt[:,2].T, obj.hax(2)/2 - obj.cent*Bt[:,2]]
#            obj.hax[:2] = obj.hax[:2]*np.sqrt(2)
#            obj.hax /= 2
#        else:
#            if np.sum((obj.type[:3]) == 'Sph') == 3:
#                obj.type = 1
#            else:
#                if np.sum(obj.type[:3] == 'Ell') == 3:
#                    obj_type = 1
#                else:
#                    sys.exit('Unrecognized obj type\\n\\r')
#    
#    
#    obj.material = lines[18 + offset]
#    obj.density = lines[19 + offset]
#    offset += 21
#
#    return obj, clip, lines, offset
    
def euler_angs(R=None):

    # R =
	# [  cos(p1)*cos(p3)-sin(p1)*cos(p2)*sin(p3),  cos(p1)*sin(p3)+sin(p1)*cos(p2)*cos(p3), sin(p1)*sin(p2)]
	# [ -sin(p1)*cos(p3)-cos(p1)*cos(p2)*sin(p3), -sin(p1)*sin(p3)+cos(p1)*cos(p2)*cos(p3), cos(p1)*sin(p2)]
	# [                          sin(p2)*sin(p3),                         -sin(p2)*cos(p3),         cos(p2)]
    
    p2 = np.real(np.acos(R(3,3)))
    
    if np.abs(np.sin(p2)) > 1e-06:
        p1 = atan2(R(1,3) / np.sin(p2),R(2,3) / np.sin(p2))
        p3=atan2(R(3,1) / sin(p2),- R(3,2) / sin(p2))
    else:
        if sign1(cos(p2)) > 0:
    # R =
	# [  cos(p1+p3), -R(2,1), 0]
	# [ -sin(p1+p3),  R(1,1), 0]
	# [           0,       0, 1]
            p2=0
            p3=0
            p1=atan2(R(1,2),R(1,1))
        else:
    # R =
	# [  cos(p1-p3),  R(2,1), 0]
	# [ -sin(p1-p3), -R(1,1), 0]
	# [           0,       0,-1]
            p2=copy(pi)
            p3=0
            p1=atan2(- R(1,2),R(1,1))
    
    Rnew=concat([[dot(cos(p1),cos(p3)) - dot(dot(sin(p1),cos(p2)),sin(p3)),dot(cos(p1),sin(p3)) + dot(dot(sin(p1),cos(p2)),cos(p3)),dot(sin(p1),sin(p2))],[dot(- sin(p1),cos(p3)) - dot(dot(cos(p1),cos(p2)),sin(p3)),dot(- sin(p1),sin(p3)) + dot(dot(cos(p1),cos(p2)),cos(p3)),dot(cos(p1),sin(p2))],[dot(sin(p2),sin(p3)),dot(- sin(p2),cos(p3)),cos(p2)]])
# Phantom_Analytic_pp_to_ppm.m:193
    err = np.linalg.norm(Rnew - R,[],1)
    if err > 1e-05:
        sys.exit('Error computing Euler angles')
    
    ea = [p1,p2,p3]
    
    return ea

def Phantom_Polygonal_ReadPolygon(Verts):

    ddir = my_path.find_dir("phantom", "poly_bin")
    filename = 'poly{}'.format(Verts)
    with open(os.path.join(ddir, filename),'rb') as fid:
        data_array = np.fromfile(fid, dtype=np.int32, count=4)
        nv_sz1, nv_sz2, vx_sz1, vx_sz2 = data_array
        tmp = np.fromfile(fid, dtype=np.float64)
        nV = tmp[:nv_sz1*nv_sz2].reshape(nv_sz2,nv_sz1).T
        Vx = tmp[nv_sz1*nv_sz2:].reshape(vx_sz2,vx_sz1).T

    return Vx,nV
