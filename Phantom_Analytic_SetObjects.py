# Generated with SMOP  0.41
# Phantom_Analytic_SetObjects.m

# -----------------------------------------------------------------------
#   Program Name: Phantom_Analytic_SetObjects.m                                      
#   Authors:  Samit Basu, Bruno De Man, Jed Pack (GE Global Research)
#   Organization: GE Global Research and GE Healthcare, General Electric Company
#   Version:  6.0.3
#   Date:  Feb 3, 2015
#   Class: GE Confidential. General Electric Proprietary Data (c) 2012 General Electric Company
    
#   This function set the objects in an analytic phantom and passes these to C.
    
# History:
#   2012-10-12 Paul FitzGerald (GE Global Research)
#              Changed call to set_phantom_info to call wrapper rather than calling library function.
#              Renamed - was YYreexpress.m.
#              Added "Verbose" output.
# -----------------------------------------------------------------------


'''
log:
    20220323: verified by comparing std, sum for all variables used in feval
'''

'''
input: objs, 
output: diction phantObject
'''
import numpy as np
from catsim.Phantom_Polygonal_ReadPolygon import Phantom_Polygonal_ReadPolygon
from catsim.CommonTools import *
from ctypes import *
from numpy.ctypeslib import ndpointer

# NOTE: phantObject should be a list of dictionaries
def Phantom_Analytic_SetObjects(objs=None, debug=False):
    numObjects = objs['params'].shape[0]
    print('Setting {} objects in the ANALYTIC phantom.'.format(numObjects))
    # create an tmpty phantObject
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
                #d = np.insert(pars[3:5], , np.inf) ** 2 this will not work, -1 will be in the middle

                if phantObject[i]['eta'] is None:
                    phantObject[i]['eta'] = np.vstack((R[2], -R[2])).T
                else:
				    #matlab code phantObject{i}.eta=[phantObject{i}.eta R(3,:)' -R(3,:)'];
                    #breakpoint()
                    phantObject[i]['eta'] = np.hstack((phantObject[i]['eta'], R[2][:,None], -R[2][:,None])) 
                    #print("ERROR! see comment in source code")
                if phantObject[i]['s'] is None:
				    #phantObject{i}.s=[phantObject{i}.s pars(15)+pars(1:3)*R(3,:)' -pars(14)-pars(1:3)*R(3,:)'];
                    phantObject[i]['s'] = np.vstack((pars[5] + pars[:3]@R[2].T, pars[5] - pars[:3]@R[2].T))
                else:
                    phantObject[i]['s'] = np.hstack((phantObject[i]['s'], pars[5] + pars[:3]@R[2].T, pars[5] - pars[:3]@R[2].T))
                    #print("ERROR! see comment in source code")
                    #breakpoint()
            else: print("new caterogies", pars[10])
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
            #if i==9: breakpoint()
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
        #breakpoint()
        if phantObject[i]['s'] is None:
            phantObject[i]['cp'] = 0
        else:
            phantObject[i]['cp'] = len(phantObject[i]['s'])
        #phantObject[i]['cp'] = np.size(phantObject[i]['s'])
        phantObject[i]['typ'] = pars[10]
        phantObject[i]['matIndex'] = pars[15]
        #print(pars[15], phantObject[i]['matIndex'])
        phantObject[i]['R'] = R
    
    # Write some of the data in phantObject to the memory of C
    Q=[]
    qt=np.zeros((3,6), dtype=np.float64)
    #breakpoint()
    Eta= np.hstack([x['eta'] for x in phantObject if x['eta'] is not None])
    #breakpoint()
    #S = np.hstack([x['s'] for x in phantObject if x['s'] is not None])
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
    #NumCP = np.hstack([np.size(x['s']) for x in phantObject])
    K = np.hstack([x['k'] for x in phantObject  if x['k'] is not None])
    X = np.vstack([x['xo'] for x in phantObject if x['xo'] is not None])
    C = np.hstack([x['cp'] for x in phantObject if x['cp'] is not None])
    T = np.hstack([int(x['typ']) for x in phantObject if x['typ'] is not None])
    D = np.hstack([x[9] for x in objs['params'] if x[9] is not None])
    materialIndex = np.hstack([x['matIndex'] for x in phantObject if x['matIndex'] is not None])
    #Q = np.hstack([x['Q']) for x in phantObject if x['Q'].shape[1]==3])
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
    #feval('C_Phantom_Analytic_Set', numObjects, T, cumCP, NumCP, materialIndex, X, K, Q, Eta, S, D, cumCP)
    T = T.astype(np.int32)
    S = np.array(S, dtype=np.float64)
    X = X.T
    cumCP = cumCP.astype(np.int32)
    NumCP = np.array(NumCP).astype(np.int32)
    materialIndex = materialIndex.astype(np.int32)
    #K = K.astype(np.int32)
    clib = load_C_lib()
    func = clib.set_phantom_info
    func.argtypes = [c_int, ndpointer(c_int), ndpointer(c_int), ndpointer(c_int), ndpointer(c_int), ndpointer(c_double), ndpointer(c_double), ndpointer(c_double), ndpointer(c_double), ndpointer(c_double), ndpointer(c_double), c_int]
    func.restype = None
    #breakpoint()
    func(numObjects, T, cumCP, NumCP, materialIndex, X, K, Q, Eta, S, D, cumCP[-1])
    # we need to return T to cumCP back to Phantom_analytic so that we can set volumne correctly
    return phantObject, numObjects, T, cumCP, NumCP, materialIndex, X, K, Q, Eta, S, D
    
if __name__ == '__main__':
    from scipy import io as sio
    mat_in = sio.loadmat("Phantom_Analytic_SetObjects_objs_in.mat")['objs']
    #breakpoint()
    objs = {}
    objs['params'] = mat_in[0][0][0]
    objs['clip'] = mat_in[0][0][1][0]

    #phantObject = [] # DONNOT write to phant = [{}]*len
    #for i in range( len(objs['params'])):
    #    phantObject.append({})
    #breakpoint()

    mat_out = sio.loadmat("Phantom_Analytic_SetObjects_out.mat")
    #phantObject_mat = mat_out['phantObject']
    numObjects_mat = mat_out['numObjects']
    T_mat = mat_out['T']
    cumCP_mat = mat_out['cumCP']
    NumCP_mat = mat_out['NumCP']
    materialIndex_mat = mat_out['materialIndex']
    X_mat = mat_out['X']
    K_mat = mat_out['K']
    Q_mat = mat_out['Q']
    Eta_mat = mat_out['Eta']
    S_mat = mat_out['S']
    D_mat = mat_out['D']

    phantObject, numObjects, T, cumCP, NumCP, materialIndex, X, K, Q, Eta, S, D = Phantom_Analytic_SetObjects(objs, debug=False)
    #breakpoint()
    try:
        assert(np.allclose(numObjects, numObjects_mat, atol=1.E-8))
        assert(np.allclose(T, T_mat, atol=1.E-8))
        assert(np.allclose(cumCP, cumCP_mat, atol=1.E-8))
        assert(np.allclose(NumCP, NumCP_mat, atol=1.E-8))
        assert(np.allclose(materialIndex, materialIndex_mat, atol=1.E-8))
        assert np.allclose(X, X_mat, atol=1.E-8), 'X'
        assert np.allclose(K, K_mat, atol=1.E-8), 'K'
        assert np.allclose(Q, Q_mat, atol=1.E-8), "Q"
        assert(np.allclose(Eta, Eta_mat, atol=1.E-8))
        assert(np.allclose(S, S_mat, atol=1.E-8))
        assert(np.allclose(D, D_mat, atol=1.E-8))
    except AssertionError as err:
        print(err)
        breakpoint()
