# Generated with SMOP  0.41
# Phantom_Analytic_pp_to_ppm.m

    # -----------------------------------------------------------------------
#   Program Name: Phantom_Analytic_pp_to_ppm.m                                        
#   Author:  Samit Basu, Bruno De Man, Jed Pack (GE Global Research)
#   Organization: GE Global Research and GE Healthcare, General Electric Company
#   Version:  6.0.3
#   Date:  Feb 3, 2015
#   Class: GE Confidential. General Electric Proprietary Data (c) 2012 General Electric Company
    
    # Aim
#   This function converts a .pp phantom file to a .ppm file.
    
    # Inputs:
#   PhantomFileBasename  The name of the phantom file, with no extension.
#                        A file named [PhantomFileBasename '.pp'] must exist.
#   Scale                (optional) A scale factor, to resize phantom, also used to convert from cm to mm (scale=10).
#                        Default: 1
# Outputs:
#   A new file is written to [PhantomFileBasename '.ppm']
    
    # History: 
#   2012-10-12 Paul FitzGerald (GE Global Research)
#              Renamed - was translate.m. 
#              Cleaned up and added "Verbose" output. 
# -----------------------------------------------------------------------

'''
log:
    not tested yet
'''
from catsim.CommonTools import feval

# not verified yet
def Phantom_Analytic_pp_to_ppm(PhantomFileBasename, Scale=1., debug=False):

    print('Converting phantom file {}.pp to .ppm file...'.format(PhantomFileBasename))

    ppPhantomFilename = PhantomFileBasename + '.pp'
    tmpPhantomFilename = PhantomFileBasename + '.tmp'

    #TODO: need to change this to python wrapper
    feval('C_Phantom_Analytic_FORBILD_to_tmp',Scale,ppPhantomFilename,tmpPhantomFilename)

    print('Reading and deleting {}.'.format(tmpPhantomFilename))

    with open(tmpPhantomFilename, 'r') as f:
        lns = f.readlines()

    os.remove(tmpPhantomFilename)

    ppmPhantomFilename = PhantomFileBasename + '.ppm'
    print('Writing {}.'.format(ppmPhantomFilename))

    f = open(ppmPhantomFilename, 'w')
    ind = [i for i in range(len(lns)) if '35' in lns[i, 1]]
    if len(ind) < 2:
        sys.exit('material table not found')
    else:
        f.write('materialList = {')
        for i in range(1,ind(2) - 2):
            f.write('\'%s\' ', lns[i + 1])
        f.write('};\n\n')
        offset = ind[2]
    
    i=0
    while (offset < size(lns,1) - 15):
        i += 1
        tmpo,tmpc,lns,offset,clip = read_object(lns, offset)
        # TODO: first save string then flush all at once
        f.write('object.center({%d},:) = [{%f} {%f} {%f}];\n'.format(i,tmpo.cent(1),tmpo.cent(2),tmpo.cent(3)));
        f.write('object.half_axes({%d},:) = [{%f} {%f} {%f}];\n'.format(i,tmpo.hax(1),tmpo.hax(2),tmpo.hax(3)));
        f.write('object.euler_angs({%d},:) = [{%f} {%f} {%f}];\n'.format(i,tmpo.ea(1),tmpo.ea(2),tmpo.ea(3)));
        f.write('object.density({%d}) = {%f};\\n'.format(i,tmpo.density))
        f.write('object.type({%d}) = {%d};\\n'.format(i,tmpo.type))
        f.write('object.material({%d}) = {%d};\\n'.format(i,tmpo.material))
        f.write('object.axial_lims({%d},:) = [0 0];\\n'.format(i))
        f.write('object.shape(%d) = 0;\\n',i)
        f.write('object.clip{%d} = ['.format(i))
        for j in range(1,tmpc.shape[0]):
            f.write('{%f} {%f} {%f} {%f};'.format(tmpc(j,1),tmpc(j,2),tmpc(j,3),tmpc(j,4)))
        f.write('];\\n\\n')
        clip[i]=tmpc
        obj[i]=tmpo
    
    f.close()
    print('... done writing {}'.format(ppmPhantomFilename))
    
def read_object(lines=None,offset=None):

    obj.type = lines[1 + offset]
    for i in range(4):
        A[i] = float(lines[i + 4 + offset])
    
    if np.linalg.norm(A - [0,0,0,1]) > 1e-05:
        sys.exit('Unexpected transform')
    
    obj.cent = np.copy(A[:3,3].T)
    B = A[0:3, 0:3]
    obj.hax = np.sqrt(np.sum(B*B))
    B /= np.ones(3,1)*obj.hax
    B = B.T
    obj.ea = euler_angs(B)*180/pi
    obj.transform = copy(A)
    tmp=sscanf(lines(15 + offset,arange()),'%s',1)
    tmp = lines[15+offset]
    #NOTE: we can step matlab code and compare python code
    if tmp(1) == 'C':
        i=1
        while tmp(1) == 'C':
            v = list(map(float, lines[16+offset, 21:]))
            clip[i,0:3] = v
            clip[i,4] = float(lines[16 + offset, 26:])
            clip[i,:] = -clip(i)
            offset += 3
            tmp = lines[15+offset]
            i += 1
    else:
        clip=[]
    
    if np.sum((obj.type[1:3]) == 'Cyl') == 3:
        obj.type = 2
        obj.hax[3]=obj.hax(3) / 2
    else:
        if sum(obj.type[:3] == 'Cub') == 3:
            obj.type = 2
            old=0
            if old:
                Bt=copy(B)
            else:
                Bt=B.T

            clip = clip + [Bt[:,1].T, obj.hax(1)/2 + obj.cent*Bt[:,1]]
            clip = clip + [-Bt[:,1].T, obj.hax(1)/2 - obj.cent*Bt[:,1]]
            clip = clip + [Bt[:,2].T, obj.hax(2)/2 + obj.cent*Bt[:,2]]
            clip = clip + [-Bt[:,2].T, obj.hax(2)/2 - obj.cent*Bt[:,2]]
            obj.hax[:2] = obj.hax[:2]*np.sqrt(2)
            obj.hax /= 2
        else:
            if np.sum((obj.type[:3]) == 'Sph') == 3:
                obj.type = 1
            else:
                if np.sum(obj.type[:3] == 'Ell') == 3:
                    obj_type = 1
                else:
                    sys.exit('Unrecognized obj type\\n\\r')
    
    
    obj.material = lines[18 + offset]
    obj.density = lines[19 + offset]
    offset += 21

    return obj, clip, lines, offset
    
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
    
if __name__ == '__main__':
    ppfile = "FB_head.pp"
    Phantom_Analytic_pp_to_ppm(ppfile, Scale=1., debug=True)
