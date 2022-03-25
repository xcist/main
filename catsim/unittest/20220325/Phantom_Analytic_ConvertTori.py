# Generated with SMOP  0.41
# Phantom_Analytic_ConvertTori.m

    # -----------------------------------------------------------------------
#   Program Name: Phantom_Analytic_ConvertTori.m                                     
#   Author:  Samit Basu, Bruno De Man, Jed Pack (GE Global Research)
#   Organization: GE Global Research and GE Healthcare, General Electric Company
#   Version:  6.0.3
#   Date:  Feb 3, 2015
#   Class: GE Confidential. General Electric Proprietary Data (c) 2012 General Electric Company
    
    # History:
#   2012-09-19 Paul FitzGerald (GE Global Research)
#              Renamed - was convert_tori.m. 
#              Added "Verbose" output.
# -----------------------------------------------------------------------
    
'''
log:
    20220322: tested by comparing std, sum
'''

# phantObject: a list of dictionaries
def Phantom_Analytic_ConvertTori(phantObject=None,params=None,clip=None, debug=False):

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
   
    if debug: breakpoint()
    return params
    
if __name__ == '__main__':
    from scipy import io as sio
    mat_in = sio.loadmat("Phantom_Analytic_ConvertTori_phantObject_in.mat")['phantObject']
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

    
    params_in = sio.loadmat("Phantom_Analytic_ConvertTori_params_in.mat")['params']
    #params_in = params_in.astype(int)
    #clip_in = sio.loadmat("converttori_clip_in.mat")['clip'][0]
    clip_in = sio.loadmat("Phantom_Analytic_ConvertTori_clip_in.mat")['clip'][0]
    #clip_in = params_in.astype(int)
    #breakpoint()
    params = Phantom_Analytic_ConvertTori(phantObject=python_in, params=params_in, clip=clip_in, debug=False)

    mat_out = sio.loadmat("Phantom_Analytic_ConvertTori_out.mat")
    params_mat = mat_out['params']

    import numpy as np
    #breakpoint()
    try:
        assert np.allclose(params, params_mat, atol=1.E-8), 'params'
    except AssertionError as err:
        print(err)
        breakpoint()
