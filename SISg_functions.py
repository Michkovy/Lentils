import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as pch


def minmax(pos):
    mima=np.zeros(2)
    mima[0]=np.min(pos)
    mima[1]=np.max(pos)
    return mima


def SISg(img,lens,R0,GE):
    r=img-lens
    r_abs=np.sqrt((r[0]**2+r[1]**2))

#field
    k=1/2*R0/r_abs
    gsis=np.array([1/2*R0*((-r[0]**2+r[1]**2)/(r_abs**3)),-R0*((r[0]*r[1])/(r_abs**3))])
    g=[gsis[0]+GE[0],gsis[1]+GE[1]]
    g_abs=np.sqrt((g[0]**2+g[1]**2))

    mu=1/((1-k)**2-g_abs**2).flatten()
#field
    psi11= R0*((r[1]**2)/(r_abs**3))+GE[0]
    psi22= R0*((r[0]**2)/(r_abs**3))-GE[0]
    psi12= -R0*((r[0]*r[1])/(r_abs**3))+GE[1]

    psi=np.zeros([img.size,2,2])
    #for i in range(len(ix)):
    #    psi[i]=np.array([[psi11[i],psi12[i]],[psi12[i],psi22[i]]])

    
    #mu_ratios=np.array([[mu/mu[0]],[mu/mu[1]],[mu/mu[2]],[mu/mu[3]]])
    return k,g_abs,psi,mu

def draw_err_ellipse(center,axis,colour='r'):
        '''draws error ellipses to figure
        center ... 2 x n array with ellipse center coordinates
        axis   ... 2 x n array with sMaj and sMin axis lenghts
        colour ... fill colour'''
        fig=plt.gcf()
        ax=fig.gca()
        for i in range(center.shape[1]):
                ell=pch.Ellipse(center[:,i],axis[0,i],axis[1,i],color=colour,alpha=.25)
                ax.add_artist(ell)

def get_coords(QSO):
    '''filename or path'''
    #if path
    
    file=open(QSO)
    name=file.readline()[:-1]
    components=file.readline()[:-1].split('\t')
    a=np.genfromtxt(file,delimiter='\t')
    ind_g= components.index('G')

    ix=a[0,:ind_g]
    iy=a[2,:ind_g]
    igx=a[0,ind_g]
    igy=a[2,ind_g]

    iex=a[1,:ind_g]
    iey=a[3,:ind_g]
    iegx=a[1,ind_g]
    iegy=a[3,ind_g]

    file.close()
    return ix,iy,igx,igy,iex,iey,iegx,iegy,name,components

#Debug
##x= np.array([0.385, 0, -0.336, 0.948])
##y= np.array([0.317, 0, -0.750, -0.802]) #image y position
##gx= 0.742   #lense x position
##gy= -0.656   #lense y position
##R0=0.7781806812845325
##GE=np.array([[-0.08265199],[ 0.25384874]])
##
##lims=np.vstack((minmax(x),minmax(y)))
##X,Y=np.meshgrid(np.linspace(lims[0,0],lims[0,1],200),np.linspace(lims[1,0],lims[1,1],200))
##
##kf,gf,_,_=SISg(X,Y,gx,gy,R0,GE)
