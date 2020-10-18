import numpy as np
import sympy as sp
from solvers import *
from SISg_functions import *
from matplotlib.offsetbox import AnchoredText
import matplotlib.pyplot as plt
import glob as glob
plt.ion()
r0,ge_x,ge_y,bx,by = sp.symbols('r0,ge_x,ge_y,bx,by', real=True) # needed to make the functions work
keys=(r0,bx,by,ge_x,ge_y)

'''
i       image plane
x/y     intra plane dimensions
e       uncertainty
0/1/2/3 images
g       lense galaxy
'''

def lensplot(res=1000, mul=1.25):

        '''Magnification/Convergence field
        Computes kappa and gamma on whole image plane in a res x res grid
        Stretching from biggest image-lens distance * mul'''
        ##config parameters
        res=1000 # higher is slower
        mul=1.5 

        #create grid
        dist=np.sqrt(np.sum((img-lens)**2,axis=0))
        dist=np.append(np.sqrt(np.sum(Source**2)),dist)
        lim=np.max(dist)*mul
        X,Y=np.meshgrid(np.linspace(-lim+lens[0],lim+lens[0],res),np.linspace(-lim+lens[1],lim+lens[1],res))

        Z=np.stack((X.reshape(res*res,1),Y.reshape(res*res,1)),axis=1)
        Z=np.reshape(Z,(res*res,2)).T
        #compute kappa/gamma
        # ?f ... field values of point variables
        _,_,_,muf=SISg(Z,lens,R0,GE)
        muf=muf.reshape(res,res) #muf comes back as 1D vector

        '''Actual plotting routine
        Caustics ... stroked line
        Critical ... solid lines
        red      ... tangential
        blue     ... radial
        '''
        Einstein_radius= plt.Circle((lens[0], lens[1]), R0, fill=0, linewidth=.66, ec='b', ls='-.')     #radial caustic

        fig, ax = plt.subplots()

        #Plotting
        plt.contourf(X,Y,np.log10(np.abs(1/muf)), cmap='gist_gray_r')   #log of absolute magnification background
        imu_0= plt.contour(X,Y,1/muf,levels=[0],colors='r', linewidths=.66)    #1/mu=0 tangential critical line

        #Get 1/mu=0 isoline coordinates
        imu_x= np.array([])
        imu_y= np.array([])

        for i in range(np.shape(imu_0.allsegs[0])[0]):
                imu_x=np.append(imu_x, imu_0.allsegs[0][i][:,0])
                imu_y=np.append(imu_y, imu_0.allsegs[0][i][:,1])

        imu=np.vstack((imu_x,imu_y))         

        #Plotting II: The plot thickens
        ax.add_artist(Einstein_radius)

        #uncertainty ellipses
        elaxis=np.abs(epn)+np.abs(epp)
        elcent=pos+lens
        eSaxis=np.array([[np.abs(BX[1])+np.abs(BX[2])],[np.abs(BY[1])+np.abs(BY[2])]])

        draw_err_ellipse(img,2*eimg,colour='k')# observed img
        draw_err_ellipse(lens,2*elens,colour='w')# observed lens
        draw_err_ellipse(pos+lens,np.abs(epn+epp),colour='r')# computed img
        draw_err_ellipse(Source+lens,eSaxis,colour='b')# source

        plt.plot(img[0],img[1],'xk') #observed images
        plt.plot(Source[0]+lens[0],Source[1]+lens[1],'xb') #computed source
        plt.plot(pos[0]+lens[0],pos[1]+lens[1],'+r') #computed images
        plt.plot(lens[0],lens[1],'+w',alpha=1) #observed lense
        plt.legend(('Observed Images','Computed Source','Computed Image','Observed Lense'),frameon=False)
        plt.title(name)

        #Reverse critical line to caustic
        XL=ax.get_xlim()
        YL=ax.get_ylim()

        cx=-1*(R0*(imu[0]/(np.sqrt(imu[0]**2+imu[1]**2)))+sis[ge_x][0]*imu[0]+sis[ge_y][0]*imu[1]-imu[0])
        cy=-1*(R0*(imu[1]/(np.sqrt(imu[0]**2+imu[1]**2)))-sis[ge_x][0]*imu[1]+sis[ge_y][0]*imu[0]-imu[1])

        #Plotting III: The return of the text
        plt.xlabel('RA [as]')
        plt.ylabel('Dec [as]')

        KGtxt='\t'+r'$\kappa$'+'\t'+r'$\gamma$'+'\n'
        for i in range(len(k)):
                KGtxt=KGtxt+components[i]+'  '+'%.3f' % (k[i])+'  '+'%.3f' % (g[i])+'\n'
                ax.text(pos[0,i]*1.2+lens[0],pos[1,i]*1.2+lens[1], components[i], fontsize=10,color='y')
                box=AnchoredText(KGtxt,loc=4,frameon=False)
        ax.add_artist(box)
        plt.plot(cx+lens[0],cy+lens[1],linewidth=.66,ls='--',color='r')
        plt.xlim(XL)
        plt.ylim(YL)

        plt.gca().invert_xaxis()


try:
        path=input('Enter filename / path:\n\nor manual data entry press enter\n')
        ix,iy,igx,igy,iex,iey,iegx,iegy,name,components = get_coords(path)
        if len(ix)==2:
                mag=input('\nPlease enter magnitude [mag] of your Images:\n')
                emag=input('\nand what would the uncertainties be on that?:\n')
                mag=np.array(mag.split(',')).astype(float)
                emag=np.array(emag.split(',')).astype(float)
        else:
                mag=np.zeros(4)
                emag=np.zeros(4)
except:
        done=False
        parts=[]
        name=input('\n\nOK, lets get typing then\nSystem Name:\n')
        parts.append(input('\nLens data:\nRA[as], RA_error, DEC[as], DEC_error\n'))
        parts.append(input('\nFirst image data:\nName, RA[as], RA_error, DEC[as], DEC_error, mag[mag], mag_error\n'))
        parts.append(input('\nSecond image data:\nName, RA[as], RA_error, DEC[as], DEC_error, mag[mag], mag_error\n'))
        while done==False:
                parts.append(input('\nPress enter to start or enter two additional images:\nName, RA[as], RA_error, DEC[as], DEC_error, mag[mag], mag_error\n'))
                if parts[-1]=='':
                        done=True

        Lens=  parts.pop(0).split(',')
        igx=float(Lens[0])
        igy=float(Lens[2])
        iegx=float(Lens[1])
        iegy=float(Lens[3])

        components=[]
        ix=np.array([])
        iy=np.array([])
        iex=np.array([])
        iey=np.array([])
        mag=np.array([])
        emag=np.array([])
        for i in range(len(parts)-1):
                particle=parts[i].split(',')
                components.append(particle[0])
                ix=np.append(ix,particle[1]).astype(float)
                iex=np.append(iex,particle[2]).astype(float)
                iy=np.append(iy,particle[3]).astype(float)
                iey=np.append(iey,particle[4]).astype(float)
                mag=np.append(mag,particle[5]).astype(float)
                emag=np.append(emag,particle[6]).astype(float)

###Prep
lens=np.array(([igx],[igy]))
img=np.vstack((ix,iy))
#mag= np.array([15.33,16.83]) #asumes magB-magA in Twin2 solver

elens=np.array(([iegx],[iegy]))
eimg=np.vstack((iex,iey))
#emag= np.array([0.02,0.03])

eln=lens-elens
ein=img-eimg
emn=mag-emag

elp=lens+elens
eip=img+eimg
emp=mag+emag

###SIS Parameter and computed Images positions
if len(ix)==2: #Twin
        t2=Twin2(img,lens,mag)
        #sis=dict(zip(keys,t2))    
        tep=Twin2(eip,elp,emp) 
        ten=Twin2(ein,eln,emn)
        tep=np.array(tep)-np.array(t2)
        ten=np.array(t2)-np.array(ten)
        SIS=np.vstack((t2,tep,ten))

        sis={} #[SIS+g Parameter, +Uncertainty, -Uncertainty]
        for i in range(SIS.shape[1]):
                sis[keys[i]]=SIS[:,i]        
       
        pos=Source2ImgRT(sis[r0][0],sis[bx][0],sis[by][0],sis[ge_x][0],sis[ge_y][0],img-lens)
        epp=Source2ImgRT(sis[r0][0]+sis[r0][1],
                         sis[bx][0]+sis[bx][1],
                         sis[by][0]+sis[by][1],
                         sis[ge_x][0]+sis[ge_x][1],
                         sis[ge_y][0]+sis[ge_y][1],
                         eip-elp)
        epn=Source2ImgRT(sis[r0][0]-sis[r0][2],
                         sis[bx][0]-sis[bx][2],
                         sis[by][0]-sis[by][2],
                         sis[ge_x][0]-sis[ge_x][2],
                         sis[ge_y][0]-sis[ge_y][2],
                         ein-eln)
        BX=np.array(sis[bx])
        BY=np.array(sis[by])
else: #Quads
        sis= Quad2(img,lens)
        sip= Quad2(eip,elp)
        sin= Quad2(ein,eln)

        R0=np.array([sis[r0],sip[r0]-sis[r0],sis[r0]-sin[r0]])
        BX=np.array([sis[bx],sip[bx]-sis[bx],sis[bx]-sin[bx]])
        BY=np.array([sis[by],sip[by]-sis[by],sis[by]-sin[by]])
        GE_X=np.array([sis[ge_x],sip[ge_x]-sis[ge_x],sis[ge_x]-sin[ge_x]])
        GE_Y=np.array([sis[ge_y],sip[ge_y]-sis[ge_y],sis[ge_y]-sin[ge_y]])

        SIS=np.stack((R0,BX,BY,GE_X,GE_Y),axis=0)
        sis={} #[SIS+g Parameter, +Uncertainty, -Uncertainty]
        for i in range(SIS.shape[0]):
                sis[keys[i]]=SIS[i,:]     
        
        pos=Source2ImgRT(sis[r0][0],sis[bx][0],sis[by][0],sis[ge_x][0],sis[ge_y][0],img-lens)
        epp=Source2ImgRT(sis[r0][0]+sis[r0][1],
                         sis[bx][0]+sis[bx][1],
                         sis[by][0]+sis[by][1],
                         sis[ge_x][0]+sis[ge_x][1],
                         sis[ge_y][0]+sis[ge_y][1],
                         eip-elp)
        epn=Source2ImgRT(sis[r0][0]-sis[r0][2],
                         sis[bx][0]-sis[bx][2],
                         sis[by][0]-sis[by][2],
                         sis[ge_x][0]-sis[ge_x][2],
                         sis[ge_y][0]-sis[ge_y][2],
                         ein-eln)


###Get kappa and gamma
GE=np.array([[float(sis[ge_x][0])],[float(sis[ge_y][0])]])
eGEp=np.array([[float(sis[ge_x][1])],[float(sis[ge_y][1])]])
eGEn=np.array([[float(sis[ge_x][2])],[float(sis[ge_y][2])]])

R0=(float(sis[r0][0]))
eR0=np.array([(float(sis[r0][1])),(float(sis[r0][2]))])

Source=np.array([[float(sis[bx][0])],[float(sis[by][0])]])
eSp=np.array([[float(sis[bx][1])],[float(sis[by][1])]])
eSn=np.array([[float(sis[bx][2])],[float(sis[by][2])]])

k,g,_,_=SISg(img,lens,R0,GE)
kp,gp,_,_=SISg(eip,elp,R0+eR0[0],GE+eGEp)
kn,gn,_,_=SISg(ein,eln,R0-eR0[0],GE-eGEn)

ek=np.array([kp-k,k-kn])
eg=np.array([gp-g,g-gn])

epp=epp-pos
epn=pos-epn

print('\n'+name)
print(components[0:len(k)])
print('R0:'+str(R0)+' +'+str(eR0[0])+' '+str(eR0[1]))
print()
print('k: '+str(k)+'\n+  '+str(ek[0])+'\n-  '+str(ek[1]))
print()
print('g: '+str(g)+'\n+  '+str(eg[0])+'\n-  '+str(eg[1]))

if input('\nDo you want a plot with that?\n yes / no\n').lower()[0]=='y':
        lensplot()
