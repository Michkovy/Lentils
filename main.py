import numpy as np
import sympy as sp
from solvers import *
from SISg_functions import *
from matplotlib.offsetbox import AnchoredText
import matplotlib.pyplot as plt
import glob
import os
from scipy.optimize import fsolve

plt.ion()
r0,ge_x,ge_y,bx,by = sp.symbols('r0,ge_x,ge_y,bx,by', real=True) # needed to make the functions work
keys=(r0,bx,by,ge_x,ge_y)
eGE=0

##
##'''
##i       image plane
##x/y     intra plane dimensions
##e       uncertainty
##0/1/2/3 images
##g       lense galaxy
##'''

def tsurfplot(res=250, mul=1.25):
        dist=np.sqrt(np.sum((img-lens)**2,axis=0))
        dist=np.append(np.sqrt(np.sum(Source+lens)**2),dist)
        lim=np.max(dist)*mul
        X,Y=np.meshgrid(np.linspace(-lim+lens[0],lim+lens[0],res),np.linspace(-lim+lens[1],lim+lens[1],res))

        Z=np.stack((X.reshape(res*res,1),Y.reshape(res*res,1)),axis=1)
        Z=np.reshape(Z,(res*res,2)).T

        z=Z-lens
        fig, ax = plt.subplots(figsize=(10,8))

        e1=BX[0]+R0*(z[0]/(np.sqrt(z[0]**2+z[1]**2)))+GE[0]*z[0]+GE[1]*z[1]-z[0]
        e2=BY[0]+R0*(z[1]/(np.sqrt(z[0]**2+z[1]**2)))-GE[0]*z[1]+GE[1]*z[0]-z[1]
        e1=np.reshape(e1,(res,res))
        e2=np.reshape(e2,(res,res))

        plt.contourf(X,Y,e2,cmap='BrBG_r',alpha=0.5)
        plt.contourf(X,Y,e1,cmap='BrBG',alpha=0.5)
        plt.contour(X,Y,e1+e2,levels=[0],linewidths=.75)
        plt.plot((img)[0],(img)[1],'xk')
        plt.plot(lens[0],lens[1],'+k')

        plt.gca().invert_xaxis()
        plt.title(name)
        for i in range(len(k)):
                ax.text(img[0,i]*1.2,img[1,i]*1.2, components[i], fontsize=10)

        plt.xlabel('RA [as]')
        plt.ylabel('Dec [as]')


        
def lensplot(res=250, mul=1.25):

        '''Magnification/Convergence field
        Computes kappa and gamma on whole image plane in a res x res grid
        Stretching from biggest image-lens distance * mul
        '''
        ##config parameters
        #res=1000 # higher is slower
        #mul=1.5 

        #create grid
        dist=np.sqrt(np.sum((img-lens)**2,axis=0))
        dist=np.append(np.sqrt(np.sum(Source+lens)**2),dist)
        lim=np.max(dist)*mul
        X,Y=np.meshgrid(np.linspace(-lim+lens[0],lim+lens[0],res),np.linspace(-lim+lens[1],lim+lens[1],res))

        Z=np.stack((X.reshape(res*res,1),Y.reshape(res*res,1)),axis=1)
        Z=np.reshape(Z,(res*res,2)).T
        #Z=Z-lens
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
        
        fig, ax = plt.subplots(figsize=(10,8))

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
        
        if no_error==False:
            elaxis=np.abs(epn)+np.abs(epp)
            elcent=pos+lens
            eSaxis=np.array([[np.abs(BX[1])+np.abs(BX[2])],[np.abs(BY[1])+np.abs(BY[2])]])

            draw_err_ellipse(img,2*eimg,colour='k')# observed img
            draw_err_ellipse(lens,2*elens,colour='w')# observed lens
            draw_err_ellipse(pos+lens,np.abs(epn+epp),colour='r')# computed img
            draw_err_ellipse(Source+lens,eSaxis,colour='b')# source

        plt.plot(img[0],img[1],'xk') #observed images
        plt.plot(pos[0]+lens[0],pos[1]+lens[1],'+r') #computed images
        plt.plot(Source[0]+lens[0],Source[1]+lens[1],'xb') #computed source
        plt.plot(lens[0],lens[1],'+w',alpha=1) #observed lense
        plt.legend(('Observed Images','Computed Images','Computed Source','Observed Lense'),frameon=False,ncol=4, loc='best',
                   bbox_to_anchor=(1,-0.06))
        plt.title(name)

        #Reverse critical line to caustic
        XL=ax.get_xlim()
        YL=ax.get_ylim()

        c=Img2SourceRT(sis,imu)
        
        #Plotting III: The return of the text
        plt.xlabel('RA [as]')
        plt.ylabel('Dec [as]')

        KGtxt='\t'+r'$\kappa$'+'\t'+r'$\gamma$'+'\n'
        for i in range(len(k)):
                KGtxt=KGtxt+components[i]+'  '+'%.3f' % (k[i])+'  '+'%.3f' % (g[i])+'\n'
                ax.text(pos[0,i]*1.2+lens[0],pos[1,i]*1.2+lens[1], components[i], fontsize=10,color='y')
                box=AnchoredText(KGtxt,loc=4,frameon=False)
        #ax.add_artist(box)
        #plt.plot(c[0]+lens[0],c[1]+lens[1],linewidth=.66,ls='--',color='g')
        #plt.plot(-c[0]+lens[0],-c[1]+lens[1],linewidth=.66,ls='--',color='g')
        plt.xlim(XL)
        plt.ylim(YL)

        plt.gca().invert_xaxis()


        #sc=Img2SourceRT(sis,Z)
        #plt.plot(-sc[0]+lens[0],-sc[1]+lens[1],',b',alpha=.0125)
        #plt.plot(sc[0]+lens[0],sc[1]+lens[1],',b',alpha=.0125)
        return Z

def SISg(img,lens,R0,GE):
        '''Calculates convergence, shear, lens mapping and magnification'''
        r=img-lens
        r_abs=np.sqrt((r[0]**2+r[1]**2))

        k=1/2*R0/r_abs
        gsis=np.array([1/2*R0*((-r[0]**2+r[1]**2)/(r_abs**3)),-R0*((r[0]*r[1])/(r_abs**3))])
        g=[gsis[0]+GE[0],gsis[1]+GE[1]]
        g_abs=np.sqrt((g[0]**2+g[1]**2))

        mu=1/((1-k)**2-g_abs**2).flatten()

        psi11= R0*((r[1]**2)/(r_abs**3))+GE[0]
        psi22= R0*((r[0]**2)/(r_abs**3))-GE[0]
        psi12= -R0*((r[0]*r[1])/(r_abs**3))+GE[1]

        i=np.arange(0,len(img[0]))
        psi=np.zeros((len(img[0]),2,2))
        psi[i]=np.array(((psi11[i],psi12[i]),(psi12[i],psi22[i]))).T
        #mu=1/np.linalg.det(psi)
        return k,g_abs,psi,mu

def Img2SourceRT(sis, coord):
    '''Raytrace from image to source plane. Using SIS+g lens equation'''

    R0= sis[r0][0]
    GE_X= sis[ge_x][0]
    GE_Y= sis[ge_y][0]
    ix=coord[0]
    iy=coord[1]
    n=len(ix)

    sx=lambda ix,iy: -GE_X*ix - GE_Y*iy - ix*R0/np.sqrt(ix**2 + iy**2) + ix
    sy=lambda ix,iy:  GE_X*iy - GE_Y*ix - iy*R0/np.sqrt(ix**2 + iy**2) + iy
    
    sCoor=np.zeros((2,n))
    for j in range(n):
            sCoor[0,j]=sx(ix[j],iy[j])
            sCoor[1,j]=sy(ix[j],iy[j])

    return sCoor

def Source2ImgRT2(sis,coord):
    '''Raytrace from source to image plane. Using SIS+g lens equation'''
    R0= sis[r0]
    GE_X= sis[ge_x]
    GE_Y= sis[ge_y]
    BX= sis[bx]
    BY= sis[by]
    ix=coord[0]
    iy=coord[1]
    n=len(ix)

    def eqs(var):
        tx,ty=var
        e1=BX+R0*(tx/(np.sqrt(tx**2+ty**2)))+GE_X*tx+GE_Y*ty-tx
        e2=BY+R0*(ty/(np.sqrt(tx**2+ty**2)))-GE_X*ty+GE_Y*tx-ty
        return (e1,e2)

    iCoor=np.zeros_like(coord)
    for i in range(n):
        iCoor[:,i]=fsolve(eqs,coord[:,i])
    
    return iCoor

######### Main #############
def single_file(file):
        file=open(path, 'r')
        parts=file.read()
        parts=parts.split('\n')
        parts= [x for x in parts if not x.startswith('#')] #remove comment lines
        file.close()
        return parts

def batch(path):
        path=path+'/*.txt'
        files= glob.glob(path)
        return files


single=True #controls printout, surpressed if false

print('LENTILS assumes dmag=mag_B-mag_A for calculation.\n\nMagnitude differences are input using mag_A=0 & mag_B=dmag.\ndmag caluclated using the equation in the first line.\n\nIf flux ratios are used convert them to dmag; using dmag=-2.5*log(f_ba)\n\n\n')

try:
        path=input('Enter filename / path:\n\npress enter for manual data entry\n')
        if os.path.isdir(path)==True:
                files= batch(path)
                single=False

        else:
                os.path.isfile(path)
                parts= single_file(path)
               
except:
        done=False
        parts=[]
        parts.append(input('\n\nOK, lets get typing then\nSystem Name:\n'))
        parts.append(input('\nLens data:\nRA[as], RA_error, DEC[as], DEC_error\n'))
        parts.append(input('\nFirst image data:\nName, RA[as], RA_error, DEC[as], DEC_error, mag[mag], mag_error\n'))
        parts.append(input('\nSecond image data:\nName, RA[as], RA_error, DEC[as], DEC_error, mag[mag], mag_error\n'))
        while done==False:
                parts.append(input('\nPress enter to start or enter two additional images:\nName, RA[as], RA_error, DEC[as], DEC_error, mag[mag], mag_error\n'))
                if parts[-1]=='':
                        parts.pop() #remove empty last entry if entering manually
                        done=True


def main(parts):
        
        global name,components,R0,eR0,GE,k,ek,g,eg,anti_div0,no_error,img,lens,Source,eSp,epn,epp,pos,BX,BY,sis,eimg,elens,eGE

        name=parts.pop(0)
        Lens= parts.pop(0).split(',')
        igx= float(Lens[0])
        igy= float(Lens[2])
        iegx= float(Lens[1])
        iegy= float(Lens[3])

        components=[]
        ix=np.array([])
        iy=np.array([])
        iex=np.array([])
        iey=np.array([])
        mag=np.array([])
        emag=np.array([])

        for i in range(len(parts)):
                if parts[i]=='':
                        break
                particle=parts[i].split(',')
                components.append(particle[0])
                ix=np.append(ix,particle[1]).astype(float)
                iex=np.append(iex,particle[2]).astype(float)
                iy=np.append(iy,particle[3]).astype(float)
                iey=np.append(iey,particle[4]).astype(float)
                mag=np.append(mag,particle[5]).astype(float)
                emag=np.append(emag,particle[6]).astype(float)

        #Prep
        lens=np.array(([igx],[igy]))
        img=np.vstack((ix,iy))

        elens=np.array(([iegx],[iegy]))
        eimg=np.vstack((iex,iey))

        no_error= True
        eln=ein=emn=elp=eip=emp=0

        anti_div0= False


        if np.sum(elens+eimg+emag)!=0:
                eln=lens-elens
                ein=img-eimg
                emn=mag-emag

                elp=lens+elens
                eip=img+eimg
                emp=mag+emag

                no_error= False
                
        if 0 in np.sum((img-lens)**2,0):
                lens=lens+0.001
                anti_div0=True

        ###SIS Parameter and computed Images positions
        if no_error==False:
            if len(ix)==2: #Twin
                    t2=Twin2(img,lens,mag)
                    tep=Twin2(eip,elp,emp)
                    ten=Twin2(ein,eln,emn)

                    sis={}
                    for i in range(len(keys)):
                            sis[keys[i]]=t2[i]
                    sip={}
                    for i in range(len(keys)):
                            sip[keys[i]]=tep[i]
                    sin={}
                    for i in range(len(keys)):
                            sin[keys[i]]=ten[i]
                 
                    pos=Source2ImgRT2(sis,img-lens)
                    epp=Source2ImgRT2(sip,eip-elp)
                    epn=Source2ImgRT2(sin,ein-eln)

                    tep=np.array(tep)-np.array(t2) #
                    ten=np.array(t2)-np.array(ten) #
                    SIS=np.vstack((t2,tep,ten)) #

                    sis={} #[SIS+g Parameter, +Uncertainty, -Uncertainty]
                    for i in range(SIS.shape[1]):
                            sis[keys[i]]=SIS[:,i]        
                   
                    BX=np.array(sis[bx])
                    BY=np.array(sis[by])
            else: #Quads
                    sis= Quad2(img,lens)
                    sip= Quad2(eip,elp)
                    sin= Quad2(ein,eln)

                    pos=Source2ImgRT2(sis,img-lens)
                    epp=Source2ImgRT2(sip,eip-elp)
                    epn=Source2ImgRT2(sin,ein-eln)

                    R0=np.array([sis[r0],sip[r0]-sis[r0],sis[r0]-sin[r0]])
                    BX=np.array([sis[bx],sip[bx]-sis[bx],sis[bx]-sin[bx]])
                    BY=np.array([sis[by],sip[by]-sis[by],sis[by]-sin[by]])
                    GE_X=np.array([sis[ge_x],sip[ge_x]-sis[ge_x],sis[ge_x]-sin[ge_x]])
                    GE_Y=np.array([sis[ge_y],sip[ge_y]-sis[ge_y],sis[ge_y]-sin[ge_y]])

                    SIS=np.stack((R0,BX,BY,GE_X,GE_Y),axis=0)
                    sis={} #[SIS+g Parameter, +Uncertainty, -Uncertainty]
                    for i in range(SIS.shape[0]):
                            sis[keys[i]]=SIS[i,:]

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

            eGE=np.max((np.sqrt(np.sum(eGEp**2)),np.sqrt(np.sum(eGEn**2))))
            
        else:
            if len(ix)==2:
                t2=Twin(img,lens,mag)

                sis={}
                for i in range(len(keys)):
                    sis[keys[i]]=t2[i]

                pos=Source2ImgRT2(sis,img-lens)

                ten=tep=np.zeros_like(t2)
                SIS=np.vstack((t2,tep,ten))

                sis={}
                for i in range(SIS.shape[1]):
                    sis[keys[i]]=SIS[:,i]


            else:
                sis= Quad2(img,lens)
                pos=Source2ImgRT2(sis, img-lens)

                R0=np.array([sis[r0],0,0])
                BX=np.array([sis[bx],0,0])
                BY=np.array([sis[by],0,0])
                GE_X=np.array([sis[ge_x],0,0])
                GE_Y=np.array([sis[ge_y],0,0])

                SIS=np.stack((R0,BX,BY,GE_X,GE_Y),axis=0)
                sis={} #[SIS+g Parameter, +Uncertainty, -Uncertainty]
                for i in range(SIS.shape[0]):
                    sis[keys[i]]=SIS[i,:]
                
            ###Get kappa and gamma
            GE=np.array([[float(sis[ge_x][0])],[float(sis[ge_y][0])]])
            eGEp= eGEn= np.zeros_like(GE)

            R0=(float(sis[r0][0]))
            eR0=np.array([0,0])

            Source=np.array([[float(sis[bx][0])],[float(sis[by][0])]])
            eSp= eSn=np.zeros_like(Source)

            k,g,_,_=SISg(img,lens,R0,GE)

            ek= eg= epp= epn= np.zeros_like(k)

if single==True:
        main(parts)

        print('\n'+name)
        print('R0:'+str(R0)+' +'+str(eR0[0])+' '+str(eR0[1]))
        print('GE:'+str(np.sqrt(GE[0]**2+GE[1]**2)[0])+'+-'+str(eGE))
        print()
        print(components[0:len(k)])
        print('k: '+str(k)+'\n+  '+str(ek[0])+'\n-  '+str(ek[1]))
        print()
        print('g: '+str(g)+'\n+  '+str(eg[0])+'\n-  '+str(eg[1]))
        print()
        if anti_div0==True:
                print('Note:\nLens and one image overlapped, shifted lens by .001as to prevent division by 0')
        if no_error==True:
                print('Note:\nNo errors computed, since none were given')



        if input('\nDo you want a plot with that?\n yes / no\n').lower()[0]=='y':
                lensplot()
                #tsurfplot()

else:
        print('BATCH Mode!!!')
        out=open('results.txt','w')
        name_str=str('ID, Einstein radius, +Einstein radius error, -Einstein radius error,total external shear,Source x coordinate, Source x +error, Source x -error,Source y coordinate, Source y +error, Source y -error,shear x component, shear x +error, shear x -error,shear y component, shear y +error, shear y -error\nImage, convergence, convergence +error, convergence -error, shear, shear +error, shear -error\n')
        out.write(name_str)
        for i in files:
                file=open(i,'r')
                parts=file.read()
                parts=parts.split('\n')
                parts= [x for x in parts if not x.startswith('#')]
                file.close()
                main(parts)
                
                lens_str=name+','+str(R0)+','+str(eR0[0])+','+str(eR0[1])+','+str(np.sqrt(GE[0]**2+GE[1]**2))+','+str(BX[0])+','+str(BX[1])+','+str(BX[2])+','+str(BY[0])+','+str(BY[1])+','+str(BY[2])+','+str(sis[ge_x][0])+','+str(sis[ge_x][1])+','+str(sis[ge_x][2])+','+str(sis[ge_y][0])+','+str(sis[ge_y][1])+','+str(sis[ge_y][2])
                lens_str=lens_str.replace(']','')
                lens_str=lens_str.replace('[','')
                
                comp_str=str()
                if no_error==True:
                        for j in range(len(components)):
                                comp_str=comp_str+str(components[j])+','+str(k[j])+','+str(ek[j])+','+str(ek[j])+','+str(g[j])+','+str(eg[j])+','+str(eg[j])+','+'NO ERROR GIVEN'+'\n'
                else:
                        for j in range(len(components)):
                                comp_str=comp_str+str(components[j])+','+str(k[j])+','+str(ek[0,j])+','+str(ek[1,j])+','+str(g[j])+','+str(eg[0,j])+','+str(eg[1,j])+'\n'
                out.write(lens_str+'\n'+comp_str+'\n')

        out.close()
