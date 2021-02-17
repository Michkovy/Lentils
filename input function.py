import numpy as np

print('LENTILS assumes dmag=mag_B-mag_A for calculation.\n\nMagnitude differences are input using mag_A=0 & mag_B=dmag.\ndmag caluclated using the equation in the first line.\n\nIf flux ratios are used convert them to dmag; using dmag=-2.5*log(f_ba)\n\n\n\n')

try:
        path=input('Enter filename / path:\n\npress enter for manual data entry\n')
        file=open(path, 'r')
        parts=file.read()
        parts=parts.split('\n')
        parts= [x for x in parts if not x.startswith('#')] #remove comment lines
       
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
anti_div0= False

if np.sum(elens+eimg+emag)!=0:
        eln=lens-elens
        ein=img-eimg
        emn=mag-emag

        elp=lens+elens
        eip=img+eimg
        emp=mag+emag

        no_error= False
        
if 0 in np.sum(img-lens,0):
        lens=lens+0.001
        anti_div0=True

