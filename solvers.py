import numpy as np
import sympy as sp
import dill

def Quad2(img, lens):

    #Prep###############################################################################
    n=4
    coord=img-lens
    rx=coord[0]
    ry=coord[1]
    r_abs= np.sqrt(rx**2+ry**2)

    r0,ge_x,ge_y,bx,by = sp.symbols('r0,ge_x,ge_y,bx,by', real=True) #Main symbols

    # define functions that will take the role of indexed symbols
    x= sp.Function('x')
    y= sp.Function('y')
    R_abs= sp.Function('R_abs')
    i= sp.symbols('i',posititve=True, real=True)

    x_fun = lambda index: rx[index-1]
    y_fun= lambda index: ry[index-1]
    R_abs_fun= lambda index: r_abs[index-1]
    #Symbolic##################################################################
    s2 = sp.Sum((bx + r0*x(i)/R_abs(i) + ge_x*x(i) + ge_y*y(i) - x(i))**2 + (by + r0*y(i)/R_abs(i) - ge_x*y(i) + ge_y*x(i) - y(i))** 2, (i, 1, n))

    dbx= sp.diff(s2,bx).doit().replace(y,y_fun).replace(R_abs,R_abs_fun).replace(x,x_fun)
    dby= sp.diff(s2,by).doit().replace(y,y_fun).replace(R_abs,R_abs_fun).replace(x,x_fun)
    dgx= sp.diff(s2,ge_x).doit().replace(y,y_fun).replace(R_abs,R_abs_fun).replace(x,x_fun)
    dgy= sp.diff(s2,ge_y).doit().replace(y,y_fun).replace(R_abs,R_abs_fun).replace(x,x_fun)
    dr0= sp.diff(s2,r0).doit().replace(y,y_fun).replace(R_abs,R_abs_fun).replace(x,x_fun)
    sol=sp.solve((dbx,dby,dgx,dgy,dr0),(r0,ge_x,ge_y,bx,by), simplify=True)
    
    R0= sol[r0]
    GE_X= sol[ge_x]
    GE_Y= sol[ge_y]
    BX= sol[bx]
    BY= sol[by]

    #theta=Source2ImgRT(R0,BX,BY,GE_X,GE_Y,coord)

    return sol#, theta


def Twin2(img,lens,mag):
    r=img-lens
    rxa=r[0,0]
    rya=r[1,0]
    rxb=r[0,1]
    ryb=r[1,1]
    dmag=mag[1]-mag[0]
    R0,BX,BY,GE_X,GE_Y=dill.load(open('SISe_Solved','rb'))

    r0=list(R0(rxa,rya,rxb,ryb,dmag))
    r0=r0[r0!=0] #selects non zero solution
    bx=BX(rxa,rya,rxb,ryb,r0)
    by=BY(rxa,rya,rxb,ryb,r0)
    ge_x=GE_X(rxa,rya,rxb,ryb,r0)
    ge_y=GE_Y(rxa,rya,rxb,ryb,r0)

    return r0,bx,by,ge_x,ge_y
#{b1 -> -0.607617, b2 -> -1.21824, t0 -> 2.19, g1 -> 0.570545, g2 -> 0.216076}}

def Img2SourceRT(sis, coord):

    R0= sis[r0]
    GE_X= sis[ge_x]
    GE_Y= sis[ge_y]
    ix=coord[0]
    iy=coord[1]
    n=len(ix)

    sx=lambda ix,iy: -GE_X*ix - GE_Y*iy - ix*R0/np.sqrt(ix**2 + iy**2) + ix
    sy=lambda ix,iy:  GE_X*iy - GE_Y*ix - iy*R0/np.sqrt(ix**2 + iy**2) + iy

    sCoor=np.zeros((2,n))
    for j in range(n):
            theta[0,j]=sx(ix[j],iy[j])
            theta[1,j]=sy(ix[j],iy[j])

    return sCoor

def Source2ImgRT(R0,BX,BY,GE_X,GE_Y,coord):
    
    n=np.shape(coord)[1]
    theta=np.zeros((2,n))

    tx, ty= sp.symbols('tx, ty', real=True)
    e1= sp.Eq(BX+R0*(tx/(sp.sqrt(tx**2+ty**2)))+GE_X*tx+GE_Y*ty-tx,0)
    e2= sp.Eq(BY+R0*(ty/(sp.sqrt(tx**2+ty**2)))-GE_X*ty+GE_Y*tx-ty,0)

    try:
        for j in range(n):
            a=sp.nsolve((e1,e2),(tx,ty),(coord[0,j],coord[1,j]))
            theta[0,j]=a[0]
            theta[1,j]=a[1]
    except ValueError:
        print('No computed image found wihtin tolerance')
        
    return theta

#Old Code
####################################################################################################################################
def Twin(ix,iy,igx,igy,mag,solver_type='analytic'): #obsolete

    n=2
    rx=ix-igx
    ry=iy-igy
    dmag=mag[1]-mag[0]
    f_ab=1/(10**(dmag/-2.5))
    r_abs= np.sqrt(rx**2+ry**2)

    r0,ge_x,ge_y,bx,by = sp.symbols('r0,ge_x,ge_y,bx,by', real=True)

    #Symbolic##################################################################
    e1= sp.Eq(rx[0]-r0*(rx[0])/(r_abs[0])-ge_x*rx[0]-ge_y*ry[0],bx)
    e2= sp.Eq(ry[0]-r0*(ry[0])/(r_abs[0])-ge_y*rx[0]+ge_x*ry[0],by)
    e3= sp.Eq(rx[1]-r0*(rx[1])/(r_abs[1])-ge_x*rx[1]-ge_y*ry[1],bx)
    e4= sp.Eq(ry[1]-r0*(ry[1])/(r_abs[1])-ge_y*rx[1]+ge_x*ry[1],by)

    k=1/2*r0/r_abs
    gsis_x= 1/2*r0*((-rx**2+ry**2)/(r_abs**3))
    gsis_y= -r0*((rx*ry)/(r_abs**3))
    gx= gsis_x+ge_x
    gy= gsis_y+ge_y
    ga=sp.sqrt(gx[0]**2+gy[0]**2)
    gb=sp.sqrt(gx[1]**2+gy[1]**2)
    
    e5= sp.Eq(f_ab*((1-k[0])**2-ga**2),(1-k[1])**2-gb**2)
    
    if solver_type=='numeric':
        sol=sp.nsolve((e1,e2,e3,e4,e5),(r0,ge_x,ge_y,bx,by),(np.max(r_abs)*1.5,0,0,igx,igy),dict=True)
    else:
        sol=sp.solve((e1,e2,e3,e4,e5),(r0,ge_x,ge_y,bx,by),dict=True)
    
    if sol[0][r0]==0: i=1
    else: i=0

    R0= sol[i][r0]
    GE_X= sol[i][ge_x]
    GE_Y= sol[i][ge_y]
    BX= sol[i][bx]
    BY= sol[i][by]

    tx, ty= sp.symbols('tx, ty', real=True)
    e1= sp.Eq(BX+R0*(tx/(sp.sqrt(tx**2+ty**2)))+GE_X*tx+GE_Y*ty-tx,0)
    e2= sp.Eq(BY+R0*(ty/(sp.sqrt(tx**2+ty**2)))-GE_X*ty+GE_Y*tx-ty,0)

    theta=np.zeros((2,n))
    for j in range(n):
            a=sp.nsolve((e1,e2),(tx,ty),(rx[j],ry[j]))
            theta[0,j]=a[0]
            theta[1,j]=a[1]


    return sol[i], theta

def Quad(ix,iy,igx,igy):

#Prep###############################################################################
    n=4
    rx=ix-igx
    ry=iy-igy
    r_abs= np.sqrt(rx**2+ry**2)

    r0,ge_x,ge_y,bx,by = sp.symbols('r0,ge_x,ge_y,bx,by', real=True) #Main symbols

    # define functions that will take the role of indexed symbols
    x= sp.Function('x')
    y= sp.Function('y')
    R_abs= sp.Function('R_abs')
    i= sp.symbols('i',posititve=True, real=True)

    x_fun = lambda index: rx[index-1]
    y_fun= lambda index: ry[index-1]
    R_abs_fun= lambda index: r_abs[index-1]
#Symbolic##################################################################
    s2 = sp.Sum((bx + r0*x(i)/R_abs(i) + ge_x*x(i) + ge_y*y(i) - x(i))**2 + (by + r0*y(i)/R_abs(i) - ge_x*y(i) + ge_y*x(i) - y(i))** 2, (i, 1, n))

    dbx= sp.diff(s2,bx).doit().replace(y,y_fun).replace(R_abs,R_abs_fun).replace(x,x_fun)
    dby= sp.diff(s2,by).doit().replace(y,y_fun).replace(R_abs,R_abs_fun).replace(x,x_fun)
    dgx= sp.diff(s2,ge_x).doit().replace(y,y_fun).replace(R_abs,R_abs_fun).replace(x,x_fun)
    dgy= sp.diff(s2,ge_y).doit().replace(y,y_fun).replace(R_abs,R_abs_fun).replace(x,x_fun)
    dr0= sp.diff(s2,r0).doit().replace(y,y_fun).replace(R_abs,R_abs_fun).replace(x,x_fun)
    sol=sp.solve((dbx,dby,dgx,dgy,dr0),(r0,ge_x,ge_y,bx,by), simplify=True)
    
    R0= sol[r0]
    GE_X= sol[ge_x]
    GE_Y= sol[ge_y]
    BX= sol[bx]
    BY= sol[by]

    tx, ty= sp.symbols('tx, ty', real=True)
    e1= sp.Eq(BX+R0*(tx/(sp.sqrt(tx**2+ty**2)))+GE_X*tx+GE_Y*ty-tx,0)
    e2= sp.Eq(BY+R0*(ty/(sp.sqrt(tx**2+ty**2)))-GE_X*ty+GE_Y*tx-ty,0)

    theta=np.zeros((2,n))
    try:
        for i in range(n):
            a=sp.nsolve((e1,e2),(tx,ty),(rx[i],ry[i]))
            theta[0,i]=a[0]
            theta[1,i]=a[1]
    except ValueError:
        print('No theta found within tolerance')

    return sol, theta
