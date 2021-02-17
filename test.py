import numpy as np
import sympy as sp

#n=2
#rx=ix-igx
#ry=iy-igy
#dmag=mag[1]-mag[0]
#f_ab=1/(10**(dmag/-2.5))
#r_abs= np.sqrt(rx**2+ry**2)

r0,ge_x,ge_y,bx,by = sp.symbols('r0,ge_x,ge_y,bx,by', real=True)
rx1,rx2,ry1,ry2,r_abs1,r_abs2,f_ab=sp.symbols('rx1,rx2,ry1,ry2,r_abs1,r_abs2,f_ab', real=True)

#Symbolic##################################################################

dbx=sp.Eq(2*2*bx + 2*ge_x*rx1 + 2*ge_x*rx2 + 2*ge_y*ry1 + 2*ge_y*ry2 + 2*r0*rx2/r_abs2 + 2*r0*rx1/r_abs1 - 2*rx1 - 2*rx2)
dby=sp.Eq(2*2*by - 2*ge_x*ry1 - 2*ge_x*ry2 + 2*ge_y*rx1 + 2*ge_y*rx2 + 2*r0*ry2/r_abs2 + 2*r0*ry1/r_abs1 - 2*ry1 - 2*ry2)
dgx=sp.Eq(2*bx*rx1 + 2*bx*rx2 + 2*by*ry1 - 2*by*ry2 + f_ab**2*4*ge_x**3 - 2.0*f_ab**2*3*ge_x**2*r0*rx1**2/r_abs1**3 + 2.0*f_ab**2*3*ge_x**2*r0*ry1**2/r_abs1**3 + 2*f_ab**2*2*ge_x*ge_y**2 - 4*f_ab**2*2*ge_x*ge_y*r0*rx1*ry1/r_abs1**3 - 0.5*f_ab**2*2*ge_x*r0**2/r_abs1**2 + 1.5*f_ab**2*2*ge_x*r0**2*rx1**4/r_abs1**6 - 1.0*f_ab**2*2*ge_x*r0**2*rx1**2*ry1**2/r_abs1**6 + 1.5*f_ab**2*2*ge_x*r0**2*ry1**4/r_abs1**6 + 2.0*f_ab**2*2*ge_x*r0/r_abs1 - 2*f_ab**2*2*ge_x - 2.0*f_ab**2*ge_y**2*r0*rx1**2/r_abs1**3 + 2.0*f_ab**2*ge_y**2*r0*ry1**2/r_abs1**3 + 4.0*f_ab**2*ge_y*r0**2*rx1**3*ry1/r_abs1**6 - 4.0*f_ab**2*ge_y*r0**2*rx1*ry1**3/r_abs1**6 + 0.5*f_ab**2*r0**3*rx1**2/r_abs1**5 - 0.5*f_ab**2*r0**3*ry1**2/r_abs1**5 - 0.5*f_ab**2*r0**3*rx1**6/r_abs1**9 - 0.5*f_ab**2*r0**3*rx1**4*ry1**2/r_abs1**9 + 0.5*f_ab**2*r0**3*rx1**2*ry1**4/r_abs1**9 + 0.5*f_ab**2*r0**3*ry1**6/r_abs1**9 - 2.0*f_ab**2*r0**2*rx1**2/r_abs1**4 + 2.0*f_ab**2*r0**2*ry1**2/r_abs1**4 + 2.0*f_ab**2*r0*rx1**2/r_abs1**3 - 2.0*f_ab**2*r0*ry1**2/r_abs1**3 + 2*f_ab*4*ge_x**3 - 2.0*f_ab*3*ge_x**2*r0*rx2**2/r_abs2**3 + 2.0*f_ab*3*ge_x**2*r0*ry2**2/r_abs2**3 - 2.0*f_ab*3*ge_x**2*r0*rx1**2/r_abs1**3 + 2.0*f_ab*3*ge_x**2*r0*ry1**2/r_abs1**3 + 4*f_ab*2*ge_x*ge_y**2 - 4*f_ab*2*ge_x*ge_y*r0*rx2*ry2/r_abs2**3 - 4*f_ab*2*ge_x*ge_y*r0*rx1*ry1/r_abs1**3 + 0.5*f_ab*2*ge_x*r0**2/r_abs2**2 + 0.5*f_ab*2*ge_x*r0**2*rx2**4/r_abs2**6 + 1.0*f_ab*2*ge_x*r0**2*rx2**2*ry2**2/r_abs2**6 + 0.5*f_ab*2*ge_x*r0**2*ry2**4/r_abs2**6 - 0.5*f_ab*2*ge_x*r0**2/r_abs1**2 + 2.0*f_ab*2*ge_x*r0**2*rx1**2*rx2**2/(r_abs1**3*r_abs2**3) - 2.0*f_ab*2*ge_x*r0**2*rx1**2*ry2**2/(r_abs1**3*r_abs2**3) - 2.0*f_ab*2*ge_x*r0**2*rx2**2*ry1**2/(r_abs1**3*r_abs2**3) + 2.0*f_ab*2*ge_x*r0**2*ry1**2*ry2**2/(r_abs1**3*r_abs2**3) + 0.5*f_ab*2*ge_x*r0**2*rx1**4/r_abs1**6 + 1.0*f_ab*2*ge_x*r0**2*rx1**2*ry1**2/r_abs1**6 + 0.5*f_ab*2*ge_x*r0**2*ry1**4/r_abs1**6 - 2.0*f_ab*2*ge_x*r0/r_abs2 + 2.0*f_ab*2*ge_x*r0/r_abs1 - 2.0*f_ab*ge_y**2*r0*rx2**2/r_abs2**3 + 2.0*f_ab*ge_y**2*r0*ry2**2/r_abs2**3 - 2.0*f_ab*ge_y**2*r0*rx1**2/r_abs1**3 + 2.0*f_ab*ge_y**2*r0*ry1**2/r_abs1**3 + 4.0*f_ab*ge_y*r0**2*rx1**2*rx2*ry2/(r_abs1**3*r_abs2**3) + 4.0*f_ab*ge_y*r0**2*rx1*rx2**2*ry1/(r_abs1**3*r_abs2**3) - 4.0*f_ab*ge_y*r0**2*rx1*ry1*ry2**2/(r_abs1**3*r_abs2**3) - 4.0*f_ab*ge_y*r0**2*rx2*ry1**2*ry2/(r_abs1**3*r_abs2**3) + 0.5*f_ab*r0**3*rx2**2/(r_abs1**2*r_abs2**3) - 0.5*f_ab*r0**3*ry2**2/(r_abs1**2*r_abs2**3) - 0.5*f_ab*r0**3*rx1**2/(r_abs1**3*r_abs2**2) + 0.5*f_ab*r0**3*ry1**2/(r_abs1**3*r_abs2**2) - 0.5*f_ab*r0**3*rx1**2*rx2**4/(r_abs1**3*r_abs2**6) - 1.0*f_ab*r0**3*rx1**2*rx2**2*ry2**2/(r_abs1**3*r_abs2**6) - 0.5*f_ab*r0**3*rx1**2*ry2**4/(r_abs1**3*r_abs2**6) + 0.5*f_ab*r0**3*rx2**4*ry1**2/(r_abs1**3*r_abs2**6) + 1.0*f_ab*r0**3*rx2**2*ry1**2*ry2**2/(r_abs1**3*r_abs2**6) + 0.5*f_ab*r0**3*ry1**2*ry2**4/(r_abs1**3*r_abs2**6) - 0.5*f_ab*r0**3*rx1**4*rx2**2/(r_abs1**6*r_abs2**3) + 0.5*f_ab*r0**3*rx1**4*ry2**2/(r_abs1**6*r_abs2**3) - 1.0*f_ab*r0**3*rx1**2*rx2**2*ry1**2/(r_abs1**6*r_abs2**3) + 1.0*f_ab*r0**3*rx1**2*ry1**2*ry2**2/(r_abs1**6*r_abs2**3) - 0.5*f_ab*r0**3*rx2**2*ry1**4/(r_abs1**6*r_abs2**3) + 0.5*f_ab*r0**3*ry1**4*ry2**2/(r_abs1**6*r_abs2**3) - 2.0*f_ab*r0**2*rx2**2/(r_abs1*r_abs2**3) + 2.0*f_ab*r0**2*ry2**2/(r_abs1*r_abs2**3) + 2.0*f_ab*r0**2*rx1**2/(r_abs1**3*r_abs2) - 2.0*f_ab*r0**2*ry1**2/(r_abs1**3*r_abs2) + 2.0*f_ab*r0*rx2**2/r_abs2**3 - 2.0*f_ab*r0*ry2**2/r_abs2**3 - 2.0*f_ab*r0*rx1**2/r_abs1**3 + 2.0*f_ab*r0*ry1**2/r_abs1**3 + 4*ge_x**3 - 2.0*3*ge_x**2*r0*rx2**2/r_abs2**3 + 2.0*3*ge_x**2*r0*ry2**2/r_abs2**3 + 2*2*ge_x*ge_y**2 - 4*2*ge_x*ge_y*r0*rx2*ry2/r_abs2**3 + 0.5*2*ge_x*r0**2/r_abs2**2 + 1.5*2*ge_x*r0**2*rx2**4/r_abs2**6 - 1.0*2*ge_x*r0**2*rx2**2*ry2**2/r_abs2**6 + 1.5*2*ge_x*r0**2*ry2**4/r_abs2**6 - 2.0*2*ge_x*r0/r_abs2 + 2*ge_x*rx1**2 + 2*ge_x*rx2**2 + 2*ge_x*ry1**2 + 2*ge_x*ry2**2 + 2*2*ge_x - 2.0*ge_y**2*r0*rx2**2/r_abs2**3 + 2.0*ge_y**2*r0*ry2**2/r_abs2**3 + 4.0*ge_y*r0**2*rx2**3*ry2/r_abs2**6 - 4.0*ge_y*r0**2*rx2*ry2**3/r_abs2**6 - 0.5*r0**3*rx2**2/r_abs2**5 + 0.5*r0**3*ry2**2/r_abs2**5 - 0.5*r0**3*rx2**6/r_abs2**9 - 0.5*r0**3*rx2**4*ry2**2/r_abs2**9 + 0.5*r0**3*rx2**2*ry2**4/r_abs2**9 + 0.5*r0**3*ry2**6/r_abs2**9 + 2.0*r0**2*rx2**2/r_abs2**4 - 2.0*r0**2*ry2**2/r_abs2**4 + 2*r0*rx2**2/r_abs2 - 2*r0*ry2**2/r_abs2 - 2.0*r0*rx2**2/r_abs2**3 + 2.0*r0*ry2**2/r_abs2**3 + 2*r0*rx1**2/r_abs1 - 2*r0*ry1**2/r_abs1 - 2*rx1**2 - 2*rx2**2 + 2*ry1**2 + 2*ry2**2)
dgy=sp.Eq(2*bx*ry1 + 2*bx*ry2 + 2*by*rx1 + 2*by*rx2 + 2*f_ab**2*ge_x**2*2*ge_y - 4*f_ab**2*ge_x**2*r0*rx1*ry1/r_abs1**3 - 2.0*f_ab**2*ge_x*2*ge_y*r0*rx1**2/r_abs1**3 + 2.0*f_ab**2*ge_x*2*ge_y*r0*ry1**2/r_abs1**3 + 4.0*f_ab**2*ge_x*r0**2*rx1**3*ry1/r_abs1**6 - 4.0*f_ab**2*ge_x*r0**2*rx1*ry1**3/r_abs1**6 + f_ab**2*4*ge_y - 4*f_ab**2*3*ge_y*r0*rx1*ry1/r_abs1**3 - 0.5*f_ab**2*2*ge_y*r0**2/r_abs1**2 + 0.5*f_ab**2*2*ge_y*r0**2*rx1**4/r_abs1**6 + 5.0*f_ab**2*2*ge_y*r0**2*rx1**2*ry1**2/r_abs1**6 + 0.5*f_ab**2*2*ge_y*r0**2*ry1**4/r_abs1**6 + 2.0*f_ab**2*2*ge_y*r0/r_abs1 - 2*f_ab**2*2*ge_y + 1.0*f_ab**2*r0**3*rx1*ry1/r_abs1**5 - 1.0*f_ab**2*r0**3*rx1**5*ry1/r_abs1**9 - 2.0*f_ab**2*r0**3*rx1**3*ry1**3/r_abs1**9 - 1.0*f_ab**2*r0**3*rx1*ry1**5/r_abs1**9 - 4.0*f_ab**2*r0**2*rx1*ry1/r_abs1**4 + 4*f_ab**2*r0*rx1*ry1/r_abs1**3 + 4*f_ab*ge_x**2*2*ge_y - 4*f_ab*ge_x**2*r0*rx2*ry2/r_abs2**3 - 4*f_ab*ge_x**2*r0*rx1*ry1/r_abs1**3 + 2.0*f_ab*ge_x*2*ge_y*r0*rx2**2/r_abs2**3 + 2.0*f_ab*ge_x*2*ge_y*r0*ry2**2/r_abs2**3 - 2.0*f_ab*ge_x*2*ge_y*r0*rx1**2/r_abs1**3 + 2.0*f_ab*ge_x*2*ge_y*r0*ry1**2/r_abs1**3 + 4.0*f_ab*ge_x*r0**2*rx1**2*rx2*ry2/(r_abs1**3*r_abs2**3) + 4.0*f_ab*ge_x*r0**2*rx1*rx2**2*ry1/(r_abs1**3*r_abs2**3) - 4.0*f_ab*ge_x*r0**2*rx1*ry1*ry2**2/(r_abs1**3*r_abs2**3) - 4.0*f_ab*ge_x*r0**2*rx2*ry1**2*ry2/(r_abs1**3*r_abs2**3) + 2*f_ab*4*ge_y - 4*f_ab*3*ge_y*r0*rx2*ry2/r_abs2**3 - 4*f_ab*3*ge_y*r0*rx1*ry1/r_abs1**3 + 0.5*f_ab*2*ge_y*r0**2/r_abs2**2 + 0.5*f_ab*2*ge_y*r0**2*rx2**4/r_abs2**6 + 1.0*f_ab*2*ge_y*r0**2*rx2**2*ry2**2/r_abs2**6 + 0.5*f_ab*2*ge_y*r0**2*ry2**4/r_abs2**6 - 0.5*f_ab*2*ge_y*r0**2/r_abs1**2 + 8*f_ab*2*ge_y*r0**2*rx1*rx2*ry1*ry2/(r_abs1**3*r_abs2**3) + 0.5*f_ab*2*ge_y*r0**2*rx1**4/r_abs1**6 + 1.0*f_ab*2*ge_y*r0**2*rx1**2*ry1**2/r_abs1**6 + 0.5*f_ab*2*ge_y*r0**2*ry1**4/r_abs1**6 - 2.0*f_ab*2*ge_y*r0/r_abs2 + 2.0*f_ab*2*ge_y*r0/r_abs1 + 1.0*f_ab*r0**3*rx2*ry2/(r_abs1**2*r_abs2**3) - 1.0*f_ab*r0**3*rx1*ry1/(r_abs1**3*r_abs2**2) - 1.0*f_ab*r0**3*rx1*rx2**4*ry1/(r_abs1**3*r_abs2**6) - 2.0*f_ab*r0**3*rx1*rx2**2*ry1*ry2**2/(r_abs1**3*r_abs2**6) - 1.0*f_ab*r0**3*rx1*ry1*ry2**4/(r_abs1**3*r_abs2**6) - 1.0*f_ab*r0**3*rx1**4*rx2*ry2/(r_abs1**6*r_abs2**3) - 2.0*f_ab*r0**3*rx1**2*rx2*ry1**2*ry2/(r_abs1**6*r_abs2**3) - 1.0*f_ab*r0**3*rx2*ry1**4*ry2/(r_abs1**6*r_abs2**3) - 4.0*f_ab*r0**2*rx2*ry2/(r_abs1*r_abs2**3) + 4.0*f_ab*r0**2*rx1*ry1/(r_abs1**3*r_abs2) + 4*f_ab*r0*rx2*ry2/r_abs2**3 - 4*f_ab*r0*rx1*ry1/r_abs1**3 - 2*ge_x**2*2*ge_y - 4*ge_x**2*r0*rx2*ry2/r_abs2**3 + 2.0*ge_x*2*ge_y*r0*rx2**2/r_abs2**3 + 2.0*ge_x*2*ge_y*r0*ry2**2/r_abs2**3 + 4.0*ge_x*r0**2*rx2**3*ry2/r_abs2**6 - 4.0*ge_x*r0**2*rx2*ry2**3/r_abs2**6 - 4*ge_y - 4*3*ge_y*r0*rx2*ry2/r_abs2**3 + 0.5*2*ge_y*r0**2/r_abs2**2 + 0.5*2*ge_y*r0**2*rx2**4/r_abs2**6 + 5.0*2*ge_y*r0**2*rx2**2*ry2**2/r_abs2**6 + 0.5*2*ge_y*r0**2*ry2**4/r_abs2**6 - 2.0*2*ge_y*r0/r_abs2 + 2*ge_y*rx1**2 + 2*ge_y*rx2**2 + 2*ge_y*ry1**2 + 2*ge_y*ry2**2 + 2*2*ge_y - 1.0*r0**3*rx2*ry2/r_abs2**5 - 1.0*r0**3*rx2**5*ry2/r_abs2**9 - 2.0*r0**3*rx2**3*ry2**3/r_abs2**9 - 1.0*r0**3*rx2*ry2**5/r_abs2**9 + 4.0*r0**2*rx2*ry2/r_abs2**4 + 4*r0*rx2*ry2/r_abs2 - 4*r0*rx2*ry2/r_abs2**3 + 4*r0*rx1*ry1/r_abs1 - 4*rx1*ry1 - 4*rx2*ry2)
dr0=sp.Eq(2*bx*rx2/r_abs2 + 2*bx*rx1/r_abs1 - 2*by*ry2/r_abs2 + 2*by*ry1/r_abs1 - 2.0*f_ab**2*ge_x**3*rx1**2/r_abs1**3 + 2.0*f_ab**2*ge_x**3*ry1**2/r_abs1**3 + 4*f_ab**2*ge_x**2*ge_y*rx1*ry1/r_abs1**3 - 0.5*f_ab**2*ge_x**2*2*r0/r_abs1**2 + 1.5*f_ab**2*ge_x**2*2*r0*rx1**4/r_abs1**6 - 1.0*f_ab**2*ge_x**2*2*r0*rx1**2*ry1**2/r_abs1**6 + 1.5*f_ab**2*ge_x**2*2*r0*ry1**4/r_abs1**6 + 2.0*f_ab**2*ge_x**2/r_abs1 - 2.0*f_ab**2*ge_x*ge_y**2*rx1**2/r_abs1**3 + 2.0*f_ab**2*ge_x*ge_y**2*ry1**2/r_abs1**3 + 4.0*f_ab**2*ge_x*ge_y*2*r0*rx1**3*ry1/r_abs1**6 - 4.0*f_ab**2*ge_x*ge_y*2*r0*rx1*ry1**3/r_abs1**6 + 0.5*f_ab**2*ge_x*3*r0**2*rx1**2/r_abs1**5 - 0.5*f_ab**2*ge_x*3*r0**2*ry1**2/r_abs1**5 - 0.5*f_ab**2*ge_x*3*r0**2*rx1**6/r_abs1**9 - 0.5*f_ab**2*ge_x*3*r0**2*rx1**4*ry1**2/r_abs1**9 + 0.5*f_ab**2*ge_x*3*r0**2*rx1**2*ry1**4/r_abs1**9 + 0.5*f_ab**2*ge_x*3*r0**2*ry1**6/r_abs1**9 - 2.0*f_ab**2*ge_x*2*r0*rx1**2/r_abs1**4 + 2.0*f_ab**2*ge_x*2*r0*ry1**2/r_abs1**4 + 2.0*f_ab**2*ge_x*rx1**2/r_abs1**3 - 2.0*f_ab**2*ge_x*ry1**2/r_abs1**3 + 4*f_ab**2*ge_y**3*rx1*ry1/r_abs1**3 - 0.5*f_ab**2*ge_y**2*2*r0/r_abs1**2 + 0.5*f_ab**2*ge_y**2*2*r0*rx1**4/r_abs1**6 + 5.0*f_ab**2*ge_y**2*2*r0*rx1**2*ry1**2/r_abs1**6 + 0.5*f_ab**2*ge_y**2*2*r0*ry1**4/r_abs1**6 + 2.0*f_ab**2*ge_y**2/r_abs1 - 1.0*f_ab**2*ge_y*3*r0**2*rx1*ry1/r_abs1**5 - 1.0*f_ab**2*ge_y*3*r0**2*rx1**5*ry1/r_abs1**9 - 2.0*f_ab**2*ge_y*3*r0**2*rx1**3*ry1**3/r_abs1**9 - 1.0*f_ab**2*ge_y*3*r0**2*rx1*ry1**5/r_abs1**9 - 4.0*f_ab**2*ge_y*2*r0*rx1*ry1/r_abs1**4 + 4*f_ab**2*ge_y*rx1*ry1/r_abs1**3 + 0.0625*f_ab**2*4*r0**3/r_abs1**4 - 0.125*f_ab**2*4*r0**3*rx1**4/r_abs1**8 - 0.25*f_ab**2*4*r0**3*rx1**2*ry1**2/r_abs1**8 - 0.125*f_ab**2*4*r0**3*ry1**4/r_abs1**8 + 0.0625*f_ab**2*4*r0**3*rx1**8/r_abs1**12 + 0.25*f_ab**2*4*r0**3*rx1**6*ry1**2/r_abs1**12 + 0.375*f_ab**2*4*r0**3*rx1**4*ry1**4/r_abs1**12 + 0.25*f_ab**2*4*r0**3*rx1**2*ry1**6/r_abs1**12 + 0.0625*f_ab**2*4*r0**3*ry1**8/r_abs1**12 - 0.5*f_ab**2*3*r0**2/r_abs1**3 + 0.5*f_ab**2*3*r0**2*rx1**4/r_abs1**7 + 1.0*f_ab**2*3*r0**2*rx1**2*ry1**2/r_abs1**7 + 0.5*f_ab**2*3*r0**2*ry1**4/r_abs1**7 + 1.5*f_ab**2*2*r0/r_abs1**2 - 0.5*f_ab**2*2*r0*rx1**4/r_abs1**6 - 1.0*f_ab**2*2*r0*rx1**2*ry1**2/r_abs1**6 - 0.5*f_ab**2*2*r0*ry1**4/r_abs1**6 - 2.0*f_ab**2/r_abs1 + 2.0*f_ab*ge_x**3*rx2**2/r_abs2**3 + 2.0*f_ab*ge_x**3*ry2**2/r_abs2**3 - 2.0*f_ab*ge_x**3*rx1**2/r_abs1**3 + 2.0*f_ab*ge_x**3*ry1**2/r_abs1**3 + 4*f_ab*ge_x**2*ge_y*rx2*ry2/r_abs2**3 - 4*f_ab*ge_x**2*ge_y*rx1*ry1/r_abs1**3 + 0.5*f_ab*ge_x**2*2*r0/r_abs2**2 + 0.5*f_ab*ge_x**2*2*r0*rx2**4/r_abs2**6 + 1.0*f_ab*ge_x**2*2*r0*rx2**2*ry2**2/r_abs2**6 + 0.5*f_ab*ge_x**2*2*r0*ry2**4/r_abs2**6 - 0.5*f_ab*ge_x**2*2*r0/r_abs1**2 + 2.0*f_ab*ge_x**2*2*r0*rx1**2*rx2**2/(r_abs1**3*r_abs2**3) - 2.0*f_ab*ge_x**2*2*r0*rx1**2*ry2**2/(r_abs1**3*r_abs2**3) - 2.0*f_ab*ge_x**2*2*r0*rx2**2*ry1**2/(r_abs1**3*r_abs2**3) + 2.0*f_ab*ge_x**2*2*r0*ry1**2*ry2**2/(r_abs1**3*r_abs2**3) + 0.5*f_ab*ge_x**2*2*r0*rx1**4/r_abs1**6 + 1.0*f_ab*ge_x**2*2*r0*rx1**2*ry1**2/r_abs1**6 + 0.5*f_ab*ge_x**2*2*r0*ry1**4/r_abs1**6 - 2.0*f_ab*ge_x**2/r_abs2 + 2.0*f_ab*ge_x**2/r_abs1 - 2.0*f_ab*ge_x*ge_y**2*rx2**2/r_abs2**3 + 2.0*f_ab*ge_x*ge_y**2*ry2**2/r_abs2**3 - 2.0*f_ab*ge_x*ge_y**2*rx1**2/r_abs1**3 + 2.0*f_ab*ge_x*ge_y**2*ry1**2/r_abs1**3 + 4.0*f_ab*ge_x*ge_y*2*r0*rx1**2*rx2*ry2/(r_abs1**3*r_abs2**3) + 4.0*f_ab*ge_x*ge_y*2*r0*rx1*rx2**2*ry1/(r_abs1**3*r_abs2**3) - 4.0*f_ab*ge_x*ge_y*2*r0*rx1*ry1*ry2**2/(r_abs1**3*r_abs2**3) - 4.0*f_ab*ge_x*ge_y*2*r0*rx2*ry1**2*ry2/(r_abs1**3*r_abs2**3) + 0.5*f_ab*ge_x*3*r0**2*rx2**2/(r_abs1**2*r_abs2**3) - 0.5*f_ab*ge_x*3*r0**2*ry2**2/(r_abs1**2*r_abs2**3) - 0.5*f_ab*ge_x*3*r0**2*rx1**2/(r_abs1**3*r_abs2**2) + 0.5*f_ab*ge_x*3*r0**2*ry1**2/(r_abs1**3*r_abs2**2) - 0.5*f_ab*ge_x*3*r0**2*rx1**2*rx2**4/(r_abs1**3*r_abs2**6) - 1.0*f_ab*ge_x*3*r0**2*rx1**2*rx2**2*ry2**2/(r_abs1**3*r_abs2**6) - 0.5*f_ab*ge_x*3*r0**2*rx1**2*ry2**4/(r_abs1**3*r_abs2**6) + 0.5*f_ab*ge_x*3*r0**2*rx2**4*ry1**2/(r_abs1**3*r_abs2**6) + 1.0*f_ab*ge_x*3*r0**2*rx2**2*ry1**2*ry2**2/(r_abs1**3*r_abs2**6) + 0.5*f_ab*ge_x*3*r0**2*ry1**2*ry2**4/(r_abs1**3*r_abs2**6) - 0.5*f_ab*ge_x*3*r0**2*rx1**4*rx2**2/(r_abs1**6*r_abs2**3) + 0.5*f_ab*ge_x*3*r0**2*rx1**4*ry2**2/(r_abs1**6*r_abs2**3) - 1.0*f_ab*ge_x*3*r0**2*rx1**2*rx2**2*ry1**2/(r_abs1**6*r_abs2**3) + 1.0*f_ab*ge_x*3*r0**2*rx1**2*ry1**2*ry2**2/(r_abs1**6*r_abs2**3) - 0.5*f_ab*ge_x*3*r0**2*rx2**2*ry1**4/(r_abs1**6*r_abs2**3) + 0.5*f_ab*ge_x*3*r0**2*ry1**4*ry2**2/(r_abs1**6*r_abs2**3) - 2.0*f_ab*ge_x*2*r0*rx2**2/(r_abs1*r_abs2**3) + 2.0*f_ab*ge_x*2*r0*ry2**2/(r_abs1*r_abs2**3) + 2.0*f_ab*ge_x*2*r0*rx1**2/(r_abs1**3*r_abs2) - 2.0*f_ab*ge_x*2*r0*ry1**2/(r_abs1**3*r_abs2) + 2.0*f_ab*ge_x*rx2**2/r_abs2**3 - 2.0*f_ab*ge_x*ry2**2/r_abs2**3 - 2.0*f_ab*ge_x*rx1**2/r_abs1**3 + 2.0*f_ab*ge_x*ry1**2/r_abs1**3 + 4*f_ab*ge_y**3*rx2*ry2/r_abs2**3 - 4*f_ab*ge_y**3*rx1*ry1/r_abs1**3 + 0.5*f_ab*ge_y**2*2*r0/r_abs2**2 + 0.5*f_ab*ge_y**2*2*r0*rx2**4/r_abs2**6 + 1.0*f_ab*ge_y**2*2*r0*rx2**2*ry2**2/r_abs2**6 + 0.5*f_ab*ge_y**2*2*r0*ry2**4/r_abs2**6 - 0.5*f_ab*ge_y**2*2*r0/r_abs1**2 + 8*f_ab*ge_y**2*2*r0*rx1*rx2*ry1*ry2/(r_abs1**3*r_abs2**3) + 0.5*f_ab*ge_y**2*2*r0*rx1**4/r_abs1**6 + 1.0*f_ab*ge_y**2*2*r0*rx1**2*ry1**2/r_abs1**6 + 0.5*f_ab*ge_y**2*2*r0*ry1**4/r_abs1**6 - 2.0*f_ab*ge_y**2/r_abs2 + 2.0*f_ab*ge_y**2/r_abs1 + 1.0*f_ab*ge_y*3*r0**2*rx2*ry2/(r_abs1**2*r_abs2**3) - 1.0*f_ab*ge_y*3*r0**2*rx1*ry1/(r_abs1**3*r_abs2**2) - 1.0*f_ab*ge_y*3*r0**2*rx1*rx2**4*ry1/(r_abs1**3*r_abs2**6) - 2.0*f_ab*ge_y*3*r0**2*rx1*rx2**2*ry1*ry2**2/(r_abs1**3*r_abs2**6) - 1.0*f_ab*ge_y*3*r0**2*rx1*ry1*ry2**4/(r_abs1**3*r_abs2**6) - 1.0*f_ab*ge_y*3*r0**2*rx1**4*rx2*ry2/(r_abs1**6*r_abs2**3) - 2.0*f_ab*ge_y*3*r0**2*rx1**2*rx2*ry1**2*ry2/(r_abs1**6*r_abs2**3) - 1.0*f_ab*ge_y*3*r0**2*rx2*ry1**4*ry2/(r_abs1**6*r_abs2**3) - 4.0*f_ab*ge_y*2*r0*rx2*ry2/(r_abs1*r_abs2**3) + 4.0*f_ab*ge_y*2*r0*rx1*ry1/(r_abs1**3*r_abs2) + 4*f_ab*ge_y*rx2*ry2/r_abs2**3 - 4*f_ab*ge_y*rx1*ry1/r_abs1**3 - 0.125*f_ab*4*r0**3/(r_abs1**2*r_abs2**2) - 0.125*f_ab*4*r0**3*rx2**4/(r_abs1**2*r_abs2**6) - 0.25*f_ab*4*r0**3*rx2**2*ry2**2/(r_abs1**2*r_abs2**6) - 0.125*f_ab*4*r0**3*ry2**4/(r_abs1**2*r_abs2**6) + 0.125*f_ab*4*r0**3*rx1**4/(r_abs1**6*r_abs2**2) + 0.25*f_ab*4*r0**3*rx1**2*ry1**2/(r_abs1**6*r_abs2**2) + 0.125*f_ab*4*r0**3*ry1**4/(r_abs1**6*r_abs2**2) + 0.125*f_ab*4*r0**3*rx1**4*rx2**4/(r_abs1**6*r_abs2**6) + 0.25*f_ab*4*r0**3*rx1**4*rx2**2*ry2**2/(r_abs1**6*r_abs2**6) + 0.125*f_ab*4*r0**3*rx1**4*ry2**4/(r_abs1**6*r_abs2**6) + 0.25*f_ab*4*r0**3*rx1**2*rx2**4*ry1**2/(r_abs1**6*r_abs2**6) + 0.5*f_ab*4*r0**3*rx1**2*rx2**2*ry1**2*ry2**2/(r_abs1**6*r_abs2**6) + 0.25*f_ab*4*r0**3*rx1**2*ry1**2*ry2**4/(r_abs1**6*r_abs2**6) + 0.125*f_ab*4*r0**3*rx2**4*ry1**4/(r_abs1**6*r_abs2**6) + 0.25*f_ab*4*r0**3*rx2**2*ry1**4*ry2**2/(r_abs1**6*r_abs2**6) + 0.125*f_ab*4*r0**3*ry1**4*ry2**4/(r_abs1**6*r_abs2**6) + 0.5*f_ab*3*r0**2/(r_abs1*r_abs2**2) + 0.5*f_ab*3*r0**2*rx2**4/(r_abs1*r_abs2**6) + 1.0*f_ab*3*r0**2*rx2**2*ry2**2/(r_abs1*r_abs2**6) + 0.5*f_ab*3*r0**2*ry2**4/(r_abs1*r_abs2**6) + 0.5*f_ab*3*r0**2/(r_abs1**2*r_abs2) - 0.5*f_ab*3*r0**2*rx1**4/(r_abs1**6*r_abs2) - 1.0*f_ab*3*r0**2*rx1**2*ry1**2/(r_abs1**6*r_abs2) - 0.5*f_ab*3*r0**2*ry1**4/(r_abs1**6*r_abs2) - 0.5*f_ab*2*r0/r_abs2**2 - 0.5*f_ab*2*r0*rx2**4/r_abs2**6 - 1.0*f_ab*2*r0*rx2**2*ry2**2/r_abs2**6 - 0.5*f_ab*2*r0*ry2**4/r_abs2**6 - 2.0*f_ab*2*r0/(r_abs1*r_abs2) - 0.5*f_ab*2*r0/r_abs1**2 + 0.5*f_ab*2*r0*rx1**4/r_abs1**6 + 1.0*f_ab*2*r0*rx1**2*ry1**2/r_abs1**6 + 0.5*f_ab*2*r0*ry1**4/r_abs1**6 + 2.0*f_ab/r_abs2 + 2.0*f_ab/r_abs1 - 2.0*ge_x**3*rx2**2/r_abs2**3 + 2.0*ge_x**3*ry2**2/r_abs2**3 + 4*ge_x**2*ge_y*rx2*ry2/r_abs2**3 + 0.5*ge_x**2*2*r0/r_abs2**2 + 1.5*ge_x**2*2*r0*rx2**4/r_abs2**6 - 1.0*ge_x**2*2*r0*rx2**2*ry2**2/r_abs2**6 + 1.5*ge_x**2*2*r0*ry2**4/r_abs2**6 - 2.0*ge_x**2/r_abs2 + 2.0*ge_x*ge_y**2*rx2**2/r_abs2**3 + 2.0*ge_x*ge_y**2*ry2**2/r_abs2**3 + 4.0*ge_x*ge_y*2*r0*rx2**3*ry2/r_abs2**6 - 4.0*ge_x*ge_y*2*r0*rx2*ry2**3/r_abs2**6 - 0.5*ge_x*3*r0**2*rx2**2/r_abs2**5 + 0.5*ge_x*3*r0**2*ry2**2/r_abs2**5 - 0.5*ge_x*3*r0**2*rx2**6/r_abs2**9 - 0.5*ge_x*3*r0**2*rx2**4*ry2**2/r_abs2**9 + 0.5*ge_x*3*r0**2*rx2**2*ry2**4/r_abs2**9 + 0.5*ge_x*3*r0**2*ry2**6/r_abs2**9 + 2.0*ge_x*2*r0*rx2**2/r_abs2**4 - 2.0*ge_x*2*r0*ry2**2/r_abs2**4 + 2*ge_x*rx2**2/r_abs2 - 2*ge_x*ry2**2/r_abs2 - 2.0*ge_x*rx2**2/r_abs2**3 + 2.0*ge_x*ry2**2/r_abs2**3 + 2*ge_x*rx1**2/r_abs1 - 2*ge_x*ry1**2/r_abs1 - 4*ge_y**3*rx2*ry2/r_abs2**3 + 0.5*ge_y**2*2*r0/r_abs2**2 + 0.5*ge_y**2*2*r0*rx2**4/r_abs2**6 + 5.0*ge_y**2*2*r0*rx2**2*ry2**2/r_abs2**6 + 0.5*ge_y**2*2*r0*ry2**4/r_abs2**6 - 2.0*ge_y**2/r_abs2 + 1.0*ge_y*3*r0**2*rx2*ry2/r_abs2**5 - 1.0*ge_y*3*r0**2*rx2**5*ry2/r_abs2**9 - 2.0*ge_y*3*r0**2*rx2**3*ry2**3/r_abs2**9 - 1.0*ge_y*3*r0**2*rx2*ry2**5/r_abs2**9 + 4.0*ge_y*2*r0*rx2*ry2/r_abs2**4 + 4*ge_y*rx2*ry2/r_abs2 - 4*ge_y*rx2*ry2/r_abs2**3 + 4*ge_y*rx1*ry1/r_abs1 - 0.0625*4*r0**3/r_abs2**4 + 0.125*4*r0**3*rx2**4/r_abs2**8 + 0.25*4*r0**3*rx2**2*ry2**2/r_abs2**8 + 0.125*4*r0**3*ry2**4/r_abs2**8 + 0.0625*4*r0**3*rx2**8/r_abs2**12 + 0.25*4*r0**3*rx2**6*ry2**2/r_abs2**12 + 0.375*4*r0**3*rx2**4*ry2**4/r_abs2**12 + 0.25*4*r0**3*rx2**2*ry2**6/r_abs2**12 + 0.0625*4*r0**3*ry2**8/r_abs2**12 - 0.5*3*r0**2/r_abs2**3 - 0.5*3*r0**2*rx2**4/r_abs2**7 - 1.0*3*r0**2*rx2**2*ry2**2/r_abs2**7 - 0.5*3*r0**2*ry2**4/r_abs2**7 + 2*r0*rx2**2/r_abs2**2 + 2*r0*ry2**2/r_abs2**2 + 1.5*2*r0/r_abs2**2 + 0.5*2*r0*rx2**4/r_abs2**6 + 1.0*2*r0*rx2**2*ry2**2/r_abs2**6 + 0.5*2*r0*ry2**4/r_abs2**6 + 2*r0*rx1**2/r_abs1**2 + 2*r0*ry1**2/r_abs1**2 - 2*rx2**2/r_abs2 - 2*ry2**2/r_abs2 - 2.0/r_abs2 - 2*rx1**2/r_abs1 - 2*ry1**2/r_abs1)

gx=sp.solveset(dbx,ge_x).args[0]

dby=dby.subs(ge_x,gx)
gy=sp.solveset(dby,ge_y).args[0]

dgx=dgx.subs({ge_x:gx,ge_y:gy})
#ux=sp.solveset(dgx,bx)

dgy=dgy.subs({ge_x:gx,ge_y:gy})

dr0=dr0.subs({ge_x:gx,ge_y:gy})


rAx = 0
rAy = 0
rBx = -1.508
rBy = 2.068
rGx = -1.661
rGy = 1.472
img=np.array([[rAx,rAy],[rBx,rBy]])
lens=np.array([rGx,rGy])
fAB = 8.62979
coord=img-lens
r_abs=np.sqrt(np.sum(coord**2,1))


dr0=dr0.subs({rx1:coord[0,0],rx2:coord[1,0],ry1:coord[0,1],ry2:coord[1,1],r_abs1:r_abs[0],r_abs2:r_abs[1],f_ab:fAB})
R0=sp.solveset(dr0,r0).args[0]

dgx=dgx.subs({rx1:coord[0,0],rx2:coord[1,0],ry1:coord[0,1],ry2:coord[1,1],r_abs1:r_abs[0],r_abs2:r_abs[1],f_ab:fAB})
ux=sp.solveset(dgx,bx).args[0]
ux=ux.subs(r0,R0)

dgy=dgy.subs({rx1:coord[0,0],rx2:coord[1,0],ry1:coord[0,1],ry2:coord[1,1],r_abs1:r_abs[0],r_abs2:r_abs[1],f_ab:fAB})
uy=sp.solveset(dgy,by).args[0]
uy=uy.subs(r0,R0)

ööö

k1=1/2*r0/r_abs1
k2=1/2*r0/r_abs2

gsis_x1= 1/2*r0*((-rx1**2+ry1**2)/(r_abs1**3))
gsis_x2= 1/2*r0*((-rx2**2+ry2**2)/(r_abs2**3))

gsis_y1= -r0*((rx1*ry1)/(r_abs1**3))
gsis_y2= -r0*((rx2*ry2)/(r_abs2**3))

gx1= gsis_x1+ge_x
gx2= gsis_x2+ge_x

gy1= gsis_y1+ge_y
gy2= gsis_y2+ge_y

ga=sp.sqrt(gx1**2+gy1**2)
gb=sp.sqrt(gx2**2+gy2**2)

e6=sp.Eq((rx1-r0*(rx1)/(r_abs1)-ge_x*rx1-ge_y*ry1-bx)**2+
         (ry1-r0*(ry1)/(r_abs1)-ge_y*rx1+ge_x*ry1-by)**2+
         (rx2-r0*(rx2)/(r_abs2)-ge_x*rx2-ge_y*ry2-bx)**2+
         (ry2-r0*(ry2)/(r_abs2)-ge_y*rx2+ge_x*ry2-by)**2+
         (f_ab*((1-k1)**2-ga**2)-(1-k2)**2-gb**2)**2)



dbx= sp.Eq(4*bx + 3.628*ge_x - 0.672*ge_y + 1.76376155287604*r0 - 3.628)#sp.diff(e6,bx).doit()
dby= sp.Eq(4*by + 0.672*ge_x + 3.628*ge_y + 0.655614714808487*r0 + 0.672)#sp.diff(e6,by).doit()
dgx= sp.Eq(3.628*bx + 0.672*by + 92.7328554441*4*ge_x**3 + 7.20086997925168*3*ge_x**2*r0 + 185.4657108882*2*ge_x*ge_y**2 + 69.8997242195574*2*ge_x*ge_y*r0 + 7.4689272577033*2*ge_x*r0**2 + 58.0859559673521*2*ge_x*r0 - 140.7069408882*2*ge_x + 7.20086997925168*ge_y**2*r0 + 2.71391850968074*ge_y*r0**2 + 0.284560222558586*r0**3 + 2.25523852650896*r0**2 - 7.38255232413207*r0 + 1.3499)#sp.diff(e6,ge_x).doit()
dgy= sp.Eq(- 0.672*bx + 3.628*by + 185.4657108882*ge_x**2*2*ge_y + 69.8997242195574*ge_x**2*r0 + 7.20086997925168*ge_x*2*ge_y*r0 + 2.71391850968074*ge_x*r0**2 + 92.7328554441*4*ge_y**3 + 69.8997242195574*3*ge_y**2*r0 + 20.5013063814814*2*ge_y*r0**2 + 58.0859559673521*2*ge_y*r0 - 140.7069408882*2*ge_y + 2.7622608293183*r0**3 + 21.891875774249*r0**2 - 59.1824032530127*r0 + 9.084736)#sp.diff(e6,ge_y).doit()
dr0= sp.Eq(1.76376155287604*bx + 0.655614714808487*by*r0+ 7.20086997925168*ge_x**3 + 69.8997242195574*ge_x**2*ge_y + 7.4689272577033*ge_x**2*2*r0 + 58.0859559673521*ge_x**2 + 7.20086997925168*ge_x*ge_y**2 + 2.71391850968074*ge_x*ge_y*2*r0 + 0.284560222558586*ge_x*3*r0**2 + 2.25523852650896*ge_x*2*r0 - 7.38255232413207*ge_x + 69.8997242195574*ge_y**3 + 20.5013063814814*ge_y**2*2*r0 + 58.0859559673521*ge_y**2 + 2.7622608293183*ge_y*3*r0**2 + 21.891875774249*ge_y*2*r0 - 59.1824032530127*ge_y + 0.144814509638343*4*r0**3 + 2.29541049973459*3*r0**2 + 5.28900531694313*2*r0 - 52.7534504107622)#sp.diff(e6,r0).doit()
    

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



