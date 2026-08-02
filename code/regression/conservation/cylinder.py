from math import sqrt

L = 640.0e+3 #domain half width/length

rf = 256.0e+3
H = 1000.0
D = -500.0
umax = 500.0

def get_r(x,y):
    return sqrt((x-L)**2+(y-L)**2)

def velocity(x,y,t,thk,*etc):

    vx, vy = 0.0, 0.0
    if (thk > 1.0e-10):
        r_cos_theta = x - L
        r_sin_theta = y - L
        vx = umax/rf * r_cos_theta
        vy = umax/rf * r_sin_theta
    return (vx,vy)

def thickness(x,y):
    if (x-L)**2 + (y-L)**2 > rf**2:
        return 0.0
    return H

def topography(x,y):
    return D

def constfriction(x,y,t,thck,topg):
    return 1.0e+4

def acab(x,y,*etc):
    r = get_r(x,y)
    a = 2*H*umax/rf
    a = 0.3
    #if ((x-L)**2 + (y-L)**2) > rf**2:
    #    return -1000000.0
    return a

def fsrs(x,y,t,*etc):
    
    uc = 0.0
    if (t >= 10.0):
        uc = 1.0 # steady
    if (t >= 20.0):
        uc = 2.0 # retreat
    if (t >= 30.0):
        uc = 1.0 # steady
    if (t >= 40.0):
        uc = 0.0 # free advance

    return uc

