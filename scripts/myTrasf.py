import numpy as np

def coord_geog(x,y,z):
    '''
    input:  X,Y,Z
    output: Latitude (F) and Longitude (L) and h WGS84
    '''
	# Longitudine calcolabile senza iterazioni
    L = np.arctan2(y,x)
    r = (x**2 + y**2)**0.5
    # WGS84 parameters
    a = 6378137
    f = 1/298.257223563
    b = a*(1-f)
    e = (1-(b**2)/(a**2))**0.5
    # I step
    # ip h=0
    h1 = 0.0
    F  = np.arctan2(z, r*(1-e**2))
    Rn = a/((1-(e**2)*(np.sin(F))**2))**0.5
    # Max 100 iterations --> usually break after 4 iterations
    for i in xrange(0,100):
            h  = r/np.cos(F) - Rn
            if abs(h - h1) <= 0.00001: break
            F  = np.arctan2((z*(Rn+h)),(Rn*(1-e**2)+h)*r)
            Rn = a/((1-(e**2)*(np.sin(F))**2))**0.5
            h1 = h
    L_grad = L/(np.pi)*180
    F_grad = F/(np.pi)*180
    return F_grad, L_grad, h
