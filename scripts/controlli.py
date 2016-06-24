import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import glob
import myTrasf as mT
plt.close('all')


path = os.path.dirname(os.path.abspath(__file__))+'//../'
path1 = os.path.dirname(os.path.abspath(__file__))+'//..//obs_16/' #cambia il verso delle slash su ubuntu!!!!
lista=glob.glob(path+'*.sky')
lista1=glob.glob(path1+'*.o')
#per leggere tutti il nome dei file con quel suffiso in una cartella(ho messo path altrimenti sto nella cartella del map.py che non ha i file che cerco)
#glob crea una lista di stringa che contiene tutti i nomi dei file aventi il suffisso scelto che si trovano nella cartella scelta
lista.sort()

#sod=28800 e PRN=8

Xs = -24451536.6367	                       
Ys=-177983.312534	  	
Zs=-10989692.8441                     

		#-24451536.6367	-177983.312534	-10989692.8441                   
hp=350000



#Xr = np.genfromtxt(lista[i], skip_header=1, usecols=(2))                  
#Yr= np.genfromtxt(lista[i], skip_header=1, usecols=(3))
#Zr= np.genfromtxt(lista[i], skip_header=1, usecols=(4))                    
#

#if 'approx' in line[]:
Xr=-5467762.7837   
Yr= -2518810.5051       
Zr= 2103349.5793 

X1=Xr
Y1=Yr
Z1=Zr

X2=Xs
Y2=Ys
Z2=Zs


counter=0
while True:
    counter=counter+1 #conta le iterazioni
    X3=(X1+X2)/2
    Y3=(Y1+Y2)/2
    Z3=(Z1+Z2)/2
    phi3,lamda3,h3=mT.coord_geog(X3,Y3,Z3)

    if abs(h3-hp)<1000:
        break
    if h3>hp:
        X2=X3 
        Y2=Y3
        Z2=Z3
    if h3<hp:
        X1=X3 
        Y1=Y3
        Z1=Z3        


#for i in xrange(0,Xm):
#    if h<350000:
#      Xm[i+1]=(Xm[i]+Xr)/2
#      Ym[i+1]=(Ym[i]+Yr)/2
#      Zm[i+1]=(Ym[i]+Zr)/2
#   # phi1,lamda1,h1=mT.coord_geog(Xm[i+1],Ym[i+1],Zm[i+1])
#    if h>350000:
#      Xm[i+1]=(Xm[i]+Xm[i-1])/2
#      Ym[i+1]=(Ym[i]+Ym[i-1])/2
#      Zm[i+1]=(Ym[i]+Zm[i-1])/2
#    if abs(h-350000)<1000:
#        break
#        
#        
#
#
#
#
#
#
#
#
#


















#dist = sqrt((Xs-Xr)^2 + (Ys-Yr)^2 + (Zs-Zr)^2)

