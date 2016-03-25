import numpy as np
import myObs as mO
import re
################################################################################
def iono(lista_obs, sat):
    """
    input:
        - list of observation files
        - name of the sat (e.g. "G18")
    output a list with:
        - lista_sat[i][0] the sod of the i file
        - lista_sat[i][1] the Iono_vad of the i file
    """

    L1 = 1.57542e9    # Hz
    L2 = 1.22760e9    # Hz
    c  = 299792458    # m/s

    lam1=c/L1         # m
    lam2=c/L2         # m

    lista_sat = []
    for obs in lista_obs:
        f   = open(obs, "r")
        lns = f.readlines()
        f.close()
        
        sats = []
        ora  = []
        sod  = []
        l1   = []
        l2   = []

        s_sats= 0
        count = 0
        while "END OF HEADER" not in lns[count]:
            count += 1
        count += 1
        for i in xrange(count):
            if "TYPES OF OBSERV" in lns[i]:
                obs_ord = re.findall('\S+\S+',lns[i])
                obs_orb_arr = np.array(obs_ord)
                
                c1_mas = (obs_orb_arr == "C1")
                l1_mas = (obs_orb_arr == "L1")
                l2_mas = (obs_orb_arr == "L2")

                index = np.arange(len(obs_orb_arr))
                
                c1_index = index[c1_mas][0]
                l1_index = index[l1_mas][0]
                l2_index = index[l2_mas][0]
                break

        ## skip the comment inside the file
        for i in xrange(count,len(lns)):  
            if "COMMENT" in lns[i]:
                continue 
            else:
                pass  
        # START stocking the obs                           
            ln = lns[i].split()                                                 # split the line[i]
            if len(ln) == 8:                                                   # check the number of the split elements, if equal to 8 is the line with "9G27G01G07G32G..."
                ## count how many S sats there are
                if "S" in ln[-1]:
                    s_sats=0
                    for s in ln[-1]:
                        if s == "S":
                            s_sats+=1
                        else:
                            pass
                                                       
                if (len(ln[-1]) - 1) % 3 == 0:                                 # if the last element (eg 9G27G01G...), it means that n_sats < 10, skip only the first character of the element
                        n_sats = int(ln[-1][0])                                # the first character is the number of the satellites observed in that epoch e.g 9
                        for k in xrange(1,(n_sats-s_sats)*2+1,2):              # e.g. k in [1,3,5,7,9,...,13] I am skipping the last satellites because thew are always
                            
                            obs = re.findall('\S+\.\S+',lns[i+k])          # use regular expression to select only the float number
                        
                            if len(obs) <= 3: continue
                            ##
                            l1.append(float(obs[l1_index])*lam1)
                            l2.append(float(obs[l2_index])*lam2)                     # In this RINEX format L1 and L2 are the first two numbers
                            
                            sats.append(ln[-1][1:len(ln[-1])][(k)/2*3:(k)/2*3+3])      # skip only the first character (slicing)
                            ora.append(ln[3] + ":" + ln[4] + ":" + ln[5][0:2])
                            sod.append(int(ln[3])*60*60 + int(ln[4])*60 + float(ln[5][0:3]))
                              
                elif (len(ln[-1]) - 2) % 3 == 0:                                   # skip the first two characters because the sat number is > 10
                        n_sats = int(ln[-1][0:2])
                        
                        if n_sats <= 12:
                            for k in xrange(1,(n_sats-s_sats)*2+1,2):                  # if n_sats=10, than k = [1,3,5,...,19]
                                        
                                obs = re.findall('\S+\.\S+',lns[i+k])              # use regular expression to select only the float number
                                
                                if len(obs) <= 3: continue
                                ##
                                l1.append(  float(obs[l1_index])*lam1  )
                                l2.append(  float(obs[l2_index])*lam2  )                        # In this RINEX format L1 and L2 are the first two numbers
                                
                                sats.append(ln[-1][2:len(ln[-1])][(k)/2*3:(k)/2*3+3])     # skip the first two elements
                                ora.append(ln[3] + ":" + ln[4] + ":" + ln[5][0:2])
                                sod.append(int(ln[3])*60*60 + int(ln[4])*60 + float(ln[5][0:3]))
                                
                        else:                                                      # Check the number of sats, if >12 use continuation lines
                            ln_2 = lns[i+1].split()                                # splitto anche la riga dopo
                            for k in xrange(1,(n_sats-s_sats)*2+1,2):                      # if n_sats=12, than k = [1,3,5,...,23]
                                
                                obs = re.findall('\S+\.\S+',lns[i+1+k])          # use regular expression to select only the float number
                                
                                if len( obs ) <= 3: continue
                                l1.append(  float(obs[l1_index])*lam1  )
                                l2.append(  float(obs[l2_index])*lam2  )
                                        
                                ora.append(ln[3] + ":" + ln[4] + ":" + ln[5][0:2])
                                sod.append(int(ln[3])*60*60 + int(ln[4])*60 + float(ln[5][0:3]))
                                        
                                if k <= 23:
                                            sats.append(ln[-1][2:len(ln[-1])][(k)/2*3:(k)/2*3+3])     # skip the first two elements
                                else:
                                            sats.append(ln_2[0][((k+1)/2-13)*3:((k+1)/2-13)*3+3])
                                              
        sats_arr = np.asarray(sats)
        ora_arr  = np.asarray(ora)
        sod_arr  = np.asarray(sod)
        l1_arr = np.asarray(l1)
        l2_arr = np.asarray(l2)
        
        G18_VAD, G18_ora, G18_sod = mO.obs_sat(sats_arr, ora_arr, sod_arr, l1_arr, l2_arr, sat)
        lista_sat.append([G18_sod,G18_VAD, G18_ora])
        
    return lista_sat
################################################################################
def interval(obs_file):
    '''
    Function that returns the inteval of the obs.
    '''
    f   = open(obs_file, "r")
    lns = f.readlines()
    f.close()    
    flag = 0
    count = 0
    while "END OF HEADER" not in lns[count]:
        count += 1
    count += 1
        
    for s in xrange(0, count):
        if "INTERVAL" in lns[s]:
            inter = re.findall('\S+\.\S+',lns[s])
            interval = float(inter[0])  
            flag=1
            break  
    if flag == 0:
        sod = []
        for i in xrange(count,len(lns)):
            if "COMMENT" in lns[i]:
                continue
            else:
                pass
            ln = lns[i].split()
            if len(ln) == 8:
                sod.append(int(ln[3])*60*60 + int(ln[4])*60 + float(ln[5][0:3]))
        interval + sod[1] - sod[0]
    return interval
################################################################################
def location(obs_file):
    '''
    Function that reads the coordinates of the station.
    input:
        - obs_file
    output:
        -lat,lon
    '''
    f   = open(obs_file, "r")
    lns = f.readlines()
    f.close()    

    count = 0
    while "END OF HEADER" not in lns[count]:
        count += 1
    count += 1
        
    for s in xrange(0, count):
        if "Monument location:" in lns[s]:
            coord = re.findall('\S+\.\S+',lns[s])
            lat = float(coord[0]) 
            lon = float(coord[1])   
            break
    return lat,lon
###############################################################################    
def coord_geog(staz):
    
    '''
	input: x,y,z
	output: Latitude (F) and Longitude (L) WGS84
    '''

    filename1 = staz
    f1=open(filename1, "r")
    lines_1 = f1.readlines()
    f1.close()
    # lettura automatica delle coordinate approssimate del ricevitore
    for i in xrange(0,35):
            if 'COMMENT' in lines_1[i]:
                continue 
            elif 'APPROX' in lines_1[i]:
                x = float(lines_1[i][1:14])
                y = float(lines_1[i][15:29])
                z = float(lines_1[i][29:43])
    	        # Longitudine calcolabile senza iterazioni
    	        L = np.arctan2(y,x)

    	        r = (x**2 + y**2)**0.5
    	        # parametri WGS84
    	        a = 6378137
    	        f = 1/298.257223563
    	        b = a*(1-f)
    	        e = (1-(b**2)/(a**2))**0.5

       	        # I step
       	        # ip h=0
       	        F  = np.arctan2(z, r*(1-e**2))
       	        Rn = a/((1-(e**2)*(np.sin(F))**2))**0.5
       	        # Proviamo 50 iterazioni
       	        for i in xrange(0,50):
      		        h  = r/np.cos(F) - Rn
      		        F  = np.arctan2((z*(Rn+h)),(r*Rn*(1-e**2)+h))
      		        Rn = a/((1-(e**2)*(np.sin(F))**2))**0.5
      		
       	        L_grad = L/(np.pi)*180
                F_grad = F/(np.pi)*180
                break
    return F_grad, L_grad 