import numpy as np
#import simplekml

#####
def subiono(sky_file, lat_GPS, long_GPS):
        '''
        Function that returns the geodetic coordinats of the sub-ionospheric point.
        input:
            - name of the sky file
            - latitude of the GPS antenna
            - longitude of the GPS antenna
        output:
            - sod
            - sats name
            - latitude of the sIP
            - longitude of the sIP
            - elevation
            - azimuth
        '''        
        A_VAD   = np.genfromtxt(sky_file, skip_header=1, usecols=(3))
        a   = np.genfromtxt(sky_file, skip_header=1, usecols=(3))
        e   = np.genfromtxt(sky_file, skip_header=1, usecols=(4))
        ## The vadase output gives the azimut in different notation (negative value)
        mask_A = ( A_VAD < 0 )
        A_VAD[mask_A] = 360 + A_VAD[mask_A]
        
        A   = A_VAD / 180.0
        E   = np.genfromtxt(sky_file, skip_header=1, usecols=(4)) / 180.0 
        sow = np.genfromtxt(sky_file, skip_header=1, usecols=(0)) # second of the week
        # DEBUGG
        n_days = int(sow[0]) / 86400           # in order to understand which day of the week is
        sod = sow - 86400 * n_days             # second of the day
        #
        sat = np.genfromtxt(sky_file,dtype=str, skip_header=1, usecols=(2))
        
        # Convert to semicircles (SC)
        lat_GPS_sc  = lat_GPS / 180.0
        long_GPS_sc = long_GPS / 180.0
        
        G = (0.0137/(E+0.11)) - 0.022
        
        lat_sIP = lat_GPS_sc + G*np.cos(A*np.pi)
        
        lat_sIP_mask_1 = (lat_sIP < -0.416)
        lat_sIP[lat_sIP_mask_1] = -0.416
        lat_sIP_mask_2 = (lat_sIP > 0.416)
        lat_sIP[lat_sIP_mask_2] = 0.416        
        
        long_sIP = long_GPS_sc + ( G*np.sin(A*np.pi) ) / ( np.cos(lat_GPS_sc*np.pi) )
        
        return sod, sat, lat_sIP * 180.0 , long_sIP * 180.0, e, a
#####
def track(sIP, sat_name, cut_off=False, kml_file=False, txt_file=False):
    '''
    Function that returns the geodetic coordinats of a particular satellites.
    input:
        - sIP list
        - name of the satellite (e.g. 'G13')
        - kml_file (True if you want to generate a kml file)
    output:
        - sod
        - names
        - lat
        - lon
    '''    
    mask = ( sIP[1] == sat_name )
    sod   =  sIP[0][mask]
    names =  sIP[1][mask]
    lat   =  sIP[2][mask] 
    lon   =  sIP[3][mask] 
    ele   =  sIP[4][mask] 
    azi   =  sIP[5][mask] 
    
    if cut_off == True:
        mask_elev = ( ele  > 22.0)
        sow   = sow[mask_elev]
        names = names[mask_elev]
        lat   = lat[mask_elev]
        lon   = lon[mask_elev]
        ele   = ele[mask_elev]
        azi   = azi[mask_elev]
        
    if kml_file == True:
        kml=simplekml.Kml()
        for i in xrange (1,len(names)):
            kml.newpoint(name=names[i], coords=[(lon[i],lat[i])])   
        kml.save(sat_name + '.kml')
        
    if txt_file==True:
        f = open(sat_name + '.txt', 'w')
        f.write('sow' + '\t' + '\t' + 'longitude' + '\t'+ '\t' 'latitude' + '\n')
        for i in xrange (1,len(names)):
            f.write( str(sow[i]) + '\t'+ '\t'+ str(lon[i]) + '\t'+ '\t'+str(lat[i]) +'\n' )
        f.close()
        
    return sod, names, lat, lon, ele, azi
#####
def track_temp(sIP_sat, start, stop, kml_file):
    '''
    input:
        - sIP_sat, that is the output of the func track
        - start, in sow
        - stop, in sow
        - kml_file ( True or False )
    output:
        - sod 
        - names
        - lat
        - long
    '''  
    mask = ( sIP_sat[0] >= start ) & ( sIP_sat[0] <= stop )
    sod   =  sIP_sat[0][mask]
    names =  sIP_sat[1][mask]
    lat   =  sIP_sat[2][mask]
    lon   =  sIP_sat[3][mask]
    if kml_file == True:
        kml=simplekml.Kml()
        for i in xrange (1,len(names)):
            kml.newpoint(name=names[i], coords=[(lon[i],lat[i])])   
        kml.save(names[0] + 'interval.kml')
        
    return sod, names, lat, lon
#####
    