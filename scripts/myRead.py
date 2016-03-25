import numpy as np
import re
################################################################################
def iono(RINEX_file, sat):
    """
    This function read a RINEX file.
    input:
        - RINEX file
        - name of the sat (e.g. "G18")
    outputs:
        - sod array
        - L1 array
        - L2 array
    """

    L1 = 1.57542e9    # Hz
    L2 = 1.22760e9    # Hz
    c  = 299792458    # m/s

    lam1=c/L1         # m
    lam2=c/L2         # m


    f   = open(RINEX_file, "r")
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
        
    return sod_arr, l1_arr, l2_arr