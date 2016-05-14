import numpy as np
import re

obs_file = 'ahup3020_test.12o'
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
        flag = 1
        break
        
if flag == 0:
    sod = []
    for i in xrange(count,len(lns)):  
        if "COMMENT" in lns[i]:
            continue 
        else:
            pass  
    # START stocking the obs                           
        ln = lns[i].split()                         # split the line[i]
        if len(ln) == 8:                           # check the number of the split elements, if equal to 8 is the line with "9G27G01G07G32G...
            sod.append(int(ln[3])*60*60 + int(ln[4])*60 + float(ln[5][0:3])) 
sod_array = np.asarray(sod)             
interval = sod_array[1:len(sod_array)] - sod_array[0:len(sod_array)-1]  