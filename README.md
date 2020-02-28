# VARION

## HOW TO USE VARION

Download your RINEX files and copy them in the "/obs" dir. Put in the same folder also the navigation file (e.g. brdc????.??n).
Then move with the terminal into the "/script" dir and execute the VARION script.


## LINUX
        $ ./VARION.py -h      --->     to run the help
        $ ./VARION.py -staz ahup ainp -time 08:00 09:30 -sat G04 G07 G08 G10
        $ ./VARION.py -staz ahup ainp -time 08:00 09:30 -sat G04 G07 -brdc    ---> to use brdc file
        $ ./VARION.py -time 08:00 09:30 -brdc -height 450       ---> change the height of the Ionosperic layer

## OS
        $ python VARION_next.py -staz ahup ainp -time 08:00 09:30 -sat G04 G07 G08 G10

## WINDOWS
        $ C:\Python27\python.exe C:\Users\Username\Desktop\my_python_script.py -h
        $ C:\Python27\python.exe C:\Users\Username\Desktop\my_python_script.py -time 08:00 09:30 -brdc

### Requirements ###

- Python 2.7+ (other versions might work but have not been tested)
- Numpy
- Pandas

### References ###

Savastano, G. et al. Real-Time Detection of Tsunami Ionospheric Disturbances with a Stand-Alone GNSS Receiver: 
A Preliminary Feasibility Demonstration. Sci. Rep. 7, 46607; doi: 10.1038/srep46607 (2017)

### Contacts ###

michela.ravanelli@uniroma1.it or giorgio.savastano@gmail.com