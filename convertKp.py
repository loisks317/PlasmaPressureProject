# convertKp.py
#
# take the Kyoto website output
# copy it as a text file
# read it in here into two arrays
# one is dates
# one is kp extended to 1 minute throughout the day
#
# LKS, SASNA January
#
import numpy as np
import pickle
#
# read in file line by line
file= open('KpFeb2013_Apr2015', 'rb')
Kps= {}
for line in file:
    try:
        date=str(int(line[0:8]))
        tempK=[]
        for i in range(8):
            tempK.extend([int(line[7+(i*2):8+(i*2)])]*180)
        Kps[date]=tempK
    except(ValueError):
        pass

# now pickle it
pickle.dump(Kps, open('KpFeb2013_Apr2015.p', 'wb'))
