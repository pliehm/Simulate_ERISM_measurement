# script to check how often a value is repeated in a row


# file, list of thicknesses
f = "simulated_thickness.txt"


import numpy as np


a = np.loadtxt(f)
counter = 0
counter_final = 0
thickness_stays = 0
for i in range(len(a[:,1])):
    if thickness_stays == a[:,1][i]:
        counter +=1
    if thickness_stays != a[:,1][i]:
        thickness_stays = a[:,1][i]
        if counter_final<counter:
            counter_final = counter
        counter = 0
    
        



b = np.array(range(7500,8500))

c = abs(a[:,1]-b)

print "Maximum thickness deviation: ",c.max()

print "Final Counter is: ", counter_final