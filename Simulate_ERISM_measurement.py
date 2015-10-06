###############################################################
### Script to simulate an ERISM measurement to test BioCalc ###
###############################################################

############### 
# not working #
###############

###########
## INPUT ##
###########

# takes a thickness image as input

# give the range of wavelengths (similar to ERISM measurement)
wave_start = 550
wave_end = 750
wave_step = 1


# give the range of thickness for which the reflectivity should be simulated, in nanometer
thickness_start = 7500
thickness_step = 1
thickness_N = 1000




    


####################################
## Code for Sim Reflection Minima ##
####################################

## implement wedge over thickness
## makes simulation data file to fit minimas of experimental microcavity data
import openTMM as tmm
import scipy as sp
import numpy as np

# For plotting only
import matplotlib.pyplot as plt 
import time



ta=time.time()

# define a class for the layers to easily access, material, thickness, optical constants
class Layer:    
    "A layer with thickness and permittivity"
    def __init__(self, material, thickness, e_r):
        self.material= material
        self.thickness=thickness*1e-9
        self.permittivity=e_r

    # function to read the n&k values from datafiles (space separated with  lines) and convert it to permitivity
    def e_r_read(self,waves):
        p = open('MaterialsDataFiles/'+self.material + '.txt','r')
        string=p.read()
        p.close()
        e_r=[]
        w=[]
        n=[]
        k=[]
        linecounter=0
        for thisline in string.split('\n'):
            for x in thisline.split():
                if linecounter==0 and len(x)>0:
                    w.append(float(x))
                if linecounter==1 and len(x)>0:
                    n.append(float(x))
                if linecounter==2 and len(x)>0:
                    k.append(float(x))
            linecounter+=1
        
        # interpolate n&k walues 
        n_new=sp.interp(waves,w,n)
        k_new=sp.interp(waves,w,k)
        e_r_new=[]      
        # calculate epsilon from n&k for every wavelength
        for i in range(len(waves)):
            e_r_new.append(n_new[i]**2-k_new[i]**2+2j*n_new[i]*k_new[i]) # calculate the complex epsilon from n&k
        e_r.append(waves)
        e_r.append(e_r_new) 
        
        self.permittivity=e_r


### Define Layer structure ###

## write in the order how your multilayer stack should look, the light will come from the side oft the last layer
## write it like: Layers.append(Layer('material', thickness, []) 
## --> thickness in nanometers

Layers=[]
# air
Layers.append(Layer('Au',15,[]))
Layers.append(Layer('Elastomer_RT601',3000,[]))
Layers.append(Layer('SIO2',50,[]))
Layers.append(Layer('Au',10,[]))
Layers.append(Layer('Cr',0.5,[]))
# substrate


### define host media, air, and glass

hm_eps = sp.array([1.5*1.5, 1.0]) 
hm_mu = sp.array([1.0,1.0])
hm_dict = {'epsilonRelative': hm_eps, 'muRelative': hm_mu}


# function to convert Frequencies to wavelengths
# takes Freq in "Hz" and gives waves in "nm"
def Freq_to_wave(f):
    cLight=299792458.0
    wave=cLight/f/1e-9
    return wave

# function to convert wavelengths Frequencies
# INPUT: waves in "nm", OUTPUT: Freq in "Hz"  
def wave_to_Freq(wave):
    wave=wave*1e-9
    cLight=299792458.0
    Freq=cLight/wave
    return Freq
        


#############################################################
## FUNCTION:                                               ##
## calculates the reflectivity as a function of wavelength ##
#############################################################

def tmm_waves(wave_min, wave_max):
    "calculates reflection and transmission for different wavelengths \n INPUT:\n 1. Minimum wavelength (in nm) \n 3. Maximum wavelength (in nm) \n\n OUTPUT:\n wavelength variation (list), transmission (list), reflection (list)"
    f_par=[]
    R_par=[]
    T_par=[]
    waves=[]
    i=0
    welle = wave_min
    while welle <= wave_max:
        waves.append(welle)
        i +=wave_step
        welle=wave_min+i

    #for i in range(wave_max-wave_min+1):
    #   waves.append(wave_min+i)
    for x in Layers:
        x.e_r_read(waves)
    for i in range(len(waves)):
        Freq_min=wave_to_Freq(waves[i])
        Freq_max=Freq_min
        eps=sp.ones(len(Layers),dtype=complex) # its important to generate a complex array, otherwise the complex permittivity will be reduced to its real part --> you get wrong results
        mu=sp.ones(len(Layers))
        height=sp.ones(len(Layers))
        for x in range(len(Layers)):
            eps[x]=Layers[x].permittivity[1][i]
            height[x]=Layers[x].thickness
        stack = {'epsilonRelative': eps, 'muRelative': mu, 'height': height, 'hostMedia': hm_dict }
        mls= tmm.Layer(stack)
        f1, T1, R1 = mls.TRvsFreq(Freq_min, Freq_max,1, 0, 'parallel')
        f_par.append(Freq_to_wave(f1[0]))
        T_par.append(T1[0])
        R_par.append(R1[0])

    return np.array(f_par), np.array(T_par), np.array(R_par)

########################
## create empty array ##
########################

ERISM = np.zeros(((wave_end-wave_start)/wave_step +1,10,thickness_N))
ERISM_line = np.zeros(((wave_end-wave_start)/wave_step +1,thickness_N))

for i in range(thickness_N):  
    thickness = thickness_start + i*thickness_step
    print thickness
    Layers[1].thickness=thickness*1e-9  
    f, T, R = tmm_waves(wave_start,wave_end) # set wavelenght range

    for j in range(ERISM.shape[1]):
        ERISM[:,j,i] = R*10000 # just scaling with 10000 to make values more comparable to measurement

for i in range(ERISM.shape[0]):
    file_name = str(wave_start+wave_step*i) + '.txt'
    np.savetxt(file_name,ERISM[i],fmt='%d')#,header=HEADER )


# thickness = 8000
# Layers[1].thickness=thickness*1e-9 
# f, T, R = tmm_waves(wave_start,wave_end)