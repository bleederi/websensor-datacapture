#
# Copyright 2017 Jesse Nieminen
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import json
import numpy
import math
from collections import OrderedDict
import pickle
from helperfuncs import prettyprint, find_interval, convert_orientation

__author__ = "Jesse Nieminen"
__status__ = "Development"

def remove_initial(jsondata):   #Do operations on JSON formatted data
        cancelledjsondata = [] #New JSON data with the effect of initial position and orientation cancelled
        for buttonpress in jsondata:
                if(buttonpress['button'] !=  None):     #skip bad values
                        cancelleddict = {}      #Dict (JSON data) for a button press with the effect of initial position and orientation cancelled
                        for key, value in buttonpress.items():
                                if key == 'button' or key == 'frequency':
                                        cancelleddict[key] = value
                                #Cancel out effect of initial position and orientation
                                if key != 'button' and key != 'frequency':
                                        i = 0
                                        cancelleddict[key] = []
                                        while(value[i] == None):
                                                i = i + 1
                                        initialvalue = value[i]
                                        if(key == 'acceleration'):
                                                for accvalue in value:
                                                        acc_c = {key: accvalue[key] - initialvalue.get(key, 0) for key in accvalue.keys()}
                                                        cancelleddict[key].append(acc_c)
                                        if(key == 'accelerationnog'):
                                                for accnogvalue in value:
                                                        accnog_c = {key: accnogvalue[key] - initialvalue.get(key, 0) for key in accnogvalue.keys()}
                                                        cancelleddict[key].append(accnog_c)
                                        if(key == 'orientation'):
                                                for orientationvalue in value:
                                                        orientationd = {key: orientationvalue[key] - initialvalue.get(key, 0) for key in orientationvalue.keys()}
                                                        cancelleddict[key].append(orientationd)           
                                        if(key == 'rotation'):
                                                for rotationvalue in value[i:]:
                                                        rotationd = {key: rotationvalue[key] - initialvalue.get(key, 0) for key in rotationvalue.keys()}
                                                        cancelleddict[key].append(rotationd)           
                        cancelledjsondata.append(cancelleddict)
        return cancelledjsondata

def calc_derivative(jsondata): #Function to calculate derivative sequence of each sequence and add that to the JSON formatted data
        jsondata_withderivative = [] #New JSON data with the derivative sequences added
        for buttonpress in jsondata:
                if(buttonpress['button'] !=  None):     #skip bad values
                        derivativedict = {}      #Dict (JSON data) for a button press with the derivative sequences added
                        #Add derivative sequences
                        for key, value in buttonpress.items():
                                derivativedict[key] = value     #Recreate the original dict, below add the derivative sequences
                                if(key == 'acceleration' or key == 'accelerationnog' or key == 'orientation' or key == 'rotation'):
                                        derivativedict[key + '_d'] = []   #Derivative sequence
                                        for i in value:
                                                if(value.index(i) == 0): #First value special case: always 0
                                                        der = {key: 0 for key in i.keys()}
                                                else:
                                                        der = {key: i[key] - previous.get(key, 0) for key in i.keys()}
                                                previous = der
                                                derivativedict[key + '_d'].append(der)
                        jsondata_withderivative.append(derivativedict)
        return jsondata_withderivative

def calc_dac(jsondata): #Function to calculate DAC of acceleration data and add that to the JSON formatted data
        jsondata_withdac = []
        for buttonpress in jsondata:
                if(buttonpress['button'] !=  None):     #skip bad values
                        dacdict = {}      #Dict (JSON data) for a button press with the DAC sequence added
                        #Add derivative sequences
                        for key, value in buttonpress.items():
                                dacdict[key] = value     #Recreate the original dict, below add the DAC sequence
                                if(key == 'acceleration'):
                                        dacdict['dac'] = []   #DAC sequence
                                        for accvalue in value:
                                                if(value.index(accvalue) == 0): #First value special case: always 0
                                                        dac = 0
                                                else:
                                                        #Calculate Euclidean distance of current-previous
                                                        dac = math.sqrt(math.pow(accvalue['x']-previous['x'], 2) + math.pow(accvalue['y']-previous['y'], 2) + math.pow(accvalue['z']-previous['z'], 2))
                                                previous = accvalue
                                                dacdict['dac'].append(dac)
                        jsondata_withdac.append(dacdict)
        return jsondata_withdac

def calc_stats(jsondata):       #Function to obtain statistical features of the button press data
        statsdict = []     #List for the statistics of a button press
        for buttonpress in jsondata:
                if(buttonpress['button'] !=  None):     #skip bad values
                        buttonstats = {}        #Dict of all statistics related to one button press
                        for key, value in buttonpress.items():
                                buttonstats[key] = value
                                if key != 'button' and key != 'frequency':
                                        #Calculate maximums
                                        buttonstats[key + '_max'] = {}
                                        if (key == 'acceleration' or key == 'accelerationnog' or key == 'rotation' or key == 'acceleration_d' or key == 'accelerationnog_d' or key == 'rotation_d'):                                                       
                                                max_x = value[0]['x']
                                                max_y = value[0]['y']
                                                max_z = value[0]['z']
                                                for i in value: #i is a sensor reading
                                                        if(value[value.index(i)]['x'] > max_x):
                                                                max_x = value[value.index(i)]['x']
                                                        if(value[value.index(i)]['y'] > max_y):
                                                                max_y = value[value.index(i)]['y']
                                                        if(value[value.index(i)]['z'] > max_z):
                                                                max_z = value[value.index(i)]['z']
                                                maxdict = {'x':max_x, 'y': max_y, 'z': max_z}    #Dict to store max values of a single sequence
                                        #Need to handle orientation and DAC separately (different structure)
                                        if (key == 'orientation' or key == 'orientation_d'):
                                                max_alpha = value[0]['alpha']
                                                max_beta = value[0]['beta']
                                                max_gamma = value[0]['gamma']
                                                for i in value: #i is a sensor reading
                                                        if(value[value.index(i)]['alpha'] > max_alpha):
                                                                max_alpha = value[value.index(i)]['alpha']
                                                        if(value[value.index(i)]['beta'] > max_beta):
                                                                max_beta = value[value.index(i)]['beta']
                                                        if(value[value.index(i)]['gamma'] > max_gamma):
                                                                max_gamma = value[value.index(i)]['gamma']
                                                maxdict = {'alpha':max_alpha, 'beta': max_beta, 'gamma': max_gamma}    #Dict to store max values of a single sequence
                                        if (key == 'dac'):
                                                max_dac = max(value)
                                                maxdict = max_dac    #Dict to store max values of a single sequence
                                        buttonstats[key + '_max'] = maxdict
                                        #Calculate minimums
                                        buttonstats[key + '_min'] = {}
                                        if (key == 'acceleration' or key == 'accelerationnog' or key == 'rotation' or key == 'acceleration_d' or key == 'accelerationnog_d' or key == 'rotation_d'):                                                      
                                                min_x = value[0]['x']
                                                min_y = value[0]['y']
                                                min_z = value[0]['z']
                                                for i in value: #i is a sensor reading
                                                        if(value[value.index(i)]['x'] < min_x):
                                                                min_x = value[value.index(i)]['x']
                                                        if(value[value.index(i)]['y'] < min_y):
                                                                min_y = value[value.index(i)]['y']
                                                        if(value[value.index(i)]['z'] < min_z):
                                                                min_z = value[value.index(i)]['z']
                                                mindict = {'x':min_x, 'y': min_y, 'z': min_z}    #Dict to store min values of a single sequence
                                        #Need to handle orientation and DAC separately (different structure)
                                        if (key == 'orientation' or key == 'orientation_d'):
                                                min_alpha = value[0]['alpha']
                                                min_beta = value[0]['beta']
                                                min_gamma = value[0]['gamma']
                                                for i in value: #i is a sensor reading
                                                        if(value[value.index(i)]['alpha'] < min_alpha):
                                                                min_alpha = value[value.index(i)]['alpha']
                                                        if(value[value.index(i)]['beta'] < min_beta):
                                                                min_beta = value[value.index(i)]['beta']
                                                        if(value[value.index(i)]['gamma'] < min_gamma):
                                                                min_gamma = value[value.index(i)]['gamma']
                                                mindict = {'alpha':min_alpha, 'beta': min_beta, 'gamma': min_gamma}    #Dict to store min values of a single sequence
                                        if (key == 'dac'):
                                                min_dac = min(value)
                                                mindict = min_dac    #Dict to store min values of a single sequence
                                        buttonstats[key + '_min'] = mindict
                                        #Calculate averages
                                        buttonstats[key + '_avg'] = {}
                                        if (key == 'acceleration' or key == 'accelerationnog' or key == 'rotation' or key == 'acceleration_d' or key == 'accelerationnog_d' or key == 'rotation_d'):                                                      
                                                sum_x = 0
                                                sum_y = 0
                                                sum_z = 0
                                                for i in value: #i is a sensor reading
                                                                sum_x = sum_x + value[value.index(i)]['x']
                                                                sum_y = sum_y + value[value.index(i)]['y']
                                                                sum_z = sum_z + value[value.index(i)]['z']
                                                avgdict = {'x':sum_x/len(value), 'y': sum_y/len(value), 'z': sum_z/len(value)}    #Dict to store min values of a single sequence
                                        #Need to handle orientation and DAC separately (different structure)
                                        if (key == 'orientation' or key == 'orientation_d'):
                                                sum_alpha = 0
                                                sum_beta = 0
                                                sum_gamma = 0
                                                for i in value: #i is a sensor reading
                                                                sum_alpha = sum_alpha + value[value.index(i)]['alpha']
                                                                sum_beta = sum_beta + value[value.index(i)]['beta']
                                                                sum_gamma = sum_gamma + value[value.index(i)]['gamma']
                                                avgdict = {'alpha':sum_alpha/len(value), 'beta': sum_beta/len(value), 'gamma': sum_gamma/len(value)}    #Dict to store min values of a single sequence
                                        if (key == 'dac'):
                                                sum_dac = sum(value)/len(value)
                                                avgdict = sum_dac    #Dict to store avg values of a single sequence
                                        buttonstats[key + '_avg'] = avgdict
                        statsdict.append(buttonstats)
        return statsdict

def calc_total_energy(jsondata):     #Function to calculate the total energy of each sequence and add that to the JSON formatted data   
        jsondata_withenergy = []
        for buttonpress in jsondata:
                if(buttonpress['button'] !=  None):     #skip bad values
                        energydict = {}      #Dict (JSON data) for a button press with the total energy sequence added
                        for key, value in buttonpress.items():
                                energydict[key] = value     #Recreate the original dict, below add the energy
                                if key == 'button' or key == 'frequency' or key.endswith('_avg') or key.endswith('_max') or key.endswith('_min'):
                                        continue
                                if (key.startswith('acceleration') or key.startswith('rotation') or key.startswith('orientation') or key == 'dac'):
                                        energydict[key + '_e'] = {}   #Total energy
                                        if (key.startswith('acceleration') or key.startswith('rotation')):
                                                #Should do with list comprehension later
                                                e_x = math.pow(value[0]['x'], 2)
                                                e_y = math.pow(value[0]['x'], 2)
                                                e_z = math.pow(value[0]['x'], 2)
                                                for i in value: #i is a sensor reading
                                                                e_x = e_x + math.pow(value[value.index(i)]['x'], 2)
                                                                e_y = e_y + math.pow(value[value.index(i)]['y'], 2)
                                                                e_z = e_z + math.pow(value[value.index(i)]['z'], 2)
                                                e_dict = {'x':e_x, 'y': e_y, 'z': e_z} 
                                #Need to handle orientation and DAC separately (different structure)
                                        if (key.startswith('orientation')):
                                                e_alpha = 0
                                                e_beta = 0
                                                e_gamma = 0
                                                for i in value: #i is a sensor reading
                                                                e_alpha = e_alpha + math.pow(value[value.index(i)]['alpha'], 2)
                                                                e_beta = e_beta + math.pow(value[value.index(i)]['beta'], 2)
                                                                e_gamma = e_gamma + math.pow(value[value.index(i)]['gamma'], 2)
                                                e_dict = {'alpha':e_alpha, 'beta': e_beta, 'gamma': e_gamma} 
                                        if (key == 'dac'):
                                                e_dac = sum( i*i for i in value)
                                                e_dict = e_dac
                                        energydict[key + '_e'] = e_dict                            
                jsondata_withenergy.append(energydict)
        return jsondata_withenergy

def calc_fft(jsondata): #Function to calculate FFT of each sequence and add that to the data. This function returns NON-JSON-FORMATTED DATA because the FFTs are complex.
        jsondata_withfft = []
        for buttonpress in jsondata:
                if(buttonpress['button'] !=  None):     #skip bad values
                        fftdict = {}      #Dict (JSON data) for a button press with the FFTs added
                        for key, value in buttonpress.items():
                                fftdict[key] = value     #Recreate the original dict, below add the FFTs
                                if (key == 'button' or key == 'frequency' or key.endswith('_avg') or key.endswith('_max') or key.endswith('_min') or key.endswith('_e') or key.endswith('_d') or key == 'dac'):
                                        continue
                                else:
                                        fftdict[key + '_fft'] = {}
                                        if (key.startswith('acceleration') or key.startswith('rotation')):
                                                #Below make the sequences to be FFTed
                                                fft_array_x = []
                                                fft_array_y = []
                                                fft_array_z = []
                                                for i in value:
                                                        fft_array_x.append(value[value.index(i)]['x'])
                                                        fft_array_y.append(value[value.index(i)]['y'])
                                                        fft_array_z.append(value[value.index(i)]['z'])
                                                fft_x = numpy.fft.fft(fft_array_x)
                                                fft_y = numpy.fft.fft(fft_array_y)
                                                fft_z = numpy.fft.fft(fft_array_z)
                                                f_dict = {'x':fft_x, 'y': fft_y, 'z': fft_z}
                                #Need to handle orientation separately (different structure)
                                        if (key.startswith('orientation')):
                                                #Below make the sequences to be FFTed
                                                fft_array_alpha = []
                                                fft_array_beta = []
                                                fft_array_gamma = []
                                                for i in value:
                                                        fft_array_alpha.append(value[value.index(i)]['alpha'])
                                                        fft_array_beta.append(value[value.index(i)]['beta'])
                                                        fft_array_gamma.append(value[value.index(i)]['gamma'])
                                                fft_alpha = numpy.fft.fft(fft_array_alpha)
                                                fft_beta = numpy.fft.fft(fft_array_beta)
                                                fft_gamma = numpy.fft.fft(fft_array_gamma)
                                                f_dict = {'alpha':fft_alpha, 'beta': fft_beta, 'gamma': fft_gamma} 
                                        fftdict[key + '_fft'] = f_dict             
                jsondata_withfft.append(fftdict)
        return jsondata_withfft

def calc_stats_fft(jsondata):   #Function to calculate maximum, minimum, mean and energy of the FFTs
        jsondata_withstats = []
        for buttonpress in jsondata:
                if(buttonpress['button'] !=  None):     #skip bad values
                        buttonstats = {}       #Dict of all statistics related to one button press
                        for key, value in buttonpress.items():
                                buttonstats[key] = value
                                if key.endswith('_fft'):
                                        #Calculate maximums, minimums, averages and energy
                                        buttonstats[key + '_min'] = {}
                                        buttonstats[key + '_max'] = {}
                                        buttonstats[key + '_avg'] = {}
                                        buttonstats[key + '_e'] = {}
                                        if (key.startswith('acceleration') or key.startswith('rotation')):
                                                fft_x = value['x']   #Array holding the FFT
                                                fft_y = value['y']
                                                fft_z = value['z']
                                                max_x = numpy.amax(fft_x)    
                                                max_y = numpy.amax(fft_y)  
                                                max_z = numpy.amax(fft_z)   
                                                min_x = numpy.amin(fft_x)    
                                                min_y = numpy.amin(fft_y)  
                                                min_z = numpy.amin(fft_z)
                                                avg_x = numpy.average(fft_x)    
                                                avg_y = numpy.average(fft_y)  
                                                avg_z = numpy.average(fft_z)
                                                #With complex numbers, energy = F.*conj(F)
                                                e_x = (1/len(fft_x)) * numpy.sum(fft_x * numpy.conj(fft_x))
                                                e_y = (1/len(fft_y)) *numpy.sum(fft_y * numpy.conj(fft_y))
                                                e_z = (1/len(fft_z)) *numpy.sum(fft_z * numpy.conj(fft_z))
                                                maxdict = {'x':max_x, 'y': max_y, 'z': max_z}    #Dict to store max values of a single sequence
                                                mindict = {'x':min_x, 'y': min_y, 'z': min_z}    #Dict to store min values of a single sequence
                                                avgdict = {'x':avg_x, 'y': avg_y, 'z': avg_z}    #Dict to store avg values of a single sequence
                                                edict = {'x':e_x, 'y': e_y, 'z': e_z}    #Dict to store energy of a single sequence
                                        #Need to handle orientation separately (different structure)
                                        if (key.startswith('orientation')):
                                                fft_alpha = value['alpha'] #Array holding the FFT
                                                fft_beta = value['beta']
                                                fft_gamma = value['gamma']
                                                max_alpha = numpy.amax(fft_alpha)
                                                max_beta = numpy.amax(fft_beta)
                                                max_gamma = numpy.amax(fft_gamma)
                                                min_alpha = numpy.amin(fft_alpha)
                                                min_beta = numpy.amin(fft_beta)
                                                min_gamma = numpy.amin(fft_gamma)
                                                avg_alpha = numpy.average(fft_alpha)    
                                                avg_beta = numpy.average(fft_beta)  
                                                avg_gamma = numpy.average(fft_gamma)
                                                e_alpha = (1/len(fft_alpha)) * numpy.sum(fft_alpha * numpy.conj(fft_alpha))
                                                e_beta = (1/len(fft_beta)) * numpy.sum(fft_beta * numpy.conj(fft_beta))
                                                e_gamma = (1/len(fft_gamma)) * numpy.sum(fft_gamma * numpy.conj(fft_gamma))
                                                maxdict = {'alpha':max_alpha, 'beta': max_beta, 'gamma': max_gamma}    #Dict to store max values of a single sequence
                                                mindict = {'alpha':min_alpha, 'beta': min_beta, 'gamma': min_gamma}    #Dict to store max values of a single sequence
                                                avgdict = {'alpha':avg_alpha, 'beta': avg_beta, 'gamma': avg_gamma}    #Dict to store avg values of a single sequence
                                                edict = {'alpha':e_alpha, 'beta': e_beta, 'gamma': e_gamma}    #Dict to store energy of a single sequence
                                        #Append the sequences - Note, we don't need the original FFTs anymore
                                        buttonstats[key + '_max'] = maxdict 
                                        buttonstats[key + '_min'] = mindict
                                        buttonstats[key + '_avg'] = avgdict
                                        buttonstats[key + '_e'] = edict
                jsondata_withstats.append(buttonstats)
        return jsondata_withstats

def calc_interval(jsondata): #Calculates the interval of 70% energy of the acceleration sequences and adds that to the data
        jsondata_with_interval = []
        for buttonpress in jsondata:
                if(buttonpress['button'] !=  None):     #skip bad values
                        buttondictwithinterval = {}        #Dict of the original data with interval data related for a single button press
                        for key, value in buttonpress.items():
                                buttondictwithinterval[key] = value
                                #Now select keys to calculate CoE for
                                if (key == 'button' or key == 'frequency' or key.endswith('_avg') or key.endswith('_max') or key.endswith('_min') or key.endswith('_e') or key.endswith('_d') or key == 'dac' or key.startswith('orientation') or key.startswith('rotation') or key.endswith('_fft')):
                                        continue
                                else:
                                        buttondictwithinterval[key + '_interval'] = {}                                     
                                        totalEnergy = buttonpress[key + '_e']
                                        index = 1       #Does the sum begin at i=0 or i=1?
                                        CoE_x = 0
                                        CoE_y = 0
                                        CoE_z = 0
                                        if(type(value) is dict):
                                                #print(key)
                                                for i, k in value.items(): #For numpy arrays, i is a dict of numpy arrays
                                                        if (i == 'x'):
                                                                index = 0
                                                                for num in k:
                                                                        CoE_x = CoE_x + (1 / len(k)) * (index * math.pow(num, 2))
                                                                        index = index + 1
                                                        if (i == 'y'):
                                                                index = 0
                                                                for num in k:
                                                                        CoE_y = CoE_y + (1 / len(k)) * (index * math.pow(num, 2))
                                                                        index = index + 1
                                                        if (i == 'z'):
                                                                index = 0
                                                                for num in k:
                                                                        CoE_z = CoE_z + (1 / len(k)) * (index * math.pow(num, 2))
                                                                        index = index + 1
                                        elif (type(value) is list):
                                                for i in value:
                                                        CoE_x = CoE_x + ((value.index(i))* math.pow(i['x'], 2))
                                                        CoE_y = CoE_y + ((value.index(i))* math.pow(i['y'], 2))
                                                        CoE_z = CoE_z + ((value.index(i))* math.pow(i['z'], 2))
                                        else:
                                                print("ERROR")
                                                return -1
                                        CoE_x = float(CoE_x / totalEnergy['x'])
                                        CoE_y = float(CoE_y / totalEnergy['y'])
                                        CoE_z = float(CoE_z / totalEnergy['z'])
                                        ceoedict = {'x': CoE_x, 'y': CoE_y, 'z': CoE_z}
                                        #Now find shortest intervals containing 70% of the total energy
                                        interval_x = find_interval(value, ceoedict['x'], totalEnergy['x'], 'x')
                                        interval_y = find_interval(value, ceoedict['y'], totalEnergy['y'], 'y')
                                        interval_z = find_interval(value, ceoedict['z'], totalEnergy['z'], 'z')
                                        intervaldict = {'x': interval_x, 'y': interval_y, 'z': interval_z}
                                        buttondictwithinterval[key + '_interval'] = intervaldict
                jsondata_with_interval.append(buttondictwithinterval)
        return jsondata_with_interval

def remove_unnecessary(jsondata):       #Remove the original time and frequency domain sequences from the data, outputting a feature vector
        featurevector = []
        unnecessary_keys = ['acceleration', 'accelerationnog', 'orientation', 'rotation',]
        for buttonpress in jsondata:
                if(buttonpress['button'] !=  None):     #skip bad values
                        buttondict = {}        #Dict of the original data without the unnecessary data
                        for key, value in buttonpress.items():
                                if(key not in unnecessary_keys and key not in [x+'_fft' for x in unnecessary_keys] and key not in [x+'_d' for x in unnecessary_keys]):
                                        buttondict[key] = value
                featurevector.append(buttondict)
        return featurevector

filename = 'data'
datafile = open(filename,'r')
jsondata = json.load(datafile)
datafile.close()

#Process the data step by step

orientationdata = convert_orientation(jsondata)
#print(orientationdata)
noninitialdata = remove_initial(orientationdata)
#print("New dict:", noninitialdata)
derivativedata = calc_derivative(noninitialdata)
#print("New dict:", derivativedata)
datawithdac = calc_dac(derivativedata)
#print(datawithdac)
statsdict = calc_stats(datawithdac)
#print(statsdict)
energydict = calc_total_energy(statsdict)
#print(energydict)
fftdict = calc_fft(energydict)
#print(fftdict)
fftwithstatsdict = calc_stats_fft(fftdict)
#print(fftwithstatsdict)
intervaldict = calc_interval(fftwithstatsdict)
#print(intervaldict)
featurevector = remove_unnecessary(intervaldict)
#print(featurevector) 

#Now format the data to be a little bit clearer and save it to a file
#First for the feature vector data
savefilename_vector = filename + '_vector_processed'
with open(savefilename_vector, 'w') as outfile:
        for buttonpress in featurevector:    ##Loop through all the button presses
                outfile.write("New button:" + str(buttonpress['button']) + "\n")
                outfile.write("@\n")    #Use @ as splitting character
                for key in sorted(buttonpress):
                        outfile.write("'" + str(key) + "'" + ':')               
                        outfile.write(str(buttonpress[key]) + ';'  +'\n')
                outfile.write("@\n")    #Use @ as splitting character

savefilename_all = filename + '_processed'
#And then with the data that still has the original time and frequency domain data (for plotting that data)
with open(savefilename_all, 'w') as outfile:
        for buttonpress in intervaldict:    ##Loop through all the button presses
                outfile.write("New button:" + str(buttonpress['button']) + "\n")
                outfile.write("@\n")    #Use @ as splitting character
                for key in sorted(buttonpress):
                        outfile.write("'" + str(key) + "'" + ':')               
                        outfile.write(str(buttonpress[key]) + ';'  +'\n')
                outfile.write("@\n")    #Use @ as splitting character



