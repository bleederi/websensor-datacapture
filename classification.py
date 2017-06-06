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

import pandas
from pandas.tools.plotting import scatter_matrix
import matplotlib.pyplot as plt
from sklearn import model_selection
from sklearn.metrics import classification_report
from sklearn.metrics import confusion_matrix
from sklearn.metrics import accuracy_score
from sklearn.linear_model import LogisticRegression
from sklearn.tree import DecisionTreeClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.naive_bayes import GaussianNB
from sklearn.svm import SVC
import numpy
from helperfuncs import prettyprint, count_features, read_from_file, list_keys

#TODO: Name variables better, same variable names used in different places (for example 'data')

#Load dataset
datafilename = 'data_1_processed'
datafilename_vector = 'data_1_vector_processed'
with open(datafilename, 'r') as datafile:
        dataset = datafile.read()
with open(datafilename_vector, 'r') as datafile_vector:
        dataset_vector = datafile_vector.read() #Feature vector dataset


buttondata = read_from_file(dataset)
buttondata_vector = read_from_file(dataset_vector)
#prettyprint(buttondata_vector)
featurevectors = []     #Array to hold all the feature vectors
i = 0
for buttonpress in buttondata_vector:
        featurevector = [] #A single feature vector
        for key, value in sorted(buttonpress.items()):  #Feature vector has to be sorted
                if(type(value) is dict):
                        for k, v in value.items():
                                featurevector.append(v)
        featurevectors.append(featurevector)
        i = i+1
#print(featurevectors[0])
#print(len(featurevectors[0]))   #Too short?

dataset_file = 'dataset/dataset1'
#Now write the feature vectors into a dataset file (CSV format)
with open(dataset_file, 'w+') as datasetfile:
        keys = list_keys(buttondata)
        print(keys)
        datasetfile.write(str(keys).strip("[]"))
        datasetfile.write('\n')
        for vector in featurevectors:
                #Add button data (what button was pressed)
                vector.append(buttondata[featurevectors.index(vector)]['button'])
                #Add frequency data (what sensor frequency was used)
                vector.append(buttondata[featurevectors.index(vector)]['frequency'])
                datasetfile.write(str(vector).strip("[]"))
                datasetfile.write('\n')
datasetfile.close()
                

#From sensor readings (one for each reading), need to make sequences (one for each coordinate) - for feature vector, don't need to do this
dictkeys = []
listkeys = []
xyzkeys = []    #Keys that have x, y, z data
abgkeys = []    #Keys that have alpha, beta, gamma data
numpyarrkeys = []    #Keys that have numpy array data
#Still need to find a rule to separate xyz and abg keys
for key, value in buttondata[0].items():
        if(key != 'button' and key != 'frequency'):
                #Keys that have dict data
                if(type(value) is dict):
                        dictkeys.append(key)
                #Keys that have list data
                if(type(value) is list):
                        listkeys.append(key)
                #Keys that have alpha, beta, gamma data
                if('orientation' in key):
                        abgkeys.append(key)
                #Keys that have x,y,z data (all the others except orientation, dac)
                elif('dac' not in key and key != 'button' and key != 'frequency'):
                        xyzkeys.append(key)
                if('_fft' in key):
                        numpyarrkeys.append(key)
buttondata_array = []   #Holds the data for all the buttons
for buttonpress in buttondata:
        data = {}       #Dict holding the sequences for a single button press
        data['button'] = buttonpress['button']
        for key, value in buttonpress.items():
                seq = {}
                seqlist = []
                if(key == 'button' or key == 'frequency'):    #Don't need to make sequences for these
                        continue
                #Make sequences for each key
                #First handle list keys
                if key in listkeys and key in xyzkeys and key not in numpyarrkeys:
                        seq_x = []
                        seq_y = []
                        seq_z = []
                        for i in value:
                                seq_x.append(i['x'])
                                seq_y.append(i['y'])
                                seq_z.append(i['z'])
                        seq = {'x':seq_x, 'y':seq_y, 'z':seq_z}
                elif key in listkeys and key in abgkeys and key not in numpyarrkeys:
                        seq_alpha = []
                        seq_beta = []
                        seq_gamma = []
                        for i in value:
                                seq_alpha.append(i['alpha'])
                                seq_beta.append(i['beta'])
                                seq_gamma.append(i['gamma'])
                        seq = {'alpha':seq_alpha, 'beta':seq_beta, 'gamma':seq_gamma}
                elif key in listkeys and key in numpyarrkeys:
                        seqlist = value
                #Then handle dict keys
                elif key in dictkeys and key in xyzkeys:
                        seq = {'x':value['x'], 'y':value['y'], 'z':value['z']}
                elif key in dictkeys and key in abgkeys:
                        seq = {'alpha':value['alpha'], 'beta':value['beta'], 'gamma':value['gamma']}
                else:   #DAC keys
                        data[key] = value
                if(seq):
                    data[key] = seq
                elif(seqlist):
                    data[key] = seqlist
        buttondata_array.append(data)

def buttonselection():      #Condition for selecting the buttons to be plotted
    for x in buttondata_array:
        if x['button']  == 2:    #Select all buttons that fulfil this condition
            yield x

#prettyprint(buttonselection())

def plot(buttons, sameplot=False):
        index = 1
        for button in buttons:
                #Plot sequences
                fig = plt.figure(index)
                fig.suptitle('Data for button ' + str(button['button']))
                ax1 = plt.subplot(221)
                plt.plot(button['acceleration']['x'], color='r', label='accelx')
                plt.plot(button['acceleration']['y'], color='b', label='accely')
                plt.plot(button['acceleration']['z'], color='g', label='accelz')
                ax1.legend(["accx", "accy", "accz"], loc='upper center', bbox_to_anchor=(0.5, 1.10),
                        ncol=3, fancybox=True, shadow=True)
                plt.ylabel('Acceleration')

                ax2 = plt.subplot(222)
                plt.plot(button['rotation']['x'], color='r', label='rotx')
                plt.plot(button['rotation']['y'], color='b', label='roty')
                plt.plot(button['rotation']['z'], color='g', label='rotz')
                ax2.legend(["rotx", "roty", "rotz"], loc='upper center', bbox_to_anchor=(0.5, 1.10),
                        ncol=3, fancybox=True, shadow=True)
                plt.ylabel('Rotation')

                ax3 = plt.subplot(223)
                plt.plot(button['orientation']['alpha'], color='r', label='orix')
                plt.plot(button['orientation']['beta'], color='b', label='oriy')
                plt.plot(button['orientation']['gamma'], color='g', label='oriz')
                ax3.legend(["orix", "oriy", "oriz"],loc='upper center', bbox_to_anchor=(0.5, 1.10),
                        ncol=3, fancybox=True, shadow=True)
                plt.ylabel('Orientation')
                if not sameplot:   
                        index = index+1 #Plot each button press in different plot window
                manager = plt.get_current_fig_manager()
                manager.resize(*manager.window.maxsize())
                plt.savefig('button_' + str(button['button']) + '_' +str(index)  + '.svg', format='svg')
        plt.show() 

#plot(buttondata_array, True)
