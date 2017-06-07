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

import math
import ast

__author__ = "Jesse Nieminen"
__status__ = "Development"

def prettyprint(jsondata, *args): #Clearly print the data - input is all the button data and optional argument specifies keys to be printed
        for buttonpress in jsondata:    ##Loop through all the button presses
                print("\nNew button: ", buttonpress['button'], "\n")
                if(args):
                        for key in sorted(buttonpress):
                                if args[0] in key:
                                        print(key)                
                                        print(buttonpress[key])
                else:
                        for key in sorted(buttonpress):
                                print(key)                
                                print(buttonpress[key])


def find_interval(seq, coe, totalEnergy, coord):     #Find the interval of 70% energy for a sequence, centered on coe. For single coord only
        interval = []
        #Now find shortest intervals containing 70% of the total energy
        index = round(coe)
        #Calculate energy for interval i centered on index
        bound_lower = 0
        bound_upper = 0
        for i in range(0, len(seq)):
                e = 0
                bound_lower = index-i
                bound_upper = index+i
                #Check that we are not out of bounds
                if(bound_lower < 0):
                        bound_lower = 0
                if(bound_upper > len(seq)):
                        bound_upper = len(seq)
                for v in seq[bound_lower:bound_upper]:
                        e = e + math.pow(seq[seq.index(v)][coord], 2)
                if(e > 0.7 * totalEnergy):
                        interval = bound_upper - bound_lower     
                        return interval

def count_features(seq):        #Counts the number of features in the list of sequences or "feature vectors"
        nfeatures = 0        
        #prettyprint(seq)
        for key, value in seq[0].items():
                #print(key)
                if(type(value) is dict):
                        #print("Key with", len(value.keys()), "keys in value", value.keys())
                        nfeatures = nfeatures + len(value.keys())
                else:
                        #print(value)
                        #print("Key with", 1, "value")
                        nfeatures = nfeatures + 1
        return nfeatures

def list_keys(seq):        #Creates the list of keys for use in making the CSV file - biutt
        keys = []
        for key, value in seq[0].items():
                #print(key)
                if(type(value) is dict):
                        #print("Key with", len(value.keys()), "keys in value", value.keys())
                        for key2 in value.keys():                        
                                keys.append(key + '_' + key2)
        #Have to separately add these two
        keys.append('button')
        keys.append('frequency')
        return keys

def read_from_file(dataset):   #Read string data from file and convert it to Python types
        buttondata = [] #List of button presses, here 
        data_split = dataset.split('@')    #Split dataset into distinct button presses
        for buttonpress in data_split:
                datadict = {}   #Here all the data read from the dataset will be saved, one per button
                data = buttonpress.split(';')
                for seq in data:
                        seq = seq.strip()
                        if ("New button" in seq or not seq) :
                            continue
                        else:
                                #Read the data into a format that we can easily manipulate (string -> new format)
                                key = seq.split(':', 1)[0].translate({ord(c): None for c in "'"})       #Translate for removing the "'"
                                data = seq.split(':', 1)[1]
                                if not(key.endswith("fft")):
                                        data_read = ast.literal_eval(data)      #Cannot read numpy arrays
                                else:   #Handle numpy arrays separately
                                        #First find the dict keys
                                        datasplit = data.split('), ')
                                        for i in datasplit:
                                                spliti = i.split(':')
                                                key2 = spliti[0].translate({ord(c): None for c in "{'"})
                                                data_read = ast.literal_eval(spliti[1].translate({ord(c): None for c in "()}'array\n "}))
                                datadict[key] = data_read
                if datadict:
                        buttondata.append(datadict)
        return buttondata

        #prettyprint(buttondata)

def convert_orientation(jsondata):      #Convert orientation matrix to Euler angles
        orientationjsondata = []        #New JSON data with orientation converted to Euler angles
        for buttonpress in jsondata:
                if(buttonpress['button'] !=  None):     #skip bad values
                        orientationdict = {}      #Dict (JSON data) for a button press with orientation converted
                        for key, value in buttonpress.items():
                                if(key != 'orientation'):
                                        orientationdict[key] = value
                                else:   #Here convert orientation
                                        orientationdict[key] = []
                                        for orimatrix in value:  #Each orientation matrix
                                                tempdict = {}   #Temp dict to store Euler angle values
                                                r11 = orimatrix['0']
                                                r21 = orimatrix['4']
                                                r31 = orimatrix['8']
                                                r32 = orimatrix['9']
                                                r33 = orimatrix['10']
                                                betadivisor = math.sqrt(math.pow(r32,2) + math.pow(r33,2))
                                                if(r11 != 0 and r33 != 0 and betadivisor != 0): #Can't divide by zero
                                                        alpha = math.atan2(r21, r11)
                                                        beta = math.atan2(-r31, (math.sqrt(math.pow(r32,2) + math.pow(r33,2))))
                                                        gamma = math.atan2(r32, r33)
                                                tempdict = {'alpha': alpha, 'beta': beta, 'gamma': gamma}
                                                orientationdict['orientation'].append(tempdict)
                orientationjsondata.append(orientationdict)
        return orientationjsondata
