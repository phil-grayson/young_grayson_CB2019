# -*- coding: utf-8 -*-
"""
Created on Wed Jul 12 14:46:43 2017

@author: Phil
"""


import os
import csv
import sys

allGenes = {}
fileout=(str(sys.argv[1].split('.')[0])+'LO_singles.txt')
fout= open(fileout, 'w') 

# start with the LO intersect file containing each transcript ID and galGal gene name associated with it
with open(sys.argv[1]) as handle: 
    reader=csv.reader(handle,delimiter='\t')
    for strLine in reader:
        if strLine[0] not in allGenes.keys(): #if the transcript isn't in the dictionary, add it with the galGal name as a list of 1
            allGenes[strLine[0]]=[strLine[1]]
        else:
            allGenes[strLine[0]].append(strLine[1]) # if the transcript is in the dictionary, append the galGal name to the associated list
        
flipped = {} #going to create a new dictionary from the old one

for dro, gal in allGenes.items(): #for key, items in original dict
    if len(gal) == 1: # if the length of the item list is 1 (i.e., only one galGal name with this transcript)
        if gal[0] not in flipped: # if the galGal name isn't already in the flipped dict
            flipped[gal[0]] = [dro] # create a key for the galGal name and add the transcript ID as a list of 1
        else: # if the galGal name already exists, add the transcript to the item list
            flipped[gal[0]].append(dro)

for dro, gal in allGenes.items(): #this time, if length of item list is not equal to 1
    if len(gal) != 1:
        for gene in gal: # for every galGal name in the item list
            if gene in flipped: # if galGal name is in flipped
                flipped[gene].append(dro) # append the transcript to the item list.  no need for else, since this is only to append additional transcript IDs to galGal names that have only appeared once up until this point.
        
for gene in flipped: # now we go galGal gene by galGal gene through flipped and if there is only one transcript associated with it ...
    if len(flipped[gene]) == 1:
        fout.write(str(flipped[gene][0]) + "\t" + str(gene) + "\n") # ... print that out

fout.close()
    
