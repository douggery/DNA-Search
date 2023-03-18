"""

This script is meant to investigate the 
effect of mutating genomic data with and
without site specificity.



___________________
TODO
___________________
1) have each muation event accompany a fission event
2) fit the search space vs time
3) fit the search space, mutation rate, and number of gens
"""
#from Gen_Random_DNA import Rand_DNA

#Rand_DNA('exDNA',99)

import numpy as np
from mutator import mutator
import matplotlib.pyplot as plt
import time

start_time=time.time()

def read_seq(inputfile): 
	with open(inputfile, "r") as f: 
		seq = f.read() 
	seq = seq.replace("\n", "") 
	seq = seq.replace("\r", "") 
	return seq 

def translate(seq):
	table = { 
		'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
		'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
		'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
		'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
		'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
		'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
		'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
		'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
		'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
		'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
		'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
		'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
		'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
		'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
		'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
		'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
	}
	protein =""
	if len(seq)%3 == 0:
		for i in range(0, len(seq), 3):
			codon = seq[i:i + 3]
			protein+= table[codon]
	return protein

inputfile="exDNA.txt"
f=open(inputfile,"r")
seq=f.read()

seq=seq.replace("\n","")
seq=seq.replace("\r","")

prt=read_seq(inputfile)


def check_if_unique(seq,unique_sequences):
	if seq not in unique_sequences:
		unique_sequences.append(seq)
	return unique_sequences

def mutate_and_check(mut_rate,seq,num_of_gens):
	unique_sequences=[]
	num_seq=[]
	for i in range(num_of_gens):
		unique_sequences=check_if_unique(seq,unique_sequences)
		seq=mutator(mut_rate,seq)
		num_seq.append(len(unique_sequences)/(4**len(seq)))
	print("Mut Rate:", mut_rate, "		num of gens:", num_of_gens, "	num of uniques", len(unique_sequences))
	return num_seq

num_of_gens=1e4

one=mutate_and_check(1e-1,prt,int(num_of_gens))
two=mutate_and_check(1e-2,prt,int(num_of_gens))
three=mutate_and_check(1e-3,prt,int(num_of_gens))
four=mutate_and_check(1e-4,prt,int(num_of_gens))
five=mutate_and_check(1e-5,prt,int(num_of_gens))

x=[k*(20/60)/24 for k in range(len(one))]

def find_slope(array):
	a=np.polyfit(x,array,1)
	return a[0]

def doubler(array):
	new_array=[]
	for element in array:
		new_array.append(element)
	for element in array:
		new_array.append(element)
	return new_array

print('------ %s seconds -------' % np.round((time.time()-start_time),decimals=3))

plt.plot(x,one)
plt.plot(x,two)
plt.plot(x,three)
plt.plot(x,four)
plt.plot(x,five)
plt.plot()
plt.yscale('log')
plt.xscale('log')
labels=('μ=1e-1','μ=1e-2','μ=1e-3','μ=1e-4','μ=1e-5')
plt.legend(labels)
plt.title('Number of unique genetic sequences vs \n generation for various mutation rates \n Number of bp is %s' % len(seq))
plt.ylabel('Fractional sequence space searched')
plt.xlabel('Number of days @ 20 min replication time')

plt.text(1*0.8,1e-2, 'One Day',
         rotation=90,
         horizontalalignment='center',
         verticalalignment='top',
         multialignment='center')
plt.text(30*0.8,1e-2, 'One Month',
         rotation=90,
         horizontalalignment='center',
         verticalalignment='top',
         multialignment='center')
plt.text(180*0.8,1e-2, 'Half Year',
         rotation=90,
         horizontalalignment='center',
         verticalalignment='top',
         multialignment='center')
plt.text(365*0.8,1e-2, 'One Year',
         rotation=90,
         horizontalalignment='center',
         verticalalignment='top',
         multialignment='center')

plt.axvline(x=1,linestyle=':',color='r')
plt.axvline(x=30,linestyle=':',color='r')
plt.axvline(x=180,linestyle=':',color='r')
plt.axvline(x=365,linestyle=':',color='r')
plt.show()
