# -*- coding: UTF-8 -*-

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
import matplotlib.pyplot as plt
from matplotlib.pyplot import colorbar
from matplotlib.colors import LogNorm
import time

start_time=time.time()

def Rand_DNA(filename,length):
	a=[]
	for x in range(0,length):
		a.append(np.random.choice([1,2,3,4]))
	print(a)
	np.savetxt(filename+str('.txt'),a,fmt='%s')

Rand_DNA('exDNA',18)

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
print(len(prt))
# prt=np.array(prt)


def check_if_unique(set_of_seq,unique_sequences):
	for seq in set_of_seq:
		if seq not in unique_sequences:
			unique_sequences.append(seq)
	return unique_sequences

def mutator(mut_rate,sequence):
	seq=""
	for x in sequence:
		new=0
		if np.random.rand()<mut_rate:
			new=np.random.choice(["1","2","3","4"])
		else:
			new=x
		seq+=new
	return seq

def doubler(array):
	new_array=[]
	if len(array)<population_cap:
		new_array=array+array
	else:
		new_array=array
		print('the population was culled')
	return new_array

def add_to_charts(charts,set_of_seq):
	for useq in set_of_seq:
		for a in range(len(useq)):
			if useq[a]=='1':
				charts[0][a]+=1
			if useq[a]=='2':
				charts[1][a]+=1
			if useq[a]=='3':
				charts[2][a]+=1
			if useq[a]=='4':
				charts[3][a]+=1
	charts[0][:]=[a/max(charts[0][:]) for a in charts[0][:]]
	charts[1][:]=[a/max(charts[1][:]) for a in charts[1][:]]
	charts[2][:]=[a/max(charts[2][:]) for a in charts[2][:]]
	charts[3][:]=[a/max(charts[3][:]) for a in charts[3][:]]

def mutate_and_double_and_check(mut_rate,set_of_seq,num_of_gens):
	unique_sequences=[]
	num_seq=[]
	frac_num_seq=[]
	counter=0
	set_of_seq=[set_of_seq]+[set_of_seq]
	while counter<num_of_gens:
		counter+=1
		set_of_seq=[mutator(mut_rate,seq) for seq in set_of_seq]
		unique_sequences=check_if_unique(set_of_seq,unique_sequences)
		set_of_seq=doubler(set_of_seq)
		frac_num_seq.append(len(unique_sequences)/(4.**len(seq)))
		num_seq.append(len(unique_sequences))

	line_opacity=0.003
	bar_charts=np.zeros((4,len(set_of_seq[0])))
	add_to_charts(bar_charts,set_of_seq)

	print(bar_charts)
	# im=ax11.imshow(bar_charts,cmap='hot',norm=LogNorm(vmin=0.01, vmax=1),origin='lower')
	
	im=ax11.imshow(bar_charts,origin='lower')
	fig.colorbar(im,ax=ax11)
	for i in range(4):
		for j in range(len(prt)):
			text = ax11.text(j, i, np.round(bar_charts[i, j],decimals=2),
				ha="center", va="center", color="w")

	sequence_numbering=np.linspace(1,len(set_of_seq[0]),len(set_of_seq[0]))
	bar_width=0.2
	opacity=0.5
	ax12.bar(sequence_numbering,bar_charts[0][:],bar_width,alpha=opacity,color='b')
	ax12.bar(sequence_numbering+bar_width,bar_charts[1][:],bar_width,alpha=opacity,color='c')
	ax12.bar(sequence_numbering+2*bar_width,bar_charts[2][:],bar_width,alpha=opacity,color='g')
	ax12.bar(sequence_numbering+3*bar_width,bar_charts[3][:],bar_width,alpha=opacity,color='r')
	return frac_num_seq

population_cap=1e3

num_of_gens=1e3

fig,(ax11,ax12,ax13)=plt.subplots(3,1)

mut_rate=1e-1
one=mutate_and_double_and_check(mut_rate,prt,int(num_of_gens))
# two=mutate_and_double_and_check(1e-2,prt,int(num_of_gens))
# three=mutate_and_double_and_check(1e-3,prt,int(num_of_gens))
# four=mutate_and_double_and_check(1e-4,prt,int(num_of_gens))
# five=mutate_and_double_and_check(1e-5,prt,int(num_of_gens))

min_per_replication=40.

x_one=[k*(min_per_replication/60.)/24. for k in range(len(one))]
# x_two=[k*(min_per_replication/60.)/24. for k in range(len(two))]
# x_three=[k*(min_per_replication/60.)/24. for k in range(len(three))]
# x_four=[k*(min_per_replication/60.)/24. for k in range(len(four))]
# x_five=[k*(min_per_replication/60.)/24. for k in range(len(five))]

def find_slope(array):
	a=np.polyfit(x,array,1)
	return a[0]

print('------ %s seconds -------' % np.round((time.time()-start_time),decimals=3))

ax11.set_title('Visualization of DNA sequence space explored')
ax11.set_ylabel('bp')
ax11.text(0-1.5,0,'A',fontsize=14,fontweight='bold')
ax11.text(0-1.5,1,'C',fontsize=14,fontweight='bold')
ax11.text(0-1.5,2,'T',fontsize=14,fontweight='bold')
ax11.text(0-1.5,3,'G',fontsize=14,fontweight='bold')
ax11.yaxis.set_visible(False)

bar_labels=('A','C','T','G')

ax12.legend(bar_labels)
ax12.set_title('Frequency of appearance of bp vs sequence position')
ax12.set_xlim(0,len(prt)+1)
ax12.set_ylabel('Frequency of appearance')
ax12.set_xlabel('bp numerical ordering')

ax13.plot(x_one,one,color='m')
# labels=(u'\u00B5=1e-1',u'\u00B5=1e-2',u'\u00B5=1e-3',u'\u00B5=1e-4',u'\u00B5=1e-5')
# plt.plot(x_two,two)
# plt.plot(x_three,three)
# plt.plot(x_four,four)
# plt.plot(x_five,five)

ax13.text(2,10,'Number of bp is %s \nNumber of gens is %s \nPopulation cap is %s \nMutation rate is %s' % (len(seq),num_of_gens,population_cap,mut_rate),bbox=dict(facecolor='green', alpha=0.25))

ax13.set_yscale('log')
ax13.set_xscale('log')
ax13.set_title('Fractional sequence space explored vs time')
ax13.set_xlabel('Number of days assuming a doubling time of %s min' %min_per_replication)
ax13.set_ylabel('Fractional sequence space explored')

plt.subplots_adjust(hspace=0.35)
# plt.axes(fig,label='Number of days @ %s min replication time' % min_per_replication)

# axs=plt.subplot('211')
# axs.set_title('bp')
# plt.axes.xlabel('Sequence position')

plt.show()

# text_height=1e-6

plt.text(1*0.8,text_height, 'One Day',
         rotation=90,
         horizontalalignment='center',
         verticalalignment='top',
         multialignment='center')
plt.text(30*0.8,text_height, 'One Month',
         rotation=90,
         horizontalalignment='center',
         verticalalignment='top',
         multialignment='center')
plt.text(180*0.8,text_height, 'Half Year',
         rotation=90,
         horizontalalignment='center',
         verticalalignment='top',
         multialignment='center')
# plt.text(365*0.8,text_height, 'One Year',
#          rotation=90,
#          horizontalalignment='center',
#          verticalalignment='top',
#          multialignment='center')

# plt.axvline(x=1,linestyle=':',color='r')
# plt.axvline(x=30,linestyle=':',color='r')
# plt.axvline(x=180,linestyle=':',color='r')
# plt.axvline(x=365,linestyle=':',color='r')
plt.show()

# output=np.column_stack((x_stitch,yref_stitch,ysig_stitch,fpcavity_stitch))
# np.savetxt(str(text) + '_assembled.csv',output, delimiter=",")