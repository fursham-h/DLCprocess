#!/usr/bin/python

import pandas as pd
import matplotlib.pyplot as plt


usage = """

*****************************************************************************
This script.....

Use Syntax
python ccut.py <option-argument> <amino acid sequence or Fasta file>

options:
opt:			argument:	Default		description:
  -o (--organisms)	  organism list  Homo_sapiens	  List of organisms (sep by comma)
  -w (--weights)	  number list  	 1 	  	  List of weight values (sep by comma)
  -s (--save)  		  NA  		 False 	  	  Option to save the compromised codon table
  -h (--help)		  NA  		 False	 	  Displays usage of the script
*****************************************************************************

"""

########################################
############ Class and Def #############
########################################

def argparse():
	# attempt to read command-line option-argument tuples and mandatory argument.
	try:
		options, rawinput = getopt.getopt(sys.argv[1:], "ho:w:s",["help","organisms=","weights=","save"])
		rawinput = str(rawinput[0]).strip()
	except getopt.GetoptError, err:
		sys.exit(usage)		# exit script and print usage if arguments to options are not provided 
	except IndexError:
		sys.exit(usage)		# exit script and print usage if command-line input is not provided

	# define necessary variables
	list_organisms = ['Homo_sapiens']
	weights = [1]
	save = False

	# parse command-line options into its appropriate variables/actions
	for opt,arg in options:
		if opt in ("-o","--organisms"):
			list_organisms = map(str,arg.split(','))
			weights = [1]*len(list_organisms)
		if opt in ("-w","--weights"):
			weights = map(float,arg.split(','))
		if opt in ("-s","--save"):
			save = True
		if opt in ("-h","--help"):
			sys.exit(usage)

	# Additional preparations and checks
	dict_codonrank,aa = pickle.load(open('definitions','r'))		#pre-load dictionaries and aa list
	dict_species = pickle.load(open('species','r'))
	if not re.search(r'.txt',rawinput):
		if not set(rawinput.upper()) <= set(aa):sys.exit('Error: Input sequence may contain illegal characters')			# check for non-amino acid characters
	if len(weights) != len(list_organisms): sys.exit('Error: weights do not match number of organisms')			# check for correct number of weights
	
	return rawinput,dict_codonrank,list_organisms,weights,save,dict_species


def processh5():
    dat = pd.read_hdf()


########################################
################ Main ##################
########################################

if __name__ == "__main__":
    argparse()
    







