#!/usr/bin/env python3

# Define arguments

import sys, getopt

directory = ''
system = ''
output = ''
try:
   opts, args = getopt.getopt(sys.argv[1:],"hd:s:n:o:",["dir=","sys=","out="])
except getopt.GetoptError as err:
   print(err)
   print('ERROR: Syntaxe should be "convert_sf.py -d <directory> -s <system_file> -o <output_name>"')
   sys.exit(2)
for opt, arg in opts:
   if opt == '-h':
      print('convert_sf.py -d <directory> -s <system_file> -o <output_name>')
      sys.exit()
   elif opt in ("-d", "--dir"):
      directory = arg
   elif opt in ("-s", "--sys"):
      system = arg
   elif opt in ("-o", "--out"):
      output = arg

##########################

import molparse as mp 

sys = mp.parse(f'{directory}/{system}')
ndx = mp.parseNDX('index_qmmm.ndx')
sys.summary()

sys['r0'].summary()
sys['r1'].summary()

indices = ndx['LINK']
print(indices)

for atom in sys.atoms:
	if atom.number in indices:
		# print(f'{atom.number} is in indices')
		# atom.summary()
		atom.set_name('H')

	atom.name = atom.symbol

mp.write(f'{directory}/{output}',sys)