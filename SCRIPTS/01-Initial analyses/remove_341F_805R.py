#!/usr/bin/env python
# -*- coding: utf-8 -*-

##########################################################################################
"""
------------------------------------------------------------------------------------------
Remove primer 341F and 805R for illumina paired-end reads. 

341F: CCTACGGGNGGCWGCAG --> R1
805R: GACTACHVGGGTATCTAATCC --> R2

Input files suhould be under the current directory in:
-R1inicial folder
-R2inicial folder

Output files will be move to:
-R1 folder
-R2 folder
or create then if they do not exists

Requirements:
-Cutadapt
------------------------------------------------------------------------------------------
"""
##########################################################################################

#--Imports
import sys
import os
import argparse
import gzip

__author__ = "Alberto Rastrojo"
__version_info__ = ('1','0','0')
__version__ = '.'.join(__version_info__)

def main():

	##########################################################################################
	#--Argument
	parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
	parser.add_argument('-c', dest = 'cores', type = str, default = 2, help = 'number of threads. Default: 2')
	parser.add_argument('-x', dest = 'debug', action='store_true', default = False, help = 'Debug. Default: False')
	parser.add_argument('-v', '--version', action='version', version=__file__ + ' v' + __version__)
	args = parser.parse_args()
	##########################################################################################
	#--Checking requirements
	checkRequirements(['cutadapt'])

	#--Checking input folder
	wd = os.getcwd()
	checkfile(f'R1inicial')
	checkfile(f'R2inicial')

	#--Creating output folders
	if not os.path.exists('R1'):
		os.mkdir('R1')
	if not os.path.exists('R2'):
		os.mkdir('R2')

	#--Removing 341F from R1 reads
	R1_files = sorted([file for file in os.listdir(f'R1inicial') if not file.startswith('.')])
	R2_files = sorted([file for file in os.listdir(f'R2inicial') if not file.startswith('.')])

	for r1, r2 in zip(R1_files, R2_files):
		r1_prefix = r1.split('.')[0]
		r2_prefix = r2.split('.')[0]

		cmd = f"""cutadapt \
		--cores={args.cores} \
		-e 0.25 \
		-g CCTACGGGNGGCWGCAG \
		-G GACTACHVGGGTATCTAATCC \
		-o ./R1/{r1_prefix}.fastq.gz \
		-p ./R2/{r2_prefix}.fastq.gz \
		--untrimmed-output ./R1/{r1_prefix}_untrimmed.fastq.gz \
		--untrimmed-paired-output ./R2/{r2_prefix}_untrimmed.fastq.gz \
		./R1inicial/{r1} \
		./R2inicial/{r2}
		"""

		if args.debug: print(cmd)
		log = runexternalcommand(cmd)

		# r1_trimmed = sum(1 for _ in gzip.open(f'./R1/{r1_prefix}.fastq.gz', 'rb'))
		# r1_untrimmed = sum(1 for _ in gzip.open(f'./R1/{r1_prefix}_untrimmed.fastq.gz', 'rb'))

		# Should be the same
		# r2_trimmed = sum(1 for _ in gzip.open(f'./R2/{r2_prefix}.fastq.gz'))
		# r2_untrimmed = sum(1 for _ in gzip.open(f'./R2/{r2_prefix}_untrimmed.fastq.gz'))

		# results = f"Inputs: {r1} ({r2})\nTrimmed: {r1_trimmed}\nUntrimmed {r1_untrimmed}"
		# print("#"*90)
		# print(results)
		# print("#"*90)

def checkRequirements(programs):

	""" Return true if program is on PATH """

	import sys
	from shutil import which

	for program in programs:
		
		if which(program) is None:
			sys.exit(f'ERROR: **{program}** were not found!\n')

	return which(program) is not None

def runexternalcommand(cmd):

	""" Facilitates the running of external commands by using subprocesses """
	import subprocess
	# from subprocess import call, Popen, communicate

	out, err = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE).communicate()

	log = ''
	if out: 
		out = out.decode()
		log += out
	if err: 
		err = err.decode()
		log += err

	return log

def checkfile(file):

	import sys

	""" File checking """
	if not os.path.exists(file):
		print ("-" * 90)
		print ("ERROR: file '" + file + "' not found")
		print ("-" * 90)
		sys.exit()

##########################################################################################
if __name__ == "__main__":

	main()

