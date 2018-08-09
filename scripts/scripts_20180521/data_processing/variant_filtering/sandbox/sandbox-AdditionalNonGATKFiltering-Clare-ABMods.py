# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

## This fixed version handles the error with 'variant' sites that aren't variant after filtering the DP/GT ... i.e. all 0/1/ or 1/1 become ./. and AC=0
## This version takes user arguments for the input values
## added 'FT' format field

### This new version filters for RGQ and GQ.

# the script takes a gatk filtered file
# there should only be variant SNPs or nonvariant sites that passed in it.

# for the genotype filtering, it is a bit hacky. After running through this, EVERY genotype will have an FT flag.
# I don't append the flags, I just make it FAILgt if it's ./. (before or after filtering), and PASS if it isn't

### tests

# new genotypes outputted in right order
# DP min and max works
# errors on unexpceted genos
# het and missing works

import sys
import gzip
import re
import time
import pyvcf

# sys argv 1 = filename
minRGQ=int(sys.argv[2]) # exclude if below 1
minGQ=int(sys.argv[3]) # exclude if below 20
minDP=int(sys.argv[4]) # 10 
maxnocall=int(sys.argv[5]) # maximum no calls in a line
maxprophet=float(sys.argv[6]) # I used 0.99 , but had tried 5/8 = 0.6, 7/10, 6/9 (6/10 gets through)

fnadd='minRGQ'+str(minRGQ) +'minGQ'+str(minGQ)+'minDP'+str(minDP)+'maxHet'+str(maxprophet)+'missing'+str(maxnocall)
print 'yourspecs', fnadd
# this is applied to my joint variant - nonvariant filtered file.

counter_nop=0
counter_p=0

in_vcf_file = sys.argv[1]
print in_vcf_file
inVCF = gzip.open(in_vcf_file, 'r')
#inVCF = open(myvcf, 'r')
name=in_vcf_file.split('.')[0]
outfname='6_bespokefixed_v4_' + str(fnadd) + '_' + str(name) + '.vcf'
errorname='ERRORs_from_' + str(name) + '.txt'
outVCF = open(outfname, 'w')
#outVCF = gzip.open('out.vcf', 'wb')
errorVCF = open(errorname, 'w')

refbases=set(['A', 'C', 'G', 'T'])
altbases=set(['A', 'C', 'G', 'T', '.'])
expected_genos=set(['0/0', '0/1', '1/1', './.'])


maxDP_D={'AW10':100,'AW11':99,'AW12':84,'AW13':97,'AW14':100,'AW15':84,'AW16':87,'AW17':95,'AW18':96,'AW19':90,'EW1':116,'EW2':113,'EW3':104,'EW4':105,'EW5':109,'EW6':102,'EW7':94,'EW8':74,'EW9':90,'EW10':92,'GS1':37,'GS2':40,'GS3':38,'GS4':38,'GS5':36,'GS6':39,'GS7':41,'GS8':36,'GS9':41,'GS10':38,'IR1':59,'IR2':57,'IR3':59,'IR4':58,'IR5':59,'IR6':57,'IR7':60,'IR8':62,'IR9':63,'IR10':60,'TM1':40,'TM2':40,'TM3':40,'TM4':35,'TM5':34,'TM6':35,'TM7':38,'TM8':37,'TM9':33,'TM10':40, 'WD1':69,'PG1':122, 'PG2':122, 'PG3':133, 'PG4':108, 'PG5':130, 'PG6':112, 'PG7':124, 'PG8':123, 'PG9':119, 'PG10':112, 'AW1':78,'AW2':125,'AW3':122,'AW4':91,'AW20':94,'AW21':79,'AW22':126,'AW23':123,'AW24':92,'PG11':113,'PG12':118,'PG13':108,'PG14':106,'PG15':110,'PG16':119,'BC1':75,'BC2':63,'BC3':64,'BC4':57,'BC5':45,'BC6':60,'BC7':57,'BC8':62,'BC9':63,'BC10':58}

## tester
#maxDP_D={'GS1':150,'GS2':20,'GS3':30,'GS4':40,'GS5':50,'GS6':60,'GS7':70,'GS8':80,'GS9':90,'GS10':100}

# gets the genotype from the genotype field e.g. 0/0:2,3:45:78 it returns 0/0
# also does a check that its 0/0, 0/1, 1/1, or ./. - nothing else should be in the file
def get_genotype(genotypefield):
	field=genotypefield.split(':')
	return field[0]

##
#chr38	36587	.	C	T	388.83	PASS	AC=2;AF=0.250;AN=8;BaseQRankSum=2.54;ClippingRankSum=0.00;DP=148;ExcessHet=3.2451;FS=0.000;InbreedingCoeff=-0.1109;MQ=60.00;MQRankSum=0.00;QD=14.40;ReadPosRankSum=1.41;SOR=0.836	GT:AD:DP:FT:GQ:PL	0/1:9,9:18:PASS:99:294,0,214	./.:18,0:18:LowGQ:48:0,48,720	0/1:5,4:9:PASS:99:136,0,153	0/0:23,0:23:PASS:60:0,60,900	./.:15,0:15:LowGQ:42:0,42,630	./.:9,0:9:LowGQ:24:0,24,360	./.:12,0:12:LowGQ:30:0,30,450	./.:14,0:14:LowGQ:39:0,39,585	0/0:18,0:18:PASS:51:0,51,765	./.:12,0:12:LowGQ:30:0,30,450

## adds the FT field to the format if needed (if already there returns it as is). FT is after the DP field.
def add_FT(currentformatfield):
#	print '***** current', currentformatfield
	newfields=''
	if 'FT' in currentformatfield:
		EXISTING='YES'
		newfields+=currentformatfield
	else:
		EXISTING='NO'
		cfields=currentformatfield.split(':')
#		print 'cfields'
		for i in cfields:
			# the FT field goes after DP field
			if i == 'GT':
				newfields += str('GT')
			elif i == 'DP':
				newfields += str(':DP:FT')
			else:
				newfields += str(':')
				newfields += str(i)
	return newfields, EXISTING
	

# the FT field goes after DP field

# existing = YES or NO, and is whether theere is an exisiting FT flag.
# oldgi is the current genotypefield e..g 0/1:10,1:35:90
# Dpp is the position of the DP field in the oldgi
# flg is the flag you want to add

## What I am actually going to do is just add FAILgt as the flag - so it doesnt get two flags.

def write_newgi(existing,oldgi,ncall,DPp,flg):
	
	prevgi=oldgi.split(':')	
	updatedgi = ''
	insertflg=''
	# If the call is ./. then the FT will be whateer flg is i.e. FAILgt
	# If the call is anything else, then FT will be PASS
	
	if ncall == './.':
		insertflg += flg
	elif ncall in expected_genos:
		insertflg += 'PASS'
	else:
		print 'ERROR on this', ncall
		sys.exit('*** ERROR exiting unexpected genotype in write_newgi function')
	
	# now make the new gi
	if existing == 'NO':
#		print '** no existing'
		# first add the new genotype
		updatedgi += ncall
		updatedgi += ':'

		# then up to and including the DP	
		updatedgi += ':'.join(prevgi[1:DPp+1])

		# then the flag
		updatedgi += ':'
		updatedgi += insertflg
		updatedgi += ':'

		# then the after DP fields
		updatedgi += ':'.join(prevgi[DPp+1:])
	
#		print 'old, new gi', oldgi, updatedgi

	# it already has an FT 	
	elif existing == 'YES':
#		print 'yes yes'
		# first add the new genotype
		updatedgi += ncall
		updatedgi += ':'

		# then up to and including the DP	
		updatedgi += ':'.join(prevgi[1:DPp+1])

		# then the flag
		updatedgi += ':'
		updatedgi += insertflg
		updatedgi += ':'

		# then the after FT fields i.e. DPp+2 onwards e..g GT:AD:DP:FT:GQ:PL i.e. GQ and PL fields in this example	
		
		updatedgi += ':'.join(prevgi[DPp+2:])
	
#		print 'old, new gi', oldgi, updatedgi
	

	else:
		sys.exit('*#*#*#*$*$#*# ERROR error exiting')
## next this

	return updatedgi

# this funciton needs to also get the genotypes AFTER DP filter	

## 0/0:13,0:13:LowGQ20;lowDP:39:0,39,453

# this function returns the new FORMAT fields, and the new genotype fields and a pass or fail for the het

def gfilter_DP(genofields, INformatfields):
#	print genofields
#	print INformatfields
	
	# get the position for the DP in the original formatfield
	formatfields=INformatfields.split(':')
	
	# get position of RGQ or GQ
	if 'RGQ' in formatfields:
 		RGQpos=[formatfields.index(x) for x in formatfields if 'RGQ' in x][0]
 	else:
 		RGQpos='NA'
	
	if 'GQ' in formatfields:
 		GQpos=[formatfields.index(x) for x in formatfields if 'GQ' in x][0]		
 	else:
 		GQpos='NA'
 
	DPpos=[formatfields.index(x) for x in formatfields if 'DP' in x][0]
	
	# add the FT to the format fields - if already there it doesnt change
	newformatfields, FTalready=add_FT(INformatfields)
#	print 'new', newformatfields
	

	#GT:AD:DP:FT:GQ:PL
	
	# i think this function really needs to work on the whole line
	counter_samples=0
	
	allcalls=[] # to check if too many missing
	FTflagadd='FAILgt'	

	newgenofields=[]
	# loop through the genotypes - check the DP, and make ./. if needed and then update and return the whole genotype field
	for gi in genofields:
	
#		print gi
		
		# get the sample so you can get the max DP for that sample from the dict

#		print 'counter', counter_samples
		samp_name=samples[counter_samples]
#		print samp_name, gi	
			
		# get the maxDP for that sample
		samplemaxDP=maxDP_D[samp_name]
		
		##### GATK is outputting incorrectly with some genotypes without full info e.g. 0/0:0,10 when the format field is GT:AD:DP:RGQ. This is problematic because it breaks my script (because i look up the DP, which isn't always there), it also means the genotype wasn't filtered by GATK e.g. with RGQ=0, GQ .. or DP.... My way to handle that is to make the call ./. vs throw the whole site, as it may only be one sample
		## I fix it by making the the sampleDP 0
				
		if len(gi.split(':')) != len(formatfields):
			sampleDP=0
		elif gi.split(':')[DPpos] == '.':
			sampleDP=0
		else:
			sampleDP=int(gi.split(':')[DPpos])		
		
		# You need to do all your filters here for the genotype - i.e. genotype becomes ./. here because of the way that I do the 
		# FT handling - i.e. it either gets a pass or fail, rather than inserting mlutiple fail flags e..g LowGQ,LowDP... i just make it FAILgt if the call is ./. before or after filtering here
		
		# New bit that is needed. This is where GATK doesnt work properly and outputs the wrong thing so there is no DP or something. These all need to be FAIL genotypes.
		#  So I just made it ./. if the format fields dont match the genotype info e.g. GT:AD:DP:FT:RGQ	./.:9,20:29:PASS:
		# basically this next line says if the length format fields doesnt match the length of genotype fields OR if the last field is empty (this has been an issue)
		
		if len(gi.split(':')) != len(formatfields) or gi.split(':')[len(gi.split(':'))-1] == '':
			newcall='./.'
			allcalls.append(newcall)
			filteredgi='./.'
			newgenofields.append(filteredgi)

		# This bit applies the filters - if it fails, change the field and genotype, and return the new genotype e.g. ./.:4,0:100:FAILgt:9			
		elif sampleDP < minDP or sampleDP > samplemaxDP:
#			print samp_name, 'FAIL', gi, sampleDP, 'too high or low, max is: ', samplemaxDP
			
			newcall='./.'
			allcalls.append(newcall)
			filteredgi=write_newgi(FTalready, gi, newcall, DPpos, FTflagadd)
			newgenofields.append(filteredgi)
			
		# assess if DP is within the range (do this rather than just an else incase of weird vcf issues)
		elif sampleDP >= minDP or sampleDP <= samplemaxDP:
			# the call passess the DP filter, so the current call is the kept call
			# n.b. this could be ./. or 0/0 ...
			
			## if it passes the DP filter, check it also passes the GQ/RGQ filter 
			# if it fails the RGQ or the GQ make it ./. If it passes, then get the new call.
			
			# invar sites have RGQ
			if RGQpos != 'NA':
				# to handle if .
				if gi.split(':')[RGQpos] == '.':
					sampleRGQ=0
				else:
					sampleRGQ=int(gi.split(':')[RGQpos])

#				print RGQpos, gi
				sampleRGQ=int(gi.split(':')[RGQpos])
#				print sampleRGQ	
				if sampleRGQ < minRGQ:
#					print '$$ FAILRGQ'
					newcall='./.'
					allcalls.append(newcall)
					filteredgi=write_newgi(FTalready, gi, newcall, DPpos, FTflagadd)
					newgenofields.append(filteredgi)
				else:
#					print '** PASS RGQ'
					gcall=gi.split(':')[0]
					newcall=gcall
					allcalls.append(newcall)
					filteredgi=write_newgi(FTalready, gi, newcall, DPpos, FTflagadd)
					newgenofields.append(filteredgi)						
			
			# var sites have GQ
			elif GQpos != 'NA':
				
				# to handle if .
				if gi.split(':')[GQpos] == '.':
					sampleGQ=0
				else:
					sampleGQ=int(gi.split(':')[GQpos])
			
				if sampleGQ < minGQ:
#					print 'FAILGQ', gi
					newcall='./.'
					allcalls.append(newcall)
					filteredgi=write_newgi(FTalready, gi, newcall, DPpos, FTflagadd)
					newgenofields.append(filteredgi)
#					print filteredgi
				else:
					gcall=gi.split(':')[0]
#					print '** PASS GQ'
					newcall=gcall
					allcalls.append(newcall)
					filteredgi=write_newgi(FTalready, gi, newcall, DPpos, FTflagadd)
					newgenofields.append(filteredgi)			
			
			else:
				print genofields
				sys.exit('**** ERROR exiting - something wrong with the new G filter loop')		
#			gcall=gi.split(':')[0]
#			newcall=gcall
#			allcalls.append(newcall)
#			filteredgi=write_newgi(FTalready, gi, newcall, DPpos, FTflagadd)
#			newgenofields.append(filteredgi)
			
		else:
			sys.exit('*** ERROR exiting - there is something wrong with your DP')
		# increment the counter
		counter_samples+=1
	
	# increment this if it fails the maxprop het or maxnocall
	calltest=0

	REFct = allcalls.count('0/1') + (allcalls.count('0/0') *2)
	ALTct = allcalls.count('0/1') + (allcalls.count('1/1') *2)
	
	# this has to come first because instances where all ./.
	if len(allcalls) != len(genofields):
		sys.exit('*** ERROR exiting - length of old genotypes and new genotypes isnt same')	
	elif allcalls.count('./.') > maxnocall:
#		print 'nocall fail', allcalls
		calltest+=1
	else:
		if allcalls.count('0/1') == 0:
			prop_het=0
		else:
			prop_het = (allcalls.count('0/1') + 0.0) / (allcalls.count('0/0') + allcalls.count('1/1') + allcalls.count('0/1') )
	
		if prop_het > maxprophet:
#			print 'hetfail', allcalls
			calltest+=1
		else:
			pass	
	return calltest, newformatfields, newgenofields, REFct, ALTct
	


## Get the sample names from the vcf
samples=[]
for line in inVCF:
	if line.startswith('##'):
		pass
	else:
		for i in line.split()[9:]: samples.append(i)
		break

# this sets the file back to the beginning
inVCF.seek(0)

# this adds the filters i apply later into the header so gatk can use the file (hopefully)
for line0 in inVCF:
	if line0.startswith('#'):
		if line0.startswith('##FORMAT'):
			outVCF.write('##FILTER=<ID=FAILgt,Description="Low quality">\n') 
			outVCF.write('##FORMAT=<ID=FT,Number=.,Type=String,Description="Genotype-level filter">\n')
			outVCF.write(line0)
			break
		else:
			outVCF.write(line0)

# this then prints out the rest of the header lines	
for line0 in inVCF:
#	print '&\t', line0
	if line0.startswith('#'): 
		outVCF.write(line0); continue

### For all other lines:

	line=line0.strip().split('\t')

	#CHROM0	POS1	ID2	REF3	ALT4	QUAL5	FILTER	INFO	FORMAT	EW1	EW10	EW2	EW3	EW4	EW5	EW6	EW7	EW8	EW9
	myref=line[3]
	myalt=line[4]
	myqual=line[5]
	myfilter=line[6]
	myinfo=line[7]
	myformat=line[8]
	mygenoinfo=line[9:]
	
	# Check for site level errors
	# Write to an error file these unexpected things
	
	if myref not in refbases:
		errorVCF.write('# ref not AGCT\n')
		errorVCF.write(line0)
		counter_nop+=1
	elif myalt not in altbases:
		errorVCF.write('# alt not . or AGCT\n')
		errorVCF.write(line0)
		counter_nop+=1
	elif myqual == '.':
		errorVCF.write('# QUAL = .\n')
		errorVCF.write(line0)
		counter_nop+=1
	elif myfilter != 'PASS':
		errorVCF.write('# FILTER is not PASS\n')
		errorVCF.write(line0)
		counter_nop+=1	
	elif 'DP' not in myinfo or 'AN' not in myinfo:
		errorVCF.write('# site info has no DP or AN\n')
		errorVCF.write(line0)
		counter_nop+=1	
	elif 'GT' not in myformat or 'AD' not in myformat or 'DP' not in myformat or 'GQ' not in myformat:
		errorVCF.write('# GT, AD, DP or GQ (which also does RGQ) not in FORMAT field\n')
		errorVCF.write(line0)
		counter_nop+=1
	else:

#		print '>> starting geno filter << '	
		# this function redoes the format to add the FT flag
		# it also filteres the genotypes based on DP
		# finally after applying the DP filter to it, count the number of ./. to see if you want to filter it

		
		call_test, newFORMAT, newGENOS, ctREFalleles, ctALTalleles = gfilter_DP(mygenoinfo, myformat)
#		print 'calltest return', call_test,newFORMAT, newGENOS, ctREFalleles, ctALTalleles		
		## Get the INFO updated
		myinfo=line[7].split(';')
		f=dict(s.split('=') for s in myinfo)
	

		if call_test == 0:
			# have to do this here otherwise errors with zero divisions for lines that fail
			
#			print 'in call test 0'
#			print ctALTalleles,ctREFalleles,
			f['AC']=ctALTalleles
			f['AN']=ctREFalleles+ctALTalleles
			f['AF']=round(float(ctALTalleles)/(float(ctREFalleles)+float(ctALTalleles)), 4)
			# need to update this - should be . for a non variant
			if f['AC'] == 0:
				f['AC'] = '.'
			if f['AF'] == 0:
				f['AF'] = '.'
			
			mynewline=''

			## This is the bit that changed this time around.
			# chr38	4162	.	T	.	30.98	PASS	AC=.;AF=.;AN=16;DP=254;InbreedingCoeff=-0.0527	GT:AD:DP:FT:RGQ	0/0:29,0:29:PASS:63	0/0:17,0:17:PASS:39	0/0:27,0:27:PASS:60	0/0:30,0:30:PASS:69	./.:21,5:26:LowRGQ0:0	0/0:27,0:27:PASS:63	0/0:30,0:30:PASS:78	./.:26,3:29:LowRGQ0:0	0/0:20,1:21:PASS:28	0/0:18,0:18:PASS:45
			# check if the variable site is still variable after filtering.
			#CHROM0[0]	POS1[1]	ID2	REF3	ALT4	QUAL5	FILTER6	INFO	FORMAT
			if ctALTalleles == 0 and myalt != '.':
#				print 'handling ALT that is no longer alt'
#				print 'old', line
#				print 'new', newGENOS
				newALT='.'
				mynewline += '\t'.join(line[:4])
				mynewline += '\t'
				mynewline += newALT
				mynewline += '\t'
				mynewline += '\t'.join(line[5:7])
				mynewline += '\t'
#				print 'newline', mynewline, newFORMAT
			else:
				# add everything up to INFO is the same (up to but not includign INFO)
				mynewline += '\t'.join(line[:7])
				mynewline += '\t'				


			# add the INFO line
			mynewline += ';'.join('{0}={1}'.format(key, val) for key, val in sorted(f.items()))
			mynewline += '\t'
			mynewline += newFORMAT
			mynewline += '\t'
			mynewline += '\t'.join(newGENOS)
			mynewline += '\n'
			
#			print 'oldline', line0
#			print 'newline', mynewline
#			print mynewline
			outVCF.write(mynewline)
#			print 'printing out'
			counter_p+=1			
		else:
#			print 'notprint'
#			print line
			counter_nop+=1

#			print 'NOT writing file - SKIPPING line - failed het or maxnocall'
						

errorVCF.write('## This many lines failed a filter and werent printed' + '\t' + str(counter_nop) + '\n')
errorVCF.write('## This many lines PASSED and were printed' + '\t' + str(counter_p) + '\n')
errorVCF.close()
outVCF.close()

print 'done'
##	GS1	GS10	GS2	GS3	GS4	GS5	GS6	GS7	GS8	GS9