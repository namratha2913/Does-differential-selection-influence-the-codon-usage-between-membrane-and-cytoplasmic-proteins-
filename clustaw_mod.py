# -*- coding: utf-8 -*-
"""
Created on Sat Jan 13 21:08:28 2024

@author: Legion
"""
import os
import sys
import csv
import re



main_fold = os.listdir('filter_membrane_align')
#print(main_fold)

for name in main_fold:
    
	file = open('filter_membrane_align/'+str(name))
    
	flag =1
	gene_count = 0
    
	final_lines = []
    
	space_line = 0
	flag = 0
    
	for l in file:
			
		file_write = open('membrane_mod_align/mod_'+str(name), 'a')
		sys.stdout = file_write
		
		#print(l)
		
		if l.startswith("lcl"):
			
			#print("1")
			last = l.split()
			last_len = len(l)
			print(l[:-1])
			flag =1
			
		elif l.startswith("  ") and flag:
			#print(flag)
		
			#print(l)
			#print("2")
			l1 = l.split("                                              ")
			space_line = 1

			
			if  space_line:
				#print(l)
				replace_len =len(last[1]) 
				print(last[0][:-3]+'ZZZ'+str(l)[len(last[0]):(last_len-replace_len)-1]+(str(l)[(last_len-replace_len)-1:-1].replace(" ", 'x')))
				print(l[:-1])
				space_line =0

				flag = 0

		   
		else:
			#print("3")
			space_line = 0
			flag =0
			print(l[:-1])
		    
		file_write.close()


main_fold = os.listdir('filter_cyto_align')
#print(main_fold)

for name in main_fold:
    
	file = open('filter_cyto_align/'+str(name))
    
	flag =1
	gene_count = 0
    
	final_lines = []
    
	space_line = 0
	flag = 0
    
	for l in file:
			
		file_write = open('cyto_mod_align/mod_'+str(name), 'a')
		sys.stdout = file_write
		
		#print(l)
		
		if l.startswith("lcl"):
			
			#print("1")
			last = l.split()
			last_len = len(l)
			print(l[:-1])
			flag =1
			
		elif l.startswith("  ") and flag:
			#print(flag)
		
			#print(l)
			#print("2")
			l1 = l.split("                                              ")
			space_line = 1

			
			if  space_line:
				#print(l)
				replace_len =len(last[1]) 
				print(last[0][:-3]+'ZZZ'+str(l)[len(last[0]):(last_len-replace_len)-1]+(str(l)[(last_len-replace_len)-1:-1].replace(" ", 'x')))
				print(l[:-1])
				space_line =0

				flag = 0

		   
		else:
			#print("3")
			space_line = 0
			flag =0
			print(l[:-1])
		    
		file_write.close()		

	 
	    
	
	    





       

             
