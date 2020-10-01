'''
Isolates a sub-shower from the shower plot and writes it into a new .dot file.
'''
import re


''' 
Edit this:
'''
targetID = ["936"]
showerdepth = 3 #must be between 1 and 8

fin = open("CovidMC_int_0.dot","r")
fout = open("CovidSubshower_0b.dot","w")





''' 
Do not edit below.
'''

infectees1 = []
infectees2 = []
infectees3 = []
infectees4 = []
infectees5 = []
infectees6 = []
infectees7 = []
infectees8 = []
infectees9 = []
outtext =[]
days = []
strings = [""]
''' 
Utility to construct a search string for grep to use
'''
def MakeSearchString(items,strings_):
	searchstring = ""
	if len(items) < 1:
		searchstring = ""
	elif len(items) > 1:
		searchstring = "("
		counter = 0
		for id in items:
			if (counter == 0):
				searchstring = searchstring + str(id)
				counter = counter+1
			else:
				searchstring = searchstring + "|" + str(id)
		searchstring = searchstring + ")"
	else:
		 searchstring = str(items[0])
	strings_[0] = searchstring
		 
strings = [""]
''' 
Go through dot file and find all the infectees belonging to all the people in infector array.
If allprevious has content, this is the last generation and we make sure not to add connections
from everyone covered so far to people in the last generation.
'''
def FindInfectees(infector, infectee, outtext, allprevious):
	MakeSearchString(infector, strings)
	searchstring = strings[0]
	print searchstring
	fin.seek(0)
	for line in fin:
		line = line.rstrip()
		# search for connections
		result = re.findall(" " + searchstring + " ->", line) 
		if len(allprevious) >0:
			MakeSearchString(allprevious,strings)
			sstring = strings[0]
			result = result + re.findall("!" + sstring + " -> " + searchstring + " ", line) 
		if len(result) > 0:
			outtext.append(line + "\n")
			if re.search("616161", line): # #616161 is the color for "infects"
				lst = re.findall('-> (\d*)', line)
				infectee.append(lst[0])

''' 
The rank determins which person appears at the height of which day.
''' 
def FindRank(ids, outtext,strings):
	fin.seek(0)
	MakeSearchString(ids, strings)
	searchstring = strings[0]
	for line in fin:	
		#search for days { rank=same 0 D0 }
		result2 = re.findall("rank=same " + searchstring + " D(\d*) }", line) 
		if len(result2) > 0:
			outtext.append(line)
			days.append(int(result2[0][1]))	
		result3 = re.findall("rank=same " + searchstring + " Dh(\d*) }", line) 
		if len(result3) > 0:
			outtext.append("{ rank=same " + result3[0][0] + " D"+result3[0][1] + "}")
			days.append(int(result3[0][1]))	

''' 
Each person is represented as a 'cluster' in dot. This finds the clusters belonging to 
every person in allids.
''' 
def FindPeopleCluster(allids, outtext):
	fin.seek(0)
	searchstring = "subgraph cluster_("
	MakeSearchString(allids, strings)
	searchstring = searchstring + strings[0] + ") {"
	findclose = 0
	for line in fin:
		if (findclose is 1):
			outtext.append(line)
			if re.search("} ", line):
				findclose = 0
		elif re.search(searchstring, line):
			outtext.append(line)
			findclose = 1 

''' 
Find the entry for each day between the minimum and maximum in days array.
''' 	
def FindDays(days, outtext):
	fin.seek(0)
	'''
	D0 -> Dh0 [constraint=true style="invis"] 
	  Dh0 -> D1 [constraint=true style="invis"] 
	  D0 [ label = "Day\n0" shape=square style=filled, fillcolor="gray"label="Days" ]
	 	'''
	if (len(days) == 0):
	 	return
	firstday = days[0] 
	lastday = days[-1]
	outtext.append("subgraph days{ \n") 
	daylist = range(firstday, lastday+1)
	MakeSearchString(daylist, strings)
	searchstring = strings[0]
	for id in daylist:
		outtext.append("D" + str(id) + "-> D" +str(id+1) + "\n")
	
	for line in fin:
		if re.search(" D" + searchstring + " \[ ", line):	
			outtext.append(line)	
	outtext.append("}\n") 
	

''' 
Now go down generation by generation, identifying the infectees.
''' 	
previous = []	
lastonestart = []
FindInfectees(targetID, infectees1, outtext, previous)
if showerdepth > 1:
	FindInfectees(infectees1, infectees2, outtext, previous)
	lastonestart = infectees2
if showerdepth > 2:
	FindInfectees(infectees2, infectees3, outtext, previous)
	lastonestart = infectees3
if showerdepth > 3:	
	FindInfectees(infectees3, infectees4, outtext, previous)
	lastonestart = infectees4	
if showerdepth > 4:
	FindInfectees(infectees4, infectees5, outtext, previous)
	lastonestart = infectees5	
if showerdepth > 5:
	FindInfectees(infectees5, infectees6, outtext, previous)
	lastonestart = infectees6	
if showerdepth > 6:
	FindInfectees(infectees6, infectees7, outtext, previous)
	lastonestart = infectees7	
if showerdepth > 7:
	FindInfectees(infectees7, infectees8, outtext, previous)
	lastonestart = infectees8
previous = targetID + infectees1 + infectees2 + infectees3 + infectees4 + infectees5 + infectees6 + infectees7 +infectees8 
FindInfectees(lastonestart, infectees9, outtext, previous)
all = previous + infectees9

FindPeopleCluster(all,outtext)
FindRank(all, outtext, strings)
days.sort()
FindDays(days, outtext)

fin.close()

# remove arrows going to out of scope people
MakeSearchString(all, strings)
sstring = strings[0]
for item in outtext:
	if re.search("days",item): 
		break
	if not re.search("->", item):
		continue
	if len(re.findall("-> " + sstring, item)) < 1: #arrows going toward people not in the shower
		outtext.remove(item)
	if len(re.findall(sstring + " -> ", item)) < 1: #arrows going out from people not in the shower
		outtext.remove(item)


fout.write("digraph infection {\n")
fout.write("newrank=true; \n")
fout.write("compound=true; \n")
fout.write('node [fontname="Helvetica" fontsize=18];\n')
for li in outtext:
	fout.write(li)
fout.write("}\n")	
fout.close
		