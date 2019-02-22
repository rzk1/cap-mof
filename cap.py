import sys, os, subprocess, random
from sys import argv
from math import sqrt, floor, ceil

rTypical = 1.7 # typical distance between the cation and oxygen
r2Typical = rTypical**2


def printAtom(record):
 print("%4s%20s%20s%20s"%(record[0],record[1],record[2],record[3]))

# =============================================================
# ------------- the main program ------------------------------
# =============================================================
# read command line arguments
script, infile, inseed = argv

# seed produces fully deterministic behavior 
random.seed(inseed)

# =============================================================
# read and process atomic coordinates
# loop over the lines of the coordinate file
linenumber = 0
xyzStr = ["","","",""]
ABC = ["","",""]
cations = []
oxygens = []
others = []

with open(infile) as fp:
 for line in fp:
     
  linenumber += 1
  
  # split the line to extract element symbol and coordinates
  if linenumber==1:
   line1=int(line)
  elif linenumber==2:
   ABC[0],ABC[1],ABC[2] = line.split()
  else:
   xyzStr[0],xyzStr[1],xyzStr[2],xyzStr[3] = line.split()
   
   # save zero-image coordinates
   if xyzStr[0]=="Rb":
    cations.append(xyzStr[:])
   elif xyzStr[0]=="O":
    oxygens.append(xyzStr[:])
   else:
    others.append(xyzStr[:])

fp.close()

NCations=len(cations)
NOxygens=len(oxygens)
if NCations%24 != 0:
 print ("Wrong number of cations: ", NCations)
 sys.exit(2)

for ication in range(NCations):
 for ioxygen in range(NOxygens):
   r2=0.0
   for icart in range(3):
Impose minimum image first
    r2 += (float(cations[ication][icart])-float(oxygens[ioxygen][icart]))**2
    if (r2 < r2Typical):
     # save this oxygen into a separate array called neighbors

# print a new XYZ file
print (line1-NVacancies)
print (line2.rstrip())
for iatom in range(NTeAtoms):
 printAtom(Te_atoms[iatom])
for iatom in range(NSbAtoms):
 printAtom(Ge_atoms[iatom])
for iatom in range(NSbAtoms):
 printAtom(Sb_atoms[iatom])

