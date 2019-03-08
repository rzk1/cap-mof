import sys, os, subprocess, random, math
from sys import argv
from math import sqrt, floor, ceil
import numpy as np

dRbOTypical = 2.92 # Ionic radii: rRb(1.66) + rO(1.26)
dOOTypical = 2*1.47 # Van der Waal radius: rO(1.52)
dRbOCutoff = 1.15 * dRbOTypical
dOXCutoff = # Covalent radii: rOx + rC
dOHTypical = 1.0 # 
dRbHTypical = # Van der Waals: rRb + rH
dHHTypical = 2* # van der Waals: rH 
NOxTarget = 8

d2RbOTypical = dRbOTypical ** 2
d2OOTypical = dOOTypical ** 2
d2RbOCutoff = dRbOCutoff ** 2
d2OXCutoff = dOXCutoff ** 2
d2OHTypical = dOHTypical ** 2 
d2RbHTypical = dRbHTypical ** 2
d2HHTypical = dHHTypical ** 2

def distance(coord1, coord2, period):
 d = float(coord1) - float(coord2)
 d -= round(d/period)*period
 return d

# RZK: it would be more elegant to return a list x,y,z of floats, not just one string
def pol2cart(r, phi, theta):
    #x = str(r * np.sin(theta) * np.cos(phi))
    #y = str(r * np.sin(theta) * np.sin(phi))
    #z = str(r * np.cos(theta))
    #str_xyz = x + " " + y + " " + z
    #return str_xyz
    x = r * np.sin(theta) * np.cos(phi)
    y = r * np.sin(theta) * np.sin(phi)
    z = r * np.cos(theta)
    return (x,y,z)

# =============================================================
# ------------- the main program ------------------------------
# =============================================================
script, inputfile, action = argv
# action = 0 : add new oxygen atoms at Rb
# action = 1 : add new hydrogen atoms at O

# read and process atomic coordinates
# loop over the lines of the coordinate file
linenumber = 0
xyzStr = ["","","",""]
ABC = [0, 0, 0]
cations = []
oxygens = []
others = []
neighbors = []
addOx = [] # store all added Oxygens
addHy = [] # store all added Oxygens
convert = [0,0,0]

with open(inputfile,"r") as fp:
 for line in fp:
     
  linenumber += 1
  
  # split the line to extract element symbol and coordinates
  if linenumber == 1:
   line1 = int(line)
  elif linenumber == 2:
   ABC[0], ABC[1], ABC[2] = line.split()
   for icart in range(3):
    ABC[icart]=float(ABC[icart])
  else:
   xyzStr[0], xyzStr[1], xyzStr[2], xyzStr[3] = line.split()
   
   # save zero-image coordinates
   if xyzStr[0] == "Rb":
    cations.append(xyzStr[:])
   elif xyzStr[0] == "O":
    oxygens.append(xyzStr[:])
   else:
    others.append(xyzStr[:])

# close the input file
fp.close()

NCations = len(cations)
NOxygens = len(oxygens)
convert=[0.0,0.0,0.0]

if NCations % 24 != 0:
 print ("Wrong number of cations: ", NCations)
 sys.exit(2)

# RZK: the output file should be opened later, not here
#with open(inputfile+".aug.xyz", "w") as f:
for ication in range(NCations):
 count = 0
 neighbors.clear()
 newOx = ['O',0,0,0]
 newHy = ['H',0,0,0]

 # RZK: more comments are necessary
 # check which oxygen atoms are neighbors
 for ioxygen in range(NOxygens):
  d2=0.0
  for icart in range(1,4): # 0th col is character 
   # impose minimum image first
   # RZK: there is a bug in the following line
   #dCO = distance(float(cations[ication][icart]),float(oxygens[ioxygen][icart]),float(ABC[0]))
   dCO = distance(float(cations[ication][icart]),float(oxygens[ioxygen][icart]),float(ABC[icart-1]))
   # RZK: division by two is not necessary
   #d2 += (dCO**2) / 2
   d2 += dCO**2
  if (d2 < d2RbOCutoff): # save this oxygen into a separate array called neighbors
   neighbors.append(oxygens[ioxygen])
   count += 1

 if (action == 1):
  for ioxygen in range(len(neighbors)):
 
   for jatom in range(len(others)):
    d2=0.0
    for icart in range(1,4): # 0th col is character 
     dCO = distance(float(neighbors[ioxygen][icart]),float(others[jatom][icart]),float(ABC[icart-1]))
     # RZK: division by two is not necessary
     #d2 += (dCO**2) / 2
     d2 += dCO**2
    if (d2 < d2OXCutoff):
     # goto another oxygen atom
     break

   # this oxygen does not have "other" neighbors
   # start adding hydrogens, add one hydrogen to the first oxygen
   # all other oxygens get two hydrogens
   if (ioxygen == 0):
    nH = 1
   else:
    nH = 2

   Naccepted = 0
   Ntrials = 1
   while (Naccepted < nH):

    phi = random.uniform(0, 2*math.pi)
    theta = random.uniform(0, math.pi)
    convert[:] = pol2cart(dOHTypical, phi, theta)
    
    # add coordinates of the center 
    for icart in range(1,4):
     newHy[icart] = float(convert[icart-1]) + float(neighbors[ioxygn][icart])
     
    # check overlap with other neighbor oxygens
    badContactFound = False
    for iNeigh in range(len(neighbors)):
     if (iNeigh == ioxygen):
      continue    
     dOO2 = 0.0
     for icart in range(1,4):
      dOO = distance(float(neighbors[iNeigh][icart]),float(newHy[icart]),float(ABC[icart-1]))
      dOO2 += dOO**2

     if (dOO2 < d2OHTypical): 
      badContactFound = True
      # no need to consider other oxigen neighbors
      break

    for iNeigh in range(Ncations):
     dOO2 = 0.0
     for icart in range(1,4):
      dOO = distance(float(cations[iNeigh][icart]),float(newHy[icart]),float(ABC[icart-1]))
      dOO2 += dOO**2

     if (dOO2 < d2RbHTypical): 
      badContactFound = True
      # no need to consider other oxigen neighbors
      break

    for iNeigh in range(len(addHy)):
     dOO2 = 0.0
     for icart in range(1,4):
      dOO = distance(float(addHy[iNeigh][icart]),float(newHy[icart]),float(ABC[icart-1]))
      dOO2 += dOO**2

     if (dOO2 < d2HHTypical): 
      badContactFound = True
      # no need to consider other oxigen neighbors
      break

    if (not badContactFound):
     Naccepted += 1
     addHy.append(newHy[:])
     print("Cation %3d. Oxygen %3d. Attempts: %6d" % (ication, ioxygen, Ntrials))   

    Ntrials += 1


 if (action == 0):
  # add new oxygen neighbors if this cation does not have enough neighbors
  if (count < NOxTarget):
   Naccepted = 0
   Ntrials = 1
   while (Naccepted + count < NOxTarget):

    # generate random coordinates (trail insertion)
    # add 1 atom at a time, check if it's a good one with appropriate
    # distance between Rb and old Oxygens,then keep adding until the array is full (8-coord)
    # generate two random angles, use dRbOTypical, compute xyz
    phi = random.uniform(0, 2*math.pi)
    theta = random.uniform(0, math.pi)
    #RZK: convert = pol2cart(dRbOTypical, phi, theta).split(" ")
    convert[:] = pol2cart(dRbOTypical, phi, theta)
    
    # add coordinates back after assumption of (0,0,0) cat coord
    for icart in range(1,4):
     # RZK: is rounding really necessary?
     #newOx[icart] = round(float(convert[icart-1])+ float(cations[ication][icart]),7)
     newOx[icart] = float(convert[icart-1]) + float(cations[ication][icart])
     
    # check overlap with other neighbor oxygens
    badContactFound = False
    for iNeigh in range(len(neighbors)):
     # RZK: your code has different identations - bad style, especially in python
     dOO2 = 0.0
     for icart in range(1,4):
      # RZK: same bug dOO = distance(float(neighbors[iNeigh][icart]),float(newOx[icart]),float(ABC[0]))
      dOO = distance(float(neighbors[iNeigh][icart]),float(newOx[icart]),float(ABC[icart-1]))
      # RZK: same bug dOO2 += (dOO**2) / 2
      dOO2 += dOO**2

     # RZK: this is a serious bug - incorrect identation, breaks the logic of the code
    #if (dOO2 > d2OOTypical): 
     if (dOO2 < d2OOTypical): 
      badContactFound = True
      # no need to consider other oxigen neighbors
      break

    if (not badContactFound):
     Naccepted += 1
     neighbors.append(newOx[:])
     addOx.append(newOx[:])
     print("Cation %3d. Oxygen neighbors (old, new): %4d%4d. Attempts: %6d" % (ication,count,Naccepted,Ntrials))   

    Ntrials += 1

# RZK: xyz file is incorrect! It contains two periodic images of the same Rb atoms!
# There should be 24 Rb atoms total, not 36
print("There are %4d new oxygen atoms " % len(addOx)) # the value should be 24*6 = 144, as there are 24 Rb missing 6 bonds

f = open(inputfile+".aug.xyz", "w")
print (line1 + len(addOx), file = f)
# RZK: learn how to print out formatted numbers in python 
print ("%15.8f%15.8f%15.8f" % (ABC[0],ABC[1],ABC[2]), file = f)
for i in range(NCations):
 # RZK: format these lines properly as well
 print(*cations[i], sep = "          ", file = f)
 #print("%3s%15.8f%15.8f%15.8f" % *cations[i], file = f)
for j in range(NOxygens):
 print(*oxygens[j], sep = "          ", file = f)
if (action == 0):
 for k in range(len(addOx)):
  print(*addOx[k], sep = "         ", file = f)
else:
 for k in range(len(addHy)):
  print(*addHy[k], sep = "         ", file = f)
for l in range(len(others)):
 print(*others[l], sep = "         ", file = f) 
f.close()

