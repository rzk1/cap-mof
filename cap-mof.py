import sys, os, subprocess, random, math
from sys import argv
from math import sqrt, floor, ceil
import numpy as np

rTypical = 2.15423 # Covalent radius = [rRb(2.1) + rO(0.66)]?
r2Typical = rTypical**2
rOOTypical = 1.15 # Van der Waal radius, decrease if this causes problem 
rOO2Typical = rOOTypical ** 2
NOxTarget = 8

def distance(coord1, coord2, period):
 d = float(coord1) - float(coord2)
 d -= round(d/period)*period
 return d

def pol2cart(r, phi, theta):
    x = str(r * np.sin(theta) * np.cos(phi))
    y = str(r * np.sin(theta) * np.sin(phi))
    z = str(r * np.cos(theta))
    str_xyz = x + " " + y + " " + z
    return str_xyz

# =============================================================
# ------------- the main program ------------------------------
# =============================================================

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
convert = [0,0,0]

with open("C:/Users/chp19/OneDrive/Documents/mofs/mof-coord","r") as fp:
 for line in fp:
     
  linenumber += 1
  
  # split the line to extract element symbol and coordinates
  if linenumber == 1:
   line1 = int(line)
  elif linenumber == 2:
   ABC[0], ABC[1], ABC[2] = line.split()
  else:
   xyzStr[0], xyzStr[1], xyzStr[2], xyzStr[3] = line.split()
   
   # save zero-image coordinates
   if xyzStr[0]== "Rb":
    cations.append(xyzStr[:])
   elif xyzStr[0]== "O":
    oxygens.append(xyzStr[:])
   else:
    others.append(xyzStr[:])
fp.close()

NCations = len(cations)
NOxygens = len(oxygens)

if NCations % 12 != 0:  # There are total 36 Rb
  print ("Wrong number of cations: ", NCations)
  sys.exit(2)
with open("C:/Users/chp19/OneDrive/Documents/mofs/new-coord", "w") as f:
 for ication in range(NCations):
  count = 0
  neighbors.clear()
  newOx = ['O',0,0,0]
  for ioxygen in range(NOxygens):
    r2=0.0
    for icart in range(1,4): # 0th col is character 
     # impose minimum image first
     dCO = distance(float(cations[ication][icart]),float(oxygens[ioxygen][icart]),float(ABC[0]))
     r2 += (dCO**2) / 2
    if (r2 < r2Typical): # save this oxygen into a separate array called neighbors
     neighbors.append(oxygens[ioxygen])
     count += 1
  # add oxygen
  if (count < NOxTarget):
   Naccepted = 0
   while (Naccepted + count < NOxTarget):
    # add 1 atom at a time, check if it's a good one with appropriate
    # distance between Rb and old Oxygens,then keep adding until the array is full (8-coord)
    # generate two random angles, use rTypical, compute xyz
    phi = random.uniform(0, 2*math.pi)
    theta = random.uniform(0, math.pi)
    convert = pol2cart(rTypical, phi, theta).split(" ")
    
    # add coordinates back after assumption of (0,0,0) cat coord
    for icart in range(1,4):
     newOx[icart] = round(float(convert[icart-1])+ float(cations[ication][icart]),7)
     
    # check overlap with other neighbor oxygens
    for iNeigh in range(len(neighbors)):
       rOO2 = 0.0
       for icart in range(1,4):
          dOO = distance(float(neighbors[iNeigh][icart]),float(newOx[icart]),float(ABC[0]))
          rOO2 += (dOO**2) / 2
    if (rOO2 > rOO2Typical): 
        Naccepted += 1
        neighbors.append(newOx[:])
        addOx.append(newOx[:])
 print(len(addOx)) # the value should be 24*6 = 144, as there are 24 Rb missing 6 bonds
  #print new XYZ file
 s = "     "
 print (line1 + len(addOx),file = f)
 print (ABC[0],s, ABC[1],s, ABC[2],s, file = f)
 for i in range(NCations):
   print(*cations[i], sep = "          ", file = f)
 for j in range(NOxygens):
   print(*oxygens[j], sep = "          ", file = f)
 for k in range(len(addOx)):
   print(*addOx[k], sep = "         ", file = f)
 for l in range(len(others)):
   print(*others[l], sep = "         ", file = f) 
f.close()
