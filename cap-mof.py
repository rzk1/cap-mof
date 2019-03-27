import sys, os, subprocess, random, math
from sys import argv
from math import sqrt, floor, ceil
import numpy as np

# Typical "insertion" distance
dRbOInsertion = 2.92 # (insertion) Ionic radii: rRb(1.66) + rO(1.26)
dOHInsertion = 0.90 # (insertion) Covalent: rOx(0.63) + rH(0.32), tune from VMD
# Typical "neighbor" ("cutoff") distance 
dRbONeighbor = 1.15 * dRbOInsertion # (neighbor) From VMD
dOONeighbor = 2 * 1.47 # (neighbor) Van der Waal radius: rO(1.52), tune from VMD
dOXNeighbor = 1.6 # (neighbor) Covalent radii: rOx(0.63) + rC(0.75)  ->1.38
dRbHNeighbor = 3.00 # (neighbor) Van der Waals: rRb(3.03) + rH(1.1), tune
dHHNeighbor = 1.2 # (neighbor) this should be a bit less than H-H distance in a water molecule (1.52)

NOxTarget = 8

d2RbOInsertion = dRbOInsertion ** 2
d2OHInsertion = dOHInsertion ** 2
d2OONeighbor = dOONeighbor ** 2
d2RbONeighbor = dRbONeighbor ** 2
d2OXNeighbor = dOXNeighbor ** 2
d2RbHNeighbor = dRbHNeighbor ** 2
d2HHNeighbor = dHHNeighbor ** 2

def distance(coord1, coord2, period):
 d = float(coord1) - float(coord2)
 d -= round(d/period)*period
 return d

def pol2cart(r, phi, theta):
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
action = int(action)

# read and process atomic coordinates
# loop over the lines of the coordinate file
linenumber = 0
xyzStr = ["","","",""]
ABC = [0, 0, 0]
cations = []
oxygens = [] # Store only oxygens have at least 1 neighbor
oxygensInitial = [] # Store all oxygens 
others = []
neighbors = []
addOx = [] # store all added Oxygens
addHy = [] # store all added Hydrogens
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
    oxygensInitial.append(xyzStr[:])
   else:
    others.append(xyzStr[:])

# close the input file
fp.close()

# Split all oxygens into 2 categories: with and without neighbor(s)
# Check with cations neighbors
# If it has neighbor, save its coord and go to the next one 
for ioxygen in range(len(oxygensInitial)):
 nextOxygen = False
 for iNeigh in range(len(cations)):
  dORb2 = 0.0
  for icart in range(1,4):
   dORb = distance(float(cations[iNeigh][icart]),float(oxygensInitial[ioxygen][icart]),float(ABC[icart-1]))
   dORb2 += dORb**2
  if (dORb2 < d2RbONeighbor):
   oxygens.append(oxygensInitial[ioxygen])
   nextOxygen = True
   break 
 if(nextOxygen):
  continue
# Check with others neighbors
# If it has neighbor, save its coord and go to the next one
 for iNeigh in range(len(others)):
  dOX2 = 0.0
  for icart in range(1,4):
   dOX = distance(float(others[iNeigh][icart]),float(oxygensInitial[ioxygen][icart]),float(ABC[icart-1]))
   dOX2 += dOX**2
  if (dOX2 < d2OXNeighbor):
   oxygens.append(oxygensInitial[ioxygen])
   nextOxygen = True
   break
 if(nextOxygen):
  continue

NCations = len(cations)
NOxygens = len(oxygens) # Use oxygens that have neighbor(s)
convert=[0.0, 0.0, 0.0]

if NCations % 24 != 0:
 print ("Wrong number of cations: ", NCations)
 sys.exit(2)


for ication in range (NCations):
 count = 0
 neighbors.clear()
 newOx = ['O',0,0,0] # store newOx per cation
 newHy = ['H',0,0,0] # store newHy per cation

 # find all oxygen neighbors of this current cation
 for ioxygen in range(NOxygens):
  d2=0.0
  for icart in range(1,4): # 0th col is character
   # impose minimum image first
   dCO = distance(float(cations[ication][icart]),float(oxygens[ioxygen][icart]),float(ABC[icart-1]))
   d2 += dCO**2
  if (d2 < d2RbONeighbor): # save this oxygen into a separate array called neighbors
   neighbors.append(oxygens[ioxygen])
   count += 1 
 
 if (action == 1): # add H to Oxygens with no "other" neighbors
  noNeighborCount = 0
  for ioxygen in range(len(neighbors)): # ioxygen refers to a neigbor of our current cation
   OtherNeighborFound = False
   for jatom in range(len(others)): # check if this O has "other" neighbors besides cation
    d2=0.0
    for icart in range(1,4):
     dCO = distance(float(neighbors[ioxygen][icart]),float(others[jatom][icart]),float(ABC[icart-1]))
     d2 += dCO**2
    if (d2 < d2OXNeighbor):
     OtherNeighborFound = True
     break  # no need to consider other neighbors
   if (OtherNeighborFound):
    print("No H added on Oxygen %d, cation %d" % (ioxygen, ication))
    continue # go to another oxygen atom
    
   noNeighborCount += 1
   if (noNeighborCount == 1): # this oxygen does not have "other" neighbors
    nH = 1 # add 1 H to the first oxygen
   else:
    nH = 2 # add 2 H to the rest of oxygens
    
   Naccepted = 0
   Ntrials = 1
   while (Naccepted < nH):
    # generate random H
    # generate two random angles, use dOHInsertion, compute xyz
    phi = random.uniform(0, 2*math.pi)
    theta = random.uniform(0, math.pi)
    convert[:] = pol2cart(dOHInsertion, phi, theta)

    # add coordinates of the center
    for icart in range(1,4):
     newHy[icart] = float(convert[icart-1]) + float(neighbors[ioxygen][icart])

    # check overlap with other oxygen neighbors in cases O-H-O
    badContactFound = False
    for iNeigh in range(len(neighbors)):
     if (iNeigh == ioxygen): # exclude checking oxygen itself
      continue
     dOO2 = 0.0
     for icart in range(1,4):
      dOO = distance(float(neighbors[iNeigh][icart]),float(newHy[icart]),float(ABC[icart-1]))
      dOO2 += dOO**2
     if (dOO2 < d2OHInsertion):
      badContactFound = True
      # no need to consider other oxygen neighbors
      break

    if (not badContactFound):
     # check overlap with other cation neighbors in cases O-H-Rb
     for iNeigh in range(NCations):
      dOO2 = 0.0
      for icart in range(1,4):
       dOO = distance(float(cations[iNeigh][icart]),float(newHy[icart]),float(ABC[icart-1]))
       dOO2 += dOO**2
      if (dOO2 < d2RbHNeighbor):
       badContactFound = True
       # no need to consider other cation neighbors
       break

    if (not badContactFound):
     # check overlap with added hydogen neighbors in cases O-H-H-O
     for iNeigh in range(len(addHy)):
      dOO2 = 0.0
      for icart in range(1,4):
       dOO = distance(float(addHy[iNeigh][icart]),float(newHy[icart]),float(ABC[icart-1]))
       dOO2 += dOO**2
      if (dOO2 < d2HHNeighbor):
       badContactFound = True
       # no need to consider other hydrogen neighbors
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
    convert[:] = pol2cart(dRbOInsertion, phi, theta)
    
    # add coordinates of the center
    for icart in range(1,4):
     newOx[icart] = float(convert[icart-1]) + float(cations[ication][icart])
     
    # check overlap with other oxygens neighbor
    badContactFound = False
    for iNeigh in range(len(neighbors)):
     dOO2 = 0.0
     for icart in range(1,4):
      dOO = distance(float(neighbors[iNeigh][icart]),float(newOx[icart]),float(ABC[icart-1]))
      dOO2 += dOO**2
 
     if (dOO2 < d2OONeighbor): 
      badContactFound = True
      # no need to consider other oxygen neighbors
      break

    if (not badContactFound):
     Naccepted += 1
     neighbors.append(newOx[:])
     addOx.append(newOx[:])
     print("Cation %3d. Oxygen neighbors (old, new): %4d%4d. Attempts: %6d" % (ication,count,Naccepted,Ntrials))   

    Ntrials += 1      
print("There are %4d new oxygen atoms " % len(addOx))
print("There are %4d new hydrogen atoms " % len(addHy))

f = open(inputfile+ "%d" % (action)+".xyz", "w")
print (line1 + len(addOx) + len(addHy), file = f)
print ("%15.8f%15.8f%15.8f" % (ABC[0],ABC[1],ABC[2]), file = f)
for i in range(NCations):
 print("%3s%15.8f%15.8f%15.8f" % (cations[i][0],float(cations[i][1]),float(cations[i][2]),float(cations[i][3])), file = f)
for j in range(NOxygens):
 print("%3s%15.8f%15.8f%15.8f" % (oxygens[j][0],float(oxygens[j][1]),float(oxygens[j][2]),float(oxygens[j][3])), file = f)
for k in range(len(others)):
 print("%3s%15.8f%15.8f%15.8f" % (others[k][0],float(others[k][1]),float(others[k][2]),float(others[k][3])), file = f)
for l in range(len(addOx)):
 print("%3s%15.8f%15.8f%15.8f" % (addOx[l][0],float(addOx[l][1]),float(addOx[l][2]),float(addOx[l][3])), file = f)
for m in range(len(addHy)):
 print("%3s%15.8f%15.8f%15.8f" % (addHy[m][0],float(addHy[m][1]),float(addHy[m][2]),float(addHy[m][3])), file = f)
f.close()


