from classes import Tile, Marina, Boat
import random
import numpy as np
import display

class Simulation:

    def __init__(self, numTiles, marinaInfos, boatInfos):
        self.tiles = [Tile(100, 0, a) if 690 < a < 700 else Tile(0, 0, a) for a in range(numTiles)]
        #There is some hydrilla around cayuga inlet
        self.marinas = [Marina(0, 0, a[0], a[1]) for a in marinaInfos]
        for marina in self.marinas:
            self.tiles[marina.loc] = marina
        self.boats = []
        for a in boatInfos:
          newBoat = Boat(self.tiles[a[0]], a[1])
          self.boats.append(newBoat)
          #self.tiles[a[0]].dockBoats.append(newBoat)
        self.docks = set() 
        for boat in self.boats:
          self.docks.add(boat.dock)
        self.obsList = []
        self.jobs = ['a', 'b', 'c', 'd']

    def smallTimeStep(self, timeOfYear):
        # for time of year: number of days since May 1 (start of hydrilla season)
        # end of season is at timeOfYear = 180
        for boat in self.boats:
            if random.random() < boat.freqUse:
                marina = boat.chooseMarina(self.marinas)[0]
                boat.numWhorls += boat.dock.whorlsCaught(timeOfYear)
                whorlsLeftBehind = int(np.random.binomial(boat.numWhorls, 1 / 300))
                marina.inactive += 0.1 * whorlsLeftBehind #survival rate of implanted whorls
                boat.numWhorls -= whorlsLeftBehind
                boat.numWhorls += marina.whorlsCaught(timeOfYear)
                whorlsLeftBehind = np.random.binomial(boat.numWhorls, 1 / 300)
                boat.dock.inactive += 0.1 * whorlsLeftBehind
                boat.numWhorls -= whorlsLeftBehind
                boat.clean(0.1)  # guessed survival rate on boat
        for tile in self.tiles:
            tile.tuberGrowth(self.tiles, timeOfYear)  # add parameters for rate of growth,
            #carring capacity, and weight of inactive compared to active

    def largeTimeStep(self):
        for boat in self.boats:
            boat.whorlsCaught = 0
        for tile in self.tiles:
            tile.active += tile.inactive
            tile.inactive = 0

    def displayResults(self, year):
      display.plotLake(self.tiles, year)

m1 = [1519, 1]
m2 = [1530, 1]
m3 = [1555, 1.5]
m4 = [23, 1.3]
m5 = [1414, 0.2]
m6 = [52, 3]
m7 = [79, 2]
m8 = [94, 0.4]
m9 = [1540, 0.6]
m10 = [254, 0.7]
m11 = [350, 0.4]
m12 = [852, 1]
m13 = [590, 2.5]
m14 = [655, 1]
m15 = [698, 2.5]

mInfo = [m1, m2, m3, m4, m5, m6, m7, m8, m9, m10, m11, m12, m13, m14, m15]

bInfo = []
for a in range(1000):
  bInfo.append([int(1600 * random.random()), 1 / (1 + np.random.binomial(30, 0.25))])
for a in range(1000):
  marinaProbs = [a[1] for a in mInfo]
  marinaChoices = [a[0] for a in mInfo]
  bInfo.append([random.choices(marinaChoices, marinaProbs)[0], 1 / (1 + np.random.binomial(30, 0.25))])

sim = Simulation(1600, mInfo, bInfo)
for year in range(0, 10):
    for day in range(0, 180):
        sim.smallTimeStep(day)
    sim.largeTimeStep()
    print(f"year {year} done")
sim.displayResults(year)
display.show()