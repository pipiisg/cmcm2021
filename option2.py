import random
import math
import numpy as np
import display


class VariableTile:
  #Because we are not sure of these values, we should choose a few to vary and see how that affects the spread in the writeup
    whorlTuberRatio = 500  # how many whorls are produced per tuber
    tileCarryingCapacity = 1000  # carrying capacity assuming only active tubers
    growthRate = 0.008 #0.008 for (not quite) doubling every year (best case), 0.011 for worst case
    energyRatio = 0.05  # amount of energy (nutrients) required by inactive tubers compared to active tubers
    diffusionRate = 0.1  # chance of spread for new inactive tubers

    def __init__(self, active, inactive, location, growthRate, diffusionRate, treatEff):
        self.active = active
        self.inactive = inactive
        self.loc = location
        self.coverDensity = 0
        self.surfaceTreatment = 0  # Between 0 and 1
        self.fullTreatment = 0  # Between and 1
        self.dockBoats = []
        self.growthRate = growthRate
        self.diffusionRate = diffusionRate
        self.treatEff = treatEff

    def tuberGrowth(self, tiles, timeOfYear):
      if self.active > 0:
        logisticFactor = (1 - (
                    self.active + VariableTile.energyRatio * self.inactive) / VariableTile.tileCarryingCapacity)
        newInactive = self.growthRate * logisticFactor * self.active * timeOfYear * (180 - timeOfYear) / 90 **2
        # if (self.loc == 0):
        #   print("Logistic Factor: ", round(logisticFactor,3))
        #   print(round(self.active, 3), round(self.inactive, 3), round(newInactive, 3))
        if newInactive < 0:
            frac = random.random()
            self.inactive += int(newInactive * frac)
            self.active += int(newInactive * (1 - frac))
            if self.active < 0:
              self.active = 0
            if self.inactive < 0:
              self.inactive = 0

        else:
            self.inactive += int((1 - self.diffusionRate) * newInactive)
            if random.random() < (1 - self.diffusionRate) * newInactive % 1:
              self.inactive += 1
            rightSpread = random.random()  # fraction of diffused new inactive tubers that go to right tile
            tileRight = tiles[(self.loc + np.random.binomial(3, 0.1) + 1) % len(tiles)]
            tileRight.inactive += int(rightSpread * self.diffusionRate * newInactive)
            if random.random() < (rightSpread * self.diffusionRate * newInactive) % 1:
              tileRight.inactive += 1
            tileLeft = tiles[(self.loc - np.random.binomial(3, 0.1) - 1) % len(tiles)]
            tileLeft.inactive += int((1 - rightSpread) * self.diffusionRate * newInactive)
            if random.random() < ((1 - rightSpread) * self.diffusionRate * newInactive) % 1:
              tileLeft.inactive += 1
        if self.fullTreatment > 0:
            # How much should we remove per day based on how much treatment we use
            spread = random.random()
            self.active -= int(self.treatEff*self.fullTreatment *spread*self.active)
            self.inactive -= int(self.inactive * self.treatEff*self.fullTreatment*(1-spread))

    # Not necessary given time constraints
    # def distance(self, tile, numTiles):
    #   return abs(abs(tile.loc - numTiles / 2) - abs(self.loc - numTiles / 2))

    def updateCoverDensity(self, timeOfYear):
        timeDensity = VariableTile.whorlTuberRatio * timeOfYear * (180 - timeOfYear) / pow(90, 2) * self.active
        self.coverDensity = timeDensity * (1 - self.surfaceTreatment) * (1 - self.fullTreatment) / 1000
        # ^ need some constants for the different treatments
        #the 1000 is an estimation of the probability that any individual whorl is picked up
        # timeOfYear is between 0 and 180
        # We are only simulating 180 days out of the year
        # capacity = someConstant * self.active * timeOfYear * (180 - timeOfYear)
        # self.coverDensity += anotherConstant * (1 - self.coverDensity / capacity) * (aThirdConstant * self.active + self.coverDensity)

    def updateTreatment(self, surface, full):
        self.surfaceTreatment += surface
        self.fullTreatment += full

    def whorlsCaught(self, timeOfYear):
        self.updateCoverDensity(timeOfYear)
        n = int(VariableTile.whorlTuberRatio * self.active)
        p = self.coverDensity / (VariableTile.whorlTuberRatio * VariableTile.tileCarryingCapacity) * random.random()
        assert p <= 1 #allow to recover from error
        return np.random.binomial(n, p)


class Marina(VariableTile):
    def __init__(self, active, inactive, location, popularity, growthRate, diffusionRate, treatEff):
        VariableTile.__init__(self, active, inactive, location, growthRate, diffusionRate, treatEff)
        self.popularity = popularity
        self.dockBoats = []
        self.dockBoats = []

class Boat:
    def __init__(self, dock, freqUse, numWhorls=0):
        self.dock = dock
        self.freqUse = freqUse
        self.numWhorls = numWhorls

    def chooseMarina(self, marinas):
        marinaProbs = [marina.popularity for marina in marinas]
        return random.choices(marinas, marinaProbs)

    def clean(self, survivalRate):
        self.numWhorls = int(self.numWhorls * survivalRate)

class Simulation:

    def __init__(self, numTiles, marinaInfos, boatInfos, growthRate, diffusionRate, fractionLeftBehind, treatEff, isReality = False):
      self.fractionLeftBehind = 0
      self.tiles = []
      self.marinas = []
      self.boats = []
      if not isReality:
        self.fractionLeftBehind = fractionLeftBehind
        self.tiles = [VariableTile(10, 0, a, growthRate * getRandomFactor(), diffusionRate* getRandomFactor(), treatEff* getRandomFactor()) if 690 < a < 700 else VariableTile(0, 0, a, growthRate* getRandomFactor(), diffusionRate* getRandomFactor(), treatEff* getRandomFactor()) for a in range(numTiles)]
        #There is some hydrilla around cayuga inlet
        self.marinas = [Marina(0, 0, a[0], a[1], growthRate* getRandomFactor(), diffusionRate* getRandomFactor(), treatEff* getRandomFactor()) for a in marinaInfos]
        for marina in self.marinas:
            self.tiles[marina.loc] = marina
        self.boats = []
        for a in boatInfos:
          newBoat = Boat(self.tiles[a[0]], a[1])
          self.boats.append(newBoat)
          self.tiles[a[0]].dockBoats.append(newBoat)
      else:
        self.fractionLeftBehind = fractionLeftBehind
        self.tiles = [VariableTile(10, 0, a, growthRate, diffusionRate, treatEff) if 690 < a < 700 else VariableTile(0, 0, a, growthRate, diffusionRate, treatEff) for a in range(numTiles)]
        #There is some hydrilla around cayuga inlet
        self.marinas = [Marina(0, 0, a[0], a[1], growthRate, diffusionRate, treatEff) for a in marinaInfos]
        for marina in self.marinas:
            self.tiles[marina.loc] = marina
        self.boats = []
        for a in boatInfos:
          newBoat = Boat(self.tiles[a[0]], a[1])
          self.boats.append(newBoat)
          self.tiles[a[0]].dockBoats.append(newBoat)

    def smallTimeStep(self, timeOfYear):
        # for time of year: number of days since May 1 (start of hydrilla season)
        # end of season is at timeOfYear = 180
        for boat in self.boats:
            if random.random() < boat.freqUse:
                marina = boat.chooseMarina(self.marinas)[0]
                boat.numWhorls += boat.dock.whorlsCaught(timeOfYear)
                whorlsLeftBehind = int(np.random.binomial(boat.numWhorls, 1 / self.fractionLeftBehind))
                marina.inactive += 0.1 * whorlsLeftBehind #survival rate of implanted whorls
                boat.numWhorls -= whorlsLeftBehind
                boat.numWhorls += marina.whorlsCaught(timeOfYear)
                whorlsLeftBehind = np.random.binomial(boat.numWhorls, 1 / self.fractionLeftBehind)
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

    def displayResults(self):
      display.plotLake(self.tiles)

def minUncertainty(sim1, sim2, reality, numCheckBoats):
  errors = [[abs(sim1.tiles[i].active - sim2.tiles[i].active), i] for i in range(len(sim1.tiles))]
  errors.sort(key=lambda a:-a[0])
  for tileLoc in errors[:40 - numCheckBoats]:
    sim1.tiles[tileLoc[1]].active = reality.tiles[tileLoc[1]].active
    sim2.tiles[tileLoc[1]].active = reality.tiles[tileLoc[1]].active
  for a in range(5 * numCheckBoats): #5??
    boat = random.choice(reality.boats)
    if boat.numWhorls > 0:
      tileLoc = boat.dock.loc
      sim1.tiles[tileLoc].active = max(1, sim1.tiles[tileLoc].active)
      sim2.tiles[tileLoc].active = max(1, sim2.tiles[tileLoc].active)
def minUncertainty2(sim1, sim2, reality, numCheckBoats):
  errors = [[abs(np.log(sim1.tiles[i].active + 1) - np.log(sim2.tiles[i].active + 1)), i] for i in range(len(sim1.tiles))]
  errors.sort(key=lambda a:-a[0])
  for tileLoc in errors[:40 - numCheckBoats]:
    sim1.tiles[tileLoc[1]].active = reality.tiles[tileLoc[1]].active
    sim2.tiles[tileLoc[1]].active = reality.tiles[tileLoc[1]].active
  for a in range(5 * numCheckBoats): #5??
    boat = random.choice(reality.boats)
    if boat.numWhorls > 0:
      tileLoc = boat.dock.loc
      sim1.tiles[tileLoc].active = max(1, sim1.tiles[tileLoc].active)
      sim2.tiles[tileLoc].active = max(1, sim2.tiles[tileLoc].active)
def maxUncertainty2(sim1, sim2, reality, numCheckBoats):
  errors = [[abs(np.log(sim1.tiles[i].active + 1) - np.log(sim2.tiles[i].active + 1)), i] for i in range(len(sim1.tiles))]
  errors.sort(key=lambda a:a[0])
  for tileLoc in errors[:40 - numCheckBoats]:
    sim1.tiles[tileLoc[1]].active = reality.tiles[tileLoc[1]].active
    sim2.tiles[tileLoc[1]].active = reality.tiles[tileLoc[1]].active
  for a in range(5 * numCheckBoats): #5??
    boat = random.choice(reality.boats)
    if boat.numWhorls > 0:
      tileLoc = boat.dock.loc
      sim1.tiles[tileLoc].active = max(1, sim1.tiles[tileLoc].active)
      sim2.tiles[tileLoc].active = max(1, sim2.tiles[tileLoc].active)

def seeMaxHydrilla(sim1, sim2, reality, numCheckBoats):
    errors = [[sim1.tiles[i].active + sim2.tiles[i].active, i] for i in range(len(sim1.tiles))]
    errors.sort(key=lambda a:-a[0])
    for tileLoc in errors[:40 - numCheckBoats]:
        sim1.tiles[tileLoc[1]].active = reality.tiles[tileLoc[1]].active
        sim2.tiles[tileLoc[1]].active = reality.tiles[tileLoc[1]].active
    for a in range(5 * numCheckBoats): #5??
        boat = random.choice(reality.boats)
        if boat.numWhorls > 0:
            tileLoc = boat.dock.loc
            sim1.tiles[tileLoc].active = max(1, sim1.tiles[tileLoc].active)
            sim2.tiles[tileLoc].active = max(1, sim2.tiles[tileLoc].active)
def randomChoices(sim1, sim2, reality, numCheckBoats):
    errors = [[np.random.rand(), i] for i in range(len(sim1.tiles))]
    errors.sort(key=lambda a:-a[0])
    for tileLoc in errors[:40 - numCheckBoats]:
        sim1.tiles[tileLoc[1]].active = reality.tiles[tileLoc[1]].active
        sim2.tiles[tileLoc[1]].active = reality.tiles[tileLoc[1]].active
    for a in range(5 * numCheckBoats): #5??
        boat = random.choice(reality.boats)
        if boat.numWhorls > 0:
            tileLoc = boat.dock.loc
            sim1.tiles[tileLoc].active = max(1, sim1.tiles[tileLoc].active)
            sim2.tiles[tileLoc].active = max(1, sim2.tiles[tileLoc].active)
def maxUncertainty(sim1, sim2, reality, numCheckBoats):
    errors = [[np.random.rand(), i] for i in range(len(sim1.tiles))]
    errors.sort(key=lambda a:a[0])
    for tileLoc in errors[:40 - numCheckBoats]:
        sim1.tiles[tileLoc[1]].active = reality.tiles[tileLoc[1]].active
        sim2.tiles[tileLoc[1]].active = reality.tiles[tileLoc[1]].active
    for a in range(5 * numCheckBoats): #5??
        boat = random.choice(reality.boats)
        if boat.numWhorls > 0:
            tileLoc = boat.dock.loc
            sim1.tiles[tileLoc].active = max(1, sim1.tiles[tileLoc].active)
            sim2.tiles[tileLoc].active = max(1, sim2.tiles[tileLoc].active)


def getTotalError(sim1, sim2, reality):
    sims = np.array([(sim1.tiles[i].inactive + sim2.tiles[i].inactive) for i in range(len(sim1.tiles))])
    realities = np.array([(reality.tiles[i].active) for i in range(len(sim1.tiles))])
    p = np.polyfit(sims, realities, 1)
    print(p)
    sumRes = 0
    for a in range(len(sims)):
        diff = (np.log(sims[a] * p[0] + p[1] + 1) - np.log(realities[a] + 1))**2
        sumRes += diff
    return sumRes / len(sims) / np.var(np.log(realities + 1))

def minUncertainty3(sim1, sim2, reality, numCheckBoats):
  errors = [[abs(sim1.tiles[i].active + sim2.tiles[i].active - 2 * reality.tiles[i].active), i] for i in range(len(sim1.tiles))]
  errors.sort(key=lambda a:-a[0])
  for tileLoc in errors[:40 - numCheckBoats]:
    sim1.tiles[tileLoc[1]].active = reality.tiles[tileLoc[1]].active
    sim2.tiles[tileLoc[1]].active = reality.tiles[tileLoc[1]].active
  for a in range(5 * numCheckBoats): #5??
    boat = random.choice(reality.boats)
    if boat.numWhorls > 0:
      tileLoc = boat.dock.loc
      sim1.tiles[tileLoc].active = max(1, sim1.tiles[tileLoc].active)
      sim2.tiles[tileLoc].active = max(1, sim2.tiles[tileLoc].active)

def getRandomFactor():
  return 1 + (np.random.normal() / 10)

class BadTile:
  carryingCapacity = 500
  def __init__(self, active, inactive, location, growthRate, diffusionRate, treatEff, sproutConst):
    self.active = active
    self.inactive = inactive
    self.loc = location
    self.growthRate =  growthRate
    self.diffusionRate = diffusionRate
    self.treatEff = treatEff
    self.sproutConst = sproutConst
  
  def grow(self, timeOfYear, treatment, tiles):
    sproutProp = self.sproutConst * timeOfYear * (90 - timeOfYear) / 45**2
    if sproutProp > 0:
      self.active += sproutProp * self.inactive
      self.inactive *= 1 - sproutProp
    currentCarryingCapacity = BadTile.carryingCapacity * np.exp(-self.treatEff * treatment) * timeOfYear * (180 - timeOfYear) / 90**2 + 1
    logisticFactor = self.growthRate * (1 - self.active / currentCarryingCapacity)
    if logisticFactor > 0:
      self.inactive += self.active * logisticFactor * (1 - self.diffusionRate)
      tiles[(self.loc + 1 + np.random.binomial(3, 0.1)) % len(tiles)].inactive += self.active * logisticFactor * self.diffusionRate / 2
      tiles[(self.loc - 1 - np.random.binomial(3, 0.1)) % len(tiles)].inactive += self.active * logisticFactor * self.diffusionRate / 2
    else:
      self.active += self.active * logisticFactor
      if self.active < 0:
        self.active = 0

class BadSimulation:
  def __init__(self, numTiles, growthRate, diffusionRate, treatEff, sproutConst):
    self.tiles = [BadTile(0, 10, a, growthRate, diffusionRate, treatEff, sproutConst) if 690 < a < 700 else BadTile(0, 0, a, growthRate, diffusionRate, treatEff, sproutConst) for a in range(numTiles)]
  def smallTimeStep(self, timeOfYear):
    for tile in self.tiles:
      tile.grow(timeOfYear, 0, self.tiles)
  def largeTimeStep(self):
    for tile in self.tiles:
      tile.inactive += tile.active
      tile.active = 0
  def displayResults(self):
      display.plotLake2(self.tiles)
bestSim = 0
worstSim = 0
reality = 0
def setup():
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
    #boats launched from docks
    for a in range(1000):
        bInfo.append([int(1600 * random.random()), 1 / (1 + np.random.binomial(30, 0.25))])
    #boats launched from marinas:
    for a in range(1000):
        marinaProbs = [a[1] for a in mInfo]
        marinaChoices = [a[0] for a in mInfo]
        bInfo.append([random.choices(marinaChoices, marinaProbs)[0], 1 / (1 + np.random.binomial(30, 0.25))])
    global bestSim
    global worstSim
    global reality
    bestSim = BadSimulation(1600, 0.001, 0.1, 0.2, 0.01)
    worstSim = BadSimulation(1600, 0.01, 0.1, 0.2, 0.05)
    reality = Simulation(1600, mInfo, bInfo, 0.0095, 0.1, 400, 0.35, True)

minErrors = []
seeErrors = []
maxErrors = []
boatErrors = []
randomErrors = []
for a in range(5):
    setup()
    for year in range(0, 5):
        for day in range(0, 180):
            worstSim.smallTimeStep(day)
            bestSim.smallTimeStep(day)
            reality.smallTimeStep(day)
            if day % 7 == 0:
                minUncertainty3(bestSim, worstSim, reality,5)
        reality.largeTimeStep()
        worstSim.largeTimeStep()
        bestSim.largeTimeStep()
        print(f"year {year} done")
    maxErrors.append(getTotalError(bestSim, worstSim, reality))
print(maxErrors)
