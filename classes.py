import random
import math
import numpy as np


class Tile:
  #Because we are not sure of these values, we should choose a few to vary and see how that affects the spread in the writeup
    whorlTuberRatio = 500  # how many whorls are produced per tuber
    tileCarryingCapacity = 1000  # carrying capacity assuming only active tubers
    growthRate = 0.011 #0.008 for (not quite) doubling every year (best case), 0.011 for worst case
    energyRatio = 0.05  # amount of energy (nutrients) required by inactive tubers compared to active tubers
    diffusionRate = 0.1  # chance of spread for new inactive tubers

    def __init__(self, active, inactive, location):
        self.active = active
        self.inactive = inactive
        self.loc = location
        self.coverDensity = 0
        self.surfaceTreatment = 0  # Between 0 and 1
        self.fullTreatment = 0  # Between and 1

    def tuberGrowth(self, tiles, timeOfYear):
      if self.active > 0:
        logisticFactor = (1 - (
                    self.active + Tile.energyRatio * self.inactive) / Tile.tileCarryingCapacity)
        newInactive = Tile.growthRate * logisticFactor * self.active * timeOfYear * (180 - timeOfYear) / 90 **2
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
            self.inactive += int((1 - Tile.diffusionRate) * newInactive)
            if random.random() < (1 - Tile.diffusionRate) * newInactive % 1:
              self.inactive += 1
            rightSpread = random.random()  # fraction of diffused new inactive tubers that go to right tile
            tileRight = tiles[(self.loc + np.random.binomial(3, 0.1) + 1) % len(tiles)]
            tileRight.inactive += int(rightSpread * Tile.diffusionRate * newInactive)
            if random.random() < (rightSpread * Tile.diffusionRate * newInactive) % 1:
              tileRight.inactive += 1
            tileLeft = tiles[(self.loc - np.random.binomial(3, 0.1) - 1) % len(tiles)]
            tileLeft.inactive += int((1 - rightSpread) * Tile.diffusionRate * newInactive)
            if random.random() < ((1 - rightSpread) * Tile.diffusionRate * newInactive) % 1:
              tileLeft.inactive += 1
        if self.fullTreatment > 0:
            effectiveness = 0.2 
            # How much should we remove per day based on how much treatment we use
            spread = random.random()
            self.active -= effectiveness*self.fullTreatment*spread
            self.inactive -= effectiveness*self.fullTreatment*(1-spread)

    # Not necessary given time constraints
    # def distance(self, tile, numTiles):
    #   return abs(abs(tile.loc - numTiles / 2) - abs(self.loc - numTiles / 2))

    def updateCoverDensity(self, timeOfYear):
        timeDensity = Tile.whorlTuberRatio * timeOfYear * (180 - timeOfYear) / pow(90, 2) * self.active
        self.coverDensity = timeDensity * (1 - self.surfaceTreatment) * (1 - self.fullTreatment) / 1000
        # ^ need some constants for the different treatments
        #the 1000 is an estimation of the probability that any individual whorl is picked up
        # timeOfYear is between 0 and 180
        # We are only simulating 180 days out of the year
        # capacity = someConstant * self.active * timeOfYear * (180 - timeOfYear)
        # self.coverDensity += anotherConstant * (1 - self.coverDensity / capacity) * (aThirdConstant * self.active + self.coverDensity)

    def updateTreatment(self, something, something_else):
        self.surfaceTreatment += something
        self.fullTreatment += something_else

    def whorlsCaught(self, timeOfYear):
        self.updateCoverDensity(timeOfYear)
        n = int(Tile.whorlTuberRatio * self.active)
        p = self.coverDensity / (Tile.whorlTuberRatio * Tile.tileCarryingCapacity) * random.random()
        assert p <= 1 #allow to recover from error
        return np.random.binomial(n, p)
3

class Marina(Tile):
    def __init__(self, active, inactive, location, popularity):
        Tile.__init__(self, active, inactive, location)
        self.popularity = popularity
        self.dockBoats = []

class PrivateDock(Tile):
    def __init__(self, active, inactive, location):
      Tile.__init__(self, active, inactive, location)
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

