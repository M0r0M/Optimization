import random
import Solution

class TabuSearch(object):
    """ object implementation the Multi-Objective Tabu Search """

    def __init__(self,datum,maxIter,stepSize,diversify=25,intensify=15,reduce=45,STMlen=20,debug=False):
        
        self.STM = list() # short term memory
        self.MTM = list() # medium term memory
        self.LTM = list() # long term memory
        self.IM = list() # intensification memory
        
        self.history = list() # store the solutions in order of visit (can store a same solution multiple times !)
        
        self.STMlen = STMlen # maximum size of the short term memory
        
        self.maxIter = maxIter # maximum iterations of the algorithm
        
        self.stepSize = stepSize # assert len(stepSize) == len(datum.DV)
        
        self.diversify = diversify
        self.intensify = intensify
        self.reduce = reduce
        
        self.iIter = 0
        self.iLocal = 0
        
        self.debug = debug
        
        self.optimize()
        
    def optimize(self):
        """ perform the Tabu Search as implemented in
        'The development of a multi-objective Tabu Search 
        algorithm for continuous optimisation problems'
        by
        D.M. Jaeggi et al. (2008)
        """
        
        self.history.append(datum)
        self.nextMoveHookeJeeves = True
        
        while self.iIter < self.maxIter: # stopping Criteria Not Met
            if self.nextMoveHookeJeeves:
                self.HookeJeevesMove()
            else:
                success = self.PatternMove()
            if success:
                self.nextMoveHookeJeeves = True
            else:
                self.HookeJeevesMove()
        self.UpdateMemories()
        if self.iLocal == self.diversify:
            if self.debug:
                print "diversifying !"
            self.DiversifyMove()
            self.nextMoveHookeJeeves = True
            self.UpdateMemories()
        elif self.iLocal == self.intensify:
            if self.debug:
                print "intensifying !"
            self.IntensifyMove()
            self.nextMoveHookeJeeves = True
            self.UpdateMemories()
        elif self.iLocal == self.reduce:
            if self.debug:
                print "reducing !"
            self.ReduceMove()
            self.nextMoveHookeJeeves = True
            self.UpdateMemories()      
    
    def removeDominatedPoints(self,container):
    	""" remove the dominated points 
    	from the container
    	container: list fo explorated points """
    	
        for pointA in container:
    		for pointB in container:
    			if pointA != pointB:
    				if pointA.dominate(pointB):
                        container.remove(pointB)

    def isTabu(self,point):
    	""" determine if a point is
    	a tabu point in the Short Term Memory
    	return:
    	True if point is tabu
    	False if point not tabu """

    	return (point in self.STM)

    def addIfNotDominated(self,newPoint,container):
    	""" add the newPoint in the container
    	if it is not dominated by any point
    	already in the container
    	return:
    	True if newPoint is added to container
    	False if newPoint not added to container """
        
    	for existingPoint in container:
    		if existingPoint.dominate(newPoint):
    			return False
    	container.append(newPoint)
    	return True

    def HookeJeevesMove(self):
    	
        currentPoint = self.history[self.iIter]
    	bestPoints = list()
        for i in range(len(currentPoint.DV)):
            newDV = list()
            for j in range(len(currentPoint.DV)):
                if i == j:
                    newDV.append(currentPoint.DV[j] + self.stepSize[j])
                else:
                    newDV.append(currentPoint.DV[j])
            newPoint = Solution(Solution.getName(),newDV)
            if (not self.isTabu(newPoint)) and (not newPoint.isValid):
                self.addIfNotDominated(newPoint,bestPoints)
                self.removeDominatedPoints(bestPoints)
            newDV = list()
            for j in range(len(currentPoint.DV)):
                if i == j:
                    newDV.append(currentPoint.DV[j] - self.stepSize[j])
                else:
                    newDV.append(currentPoint.DV[j])
            newPoint = Solution(Solution.getName(),newDV)
            if (not self.isTabu(newPoint)) and (not newPoint.isValid):
                self.addIfNotDominated(newPoint,bestPoints)
                self.removeDominatedPoints(bestPoints)
    	self.iIter += 1
    	self.iLocal += 1
    	nextPoint = self.selectRandom(bestPoints)
    	self.history.append(nextPoint)
    	bestPoints.remove(nextPoint)
    	for remainingPoints in bestPoints:
    		self.addIfNotDominated(remainingPoint,self.IM)
    	self.removeDominatedPoints(self.IM)
    	self.nextMoveHookeJeeves = False

    def PatternMove(self):
        
    	currentPoint = self.history[self.iIter]
    	lastPoint = self.history[self.iIter-1]
        # find out last move
        lastMove = list()
        for i in range(len(currentPoint.DV)):
            lastMove.append(currentPoint.DV[i] - lastPoint.DV[i])
        # recreate last move to current design vector
        nextDV = list()
        for i in range(len(currentPoint.DV)):
            nextDV.append(currentPoint.DV[i] + lastMove[i])
        
    	newPoint = Solution(Solution.getName(),nextDV)
    	newPoint.evaluate()
        
    	if newPoint.dominate(currentPoint):
    		self.iIter += 1
    		self.iLocal += 1
    		self.history.append(newPoint)
    		return True
    	return False
    
    def UpdateMemories(self):
        
    	currentPoint = self.history[self.iIter]
    	self.STM.append(currentPoint)
        while len(self.STM) > self.STMlen:
            self.STM.pop(0)
        success = self.addIfNotDominated(currentPoint,self.MTM)
        if success:
            self.removeDominatedPoints(self.MTM)
            self.iLocal = 0
        self.LTM.append(currentPoint)

    def DiversifyMove(self):
        # not functionnal
        regionCount = dict.fromkeys(listRegion)
        for point in self.LTM:
            for region in listRegion:
                if point in region:
                    regionCount[region] += 1
        newRegion = min(regionCount)
        newPoint = selectRandomPointFromRegion(newRegion)
        if self.isTabu(newPoint) or not newPoint.isValid:
            self.DiversifyMove()
        self.iIter += 1
        self.iLocal += 1
        self.history.append(newPoint)

    def IntensifyMove(self):
        
        newPoint = self.selectRandom(self.IM)
        if self.isTabu(newPoint):
            self.IntensifyMove()
        self.iIter += 1
        self.iLocal += 1
        self.history.append(newPoint)
        self.IM.remove(newPoint)

    def ReduceMove(self):
        
        newPoint = self.selectRandom(self.MTM)
        if self.isTabu(newPoint):
            self.ReduceMove()
        self.stepSizeReduction()
        self.iIter += 1
        self.iLocal = 0
        self.histroy.append(newPoint)
    
    def stepSizeReduction(self,coeff=2.0):
        """ reduce the step of each design vector of a factor of coeff """
        
        newStepSize = list()
        for DV in stepSize:
            newStepSize.append(DV/coeff)
            if self.debug:
                print "old step size: "+ DV +" new step size: "+ DV/coeff
        self.stepSize = newStepSize
    
    def selectRandom(container):
        """ randomly select a solution in the container """
        
        i = random.randint(0,len(container)-1)
        return container[i]
        
    def writePareto(self):
        
        content = str()
        file = open('result.txt','a')
        pareto = self.optimum()
        content += "iteration " + str(self.iIter) + "\n"
        content += "pareto front:\n"
        for element in self.MTM:
            content += element.description
        file.write(content)
        file.close()
            
if __name__ == "__main__":
    
    