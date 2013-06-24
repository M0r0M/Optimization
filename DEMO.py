import random, copy
import Solution

class DEMO(object):
    """ implementation of the Differential Evolution Multi-Objective """
    
    def __init__(self,InitialPopulation,iteration,crossover,differenceFactor,debug=False):
        
        self.P = InitialPopulation # list of the population
        self.lenInitPop = len(InitialPopulation) # lenght of the population
        self.iMax = iteration # number of iteration max
        self.i = 0 # number of iteration performed
        self.C = crossover # crossover parameter, in [0,1]
        self.F = differenceFactor # mate difference factor, in [0,1]
        
        self.debug = debug # debug flag
        
        self.optimize()
        
    def optimize(self):
        
        for member in self.P:
            member.evaluate()
        while self.i < len(self.P)*self.iMax:
            i = random.randint(0,len(self.P)-1)
            parent = self.P[i]
            child  = self.createChild(parent)
            child.evaluate()
            if child.dominate(parent):
                if self.debug:
                    print "child dominates !"
                    print "child:  " + child.name + "Objectives: " + str(child.objectives)
                    print "parent: " + self.P[i].name + "Objectives: " + str(child.objectives)
                self.replace(parent,child)
            elif parent.dominate(child):
                if self.debug:
                    print "parent dominates !"
                    print "child:  " + child.name + " L/D: " + str(child.liftOverDragMean) + " sigma: " + str(child.liftOverDragStdDev)
                    print "parent: " + self.P[i].name + " L/D: " + str(self.P[i].liftOverDragMean) + " sigma: " + str(self.P[i].liftOverDragStdDev)
                pass
            elif parent.objectives == child.objectives:
                if self.debug:
                    print "parent and child have the same characterisitcs !"
                    print "child:  " + child.name + " L/D: " + str(child.liftOverDragMean) + " sigma: " + str(child.liftOverDragStdDev)
                    print "parent: " + self.P[i].name + " L/D: " + str(self.P[i].liftOverDragMean) + " sigma: " + str(self.P[i].liftOverDragStdDev)
                pass
            else:
                if self.debug:
                    print "no domination !"
                    print "child:  " + child.name + " L/D: " + str(child.liftOverDragMean) + " sigma: " + str(child.liftOverDragStdDev)
                    print "parent: " + self.P[i].name + " L/D: " + str(self.P[i].liftOverDragMean) + " sigma: " + str(self.P[i].liftOverDragStdDev)
                self.P.append(child)
            self.i += 1
            if self.i % self.iMax == 0:
                while len(self.P) > self.lenInitPop:
                    if self.debug:
                        print "truncate !"
                    self.truncateCrowding()
            if self.debug:
                print self.i
            self.writePareto()
            if self.debug: 
                print "list status: "
                for members in self.P:
                    members.printChar()
    
    def createChild(self,parent):
        """ create a candidate child """
        
        i = self.P.index(parent)
        a = i
        b = i
        c = i
        while i == a or i == b or i == c or a == b or a == c or b == c:
            a = random.randint(0,len(self.P)-1)
            b = random.randint(0,len(self.P)-1)
            c = random.randint(0,len(self.P)-1)
        mateDV = self.mate(self.P[a],self.P[b],self.P[c])
        if self.debug:
            for DV in mateDV:
                print str(DV)
        childDV = self.crossover(self.P[i],DV)
        return Solution(Solution.getName(),childDV)        
        
    def mate(self,Pa,Pb,Pc):
        """ create a mate from the members Pa, Pb, and Pc
        output:
        MateDV: the new list of design vectors """
        
        PbDV = copy.deepcopy(Pb.DV)
        PcDV = copy.deepcopy(Pc.DV)
        MateDV = copy.deepcopy(Pa.DV)
        for i in range(len(MateDV)):
            MateDV[i] += self.F * (PbDV[i].x-PcDV[i].x)   
        return MateCtrlPoints
        
    def crossover(self,Pi,MateDV):
        """ perform a crossover between the member Pi and 
        the list of control points of design vector MateDV
        output:
        ChildDV: the new list of design vector """
        
        ChildDV = list()
        PiDV = copy.deepcopy(Pi.DV)
        for i in range(len(PiDV)):
            if random.random() > self.C:
                ChildDV.append(MateDV[i])
            else:
                ChildDV.append(PiDV[i])
        return ChildDV
    
    def replace(self,Pi,Pj):
        """ replace the member Pi by 
        the member Pj in the list of the population"""
        
        newPop = list()
        for member in self.P:
            if member != Pi:
                newPop.append(member)
        newPop.append(Pj)
        self.P = list(newPop)
    
    def sort(listSolution,criteria):
        """ sorts a list of Solutions in respect
        to a certain criteria, returns the list unchanged
        if the criteria not an objective function """
        
        sortedList = copy.deepcopy(listSolution)
        if criteria == 'distance':
            criteria = len(listSolution[0].objectives)
        swap = True
        while swap:
            swap = False
            for j in range(len(sortedList)-1):
                if criteria < len(listSolution[0].objectives):
                    if sortedList[j].objectives[criteria] < sortedList[j+1].objectives[criteria]:
                        temp = sortedList[j]
                        sortedList[j] = sortedList[j+1]
                        sortedList[j+1] = temp
                        swap = True
                elif criteria == len(listSolution[0].objectives):
                    if sortedList[j].distance < sortedList[j+1].distance:
                        temp = sortedList[j]
                        sortedList[j] = sortedList[j+1]
                        sortedList[j+1] = temp
                        swap = True
        return sortedList
      
    def crowdingDistance(self,listSolution):
        """ crowding distance of a list of
        FFD_Solutions as defined in the NSGA-II paper,
        returns the list ordered in decreasing distances """
        
        for sol in listSolution:
            sol.crowdingDist = 0
        sortedList = copy.deepcopy(listSolution)
        for criteria in range(len(listSolution[0].objectives)):
            sortedList = sort(sortedList,criteria)
            sortedList[0].distance = 100
            sortedList[len(sortedList)-1].distance = 100
            for i in range(1,len(listFFD_Sol)-1):
                sortedList[i].distance += (sortedList[i-1].objectives[criteria]-sortedList[i+1].objectives[criteria])/(sortedList[0].objectives[criteria]-sortedList[len(listFFD_Sol)-1].objectives[criteria])
        return sort(sortedList,"distance")
        
    def truncateDistance(self,listToTruncate=[]):
        
        listToTruncate = self.crowdingDistance(listToTruncate)
        minDistance = listToTruncate[0].distance
        MinDistanceElement = listToTruncate[0]
        for element in listToTruncate:
            if element.distance < MinDistanceElement.distance:
                MinDistanceElement = element
                minDistance = element.distance
        listToTruncate.remove(MinDistanceElement)
        return listToTruncate
        
    def truncateCrowding(self):
        """ truncation of the population 
        maximizing the sum of the distances between members """
        
        nonDominatedFronts = self.paretoRanking()
        truncatedPop = list()
        i = 0
        while (self.lenInitPop > len(truncatedPop)) and (i < len(nonDominatedFronts)):
            while (len(truncatedPop) + len(nonDominatedFronts[i])) > self.lenInitPop:
                nonDominatedFronts[i] = self.truncateDistance(nonDominatedFronts[i]) 
            truncatedPop += nonDominatedFronts[i]
            i += 1
        self.P = copy.deepcopy(truncatedPop)
        
    def optimum(self,listPoints=[]):
        """ retrun the current pareto front """
        
        if listPoints == []:
            listPoints = self.P
        paretoFront = [listPoints[0]]
        for member in listPoints:
            dominatePareto = True
            for paretoPoint in paretoFront:
                if member.dominate(paretoPoint):
                    paretoFront.remove(paretoPoint)
                else:
                    dominatePareto = False
            if dominatePareto:
                paretoFront = [member]
            else:
                dominated = False
                for paretoPoint in paretoFront:
                    if paretoPoint.dominate(member):
                        dominated = True
                if not dominated and (member not in paretoFront):
                    paretoFront.append(member)
        return paretoFront
    
    def paretoRanking(self):
        """ return a list of lists of
        non-dominated sets """
        rankingList = list()
        currentState = copy.deepcopy(self.P)
        while currentState != []:
            nextFront = self.optimum(currentState)
            rankingList.append(nextFront)
            for point in nextFront:
                currentState.remove(point)
        return rankingList

    def writePareto(self):
        
        content = str()
        file = open('result.txt','a')
        pareto = self.optimum()
        content += "iteration " + str(self.i) + "\n"
        content += "pareto front:\n"
        for element in pareto:
            content += element.description
        file.write(content)
        file.close()

if __name__ == "__main__":
    