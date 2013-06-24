""" The model class for a multi-objective optimization solution """

import time

def getName(style='time'):
    if style == 'time':
        return str(int(time.time()*1000000))

class Solution(object):
    """ Generic solution, nominatory methods are described here """
    
    def __init__(self,name,designVector,boundaries=None):
        
        self.name = name
        self.DV = designVector
        self.boundaries = boundaries # boundaries a list of tuples or list of 2 elements: a min and a max
        
        self.crowdingDist = 0 # distance between the objective functions in a set of solutions
        
        self.objectives = list() # objectives to maximize have to be negative
     
    def evaluate(self):
        """ implementation of the evaluation of the solution """
        # TO DO
    
    @property
    def description(self):
        """ the descriptor of the solution """
        print "id: " + self.name + " Design Vector: " + str(self.DV) + " Objectives: " + str(self.objectives) + "\n"
    
    @property
    def isValid(self):
        """ check if the design vector is in the boundaries """
        if self.boundaries == None:
            return True
        for i in range(len(self.DV)-1):
            if self.DV[i] < self.boundaries[i][0] or self.DV[i] > self.boundaries[i][1]:
                return False
        return True
        
    def dominate(self,Solution):
        """ check if the current solution dominates the solution Solution
        returns
        True if self dominates Solution
        False if not (either it is non-dominated or dominated) """
        
        dominates = True
        for i in range(len(self.objectives)):
            if self.objectives[i] > Solution.objectives[i]:
                dominates = False