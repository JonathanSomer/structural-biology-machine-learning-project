from enum import Enum
import os

class Rank(Enum):
    High = 1
    Medium = 2
    Acceptable = 3
    Incorrect = 4

class CapriScorer(object):
    
    # singleton implementation:
    _instance = None

    def __new__(_class):
    	if CapriScorer._instance == None:
    		CapriScorer._instance = object.__new__(CapriScorer)
    	return CapriScorer._instance

    def score():
    	pass
