import copy
from itertools import combinations
import functools
import sys


# assert len(sys.argv)==2, "Call script: python3 lattice.py number_best"

# number_best = int(sys.argv[1])

# print("Compute for the first best:" , number_best)


class UniversePower(object):
    '''Create a Universe of Powerset composed of all subsets of a set
        Returns a lattice instance.
    '''
    def __init__(self, baseElementsIterable: list):

        self.uElements = self.powerset_combinations(baseElementsIterable)

        assert all(self.uElements.count(elem) == 1 for elem in self.uElements), "uElements have duplications!"

        # self.idxToElements = {idx: uElem for idx, uElem in enumerate(self.uElements)}


    @property
    def size(self):
        return len(self.uElements)

    @property
    def indexes(self):
        return list(range(len(self.uElements)))

    @property
    def upperElement(self):
        return self.uElements[-1]

    @property
    def lowerElement(self):
        return self.uElements[0]

    @staticmethod
    def diff_size(elem1, elem2):
        elem = elem1.difference(elem2)
        return len(elem)

    def joinElements(self, *argvElements):
        assert all(el in self.uElements for el in argvElements )
        # compute lower upper bound
        lub = set.union(*argvElements)
        return lub

    def joinIndexes(self, *argvIndexes):
        assert all(0<= idx < len(self.uElements) for idx in argvIndexes )
        argvElements = [self.uElements[idx] for idx in argvIndexes]
        lub = self.joinElements(*argvElements)
        assert lub in self.uElements
        lub_index = self.uElements.index(lub)
        return lub_index


    def meetElements(self, *argvElements):
        assert all(el in self.uElements for el in argvElements )
        # compute greatest lower bound
        glb = set.intersection(*argvElements)
        return glb

    def meetIndexes(self, *argvIndexes):
        assert all(0<= idx < len(self.uElements) for idx in argvIndexes )
        argvElements = [self.uElements[idx] for idx in argvIndexes]
        # greatest lower bound
        glb = self.meetElements(*argvElements)
        assert glb in self.uElements
        glb_index = self.uElements.index(glb)
        return glb_index


    @staticmethod
    def powerset_combinations(iterable):
        "powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"
        s = list(iterable)
        lst = [ list( map(lambda x: set(x), list(combinations(s, r)) ) )  for r in range(len(s)+1) ]
        result = functools.reduce(lambda a, b: a+b, lst)
        return result

        # return functools.reduce( lambda a, b: a+b, lst)


    def getIndex(self,uElement):
        assert uElement in self.uElements, uElement
        return self.uElements.index(uElement)

    def getIndexes(self, uElements):
        assert all( uElem in self.uElements for uElem in uElements)
        return sorted([self.uElements.index(uElem) for uElem in uElements], reverse=False)


    def getElement(self, index):
        assert -len(self.uElements) <= index < len(self.uElements)
        return self.uElements[index]

    def getElements(self, indexes):
        assert all( -len(self.uElements) <= index < len(self.uElements) for index in indexes)
        return [self.uElements[index] for index in sorted(indexes, reverse=False)]


class PoSet(object):
    def __init__(self,  Universe, 
                        uElementsIndexes,
                        idxToWeights={ }) :
        '''Partially Ordered Set -> Poset 
        (Right/Left lattice if for all subset the LUB/GLB
        exists and also belongs to the Poset)

        Keyword arguments:
        Universe   -- the space where the uElements live
        uElementsIndexes  -- list. The indexes corresponding to the Poset from the Universe.
        idxToWeightsMap   -- list. A dict mapping indexes to weights associated with uElements 
        
        Returns a Poset instance.
        '''
        if len(idxToWeights) > 0:
            assert len(uElementsIndexes)==len(idxToWeights) and sorted(list(idxToWeights.keys())) == sorted(uElementsIndexes), "Elements indexes and idxToWeight should match"
        
        self.idxToWeights = copy.deepcopy(idxToWeights)

        self.Universe = Universe
        assert all(0<= x < len(Universe.uElements) for x in uElementsIndexes)

        self.Indexes = sorted(uElementsIndexes)
        self.addedIndexes = [ ]
        self.allIndexes = self.Indexes + self.addedIndexes

        # self.upperElement, self.lowerElement = None, None
        self.upperElement, self.upperElementIndex  = self.computeUpperElement()
        self.lowerElement, self.lowerElementIndex  = self.computeLowerElement()

    @property
    def idxToWeightsNormalized(self):
    # ''' To be used wisely.!! '''
        weights_ = max( list(self.idxToWeights.values()) )
        idxToWeightsNormalized = {idx: weight / weights_ for idx, weight in self.idxToWeights.items()}
        return idxToWeightsNormalized

    @property
    def uElementsToWeights(self):
        uElemToWeights = [(self.getElement(idx), weight) for idx, weight in self.idxToWeights.items()]
        return uElemToWeights

    @classmethod
    def fromElementsToWeights(cls, Universe, pairsElementsToWeights):
        # pairsElemsToWeights = sorted(pairsElemsToWeights, key=lambda p: p[])
        idxToWeights = { Universe.getIndex(el): weight for el, weight in pairsElementsToWeights }
        idxToWeights = dict(sorted(idxToWeights.items(), key=lambda p:p[0]))
        idxElements = sorted([idx for idx, weight in idxToWeights.items()])
        poset_object = cls(Universe, idxElements, idxToWeights)
        return poset_object
    # @property
    # def allIndexes(self):
    #     return sorted(self.Indexes + self.addedIndexes)

    def computeUpperElement(self):
        self.upperElementIndex     = self.Universe.joinIndexes( * self.allIndexes )
        self.upperElement = self.Universe.getElement(self.upperElementIndex)
        return self.upperElement, self.upperElementIndex

    def computeLowerElement(self):
        self.lowerElementIndex     = self.Universe.meetIndexes( *self.allIndexes )
        self.lowerElement = self.Universe.getElement(self.lowerElementIndex)
        return self.lowerElement, self.lowerElementIndex

    @property
    def isRightLattice(self):
        upperEl, upperInd = self.computeUpperElement() 
        return upperInd in self.allIndexes

    @property
    def isLeftLattice(self):
        lowerEl, lowerInd = self.computeLowerElement() 
        return lowerInd in (self.Indexes + self.addedIndexes)


    def getIndex(self,uElement):
        index = self.Universe.getIndex(uElement)
        assert index in (self.Indexes + self.addedIndexes)
        return index

    ''' !!! Next three functions can be expensive if used too often!!!'''
    @property
    def uElements(self):
        return self.Universe.getElements(self.Indexes)

    @property
    def addedElements(self):
        return self.Universe.getElements(self.addedIndexes)

    @property
    def allElements(self):
        return self.Universe.getElements(self.allIndexes)
    '''  !!!! To be used wisely !!!'''

    def getElement(self, index):
        assert index in self.allIndexes
        uElement = self.Universe.getElement(index)
        return uElement

    def addIndex(self, indexElement):
        assert 0 <= indexElement < len(self.Universe.uElements), str(indexElement) + " not in Universe indexes"
        assert indexElement not in self.Indexes, str(indexElement) + " already in " + str(self.Indexes)
        assert indexElement not in self.addedIndexes, str(indexElement) + " already added in " + str(self.addedIndexes)
        
        self.addedIndexes.append(indexElement)
        self.allIndexes = sorted(self.Indexes + self.addedIndexes)

    def addElement(self, uElement):
        assert uElement in self.Universe.uElements, "Element not in universe: " + str(uElement)
        indexElement = self.Universe.getIndex(uElement)
        assert indexElement not in self.Indexes, "Index already in Indexes: " + str(indexElement)
        
        if indexElement in self.addedIndexes:
            print( "Index already in addedIndexes: " + str(indexElement))
            return

        self.addedIndexes.append(indexElement)
        self.allIndexes = sorted(self.Indexes + self.addedIndexes)

        self.computeLowerElement()
        self.computeUpperElement()

    def addElements(self, uElements: list):
        for uElem in uElements:
            self.addElement(uElem)

    def mergeElementWeight(self, uElement, weight=None):
    # def addElement(self, uElement):
        assert uElement in self.Universe.uElements, "Element not in universe: " + str(uElement)
        indexElement = self.Universe.getIndex(uElement)
        # assert indexElement not in self.Indexes, "Index already in Indexes: " + str(indexElement)
        
        if indexElement in self.allIndexes:
            if weight is not None:
                self.idxToWeights[indexElement] = weight
            # print( "Index already in addedIndexes: " + str(indexElement))
            return
        else:
            if weight is not None:
                self.idxToWeights[indexElement] = weight

            self.addedIndexes.append(indexElement)
            self.allIndexes = sorted(self.Indexes + self.addedIndexes)
            # self.idxToWeights[indexElement] = weight

            self.computeLowerElement()
            self.computeUpperElement()


    def removeIndex(self, uElemIndex):
        '''Only added indexes can be removed'''
        assert uElemIndex in self.allIndexes, "Not in all Indexes: " + str(uElemIndex)
        if uElemIndex in self.Indexes:
            self.Indexes.remove(uElemIndex)
        if uElemIndex in self.addedIndexes:
            self.addedIndexes.remove(uElemIndex)
        self.allIndexes.remove(uElemIndex)

        if uElemIndex in self.idxToWeights.keys():
            del self.idxToWeights[uElemIndex]

        self.computeLowerElement()
        self.computeUpperElement()

    def removeElement(self, uElement):
        assert uElement in self.addedElements, "Not in added Elements: " + str(uElement) 
        elemIndex = self.getIndex(uElement)
        self.removeIndex(elemIndex)

# =============================================================

# #==========================================================


# univ01 = UniversePower([1,'ab','cd','ef','gh'])


# univ = UniversePower([1,2,3,4,5])

# indexes = univ.getIndexes([ set([2,3,4]) , set([2,5]), set() ])
# # elems = univ.uElements

# pairsElemsToWeights = [ ( set([2,3,4]), 22.22), (set([2,5]), 10.10 ), ( set(), 15.15) ]

# poset2 = PoSet.fromElementsToWeights(univ, pairsElemsToWeights)
# poset2.idxToWeightsNormalized



# idxToWeights = {5: 44.44, 8: 33.33, 25: 30.30}


# poset = PoSet(univ, list(idxToWeights.keys()), idxToWeights)

# poset.idxToWeightsNormalized


# result=combine((poset, 2), (poset2, 2))

# result.idxToWeightsNormalized
# result.idxToWeights

