import copy
from itertools import combinations
import functools
import sys


assert len(sys.argv)==2, "Call script: python3 lattice.py number_best"

number_best = int(sys.argv[1])

print("Compute for the first best:" , number_best)


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


def combine_4( pos1_nmb, pos2_nmb, pos3_nmb, pos4_nmb, df_gr):
    columns = list(df_gr.columns)
    columns.remove('id')
    columns.remove('idx_elem')
    columns.remove('auc')
    columns.remove('method')

    cols_features = columns
    columns = ['idx_elem','idx_elem1','idx_elem2','idx_elem3','idx_elem4', 'auc1', 'auc2', 'auc3', 'auc4', 
              'diff_size1','diff_size2','diff_size3','diff_size4','size_elem', 
              'weight1','weight2','weight3','weight4','weight_sum', 'weight_min'] + cols_features
    df_res = df = pd.DataFrame(columns=columns)

    print(df_res.columns)

    # df_res.drop('id', axis=1, inplace=True)
    # [(pos1, first_nmb1), (pos2, first_nmb2),...]
    pos1, nmb1 = pos1_nmb[0], pos1_nmb[1]
    sortedIdxToWeightsNormalized1 = sorted(pos1.idxToWeightsNormalized.items(), key = lambda x: x[1], reverse=True)[:nmb1]

    print(sortedIdxToWeightsNormalized1)

    pos2, nmb2 = pos2_nmb[0], pos2_nmb[1]
    sortedIdxToWeightsNormalized2 = sorted(pos2.idxToWeightsNormalized.items(), key = lambda x: x[1], reverse=True)[:nmb2]

    print(sortedIdxToWeightsNormalized2)

    pos3, nmb3 = pos3_nmb[0], pos3_nmb[1]
    sortedIdxToWeightsNormalized3 = sorted(pos3.idxToWeightsNormalized.items(), 
                                           key = lambda x: x[1], reverse=True)[:nmb3]

    print(sortedIdxToWeightsNormalized3)

    pos4, nmb4 = pos4_nmb[0], pos4_nmb[1]
    sortedIdxToWeightsNormalized4 = sorted(pos4.idxToWeightsNormalized.items(), 
                                           key = lambda x: x[1], reverse=True)[:nmb4]

    print(sortedIdxToWeightsNormalized4)



    universe = pos1.Universe
    idxToWeightPair   = { }
    idxToWeightResult = { }
    idx_df = 0



    for idx4, weight4 in sortedIdxToWeightsNormalized4:
        elem4 = universe.getElement(idx4)
        auc4 = pos4.idxToWeights[idx4]


        for idx1, weight1 in sortedIdxToWeightsNormalized1:
            elem1 = universe.getElement(idx1)
            auc1 = pos1.idxToWeights[idx1]

            for idx2, weight2 in sortedIdxToWeightsNormalized2:
                elem2 = universe.getElement(idx2)
                auc2 = pos2.idxToWeights[idx2]

                for idx3, weight3 in sortedIdxToWeightsNormalized3:
                    auc3 = pos3.idxToWeights[idx3]
                    elem3 = universe.getElement(idx3)
                    idx_glb = universe.joinIndexes(idx1, idx2, idx3, idx4)
                    elem_glb = universe.getElement(idx_glb)
                    weight_sum =  weight1 + weight2 + weight3 + weight4
                    idxToWeightResult[idx_glb] = max( idxToWeightResult.get(idx_glb, 0.0) , weight1 + weight2 + weight3+weight4)

                    diff_size1 = universe.diff_size(elem_glb, elem1)
                    diff_size2 = universe.diff_size(elem_glb, elem2)
                    diff_size3 = universe.diff_size(elem_glb, elem3)
                    diff_size4 = universe.diff_size(elem_glb, elem4)

                    # if idx_glb in idxToWeightResult.keys():
                    idxToWeightPair[idx_glb] = (weight1, weight2, weight3, weight4)
                    # idxToWeightResult[idx_glb] = max( idxToWeightResult.get(idx_glb, 0.0) , weight1 + weight2)
                    # print('-------------')
                    # print(elem_glb)
                    # print(idx_glb)
                    # print(auc1)
                    # print(auc2)
                    # print(diff_size1)
                    # print(diff_size2)
                    # print(idx1, weight1)
                    # print(idx2, weight2)   
                    # print(weight_sum)    

                    row_dict = {'idx_elem': idx_glb,
                                    'idx_elem1': idx1,
                                    'idx_elem2': idx2,
                                    'idx_elem3': idx3,
                                    'idx_elem4': idx4,
                                    'auc1': auc1, 
                                    'auc2': auc2, 
                                    'auc3': auc3,
                                    'auc4': auc4, 
                                    'diff_size1':diff_size1,
                                    'diff_size2':diff_size2, 
                                    'diff_size3':diff_size3,
                                    'diff_size4':diff_size4, 
                                    'size_elem': len(elem_glb),
                                    'weight1':weight1,
                                    'weight2':weight2,
                                    'weight3':weight3,
                                    'weight4':weight4,
                                    'weight_sum':weight_sum,
                                    'weight_min': min(weight1,weight2,weight3,weight4)}
                    
                    for col_feat in cols_features:
                        row_dict[col_feat] = 1 if col_feat in elem_glb else 0
                        
                    idx_df += 1
                    df_res.loc[idx_df] = row_dict
            

    print(idx_df)
    posRes = PoSet(universe, list( idxToWeightResult.keys() ) , idxToWeightResult)
    return df_res



def combine_3( pos1_nmb, pos2_nmb, pos3_nmb, df_gr):
    columns = list(df_gr.columns)
    columns.remove('id')
    columns.remove('idx_elem')
    columns.remove('auc')
    columns.remove('method')

    cols_features = columns
    columns = ['idx_elem','idx_elem1','idx_elem2','idx_elem3','auc1', 'auc2', 'auc3', 
              'diff_size1','diff_size2','diff_size3','size_elem', 'weight1','weight2','weight3','weight_sum', 'weight_min'] + cols_features
    df_res = df = pd.DataFrame(columns=columns)

    print(df_res.columns)

    # df_res.drop('id', axis=1, inplace=True)
    # [(pos1, first_nmb1), (pos2, first_nmb2),...]
    pos1, nmb1 = pos1_nmb[0], pos1_nmb[1]
    sortedIdxToWeightsNormalized1 = sorted(pos1.idxToWeightsNormalized.items(), key = lambda x: x[1], reverse=True)[:nmb1]

    print(sortedIdxToWeightsNormalized1)

    pos2, nmb2 = pos2_nmb[0], pos2_nmb[1]
    sortedIdxToWeightsNormalized2 = sorted(pos2.idxToWeightsNormalized.items(), key = lambda x: x[1], reverse=True)[:nmb2]

    print(sortedIdxToWeightsNormalized2)

    pos3, nmb3 = pos3_nmb[0], pos3_nmb[1]
    sortedIdxToWeightsNormalized3 = sorted(pos3.idxToWeightsNormalized.items(), 
                                           key = lambda x: x[1], reverse=True)[:nmb3]

    print(sortedIdxToWeightsNormalized3)


    universe = pos1.Universe
    idxToWeightPair   = { }
    idxToWeightResult = { }
    idx_df = 0

    for idx1, weight1 in sortedIdxToWeightsNormalized1:
        elem1 = universe.getElement(idx1)
        auc1 = pos1.idxToWeights[idx1]

        for idx2, weight2 in sortedIdxToWeightsNormalized2:
            elem2 = universe.getElement(idx2)
            auc2 = pos2.idxToWeights[idx2]

            for idx3, weight3 in sortedIdxToWeightsNormalized3:
                auc3 = pos3.idxToWeights[idx3]
                elem3 = universe.getElement(idx3)
                idx_glb = universe.joinIndexes(idx1, idx2, idx3)
                elem_glb = universe.getElement(idx_glb)
                weight_sum =  weight1 + weight2 + weight3
                idxToWeightResult[idx_glb] = max( idxToWeightResult.get(idx_glb, 0.0) , weight1 + weight2 + weight3)

                diff_size1 = universe.diff_size(elem_glb, elem1)
                diff_size2 = universe.diff_size(elem_glb, elem2)
                diff_size3 = universe.diff_size(elem_glb, elem3)

                # if idx_glb in idxToWeightResult.keys():
                idxToWeightPair[idx_glb] = (weight1, weight2, weight3)
                # idxToWeightResult[idx_glb] = max( idxToWeightResult.get(idx_glb, 0.0) , weight1 + weight2)
                print('-------------')
                # print(elem_glb)
                # print(idx_glb)
                # print(auc1)
                # print(auc2)
                # print(diff_size1)
                # print(diff_size2)
                # print(idx1, weight1)
                # print(idx2, weight2)   
                # print(weight_sum)    

                row_dict = {'idx_elem': idx_glb,
                                'idx_elem1': idx1,
                                'idx_elem2': idx2,
                                'idx_elem3': idx3,
                                'auc1': auc1, 
                                'auc2': auc2, 
                                'auc3': auc3, 
                                'diff_size1':diff_size1,
                                'diff_size2':diff_size2, 
                                'diff_size3':diff_size3, 
                                'size_elem': len(elem_glb),
                                'weight1':weight1,
                                'weight2':weight2,
                                'weight3':weight3,
                                'weight_sum':weight_sum,
                                'weight_min': min(weight1,weight2,weight3)}
                
                for col_feat in cols_features:
                    row_dict[col_feat] = 1 if col_feat in elem_glb else 0
                    
                idx_df += 1
                df_res.loc[idx_df] = row_dict
            

    print(idx_df)
    posRes = PoSet(universe, list( idxToWeightResult.keys() ) , idxToWeightResult)
    return df_res



def combine_2( pos1_nmb, pos2_nmb, df_gr):
    columns = list(df_gr.columns)
    columns.remove('id')
    columns.remove('idx_elem')
    columns.remove('auc')
    columns.remove('method')

    cols_features = columns
    columns = ['idx_elem','idx_elem1','idx_elem2','auc1', 'auc2', 
              'diff_size1','diff_size2','size_elem', 'weight1','weight2','weight_sum','weight_min'] + cols_features
    df_res = df = pd.DataFrame(columns=columns)

    print(df_res.columns)

    # df_res.drop('id', axis=1, inplace=True)
    # [(pos1, first_nmb1), (pos2, first_nmb2),...]
    pos1, nmb1 = pos1_nmb[0], pos1_nmb[1]
    sortedIdxToWeightsNormalized1 = sorted(pos1.idxToWeightsNormalized.items(), key = lambda x: x[1], reverse=True)[:nmb1]

    print(sortedIdxToWeightsNormalized1)

    pos2, nmb2 = pos2_nmb[0], pos2_nmb[1]
    sortedIdxToWeightsNormalized2 = sorted(pos2.idxToWeightsNormalized.items(), key = lambda x: x[1], reverse=True)[:nmb2]

    print(sortedIdxToWeightsNormalized2)

    universe = pos1.Universe
    idxToWeightPair   = { }
    idxToWeightResult = { }
    idx_df = 0

    for idx1, weight1 in sortedIdxToWeightsNormalized1:
        elem1 = universe.getElement(idx1)
        auc1 = pos1.idxToWeights[idx1]
        for idx2, weight2 in sortedIdxToWeightsNormalized2:
            auc2 = pos2.idxToWeights[idx2]
            elem2 = universe.getElement(idx2)
            idx_glb = universe.joinIndexes(idx1, idx2)
            elem_glb = universe.getElement(idx_glb)
            weight_sum =  weight1 + weight2
            idxToWeightResult[idx_glb] = max( idxToWeightResult.get(idx_glb, 0.0) , weight1 + weight2)

            diff_size1 = universe.diff_size(elem_glb, elem1)
            diff_size2 = universe.diff_size(elem_glb, elem2)

            # if idx_glb in idxToWeightResult.keys():
            idxToWeightPair[idx_glb] = (weight1, weight2)
            # idxToWeightResult[idx_glb] = max( idxToWeightResult.get(idx_glb, 0.0) , weight1 + weight2)
            print('-------------')
            # print(elem_glb)
            # print(idx_glb)
            # print(auc1)
            # print(auc2)
            # print(diff_size1)
            # print(diff_size2)
            # print(idx1, weight1)
            # print(idx2, weight2)   
            # print(weight_sum)    

            row_dict = {'idx_elem': idx_glb,
                            'idx_elem1': idx1,
                            'idx_elem2': idx2,
                            'auc1': auc1, 
                            'auc2': auc2, 
                            'diff_size1':diff_size1,
                            'diff_size2':diff_size2, 
                            'size_elem': len(elem_glb),
                            'weight1':weight1,
                            'weight2':weight2,
                            'weight_sum':weight_sum,
                            'weight_min': min(weight1,weight2)}
            
            for col_feat in cols_features:
                row_dict[col_feat] = 1 if col_feat in elem_glb else 0
                 
            idx_df += 1
            df_res.loc[idx_df] = row_dict
            

    print(idx_df)
    posRes = PoSet(universe, list( idxToWeightResult.keys() ) , idxToWeightResult)
    return df_res




def load_posetgr(univ, df_gr):
    poset = PoSet(univ, [0], {0:0.0})

    columns = list(df_gr.columns)
    indexes_df = [ ]
    for idx, row in df_gr.iterrows():
        # print(row)
        auc = row['auc']
        element = [columns[idx] for idx in range(3, len(row)) if  row[idx]==1.0]
        element = set(element)
        element_idx = univ.getIndex(element)
        indexes_df.append(element_idx)
        poset.mergeElementWeight(element, weight=auc)
    poset.removeIndex(0)

    df_gr['idx_elem'] = indexes_df
    return poset, df_gr


import pandas as pd
import numpy as np

method = 'xb'

def load_group_best_results(method):
    df_gr1 = pd.read_csv('csv_best_original/best_json_group_1_10k_top100.csv')
    df_gr1 = df_gr1[df_gr1['method']==method]
    df_gr1 = df_gr1.fillna(0.0)

    df_gr2 = pd.read_csv('csv_best_original/best_json_group_2_10k_top100.csv')
    df_gr2 = df_gr2[df_gr2['method']==method]
    df_gr2 = df_gr2.fillna(0.0)


    df_gr3 = pd.read_csv('csv_best_original/best_json_group_3_10k_top100.csv')
    df_gr3 = df_gr3[df_gr3['method']==method]
    df_gr3 = df_gr3.fillna(0.0)

    df_gr4 = pd.read_csv('csv_best_original/best_json_group_4_10k_top100.csv')
    df_gr4 = df_gr4[df_gr4['method']==method]
    df_gr4 = df_gr4.fillna(0.0)

    return df_gr1, df_gr2, df_gr3, df_gr4


method = 'rf'

def run_by_method(method, number_best):
    df_gr1, df_gr2, df_gr3, df_gr4 = load_group_best_results(method)

    univ02 = UniversePower(df_gr4.columns[3:])

    poset1, df_gr1 = load_posetgr(univ02, df_gr1)
    poset2, df_gr2 = load_posetgr(univ02, df_gr2)
    poset3, df_gr3 = load_posetgr(univ02, df_gr3)
    poset4, df_gr4 = load_posetgr(univ02, df_gr4)

# df = pd.DataFrame(columns=['a','b'])
# df.loc[1] = {'a':3, 'b': 4}



    # df_result=combine_2((poset1, 50), (poset2, 50), df_gr3)
    # df_result.to_csv('result_'+ method + '_gr1gr2.csv',index_label='id')

    # df_result=combine_2((poset1, 50), (poset4, 50), df_gr1)
    # df_result.to_csv('result_'+ method + '_gr1gr4.csv',index_label='id')

    # df_result=combine_3((poset1, 30), (poset2, 30), (poset4, 30), df_gr1)
    # df_result.to_csv('result_'+ method + '_gr1gr2gr4.csv',index_label='id')

    df_result=combine_4((poset1, number_best), 
                        (poset2, number_best), 
                        (poset3, number_best), 
                        (poset4, number_best), 
                        df_gr1)

    df_result.to_csv('result_'+ method +'_'+ str(number_best) +'_gr1gr2gr3gr4.csv',index_label='id')

run_by_method('lr', number_best)
run_by_method('xb', number_best)
run_by_method('rf', number_best)

# df_gr1.to_csv('df_'+ method +'_gr1.csv',index=None)
# df_gr2.to_csv('df_'+ method +'_gr2.csv',index=None)
# df_gr3.to_csv('df_'+ method +'_gr3.csv',index=None)
# df_gr4.to_csv('df_'+ method +'_gr4.csv',index=None)



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

