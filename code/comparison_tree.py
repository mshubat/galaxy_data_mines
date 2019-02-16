from treelib import Node, Tree
from tree_builder import TreeBuilder
from data_controller import DataController
from astropy.table import Table, Column
import numpy as np


class ComparisonTree:

    def __init__(self, *, run_mode):
        self.builder = TreeBuilder()
        self.run_mode = run_mode

        if run_mode:
            print("Run Mode: Enabled. Converting numeric\
             id's to SIMBAD entries.")
            dictandtree = self.builder.buildTree("simbad_raw.csv")
            self.tree = dictandtree[0]
            self.dict = dictandtree[1]
        else:
            self.tree = self.builder.buildTestTree()

    def showTree(self):
        '''
        Allows ComparisonTree user to call Tree show() function more easily.
        '''
        self.tree.show()

    def areDirectMatch(self, Type_N, Type_S):
        return Type_N == Type_S

    def areCandidateMatch(self, *, Type_N, Type_S):

        if (Type_N + "?" == Type_S):
            return True

    def areSiblings(self, firstNodeId, secNodeId):
        '''
        This function will return True if two nodes are siblings
        in the tree provided. It will return False otherwise.
        '''
        if self.run_mode:
            firstNodeId = self.dict[firstNodeId]
            secNodeId = self.dict[secNodeId]

        siblings = self.tree.siblings(firstNodeId)

        for s in siblings:
            if s.identifier == secNodeId:
                return True
        return False

    def isParentOf(self, firstNodeId, secNodeId):
        '''
        This method will return True if the first node is the parent of the
        second. It will return False otherwise.
        '''
        if self.run_mode:
            firstNodeId = self.dict[firstNodeId]
            secNodeId = self.dict[secNodeId]

        parent = self.tree.parent(secNodeId)

        # if second node has no parent
        if (parent == None):
            return False

        return parent.identifier == firstNodeId

    def isOfType(self, firstNodeId, secNodeId):
        '''
        This method will return True if the first node is a descendent of the
        second (ie. if it is the the same type). It will return False otherwise.
        '''
        if self.run_mode:
            firstNodeId = self.dict[firstNodeId]
            secNodeId = self.dict[secNodeId]

        subtree = self.tree.subtree(secNodeId)

        return subtree.contains(firstNodeId)

    def shareCommonAncestor(self, firstNodeId, secNodeId):
        '''
        This method will return True if the nodes corresponding to the given
        nodeIDs share a common ancestor. Return False othewise.
        '''
        if self.run_mode:
            firstNodeId = self.dict[firstNodeId]
            secNodeId = self.dict[secNodeId]

        first_node = self.tree.get_node(firstNodeId)
        sec_node = self.tree.get_node(secNodeId)

        # If the nodes are at the top level they do not share a common
        # ancestor. ie. Star and Galaxy. and their parent is root.
        if (self.tree.parent(firstNodeId).is_root() or
                self.tree.parent(secNodeId).is_root()):
            return False
        else:
            return True

        # Otherwise traverse up to look for common ancestor.
        current_node_id = firstNodeId
        p = self.tree.parent(current_node_id)
        while (p and not p.identifier == "root"):
            if (self.isOfType(secNodeId, p.identifier)):
                return True
            else:
                p = self.tree.parent(p.identifier)
        return False

    def compare_objects(self, t):
        '''
        Once you have the matched table is either generated or imported
        we need to do the comparisons between NED and SIMBAD objects.

        t - the combined table fetched or generated.
        '''
        tsize = len(t)

        cols = []
        cols.append(Column(name='exactMatch', length=tsize, dtype=bool))
        cols.append(Column(name='candidateMatch', length=tsize, dtype=bool))
        cols.append(Column(name='areSiblings', length=tsize, dtype=bool))
        # cols.append(Column(name='shareCommonAncestor', length=tsize, dtype=bool))

        # Insert NED Analogues right after the Name_N column.
        t.add_column(Column(name="Type_N_Analogue", length=tsize, dtype=np.dtype(('U', 8))),
                     index=t.index_column("Type_N")+1)
        t.add_column(Column(name="Type_S_cond", length=tsize, dtype=np.dtype(('U', 8))),
                     index=t.index_column("Type_S")+1)
        t.add_columns(cols)

        for i in range(0, tsize):
            ned_analogue = DataController.ned_to_simbad(t["Type_N"][i])
            t["Type_N_Analogue"][i] = ned_analogue
            t["Type_S_cond"][i] = DataController.simbad_long_to_small(t["Type_S"][i])

            if (self.areDirectMatch(t["Type_N_Analogue"][i], t["Type_S"][i])):
                t["exactMatch"][i] = True
                print("match i={} - N: {} S: {}".format(i,
                                                        ned_analogue,
                                                        t["Type_S"][i]))
            elif (self.areCandidateMatch(Type_N=t["Type_N_Analogue"][i], Type_S=t["Type_S"][i])):
                t["candidateMatch"][i] = True
                print("candidateMatch i={} - N: {} S: {}".format(i,
                                                                 ned_analogue,
                                                                 t["Type_S"][i]))
            else:
                print("non-match i={} - N: {} S: {}".format(i,
                                                            ned_analogue,
                                                            t["Type_S"][i]))

            '''elif (self.areSiblings(t["Type_N_Analogue"][i], t["Type_S"][i])):
                t["areSiblings"][i] = True
                print("areSiblings i={} - N: {} S: {}".format(i,
                                                              ned_analogue,
                                                              t["Type_S"][i]))'''
