from treelib import Node, Tree
from TreeBuilder import TreeBuilder


class ComparisonTree:

    def __init__(self):
        self.builder = TreeBuilder()
        self.tree = self.builder.buildTree("simbad_raw.csv")

    # def buildTree(self):
        '''
        This function will soon encorporate reading SIMBAD hierarchy from file
        into tree structure. Then it will focus on
        '''
        #tree = Tree()

        '''
        # create a node with name & identifier respectively
        tree.create_node("A", "a")  # root node
        tree.create_node("B", "b", parent="a")
        tree.create_node("C", "c", parent="a")
        tree.create_node("X", "x", parent="a")
        tree.create_node("Y", "y", parent="a")
        tree.create_node("Z", "z", parent="a")

        tree.create_node("D", "d", parent="b")
        tree.create_node("E", "e", parent="b")
        tree.create_node("F", "f", parent="c")
        tree.create_node("G", "g", parent="c")
        '''
        # return tree

    def showTree(self):
        '''
        Allows class user to call Tree show() function more easily.
        '''
        self.tree.show()

    def areSiblings(self, firstNodeId, secNodeId):
        '''
        This function will return True if two nodes are siblings
        in the tree provided. It will return False otherwise.
        '''
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
        subtree = self.tree.subtree(secNodeId)

        return subtree.contains(firstNodeId)

    def test(self, firstNodeId, secNodeId):
        return areSiblings(firstNodeId, secNodeId)

    def shareCommonAncestor(firstNodeId, secNodeId):
        '''
        This method will return True if the nodes corresponding to the given
        nodeIDs share a common ancestor. Return False othewise.
        '''
        if (areSiblings(firstNodeId, secNodeId)):
            # If the nodes are at the top level they are not actually related.
            # ie. Star and Galaxy.
            if ((is_root(firstNodeId) or is_root(secNodeId))):
                return False
            else:
                return True
        elif (isOfType(firstNodeId, secNodeId) or isOfType(secNodeId, firstNodeId)):
            return True
        else:
            current_node_id = firstNodeId
            p = parent(current_node_id)
            while (p):
                if (isOfType(secNodeId, p.identifier)):
                    return True
            return false
