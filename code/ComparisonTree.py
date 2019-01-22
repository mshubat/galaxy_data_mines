from treelib import Node, Tree
from TreeBuilder import TreeBuilder


class ComparisonTree:

    def __init__(self):
        self.builder = TreeBuilder()
        #self.tree = self.builder.buildTree("simbad_raw.csv")
        self.tree = self.builder.buildTestTree()

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

    def shareCommonAncestor(self, firstNodeId, secNodeId):
        '''
        This method will return True if the nodes corresponding to the given
        nodeIDs share a common ancestor. Return False othewise.
        '''
        first_node = self.tree.get_node(firstNodeId)
        sec_node = self.tree.get_node(secNodeId)

        if (self.areSiblings(firstNodeId, secNodeId)):
            # If the nodes are at the top level they do not share a common
            # ancestor. ie. Star and Galaxy.
            if ((first_node.is_root() or sec_node.is_root())):
                return False
            else:
                return True
        # If one of the nodes is of the type of the other
        elif (self.isOfType(firstNodeId, secNodeId) or self.isOfType(secNodeId, firstNodeId)):
            # they share a common ancestor only if both are not root.
            if (first_node.is_root() or sec_node.is_root()):
                return False
            return True
        else:
            current_node_id = firstNodeId
            p = self.tree.parent(current_node_id)
            while (p):
                # print(p)
                if (self.isOfType(secNodeId, p.identifier)):
                    return True
            return False
