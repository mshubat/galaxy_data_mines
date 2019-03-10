from treelib import Node, Tree
from astropy.table import Table
import numpy as np
import string
import os


class TreeBuilder:

    def buildTestTree(self):
        '''
        This method builds & returns a simple tree used for testing.
        '''
        tree = Tree()

        # create a node with name & identifier respectively
        # 0th level
        tree.create_node("ROOT", "root")  # root node

        # level 1
        tree.create_node("A", "a", parent="root")
        tree.create_node("T", "t", parent="root")

        # level 2
        tree.create_node("B", "b", parent="a")
        tree.create_node("C", "c", parent="a")
        tree.create_node("X", "x", parent="a")
        tree.create_node("Y", "y", parent="a")
        tree.create_node("Z", "z", parent="a")

        # level 3
        tree.create_node("D", "d", parent="b")
        tree.create_node("E", "e", parent="b")
        tree.create_node("F", "f", parent="c")
        tree.create_node("G", "g", parent="c")

        # level 4
        tree.create_node("H", "h", parent="d")
        # tree.create_node("I", "i", parent="d")

        return tree

    def buildTree(self, infile):
        '''
        Read SIMBAD class structure and put it into the treelib tree structure

        Returns both a tree and the dictionary needed to translate between nodeIDs
        and tags
        '''
        new_tree = Tree()
        # create classless root node Since SIMBAD classes
        # are really separate trees.
        new_tree.create_node("ROOT", "root")
        tree_dict = {}
        dat = Table.read(os.path.join(os.path.dirname(__file__), infile), format='csv')
        numeric_id = dat['id']
        tag = dat['tag']
        for i, num_id in enumerate(numeric_id):
            #        id_str = num_id.replace(".","")
            new_tree = TreeBuilder.add_to_tree(num_id, tag[i], new_tree)
            tree_dict[tag[i]] = num_id

        return(new_tree, tree_dict)

    def add_to_tree(id, tag, tree):
        '''
        Add an entry to treelib tree, using ab.cd.ef.gh as identifier
        where the parent of ab.cd.ef.gh is ab.cd.eg.00
        '''
        parent = TreeBuilder.get_parent(id)
        if parent == '00.00.00.00':
            tree.create_node(tag, id, parent=tree.root)  # add to root
        elif tree.contains(parent):  # add as child to parent
            tree.create_node(tag, id, parent=parent)
        else:  # shouldn't get here, SIMBAD tree is in order
            print('create parent first')
        return(tree)

    def get_parent(str_id):
        '''
        Find the first occurrence of 00 in string,
        replace previous 2 and all following digits with zeros
        '''
        double_zero = str_id.find("00")
        if double_zero == -1:  # not found
            parent = str_id[:-3]
        elif double_zero < 4:
            parent = '00'
        else:
            parent = str_id[:double_zero-4]

        # pad out to length 11
        while len(parent) < 11:
            parent = parent+'.00'
        return(parent)

    def match_tag_to_nid(tree, tree_dict):
        '''
        Matches identifier and tag of each node in the tree.
        '''
