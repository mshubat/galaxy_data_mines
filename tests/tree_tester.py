# Tree Tester Class
#
# Class is used to unit test code and methods. Prior to data integration.
#
#

from galaxy_data_mines.comparison_tree import ComparisonTree
from galaxy_data_mines.data_controller import DataController
from treelib import Node, Tree, exceptions
from astropy.table import Table


def siblingMethodTest(ct):
    print("Sibling Method Tests")
    print("---------------------")
    print("f and x are siblings: {}".format(ct.areSiblings("f", "x")))
    print("a and a are siblings: {}".format(ct.areSiblings("a", "a")))
    print("c and x are siblings: {}".format(ct.areSiblings("c", "x")))
    print("f and g are siblings: {}".format(ct.areSiblings("f", "g")))
    print("f and d are siblings: {}".format(ct.areSiblings("f", "d")))
    print()


def parentMethodTest(ct):
    print("Parent Method Tests")
    print("---------------------")
    print("a is parent of b: {}".format(ct.isParentOf("a", "b")))
    print("b is parent of a: {}".format(ct.isParentOf("b", "a")))
    print("b is parent of b: {}".format(ct.isParentOf("b", "b")))
    print("b is parent of d: {}".format(ct.isParentOf("b", "d")))
    print("d is parent of h: {}".format(ct.isParentOf("d", "h")))
    print()


def ofTypeMatchTest(ct):
    print("ofType Match between d and b: {}".format(ct.ofTypeMatch("d", "b")))
    print("ofType Match between c and a: {}".format(ct.ofTypeMatch("c", "a")))
    print("ofType Match between d and a: {}".format(ct.ofTypeMatch("d", "a")))
    print("ofType Match between a and a: {}".format(ct.ofTypeMatch("a", "a")))
    print("ofType Match between b and c: {}".format(ct.ofTypeMatch("b", "c")))
    print("ofType Match between d and c: {}".format(ct.ofTypeMatch("d", "c")))
    print("ofType Match between d and d: {}".format(ct.ofTypeMatch("d", "d")))
    print("ofType Match between c and x: {}".format(ct.ofTypeMatch("c", "x")))
    print("ofType Match between d and h: {}".format(ct.ofTypeMatch("d", "h")))
    print("ofType Match between h and d: {}".format(ct.ofTypeMatch("h", "d")))
    print("ofType Match between a and b: {}".format(ct.ofTypeMatch("a", "b")))
    print("ofType Match between b and a: {}".format(ct.ofTypeMatch("b", "a")))
    print("ofType Match between a and h: {}".format(ct.ofTypeMatch("a", "h")))
    print("ofType Match between a and d: {}".format(ct.ofTypeMatch("a", "d")))
    print("ofType Match between a and x: {}".format(ct.ofTypeMatch("a", "x")))
    print("ofType Match between b and h: {}".format(ct.ofTypeMatch("b", "h")))
    print("ofType Match between b and f: {}".format(ct.ofTypeMatch("b", "f")))


def isDescendantOfMethodTest(ct):
    print("isDescendantOf Method Tests")
    print("---------------------")
    print("d is descendant of b: {}".format(ct.isDescendantOf("d", "b")))
    print("c is descendant of a: {}".format(ct.isDescendantOf("c", "a")))
    print("d is descendant of a: {}".format(ct.isDescendantOf("d", "a")))
    print("a is descendant of a: {}".format(ct.isDescendantOf("a", "a")))
    print("b is descendant of c: {}".format(ct.isDescendantOf("b", "c")))
    print("d is descendant of c: {}".format(ct.isDescendantOf("d", "c")))
    print()


def shareCommonAncestorMethodTest(ct):
    print("shareCommonAncestor Tests")
    print("-------------------------")
    print("a and b share a common ancestor: {}".format(
        ct.shareCommonAncestor("a", "b")))
    print("b and d share a common ancestor: {}".format(
        ct.shareCommonAncestor("b", "d")))
    print("a and a share a common ancestor: {}".format(
        ct.shareCommonAncestor("a", "a")))
    print("x and y share a common ancestor: {}".format(
        ct.shareCommonAncestor("x", "y")))
    print("f and g share a common ancestor: {}".format(
        ct.shareCommonAncestor("f", "g")))
    print("c and f share a common ancestor: {}".format(
        ct.shareCommonAncestor("c", "f")))
    print("d and c share a common ancestor: {}".format(
        ct.shareCommonAncestor("d", "c")))
    print("a and t share a common ancestor: {}".format(
        ct.shareCommonAncestor("a", "t")))
    print("b and t share a common ancestor: {}".format(
        ct.shareCommonAncestor("b", "t")))
    print("t and h share a common ancestor: {}".format(
        ct.shareCommonAncestor("t", "h")))
    print("h and t share a common ancestor: {}".format(
        ct.shareCommonAncestor("h", "t")))
    print("h and h share a common ancestor: {}".format(
        ct.shareCommonAncestor("h", "h")))
    print("h and z share a common ancestor: {}".format(
        ct.shareCommonAncestor("h", "z")))


if __name__ == "__main__":

    # Create and Show the ComparisonTree.
    ct = ComparisonTree(run_mode=False)
    ct.showtree()

    # ct.tree.save2file("savedTree")
    # ct.tree.save2file("savedTree", nid="10.00.00.00")

    # str_test = ct.tree.to_json()
    # dic_test = ct.tree.to_dict()

    # print(dic_test)
    # print(dic_test["ROOT"])
    # print(dic_test["ROOT"]["children"])

   # Run tests
    try:
        ofTypeMatchTest(ct)
        # siblingMethodTest(ct)
        # parentMethodTest(ct)
        # isDescendantOfMethodTest(ct)
        # shareCommonAncestorMethodTest(ct)
    except exceptions.NodeIDAbsentError:
        print(print("Oh no. It looks like that node is not in the tree."))
