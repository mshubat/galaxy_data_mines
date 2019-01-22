from treelib import Node, Tree, exceptions
from ComparisonTree import ComparisonTree


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
    print()


def isOfTypeMethodTest(ct):
    print("isOfType Method Tests")
    print("---------------------")
    print("d is of type b: {}".format(ct.isOfType("d", "b")))
    print("c is of type a: {}".format(ct.isOfType("c", "a")))
    print("d is of type a: {}".format(ct.isOfType("d", "a")))
    print("a is of type a: {}".format(ct.isOfType("a", "a")))
    print("b is of type c: {}".format(ct.isOfType("b", "c")))
    print("d is of type c: {}".format(ct.isOfType("d", "c")))
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


if __name__ == "__main__":

    # Create the ComparisonTree and show the tree

    ct = ComparisonTree()
    ct.showTree()

    # Run tests
    try:
        siblingMethodTest(ct)
        parentMethodTest(ct)
        isOfTypeMethodTest(ct)
        shareCommonAncestorMethodTest(ct)
    except exceptions.NodeIDAbsentError:
        print(print("Oh no. It looks like that node is not in the tree."))
