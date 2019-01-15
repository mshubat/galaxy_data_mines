from ComparisonTree import ComparisonTree


if __name__ == "__main__":
    # Create the ComparisonTree and show the tree
    ct = ComparisonTree()
    ct.showTree()

    print("Sibling Method Tests")
    print("---------------------")
    print("f and x are siblings: {}".format(ct.areSiblings("f", "x")))
    print("a and a are siblings: {}".format(ct.areSiblings("a", "a")))
    print("c and x are siblings: {}".format(ct.areSiblings("c", "x")))
    print("f and g are siblings: {}".format(ct.areSiblings("f", "g")))
    print()

    print("Parent Method Tests")
    print("---------------------")
    print("a is parent of b: {}".format(ct.isParentOf("a", "b")))
    print("b is parent of a: {}".format(ct.isParentOf("b", "a")))
    print("b is parent of b: {}".format(ct.isParentOf("b", "b")))
    print()

    print("isOfType Method Tests")
    print("---------------------")
    print("d is of type b: {}".format(ct.isOfType("d", "b")))
    print("c is of type a: {}".format(ct.isOfType("c", "a")))
    print("d is of type a: {}".format(ct.isOfType("d", "a")))
    print("a is of type a: {}".format(ct.isOfType("a", "a")))
    print("b is of type a: {}".format(ct.isOfType("b", "c")))
    print()

    print("shareCommonAncestor Tests")
    print("-------------------------")

    #print(ct.areSiblings("a", "a"))

    #print(ct.isOfType("d", "b"))

    #print(ct.test("d", "b"))

#classtree.isOfType("f", "c")
