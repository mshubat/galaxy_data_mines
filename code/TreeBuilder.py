from treelib import Node, Tree


class TreeBuilder:

    """read SIMBAD class structure and put it into the treelib tree structure

    returns both a tree and the dictionary needed to translate between nodeIDs
    and tags"""

    def buildTree():
        new_tree = Tree()
        tree_dict = {}
        dat = Table.read(infile, format='csv')
        numeric_id = dat['id']
        tag = dat['tag']
        for i, num_id in enumerate(numeric_id):
            #        id_str = num_id.replace(".","")
            new_tree = add_to_tree(num_id, tag[i], new_tree)
            tree_dict[tag[i]] = num_id
        return(new_tree, tree_dict)

    """add an entry to treelib tree, using ab.cd.ef.gh as identifier
       where the parent of ab.cd.ef.gh is ab.cd.eg.00"""

    def add_to_tree(id, tag, tree):
        parent = get_parent(id)
        if parent == '00.00.00.00':
            tree.create_node(tag, id, parent=tree.root)  # add to root
        elif tree.contains(parent):  # add as child to parent
            tree.create_node(tag, id, parent=parent)
        else:  # shouldn't get here, SIMBAD tree is in order
            print('create parent first')
        return(tree)

    """find the first occurrence of 00 in string,
       replace previous 2 and all following digits with zeros"""

    def get_parent(str_id):
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
