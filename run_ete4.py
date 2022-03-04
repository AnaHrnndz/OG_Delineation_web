from ete4 import Tree
from ete4.smartview import TreeStyle, NodeStyle, TreeLayout
from layout_web import get_layout_leafname, get_layout_evoltype, get_layout_lca_rects, collapse_og

#t = Tree('/data/projects/find_ogs/delineation_og_test/post_14-3-3.nw', format =1)

#t = Tree('(A.A:1,(B.:1,(E:1,D:1)Internal_1:0.5)Internal_2:0.5)Root;', format = 1)
t = Tree('((9606.A:1, 9598.B:1)D:1,10090.C:1);', format = 1)
#t = Tree('/data/projects/find_ogs/delineation_og_test/post_14-3-3.nw', format = 1)


layouts = [
    TreeLayout('Scientific name', ns = get_layout_leafname(), aligned_faces = True), 
    TreeLayout('Evolution events',  ns = get_layout_evoltype()),
    TreeLayout('Last common ancestor',  ns = get_layout_lca_rects(t),  aligned_faces = True),
    TreeLayout('Collapse OG',  ns = collapse_og(),  aligned_faces = True)
        ]

for l in layouts:
    if l.name in ['Last common ancestor']:
        l.active = False
    



# Launch smartview explorer
t.dist = 0.01

name = 'upload_tree'


# Trigger the interactive web server
t.explore(tree_name = name, layouts = layouts, show_leaf_name = False)
