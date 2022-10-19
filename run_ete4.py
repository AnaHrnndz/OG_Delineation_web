from ete4 import Tree, SeqGroup
from ete4.smartview import TreeStyle, NodeStyle, TreeLayout
from ete4.smartview.renderer.layouts.seq_layouts import LayoutAlignment
from layout_web import get_layout_leafname, get_layout_evoltype, get_layout_lca_rects, collapse_og, seqs_out_og, collapse_rank, get_species_overlap, get_pnames, get_ogs, get_doms, prune_tree #sequencebouncer
from layout_possvm import get_pnames_possvm, get_og_euk, get_og_met, get_og_possvm, get_og_ref, get_pfams


t = Tree('((9606.A:1, 9598.B:1)0.5:1,10090.C:1);', format = 0)

#t = Tree(sys.argv[1], format = 0)

#raw_alg = SeqGroup(sys.argv[2])




layouts = [
    TreeLayout('Scientific name', ns = get_layout_leafname(), aligned_faces = True), 
    TreeLayout('Evolution events',  ns = get_layout_evoltype()),
    TreeLayout('Last common ancestor',  ns = get_layout_lca_rects(t),  aligned_faces = True),
    TreeLayout('Collapse OG',  ns = collapse_og(),  aligned_faces = True),
    TreeLayout('Seqs out OGS', ns = seqs_out_og()),
    TreeLayout('Collapse order rank', ns = collapse_rank()),
    TreeLayout('Species Overlap', ns = get_species_overlap()),
    TreeLayout('Pref name', ns = get_pnames(), aligned_faces = True),
    TreeLayout('EGGNOG', ns = get_ogs(), aligned_faces = True),
    TreeLayout('PFAM Domains', ns = get_doms() ),
    TreeLayout('Prune tree', ns = prune_tree()),
    #TreeLayout('SequenceBouncer', ns=sequencebouncer())


    TreeLayout('possvm_pname', ns = get_pnames_possvm(), aligned_faces = True),
    TreeLayout('possvm_OG_Euk', ns = get_og_euk(), aligned_faces = True),
    TreeLayout('possm_OG_Metazoa', ns = get_og_met(), aligned_faces = True),
    TreeLayout('possvm_OG', ns = get_og_possvm(), aligned_faces = True),
    TreeLayout('possvm_RefOG', ns = get_og_ref(), aligned_faces = True),
    TreeLayout('possvm_PFAM', ns = get_pfams(), aligned_faces = True)        
    #LayoutAlignment(sys.argv[2], width=500)
        ]



for l in layouts:
    if l.name in ['Last common ancestor', 'Seqs out OGS', 'Alignment', 'Pref name', 'EGGNOG', 'PFAM Domains', 'possvm_PFAM', 'Prune tree']:
        l.active = False
    


# Launch smartview explorer
t.dist = 0.05

name = 'upload_tree'


# Trigger the interactive web server
t.explore(tree_name = name, layouts = layouts, show_leaf_name = False , host = '138.4.138.141' )
