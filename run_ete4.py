#!/usr/bin/env python3

import sys
from ete4 import Tree, SeqGroup, PhyloTree
from ete4.smartview import TreeStyle, NodeStyle, TreeLayout
from ete4.smartview.renderer.layouts.seq_layouts import LayoutAlignment

sys.path.append('/data/projects/og_delineation_web/layouts')
import layouts_emapper
import layout_cogs

t = Tree('((9606.A:1, 9598.B:1)0.5:1,10090.C:1);', parser = 0)


if len(sys.argv) >1:
    t = PhyloTree(open(sys.argv[1]), parser = 0)
    if sys.argv[2] == 'gtdb':
        import layout_web_gtdb as layouts_ogd
    if sys.argv[2] == 'ncbi':
        import layouts_ogd
else:
    import layouts_ogd


        
layouts = [
    # Layout for OGD annotations
    TreeLayout('Scientific name', ns = layouts_ogd.get_layout_leafname(), aligned_faces = True),
    TreeLayout('Evolution events',  ns = layouts_ogd.get_layout_evoltype()),
    TreeLayout('Last common ancestor',  ns = layouts_ogd.get_layout_lca_rects(t),  aligned_faces = True),
    TreeLayout('Collapse OG',  ns = layouts_ogd.collapse_og(),  aligned_faces = True),
    TreeLayout('Seqs out OGS', ns = layouts_ogd.seqs_out_og()),
    TreeLayout('Collapse order rank', ns = layouts_ogd.collapse_rank()),
    TreeLayout('Species Overlap', ns = layouts_ogd.get_species_overlap()),
    TreeLayout('Long branches', ns = layouts_ogd.long_branch_outlier()),
    TreeLayout('Taxonomic Outliers', ns = layouts_ogd.taxo_out()),
    TreeLayout('Brackground OG', ns = layouts_ogd.background_og()),
    TreeLayout('Brackground mOG', ns = layouts_ogd.background_mog(),aligned_faces = True ),
    

    #Layout for emapper annotations
    TreeLayout('emapper_pref_name', ns = layouts_emapper.get_emapper_pref_name(), aligned_faces = True),
    TreeLayout('eggnog_OGs', ns = layouts_emapper.get_eggnog_OG(), aligned_faces = True),
    TreeLayout('emapper_kegg_ko', ns = layouts_emapper.get_emapper_kegg_ko(), aligned_faces = True),
    TreeLayout('PFAM_arq', ns = layouts_emapper.get_pfams(), aligned_faces = True ),


    #Layout for COG benkchmark
    TreeLayout('original_COG', ns = layout_cogs.get_COG(), aligned_faces = True)



    #TreeLayout('Prune tree', ns = layout_web.prune_tree()),
 
    ]


# Turn Off layouts. By default all layout are activated
for l in layouts:
    if l.name in ['Seqs out OGS', 'Collapse', 'Collapse order rank','Long branches', 'Taxonomic Outliers'
                    'emapper_pref_name', 'eggnog_OGs', 'PFAM_arq' ]:
        l.active = False



# Launch smartview explorer
t.dist = 0.05

name = 'upload_tree'

props_popup = ['node_is_og', 'dist', 'species_losses', 'node_create_og', 'lca_node_name', 'len_leaves_in', 
    'taxid', 'sci_name', 'lineage', 'lca_dup', 'inparalogs_rate', 'ch1_name', 'dup_node_name', 'is_root', 
    'total_leaves', 'dups_up', 'ogs_up', 'common_name', 'dups_down', 'so_score_dup','ogs_down', 'score1', 
    'rank', 'dup_lineage', 'lca_node', 'len_leaves_out', 'species_losses_percentage', 'name', 'ch2_name', 
    'score2', 'sp_out', 'so_score', 'leaves_out','dup_score', 'overlap', 'evoltype_2', 'mOG', 'len_sp_in', 
    'best_tax', 'node_is_mog', 'recover_seqs', 'recover_in', 'Preferred_name', 'Preferred_name_counter',
    'eggNOG_OGs_counter', 'eggNOG_OGs', 'long_branch_outlier']



# Trigger the interactive web server
t.explore(name = name, layouts = layouts, show_leaf_name = False , include_props = props_popup, keep_server=True , host = '138.4.138.141', open_browser =False, port = 5000)