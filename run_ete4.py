#!/usr/bin/env python3

import sys
from ete4 import Tree, SeqGroup, PhyloTree
from ete4.smartview import TreeStyle, NodeStyle, TreeLayout
from ete4.smartview.renderer.layouts.seq_layouts import LayoutAlignment

sys.path.append('/data/projects/og_delineation_web/layouts')
import layouts_emapper

t = Tree('((9606.A:1, 9598.B:1)0.5:1,10090.C:1);', parser = 0)


if len(sys.argv) >1:
    t = PhyloTree(open(sys.argv[1]))
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
    TreeLayout('PFAM_arq', ns = layouts_emapper.get_pfams(), aligned_faces = True )


    #TreeLayout('Prune tree', ns = layout_web.prune_tree()),
 
    ]


# Turn Off layouts. By default all layout are activated
for l in layouts:
    if l.name in ['Seqs out OGS', 'Alignment', 'EGGNOG', 'PFAM Domains', 'Long branches','possvm_PFAM', 
                    'Prune tree', 'Collapse order rank', 'possm_OG_Metazoa', 'Taxonomic Outliers']:
        l.active = False



# Launch smartview explorer
t.dist = 0.05

name = 'upload_tree'

props_popup = ['taxo_outlier', 'mOG','node_create_og', 'inparalogs_rate', 'BEST_RESP_NUM', 'evoltype_2', 
                'so_cell_org', 'outliers_node', 'dup_lineage', 'score1', 'support', 'species_losses_percentage',
                'rank', 'OGD_annot', 'sci_name', 'so_arq', 'dup_score', 'species_losses', 'node_is_og', 'lca_dup', 
                'so_score_dup', 'is_root', 'common_name', 'sp_out', 'so_bact', 'BEST_REP_NODE', 'create_mOG',
                'sp_loss_perc', 'num_lineage_losses', 'BEST_REP', 'score2', '_mems_og', 'so_euk', 'lineage', 'taxid', 
                'outliers_tax', 'dist', 'so_score', 'total_leaves', 'name', 'min_sp_overlap',
                'lca_node', 'len_sp_in', 'len_leaves_in', 'recover_seqs', 'seq_out_og', 'outliers_reftree', 
                'lca_node_name', 'best_tax', '2', '2157', '2759' , 'common_name', 'best_lost_1', 'best_lost_1']

# Trigger the interactive web server
t.explore(name = name, layouts = layouts, show_leaf_name = False , include_props = props_popup, keep_server=True, host = '138.4.138.141' ,open_browser =False)