import sys
from ete4 import Tree, SeqGroup
from ete4.smartview import TreeStyle, NodeStyle, TreeLayout
from ete4.smartview.renderer.layouts.seq_layouts import LayoutAlignment

sys.path.append('/data/projects/og_delineation_web/layouts')
import layout_possvm
import layout_web
# from layout_web import get_layout_leafname, get_layout_evoltype, get_layout_lca_rects, collapse_og, seqs_out_og, collapse_rank, get_species_overlap, get_pnames, get_ogs, get_doms, prune_tree #sequencebouncer

# from layout_possvm import get_pnames_possvm, get_og_euk, get_og_met, get_og_possvm, get_og_ref, get_pfams, background_color_ogs, recover_seqs


t = Tree('((9606.A:1, 9598.B:1)0.5:1,10090.C:1);', parser = 0)



if len(sys.argv) >1:
    t = Tree(open(sys.argv[1]))


#raw_alg = SeqGroup(sys.argv[2])



layouts = [
    TreeLayout('Scientific name', ns = layout_web.get_layout_leafname(), aligned_faces = True),
    TreeLayout('Evolution events',  ns = layout_web.get_layout_evoltype()),
    TreeLayout('Last common ancestor',  ns = layout_web.get_layout_lca_rects(t),  aligned_faces = True),
    TreeLayout('Collapse OG',  ns = layout_web.collapse_og(),  aligned_faces = True),
    TreeLayout('Seqs out OGS', ns = layout_web.seqs_out_og()),
    TreeLayout('Collapse order rank', ns = layout_web.collapse_rank()),
    TreeLayout('Species Overlap', ns = layout_web.get_species_overlap()),
    TreeLayout('Pref name', ns = layout_web.get_pnames(), aligned_faces = True),
    TreeLayout('EGGNOG', ns = layout_web.get_ogs(), aligned_faces = True),
    TreeLayout('PFAM Domains', ns = layout_web.get_doms(), aligned_faces = True ),
    TreeLayout('Long branches', ns = layout_web.long_branch_outlier()),
    TreeLayout('Taxonomic Outliers', ns = layout_web.taxo_out()),
    TreeLayout('Brackground OG', ns = layout_web.background_og()),
    TreeLayout('Brackground mOG', ns = layout_web.background_mog(),aligned_faces = True )
    #TreeLayout('background_mog_int', ns = layout_web.background_mog_int(), aligned_faces = True)
    #TreeLayout('Prune tree', ns = layout_web.prune_tree()),

    #TreeLayout('SequenceBouncer', ns=sequencebouncer())


    # TreeLayout('possvm_pname', ns = layout_possvm.get_pnames_possvm(), aligned_faces = True),
    # TreeLayout('possvm_OG_Euk', ns = layout_possvm.get_og_euk(), aligned_faces = True),
    # TreeLayout('possm_OG_Metazoa', ns = layout_possvm.get_og_met(), aligned_faces = True),
    # TreeLayout('possvm_OG', ns = layout_possvm.get_og_possvm(), aligned_faces = True),
    # TreeLayout('possvm_RefOG', ns = layout_possvm.get_og_ref(), aligned_faces = True),

    # TreeLayout('possvm_PFAM', ns = get_pfams(), aligned_faces = True),
    # TreeLayout('background_refogs', ns = background_color_ogs()),
    # TreeLayout('recover', ns = recover_seqs())
    #LayoutAlignment(sys.argv[2], width=500)
        ]



for l in layouts:
    if l.name in ['Seqs out OGS', 'Alignment', 'EGGNOG', 'PFAM Domains', 'Long branches','possvm_PFAM', 'Prune tree', 'Collapse order rank', 'possm_OG_Metazoa', 'Taxonomic Outliers']:
        l.active = False



# Launch smartview explorer
t.dist = 0.05

name = 'upload_tree'

props_popup = ['taxo_outlier', 'mOG','node_create_og', 'inparalogs_rate', 'BEST_RESP_NUM', 'evoltype_2', 'so_cell_org', 'outliers_node', 'dup_lineage', 'score1', 'support', 'species_losses_percentage',
'rank', 'OGD_annot', 'sci_name', 'so_arq', 'dup_score', 'species_losses', 'node_is_og', 'lca_dup', 'so_score_dup', 'is_root', 'common_name', 'sp_out', 'so_bact', 'BEST_REP_NODE', 'create_mOG',
'sp_loss_perc', 'num_lineage_losses', 'BEST_REP', 'score2', '_mems_og', 'so_euk', 'lineage', 'taxid', 'outliers_tax', 'dist', 'so_score', 'total_leaves', 'name', 'min_sp_overlap',
'lca_node', 'len_sp_in', 'len_leaves_in', 'recover_seqs', 'seq_out_og', 'outliers_reftree', 'lca_node_name', 'best_tax', '2', '2157', '2759' , 'common_name', 'best_lost_1', 'best_lost_1']

# Trigger the interactive web server
t.explore(name = name, layouts = layouts, show_leaf_name = False , include_props = props_popup, keep_server=True, host = '138.4.138.141' ,open_browser =False)
