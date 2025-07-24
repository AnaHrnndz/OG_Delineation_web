import sys
from ete4 import PhyloTree

sys.path.append('/data/projects/og_delineation/')
import og_delineation
sys.path.append('/data/projects/og_delineation/ogd')
import tree_setup
sys.path.append('/data/projects/og_delineation_web/ogd_web')
import util
from collections import defaultdict



def run_upload(tree, name_tree, reftree, user_counter, user_taxo, taxonomy_type, rooting, UPLOAD_FOLDER):

    """
        If tree is already annotated with OGD, just load  info from tree
        if tree is not annotated,:
            1. Load or create:
                taxonomy db
                species counter
                reference species tree
            2. Run tree_setup from og_delineation
    """

    t = PhyloTree(tree, parser = 0)

    # If tree is annotated with OGD
    if "OGD_annot" in t.props:

        print('Upload annotated OGD-tree.....')
        t_nw, general_results,  taxo_stats,  parameters, SPTOTAL, t_prune_nw, total_mems_in_tree = load_annot_tree(t)
        taxonomy_db=None
        taxonomy_counter=None
        return t_nw, general_results,  taxo_stats,  parameters, SPTOTAL, taxonomy_db, taxonomy_counter, total_mems_in_tree

    # Run tree_setup from og_delineation
    else:
        print('Upload raw tree.....')

        t.set_species_naming_function(util.parse_taxid)
        total_mems_in_tree = t.leaf_names()
        taxonomy_db = og_delineation.load_taxonomy(taxonomy = taxonomy_type, user_taxonomy= user_taxo)
        reftree = og_delineation.load_reftree(rtree= reftree, t = t, taxonomy_db = taxonomy_db)
        taxonomy_counter = og_delineation.load_taxonomy_counter(reftree=reftree, user_taxonomy_counter = user_counter)

        t.write(outfile=UPLOAD_FOLDER+'/user_data/ori_tree.nw', format_root_node = True)
        path2tree = UPLOAD_FOLDER+'/user_data/ori_tree.nw'
        path_out = UPLOAD_FOLDER+'/user_data/'

        # 2. Tree setup (Pre-analysis):  resolve polytomies, rooting, ncbi annotation, etc
        t , sp_set, total_mems_in_tree, SPTOTAL, user_props = tree_setup.run_setup(t, name_tree, taxonomy_db, rooting, path_out, path2tree, '.')


        # Clean properties
        t, all_props = util.run_clean_properties(t)
        t_nw = util.get_newick(t, all_props)


        general_results = defaultdict()
        taxo_stats = defaultdict()
        parameters = defaultdict()

        return t_nw, general_results,  taxo_stats,  parameters, SPTOTAL, taxonomy_db, taxonomy_counter, total_mems_in_tree




def load_annot_tree(t):

    t.set_species_naming_function(util.parse_taxid)
    SPTOTAL = set(t.species)

    total_mems_in_tree = set(t.leaf_names())

    parameters = defaultdict()
    parameters_info = t.props.get('parameters').split('|')
    for p in parameters_info:
        p_name, p_value = p.split('@')
        parameters[p_name] = p_value

    general_results = defaultdict()
    general_result_info = t.props.get('general_result').split('|')

    for r in general_result_info:
        r_name, r_value = r.split('@')
        general_results[r_name] = r_value


    stats_taxo = defaultdict(dict)
    taxo_stats_info = t.props.get('taxlev2ogs').split('@')

    if len(taxo_stats_info) == 1:
        sci_name_taxid = t.props.get('lca_node_name')
        stats_taxo[sci_name_taxid]["num_ogs"] = 1
        stats_taxo[sci_name_taxid]["num_mems"] = len(t)

    elif len(taxo_stats_info) > 1:
        for taxlev in taxo_stats_info:
            taxlev_info = taxlev.split('|')
            sci_name_taxid = taxlev_info[0]
            ogs = taxlev_info[1].split('_')
            num_mems = taxlev_info[2]

            stats_taxo[sci_name_taxid]["num_ogs"] = len(ogs)
            stats_taxo[sci_name_taxid]["num_mems"] = num_mems

    if len(taxo_stats_info ) >1:
        for taxo in taxo_stats_info:
            info = taxo.split('|')
            taxid = int(info[0].split('_')[-1])
            ogs = info[1].split('_')


    #Prune tree, First I need to copy the tree
    dup_tree = t.copy("deepcopy")

    seqs_out = set()
    for l in t.search_nodes(seq_out_og='true'):
        seqs_out.add(l.name)

    total_mems_in_ogs = total_mems_in_tree.difference(seqs_out)
    t_prune= util.prune_tree(dup_tree, total_mems_in_ogs)

    all_props = util.get_props(t)

    t_nw = util.get_newick(t, all_props)
    t_prune_nw = util.get_newick(t_prune, all_props)


    return  t_nw, general_results,  stats_taxo,  parameters, SPTOTAL,  t_prune_nw, total_mems_in_tree


