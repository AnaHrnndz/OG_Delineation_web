import sys
import subprocess
import re
import argparse
from collections import defaultdict
from ete4 import PhyloTree


####    FUNCTIONS FOR OG DELINEATION WEB    ####


def create_aln_dir(UPLOAD_FOLDER):
    try:
        path = UPLOAD_FOLDER+'/user_data'
        subprocess.run("rm -r %s" % path, shell = True)
        subprocess.run("mkdir -p %s" % path, shell = True)

    except OSError:
        print ("Creation of the directory user data folder failed")
    else:
        print ("Successfully created the directory user data folder ")


def parse_taxid(node):

    #TODO: add argument for split gen name
    return node.name.split('.')[0]

def clean_folder(dir, aln_path):

    """
        Remove files in tmp dir created for web analysis
    """

    filelist = glob.glob(os.path.join(dir, "*"))
    for f in filelist:

        if f != aln_path:
            os.remove(f)

# def load_annotated_tree(t, name_tree):
    # sp_set = set(t.species)
    # total_mems_in_tree = set(t.leaf_names())
    # NUM_TOTAL_SP = len(sp_set)


    # return t, sp_set, total_mems_in_tree, NUM_TOTAL_SP

def get_newick(t, all_props):

    """
        Return tree in newick format with annotations
    """

    t = t.write(props=all_props, format_root_node=True)
    return t

def get_props(t):

    all_props = set()
    for n in t.traverse():
        all_props.update(n.props.keys())

    return all_props

def get_taxlevel2ogs(t, taxid_dups_og, total_mems_in_ogs, taxonomy_db):

    """
        Table needed in web app
        For each taxlevel that create ogs, keep:
            ogs names, all mems in that ogs, ogs below that taxlevel, sci_name
    """

    taxlev2ogs = defaultdict(dict)


    for taxid in taxid_dups_og:
        taxlev2ogs[taxid]['ogs_names'] = set()
        taxlev2ogs[taxid]['mems'] = set()
        taxlev2ogs[taxid]['ogs_down_included'] = set()
        taxlev2ogs[taxid]["sci_name"] = taxonomy_db.get_taxid_translator([int(taxid)])[int(taxid)]

        for node in t.traverse('preorder'):
            #Keep ogs at that taxlevel
            if node.props.get('node_is_og') and taxid == node.props.get("lca_dup"):

                taxlev2ogs[taxid]['ogs_names'].add(node.name)

                #taxlev2ogs[taxid]['mems'].update(set(node.props.get('_leaves_in_nodes').split('|')))
                taxlev2ogs[taxid]['mems'].update(get_members(node, taxid).split('|'))

                if node.props.get('_ogs_down'):
                    taxlev2ogs[taxid]['ogs_down_included'].update(node.props.get('_ogs_down'))

            #Keep ogs below taxlevel
            elif node.props.get('node_create_og') and taxid != node.props.get("lca_node") and taxid in node.props.get("lineage"):

                candidates_nodes = node.search_nodes(node_is_og="True")
                save_candidates = list()
                for candidate in candidates_nodes:
                    if isinstance(candidate.props.get('lineage'), list) :
                        lin_cadidate = candidate.props.get('lineage')
                    else:
                        lin_cadidate = candidate.props.get('lineage').split('|')

                    if taxid in  lin_cadidate  and candidate.props.get('name') not in taxlev2ogs[taxid]['ogs_down_included']:
                        save_candidates.append(candidate.props.get('name'))

                if len(save_candidates) >0:
                    taxlev2ogs[taxid]['ogs_names'].add(node.props.get('name'))
                    #taxlev2ogs[taxid]['mems'].update(node.props.get('_leaves_in_nodes'))
                    taxlev2ogs[taxid]['mems'].update(get_members(node, taxid).split('|'))
                    taxlev2ogs[taxid]['ogs_down_included'].update(node.props.get('_ogs_down', '-'))



    return(taxlev2ogs)

def get_members(node, taxa):

    all_leafs = node.props.get('_leaves_in_nodes')

    mems_set = set()
    for l in all_leafs:
        if str(taxa) in l.props.get('lineage').split('|'):
            mems_set.add(l.name)

    mems = '|'.join(list(mems_set))

    return mems

def prune_tree(t, total_mems_in_ogs):

    """
        Prune tree to keep only mems that belong to at least one OG
        We need this version of the tree to visualization
    """

    t.prune(list(total_mems_in_ogs))
    return t

def run_clean_properties(t):

    """
        Clean problematic characters from some properties
        Call clean_string()
    """

    # clean strings
    all_props = set()

    #if tree came from web server is  str format,
    if isinstance(t, str):
        t = PhyloTree(t)

    for n in t.traverse():
        for string in ("sci_name", "lca_node_name", "common_name"):
            prop = n.props.get(string)
            if prop:
                n.props[string] = clean_string(prop)

        all_props.update(set(n.props.keys()))


    return t, all_props

def clean_string(string):


    """
        Remove problematic characters for newick format
    """

    clean_string = re.sub(r"'|\[|\]|\=|\.|\-|\:", "", string)
    return clean_string


def dict2args(args_dict):

    args = argparse.Namespace()

    args.user_taxonomy_counter = args_dict["user_counter"]
    args.outliers_node = args_dict["outliers_node"]
    args.outliers_reftree = args_dict["outliers_reftree"]
    args.inherit_out = args_dict["inherit_outliers"]
    args.so_bact = args_dict["so_bact"]
    args.so_euk = args_dict["so_euk"]
    args.so_arq = args_dict["so_arq"]
    args.so_cell_org = args_dict["so_cell_org"]
    args.sp_loss_perc = args_dict["sp_loss_perc"]
    return args