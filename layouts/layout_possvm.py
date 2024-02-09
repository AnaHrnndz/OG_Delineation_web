from ete4.smartview import NodeStyle
from ete4 import NCBITaxa
from ete4.smartview import TreeStyle
from ete4.smartview  import RectFace, CircleFace, SeqMotifFace, TextFace, OutlineFace
from ete4.utils import random_color
from collections import Counter, OrderedDict, defaultdict
import json
import random


#fpI = json.load(open('/data/projects/og_delineation/benchmark/possvm/scripts/prds_fpI.json'))
# refog_sp = set()
# with open('/data/projects/og_delineation/benchmark/possvm/original_data/refog_sp_sci_name.txt') as f:
    # for line in f:
        # info = line.strip().split('\t')
        # refog_sp.add(info[1])


# best_myogs = defaultdict()
# with open('/data/projects/og_delineation/benchmark/possvm/results_ogd/species_losses_08/prds/best_ogd_fscore.tsv') as f:
    # for line in f:
        # info = line.strip().split('\t')
        # best_myogs[info[1]] = info[0]



# recovery_info = json.load(open('/data/projects/og_delineation/benchmark/possvm/results_ogd/species_losses_08/prds/aln_hmm/recovery_info.json'))

# print(recovery_info)


def recover_seqs():
    def layout_fn(node):
        if node.is_leaf:
            if node.name in recovery_info.keys():
                for k,val in recovery_info[node.name].items():
                    og_name = k.split('_')[0]
                    if og_name in best_myogs.keys():
                        recover_face = TextFace(best_myogs[og_name])
                    else:
                        recover_face = TextFace(og_name)


                node.add_face(recover_face, column = 7, position = 'aligned' )
    layout_fn.__name__ = 'recover'
    layout_fn.contains_aligned_face = True
    return layout_fn




def background_color_ogs():
    def layout_fn(node):
        if node.name in best_myogs.keys():

            node.sm_style["bgcolor"] = '#2596be'

    layout_fn.__name__ = 'background_refogs'
    layout_fn.contains_aligned_face = True
    return layout_fn

def get_pnames_possvm():
    def layout_fn(node):
        if node.is_leaf and node.props.get('Pname'):
            face_pname = TextFace(node.props.get('Pname'))
            node.add_face(face_pname, column = 2, position = 'aligned')
        elif node.props.get('Pname'):
            face_pname = TextFace(node.props.get('Pname'))
            node.add_face(face_pname, column = 2, position = 'aligned', collapsed_only=True)

    layout_fn.__name__ = 'possvm_pname'
    layout_fn.contains_aligned_face = True
    return layout_fn


def get_og_euk():
    def layout_fn(node):
        if node.is_leaf and node.props.get('OG_Euk'):
            face_pname = TextFace(node.props.get('OG_Euk'))
            node.add_face(face_pname, column = 3, position = 'aligned')
        elif node.props.get('OG_Euk'):
            face_pname = TextFace(node.props.get('OG_Euk'))
            node.add_face(face_pname, column = 3, position = 'aligned', collapsed_only=True)

    layout_fn.__name__ = 'possvm_OG_Euk'
    layout_fn.contains_aligned_face = True
    return layout_fn


def get_og_met():
    def layout_fn(node):
        if node.is_leaf and node.props.get('OG_Metazoa'):
            face_pname = TextFace(node.props.get('OG_Metazoa'))
            node.add_face(face_pname, column = 4, position = 'aligned')
        elif node.props.get('OG_Metazoa'):
            face_pname = TextFace(node.props.get('OG_Metazoa'))
            node.add_face(face_pname, column = 4, position = 'aligned', collapsed_only=True)

    layout_fn.__name__ = 'possvm_OG_Metazoa'
    layout_fn.contains_aligned_face = True
    return layout_fn


def get_og_possvm():
    def layout_fn(node):


        # if node.is_leaf and node.props.get('Possvm_OG'):
            # color = node.props.get('Possvm_OG').split('_')[1]
            # face_pname = TextFace(node.props.get('Possvm_OG').split('_')[0], color = color)
            # node.add_face(face_pname, column = 5, position = 'aligned')
        # elif node.props.get('Possvm_OG'):
            # color = node.props.get('Possvm_OG').split('_')[1]
            # face_pname = TextFace(node.props.get('Possvm_OG').split('_')[0], color = color)
            # node.add_face(face_pname, column = 5, position = 'aligned', collapsed_only=True)

        # if node.props.get('Possvm_OG') in fpI.keys():
            # color = 'red'
        # else:
            # color = 'black'
        if node.is_leaf and node.props.get('Possvm_OG'):

            face_pname = TextFace(node.props.get('Possvm_OG'))
            node.add_face(face_pname, column = 5, position = 'aligned')
        elif node.props.get('Possvm_OG'):

            face_pname = TextFace(node.props.get('Possvm_OG'))
            node.add_face(face_pname, column = 5, position = 'aligned', collapsed_only=True)


    layout_fn.__name__ = 'possvm_OG'
    layout_fn.contains_aligned_face = True
    return layout_fn


def get_og_ref():
    def layout_fn(node):
        if node.is_leaf and node.props.get('RefOG'):
            face_pname = TextFace(node.props.get('RefOG'))
            node.add_face(face_pname, column = 6, position = 'aligned')
        elif node.props.get('RefOG'):
            face_pname = TextFace(node.props.get('RefOG'))
            node.add_face(face_pname, column = 6, position = 'aligned', collapsed_only=True)

    layout_fn.__name__ = 'possvm_RefOG'
    layout_fn.contains_aligned_face = True
    return layout_fn


def parse_pfam_doms(n):

    doms_string = n.props.get('Pfam_Doms')
    doms = []
    for d in doms_string.split('|'):
        d_info = d.split('@')
        dom = [int(d_info[1]), int(d_info[2]), "()", None, None, 'grey', 'grey' ,"arial|20|black|%s" %(d_info[0])]
        doms.append(dom)
    return(doms)


def get_pfams():
    def layout_fn(node):
        if node.is_leaf:
            doms = parse_pfam_doms(node)
            seqFace = SeqMotifFace(seq=None, motifs = doms)
            node.add_face(seqFace, column =  1, position = "aligned")
        else:
            first_node = next(node.iter_leaves())
            if first_node.name:
                doms = parse_pfam_doms(node)
                seqFace = SeqMotifFace(seq=None, motifs = doms)
                node.add_face(seqFace, column =  1, position = "aligned", collapsed_only=True)


    layout_fn.__name__ = 'possvm_PFAM'
    layout_fn.contains_aligned_face = True
    return layout_fn