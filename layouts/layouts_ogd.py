from ete4.smartview  import RectFace, TextFace
from ete4.utils import random_color
from collections import  OrderedDict, defaultdict
import random


hues = { "green": (81 / 360, 240 / 360) ,
         "purple": (240/360, 345 / 360),
         "red": (0 / 360, 10 / 360),
         "orange": (10 / 360, 80 / 360),
        }


def get_gradient(color, num, s=None, l=None):
    h0, h1 = hues[color]
    sep = (h1 - h0) / num
    return random_color(h=h0, s=s, l=l, num=num, sep=sep)

colors_bact = get_gradient("green", 50)
colors_euk = get_gradient("purple", 50)
colors_arq = get_gradient("orange", 50)

lin2colors = defaultdict()
lin2colors['Bacteria'] = '#53bd42'
lin2colors['2'] = '#53bd42'
lin2colors['Eukaryota'] = '#f200ff'
lin2colors['2759'] = '#f200ff'
lin2colors['Archaea'] = '#ff7b00'
lin2colors['2157'] = '#ff7b00'
lin2colors['cellular organisms'] = "Grey"
lin2colors['131567'] = "Grey"
lin2colors['root'] = "Red"



def get_level(node, level=0):
    if node.is_root:
        return level
    else:
        return get_level(node.up, level +1)



def get_layout_leafname():
    "name: property used to obtain leaf name"

    def summary(nodes):
        "Return a list of names summarizing the given list of nodes"
        return list(OrderedDict((first_name(node), None) for node in nodes).keys())

    def first_name(tree):
        "Return the name of the first node that has a name"

        sci_names = []
        for node in tree.traverse('preorder'):
            if node.is_leaf:
                sci_name = node.props.get('sci_name')
                sci_names.append(sci_name)

        return next(iter(sci_names))


    def layout_fn(node):
        if node.is_leaf:

            sci_name = node.props.get('sci_name')

            color = 'black'
            name_seq = node.name.split('.',1)[1]

            node.add_face(TextFace(sci_name, color = color, padding_x=2),
                column=0, position="branch_right")

            node.add_face(TextFace(name_seq, color = 'grey', padding_x=2),
                column=1, position="branch_right")


        else:
            # Collapsed face
            if node.props.get('rank') in  ['order', 'genus', 'family']:
                text = node.props.get('lca_node_name')
                node.add_face(TextFace(text, padding_x=2),
                    position="branch_right", column=1, collapsed_only=True)


            elif not node.is_root and  node.props.get('collapsed', 'false') == 'true':
                text = node.props.get('lca_node_name')
                node.add_face(TextFace(text, padding_x=2),
                    position="branch_right", column=1, collapsed_only=True)
                len_leaves = str(len(set(node.props.get('_leaves_in').split('|'))))

                node.add_face(TextFace(len_leaves, padding_x=2),
                    position="branch_right", column=1, collapsed_only=True)

                og_egg = node.props.get('og', '-')
                node.add_face(TextFace(og_egg, padding_x=2),
                    position="branch_right", column=1, collapsed_only=True)

                pname_emapper = node.props.get('pname', '-')
                node.add_face(TextFace(pname_emapper, padding_x=2),
                    position="branch_right", column=1, collapsed_only=True)

            else:
                names = summary(node.children)
                texts = names if len(names) < 6 else (names[:3] + ['...'] + names[-2:])
                for i, text in enumerate(texts):
                    node.add_face(TextFace(text, padding_x=2),
                            position="branch_right", column=1, collapsed_only=True)

    layout_fn.__name__ = 'Scientific name'
    layout_fn.contains_aligned_face = True
    return layout_fn



def get_layout_lca_rects(tree):

    def layout_fn(node):

        if node.props.get('lca_node_name'):
            color = 'Grey'
            lca = node.props.get('lca_node_name')
            
            if lca in lin2colors.keys():
                color = lin2colors[lca]
            else:
                if '2' in node.props.get('lineage').split('|'):
                    color = random.choice(colors_bact)
                    lin2colors[lca] = color
                elif '2759' in node.props.get('lineage').split('|'):
                    color = random.choice(colors_euk)
                    lin2colors[lca] = color
                elif '2157' in node.props.get('lineage').split('|'):
                    color = random.choice(colors_arq)
                    lin2colors[lca] = color

            level = get_level(node)+7
            lca_face = RectFace(15, float('inf'),
                    color = color ,
                    text = lca,
                    fgcolor = "white",
                    padding_x = 1, padding_y = 1)
            lca_face.rotate_text = True
            node.add_face(lca_face, position='aligned', column=level)
            node.add_face(lca_face, position='aligned', column=level, collapsed_only=True)

    layout_fn.__name__ = 'Last common ancestor'
    layout_fn.contains_aligned_face = True
    return layout_fn



def get_layout_evoltype():
    def layout_fn(node):

        if node.props.get('lineage') and not node.is_leaf:

            if node.props.get('evoltype_2') == 'S':
                node.sm_style["fgcolor"] = 'blue'
                node.sm_style["size"] = 2

            elif node.props.get('evoltype_2') == 'D':
                node.sm_style["fgcolor"] = 'red'
                node.sm_style["size"] = 2

            elif node.props.get('evoltype_2') == 'FD':
                node.sm_style["fgcolor"] = 'Coral'
                node.sm_style["size"] = 2

            if node.props.get('node_is_og'):
                node.sm_style['size'] = 5
                node.sm_style['shape'] = 'square'

                lca = node.props.get('lca_dup')
                if lca in lin2colors.keys():
                    color = lin2colors[lca]
                else:

                    if '2' in node.props.get('lineage').split('|'):
                        color = random.choice(colors_bact)
                        lin2colors[lca] = color
                    elif '2759' in node.props.get('lineage').split('|'):
                        color = random.choice(colors_euk)

                        lin2colors[lca] = color
                    elif '2157' in node.props.get('lineage').split('|'):
                        color = random.choice(colors_arq)
                        lin2colors[lca] = color
                    else:
                        color='grey'

                node.sm_style["fgcolor"] = color


    layout_fn.__name__ = 'Evolution events'
    return layout_fn


#Esta no se que hace
def collapse_og():
    def layout_fn(node):
        if not node.is_root and  node.props.get('collapsed', 'false') == 'true':

            node.sm_style["draw_descendants"] = False

            lca = node.props.get('lca_node_name')
            if lca in lin2colors.keys():
                color = lin2colors[lca]
            else:

                if '2' in node.props.get('lineage').split('|'):
                    color = random.choice(colors_bact)
                    lin2colors[lca] = color
                elif '2759' in node.props.get('lineage').split('|'):
                    color = random.choice(colors_euk)

                    lin2colors[lca] = color
                elif '2157' in node.props.get('lineage').split('|'):
                    color = random.choice(colors_arq)
                    lin2colors[lca] = color

            node.sm_style["outline_color"] = color


    layout_fn.__name__ = "Collapse OG"
    return layout_fn


def background_og():
    def layout_fn(node):
        if node.props.get('node_is_og'):
            lca = node.props.get('lca_node_name')
            
            if lca in lin2colors.keys():
                color = lin2colors[lca]

                node.sm_style["bgcolor"] = color


    layout_fn.__name__ = "Brackground OG"
    return layout_fn

def background_mog():
    def layout_fn(node):

        if node.is_leaf and node.props.get('mOG'):
            face = TextFace(node.props.get('mOG'))
            node.sm_style["bgcolor"] = 'blue'
            node.add_face(face, column = 2, position = 'aligned')

    layout_fn.__name__ = "Brackground mOG"
    layout_fn.contains_aligned_face = True
    return layout_fn



def seqs_out_og():
    def layout_fn(node):


        for l in node.search_nodes(seq_out_og= "true"):

            l.sm_style['fgcolor'] = 'red'
            l.sm_style['bgcolor'] = 'red'
            for anc in l.ancestors():
                anc.sm_style['hz_line_color'] = 'red'
                anc.sm_style['vt_line_color'] = 'red'

        if 'seq_out_og' in node.props.keys():
            node.sm_style['size'] = 5


    layout_fn.__name__ = "Seqs out OGS"
    return layout_fn


def taxo_out():
    def layout_fn(node):

        for l in node.search_nodes(taxo_outlier= "true"):

            l.sm_style['fgcolor'] = 'orange'
            l.sm_style['bgcolor'] = 'orange'
            for anc in l.ancestors():
                anc.sm_style['hz_line_color'] = 'orange'
                anc.sm_style['vt_line_color'] = 'orange'

        if 'taxo_outlier' in node.props.keys():
            node.sm_style['size'] = 5


    layout_fn.__name__ = "Taxonomic Outliers"
    return layout_fn


def long_branch_outlier():
    def layout_fn(node):

        for l in node.search_nodes(long_branch_outlier='True'):

            l.sm_style['fgcolor'] = "orange"
            l.sm_style['bgcolor'] = "orange"

            for anc in l.ancestors():
                anc.sm_style['hz_line_color'] = 'orange'
                anc.sm_style['vt_line_color'] = 'orange'


    layout_fn.__name__ = "Long Branches"
    return layout_fn



#esta falla
def collapse_rank():
    def layout_fn(node):
        if node.props.get('rank') in  ['order', 'genus', 'family'] :
            node.sm_style["draw_descendants"] = False

            og_nodes = node.search_nodes(node_is_og="True")
            if len(og_nodes) > 0:
                node.sm_style['outline_color'] = "Dimgray"


    layout_fn.__name__ = "Collapse order rank"
    return layout_fn






def get_species_overlap():
    def layout_fn(node):
        if node.props.get('so_score', '0.0'):
            so = round(float(node.props.get('so_score', '0.0')),3)
            so_face = TextFace(str(so), color = "green")
            node.add_face(so_face, position = "branch_bottom",  column = 0)

    layout_fn.__name__ = "Species Overlap"
    return layout_fn





# def get_ogs():
    # def layout_fn(node):
        # if node.is_leaf and node.props.get('basal_og'):
            # og_pname = TextFace(node.props.get('basal_og'))
            # node.add_face(og_pname, column = 3, position = 'aligned')
        # if node.props.get('basal_og'):
            # face_og = TextFace(node.props.get('basal_og'))
            # node.add_face(face_og, column = 3, position = 'aligned', collapsed_only=True)

    # layout_fn.__name__ = 'EGGNOG'
    # layout_fn.contains_aligned_face = True
    # return layout_fn






