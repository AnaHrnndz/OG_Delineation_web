from ete4.smartview import NodeStyle
from ete4 import NCBITaxa
from ete4.smartview import TreeStyle
from ete4.smartview  import RectFace, CircleFace, SeqMotifFace, TextFace, OutlineFace
from ete4.treeview import random_color
from collections import Counter, OrderedDict, defaultdict
import json
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
lin2colors['Eukaryota'] = '#f200ff'
lin2colors['Archaea'] = '#ff7b00'
lin2colors['cellular organisms'] = "LightGrey"



# euk_colors = {}
# colors_gradient = list()


# if len(colors_gradient) == 0:
    # colors_gradient = get_gradient('euk', 50)





def get_level(node, level=0):
    if node.is_root():
        return level
    else:
        return get_level(node.up, level + 1)


def get_layout_leafname():
    "name: property used to obtain leaf name"

    def summary(nodes):
        "Return a list of names summarizing the given list of nodes"
        return list(OrderedDict((first_name(node), None) for node in nodes).keys())

    def first_name(tree):
        "Return the name of the first node that has a name"
        
        sci_names = []
        for node in tree.traverse('preorder'):
            if node.is_leaf():
                sci_name = node.props.get('sci_name')
                sci_names.append(sci_name)

        return next(iter(sci_names))


    def layout_fn(node):
        if node.is_leaf():
           
            sci_name = node.props.get('sci_name')
            name_seq = node.name.split('.',1)[1]

            node.add_face(TextFace(sci_name, color = 'black', padding_x=2),
                column=0, position="branch_right")

            node.add_face(TextFace(name_seq, color = 'grey', padding_x=2),
                column=1, position="branch_right")

            

        else:
            # Collapsed face
            if not node.is_root() and  node.props.get('collapsed', 'false') == 'true': 
                text = node.props.get('lca_node_name')
                node.add_face(TextFace(text, padding_x=2),
                    position="branch_right", column=1, collapsed_only=True)
                len_leaves = str(len(set(node.props.get('_leaves_in').split('|'))))
                print(node.props.get('_leaves_in'))
                node.add_face(TextFace(len_leaves, padding_x=2),
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
            lca = node.props.get('lca_node_name')
            #color = node.props.get('_Lca_color')
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

           

            level = get_level(node)
            lca_face = RectFace(15, float('inf'), 
                    color = color , 
                    text = lca,
                    fgcolor = "white",
                    padding_x = 1, padding_y = 1)
            lca_face.rotate_text = True
            node.add_face(lca_face, position='aligned', column=level)
            node.add_face(lca_face, position='aligned', column=level,
                collapsed_only=True)

    

    layout_fn.__name__ = 'Last common ancestor'
    layout_fn.contains_aligned_face = True
    return layout_fn


def get_layout_evoltype():
    def layout_fn(node):

        ns = NodeStyle()

        if node.props.get('lineage') and not node.is_leaf():
            
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
                
                node.sm_style["fgcolor"] = color

            
                

    layout_fn.__name__ = 'Evolution events'
    return layout_fn


def collapse_og():
    def layout_fn(node):
        if not node.is_root() and  node.props.get('collapsed', 'false') == 'true':            
            
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