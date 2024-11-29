from ete4.smartview  import  SeqMotifFace, TextFace


def get_emapper_pref_name():
    def layout_fn(node):
        face_pname = TextFace(node.props.get('Preferred_name'))
        
        if node.is_leaf:
            node.add_face(face_pname, column = 2, position = 'aligned')
        else:
            node.add_face(face_pname, column = 2, position = 'aligned', collapsed_only=True)

    layout_fn.__name__ = 'emapper_pref_name'
    layout_fn.contains_aligned_face = True
    return layout_fn


def get_eggnog_OG():
    def layout_fn(node):
        if node.is_leaf:
            eggnog_list = node.props.get('eggNOG_OGs').split('||')
            for egg_og in eggnog_list:
                og_name, taxid, = egg_og.split('|')[0].split('@')
                
                if taxid == '2759':
                    og_euk = og_name
                    face_og_euk = TextFace(og_euk)
                    node.add_face(face_og_euk, column = 3, position = 'aligned')
                elif taxid =='2':
                    og_bact = og_name
                    face_og_euk = TextFace(og_bact)
                    node.add_face(face_og_euk, column = 4, position = 'aligned')

                elif taxid =='2157':
                    og_arq = og_name
                    face_og_euk = TextFace(og_arq)
                    node.add_face(face_og_euk, column = 5, position = 'aligned')
            
        else:
            eggnog_list = node.props.get('eggNOG_OGs_counter').split('||')
            for egg_og in eggnog_list:
                og_name, taxid, = egg_og.split('|')[0].split('@')
                
                if taxid == '2759':
                    og_euk = og_name
                    face_og_euk = TextFace(og_euk)
                    node.add_face(face_og_euk, column = 3, position = 'aligned', collapsed_only=True)
                elif taxid =='2':
                    og_bact = og_name
                    face_og_euk = TextFace(og_bact)
                    node.add_face(face_og_euk, column = 4, position = 'aligned', collapsed_only=True)

                elif taxid =='2157':
                    og_arq = og_name
                    face_og_euk = TextFace(og_arq)
                    node.add_face(face_og_euk, column = 5, position = 'aligned', collapsed_only=True)
            
    layout_fn.__name__ = 'eggnog_OGs'
    layout_fn.contains_aligned_face = True
    return layout_fn



def parse_pfam_doms(n):
    doms_string = n.props.get('dom_arq')
    doms = []
    for d in doms_string.split('||'):
        d_info = d.split('@')
        dom = [int(d_info[1]), int(d_info[2]), "()", None, None, 'grey', 'grey' ,"arial|20|black|%s" %(d_info[0])]
        doms.append(dom)
    return(doms)


def get_pfams():
    def layout_fn(node):
        if node.is_leaf:
            doms = parse_pfam_doms(node)
            seqFace = SeqMotifFace(seq=None, motifs = doms, width=500)
            node.add_face(seqFace, column =  1, position = "aligned" )
        else:
            first_node = next(node.leaves())
            if first_node.name:
                doms = parse_pfam_doms(node)
                seqFace = SeqMotifFace(seq=None, motifs = doms, width=500)
                node.add_face(seqFace, column =  1, position = "aligned", collapsed_only=True )


    layout_fn.__name__ = 'PFAM_arq'
    layout_fn.contains_aligned_face = True
    return layout_fn