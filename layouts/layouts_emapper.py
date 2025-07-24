from ete4.smartview  import  SeqMotifFace, TextFace


def get_emapper_pref_name():
    
    def layout_fn(node):    
        if node.is_leaf:
            face_pname = TextFace(node.props.get('Preferred_name'))
            #node.add_face(face_pname, column = 2, position = 'aligned')
        else:
          
            face_pname = TextFace(node.props.get('Preferred_name_counter').split('--')[0])
            #face_pname = TextFace(node.props.get('Preferred_name'))

        node.add_face(face_pname, column = 2, position = 'aligned', collapsed_only=True)    

    layout_fn.__name__ = 'emapper_pref_name'
    layout_fn.contains_aligned_face = True
    return layout_fn




def get_eggnog_OG():
    def layout_fn(node):
        if node.is_leaf:
            eggnog_str = node.props.get('eggNOG_OGs')
            if eggnog_str != None:
                eggnog_list = eggnog_str.split('||')
                for egg_og in eggnog_list:
                    og_name, taxid, = egg_og.split('|')[0].split('@')

                    if taxid in ['2759', '2', '2157']:
                        face_og = TextFace(og_name)
                        #node.add_face(face_og, column = 3, position = 'aligned')

                    
        else:
            eggnog_list = node.props.get('eggNOG_OGs_counter').split('||')
            for egg_og in eggnog_list:
                og_name, taxid, = egg_og.split('|')[0].split('@')

                if taxid in ['2759', '2', '2157']:
                    face_og = TextFace(og_name)

        node.add_face(face_og, column = 3, position = 'aligned',  collapsed_only=True)           
                
            
    layout_fn.__name__ = 'eggnog_OGs'
    layout_fn.contains_aligned_face = True
    return layout_fn



def get_emapper_kegg_ko():
   
    def layout_fn(node):   

        
     
        if node.is_leaf:
            
            face_pname = TextFace(node.props.get('KEGG_ko'))
            node.add_face(face_pname, column = 4, position = 'aligned')
        else:
            best_score = 0
            best_ko = set()
            for ko in (node.props.get('KEGG_ko_counter', str())).split('||'):
            
                if len(ko.split('--')) ==1:
                    continue
                ko_id , score = ko.split('--')
                if int(score) > best_score:
                    best_score = int(score)
                    best_ko = ko_id
            #face_pname = TextFace(node.props.get('Preferred_name_counter').split('--')[0])
            face_pname = TextFace(best_ko)
            node.add_face(face_pname, column = 4, position = 'aligned', collapsed_only=True)

    layout_fn.__name__ = 'emapper_kegg_ko'
    layout_fn.contains_aligned_face = True
    return layout_fn



def parse_pfam_doms(n):
    doms_string = n.props.get('dom_arq')
    
    try:
        doms = []
        for d in doms_string.split('||'):
            d_info = d.split('@')
            dom = [int(d_info[1]), int(d_info[2]), "()", None, None, 'grey', 'grey' ,"arial|20|black|%s" %(d_info[0])]
            doms.append(dom)
    except:
        doms = []
    return(doms)


def get_pfams():
    def layout_fn(node):
        if node.is_leaf:
            doms = parse_pfam_doms(node)
            try:
                aa = node.props.get('alignment')
                seqFace = SeqMotifFace(seq=aa, motifs = doms, width=500)
            except:
                aa = node.props.get('alignment')
                seqFace = SeqMotifFace(seq= aa, gapcolor='red')
            node.add_face(seqFace, column =  5, position = "aligned" )
        else:
            first_node = next(node.leaves())
            if first_node.name:
                doms = parse_pfam_doms(node)
                seqFace = SeqMotifFace(seq=None, motifs = doms, width=500)
                node.add_face(seqFace, column =  7, position = "aligned", collapsed_only=True )


    layout_fn.__name__ = 'PFAM_arq'
    layout_fn.contains_aligned_face = True
    return layout_fn