from flask import Flask, render_template, request, json, Response, jsonify, send_file
import requests
from werkzeug.utils import secure_filename
import os
import sys
from collections import defaultdict
from ete4 import Tree, PhyloTree, SeqGroup
from ete4.smartview import TreeStyle, NodeStyle, TreeLayout
import csv
import subprocess
from io import StringIO

sys.path.append('/data/projects/og_delineation')
from  og_delineation import run_load_tree, run_load_taxonomy, run_load_reftree, run_load_taxonomy_counter, run_preanalysis_annot_tree, get_og_info, taxlev2ogs_annotated_tree#, run_load_annotated_tree
from  og_delineation import run_preanalysis, run_outliers_dup_score, run_dups_and_ogs, run_clean_properties, get_newick,  get_taxlevel2ogs, add_ogs_up_down, get_analysis_parameters, flag_seqs_out_og, prune_tree
from  og_delineation import run_write_fastas, run_create_hmm_og, run_hmmscan, get_best_match, expand_hmm, update_taxlevel2ogs, update_og_info, update_tree, clean_tree, clean_folder, get_seq2og, write_seq2ogs
# sys.path.append('/data/projects/og_delineation/egg5_pfamA/data')
# from add_pnames_ogs import add_pnames_ogs, add_possvm_info, add_ranger_info

sys.path.append('/data/projects/og_delineation_web/bin')


UPLOAD_FOLDER = '/data/projects/og_delineation_web/tmp'
url = 'http://138.4.138.141:5000/trees'

app = Flask(__name__)
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER


current_data = defaultdict()
general_results = defaultdict()
USER_PROPS = list()

glob_og_info = defaultdict(dict)
glob_og_info_updated = defaultdict(dict)
glob_stats_taxo = defaultdict(dict)
glob_taxlev2ogs = defaultdict(dict)
glob_taxlev2ogs_updated = defaultdict(dict)
glob_content = defaultdict()


#UPLOAD TREE WITH OR WITHOUT ANNOTATATIONS
def upload_annotated_tree(tree, name_tree, reftree, taxonomy_counter, user_taxo, taxonomy_type, taxonomy_db):
    
    print('uploading tree.....')
    t , sp_set, total_mems_in_tree, SPTOTAL, CONTENT, taxid_dups, names_ogs = run_preanalysis_annot_tree(tree, name_tree)
    print('tree upload')
    
    current_data["tree_name"] = name_tree
    current_data["reftree"] = reftree
    current_data["total_species"] = SPTOTAL
    current_data["tree_cache"] = CONTENT
    current_data["taxo_counter"] = taxonomy_counter
    current_data["taxo_db"] = taxonomy_db
    current_data["taxo_type"] = taxonomy_type
    current_data["taxo_user"] = user_taxo
    current_data["total_mems"] = total_mems_in_tree
    
    
    # parameters = get_analysis_parameters(t)
    
    parameters = defaultdict()
    parameters_info = t.props.get('parameters').split('|')
    for p in parameters_info:
        p_name, p_value = p.split('@')
        parameters[p_name] = p_value

    current_data["parameters"] = parameters
    
    
    
    general_result_info = t.props.get('general_result').split('|')
    for r in general_result_info:
        r_name, r_value = r.split('@')
        general_results[r_name] = r_value

    
    # print(general_result_2)

    # general_results["Tree name"] = name_file
    # general_results["Total seqs"] = len(total_mems_in_tree)
    # general_results["Total_species"] = SPTOTAL
    # general_results["Seqs in OGs"] = len(total_mems_in_ogs)
    # general_results["Seqs out OGS"] = len(total_mems_in_tree)-len(total_mems_in_ogs)
    # general_results["Num Ogs"] = len(names_ogs)


    stats_taxo = defaultdict(dict)
    taxo_stats_info = t.props.get('taxlev2ogs').split('@')
    for taxlev in taxo_stats_info:
        taxlev_info = taxlev.split('|')
        sci_name_taxid = taxlev_info[0]
        ogs = taxlev_info[1].split('_')
        num_mems = taxlev_info[2]

        stats_taxo[sci_name_taxid]["num_ogs"] = len(ogs)
        stats_taxo[sci_name_taxid]["num_mems"] = num_mems


    
    # print('loading og info ....')
    # #ogs_info, taxid_dups_og, total_mems_in_ogs = get_og_info(t)

    
    

    # print('Loading taxlevel ogs')
    # taxlev2ogs, total_mems_in_ogs=  taxlev2ogs_annotated_tree(t, taxid_dups, taxonomy_db)
    
   
    #current_data["mems_in_ogs"] = total_mems_in_ogs
    
    #Clean properties
    t, all_props = run_clean_properties(t)

    t = get_newick(t, all_props)

    current_data['properties'] = all_props
    current_data["tree"] = t

    # glob_og_info.clear()
    # for og_name, og_info in ogs_info.items():
        # glob_og_info[og_name] = og_info

    # glob_taxlev2ogs.clear()
    # for taxid, info in taxlev2ogs.items():
        # glob_taxlev2ogs[int(taxid)] = info
    
    
    


    #stats_taxo = defaultdict(dict)
    
    # for tax, info in taxlev2ogs.items():
        # sci_name = info["sci_name"]
        # name = str(sci_name)+'_'+str(tax)
        # stats_taxo[name]["num_ogs"] = len(info["ogs_names"])
        # stats_taxo[name]["num_mems"]= len(info["mems"])
    
    glob_stats_taxo.clear()
    for k,v in stats_taxo.items():
        glob_stats_taxo[k] = v
        

    
    
    return  t, general_results,  stats_taxo,  parameters


def run_upload(tree, name_tree, reftree, user_counter, user_taxo, taxonomy_type, midpoint):
    #Load files and DBs:
    t = run_load_tree(tree= tree)
    
    taxonomy_db = run_load_taxonomy(taxonomy = taxonomy_type, user_taxonomy= user_taxo)
    reftree = run_load_reftree(rtree= reftree, t = t, taxonomy_db = taxonomy_db)
    taxonomy_counter = run_load_taxonomy_counter(reftree=reftree, user_taxonomy_counter = user_counter)

    #If tree is annotated with OGD
    if "OGD_annot" in t.props:
        tree, general_results,  taxo_stats,  parameters = upload_annotated_tree(t, name_tree, reftree, taxonomy_counter, user_taxo, taxonomy_type, taxonomy_db)
        

        return tree, general_results,  taxo_stats,  parameters

    else:
    #Pre-analysis: rooting, annotate tree, etc
        t , sp_set, total_mems_in_tree, SPTOTAL, CONTENT, user_props = run_preanalysis(t, name_tree, taxonomy_db, midpoint)
    
        current_data["tree_name"] = name_tree
        current_data["reftree"] = reftree
        current_data["total_species"] = SPTOTAL
        current_data["tree_cache"] = CONTENT
        current_data["taxo_counter"] = taxonomy_counter
        current_data["taxo_db"] = taxonomy_db
        current_data["taxo_type"] = taxonomy_type
        current_data["taxo_user"] = user_taxo
        current_data["total_mems"] = total_mems_in_tree

        for p in user_props:
            USER_PROPS.append(p)

        

        #Clean properties
        t, all_props = run_clean_properties(t)
        t = get_newick(t, all_props)
        current_data["tree"] = t

        general_results = defaultdict()
        taxo_stats = defaultdict()
        parameters = defaultdict()

        return t, general_results,  taxo_stats,  parameters



@app.route('/upload_data', methods=['GET', 'POST'])
def upload_data():    
        
   
    if bool(current_data):
        current_data.clear()
        
    # Load tree 
    if not request.files['tree'].filename == "":      
        f = request.files['tree']
        name_file = secure_filename(f.filename)     
        ori_tree = f.read().decode("utf-8")       
        success_message = "Load tree: "+ name_file
    else:
        error_message = "Upload tree"
        return render_template('index.html', error_message = error_message)


    # Load Aln
    aln = request.files["aln"]
    print(current_data)
    if aln.filename == '':
        aln_file = None
    else:
        aln_file = request.files['aln']          
        name_aln = secure_filename(aln.filename)
        current_data['aln_name'] = name_aln
        
    
        #Get fastas: Core-OGs, seqs_out
        #Create dir to compute recovery pipelin
        try:
            path = UPLOAD_FOLDER+'/user_data'
            subprocess.run("rm -r %s" % path, shell = True)
            subprocess.run("mkdir -p %s" % path, shell = True)


        except OSError:
            print ("Creation of the directory user data folder failed")
        else:
            print ("Successfully created the directory user data folder ")

        #Write fastas
        aln.save(path+'/'+name_aln)
    
        path2fasta_server = path+'/'+name_aln
        current_data['aln_path'] = path2fasta_server

   
    # Load more Info: 
    #   Taxonomy type (NCBI, GTDB) 
    #   Rooting (Midpoint, MinVar)  
    #   Reference tree
    #   User Taxonomic Counter 
    #   User Taxonomy DB
    
    taxonomy_type = request.form["taxo_type"]
    rooting = request.form["rooting"]
   # annot_tree = request.form["annotated_tree"]

    rt = request.files["rtree"]
    if rt.filename == '':    
        reftree = None
    else:
        reftree = rt.read().decode("utf-8")
    
    user_counter_file = request.files["count_taxo"]
    if user_counter_file.filename == '':    
        user_counter = None
    else:
        user_counter = json.load(user_counter_file)

    user_taxo_file = request.form["user_taxonomy_database"]
    print(user_taxo_file)
    if   user_taxo_file == 'Egg5':    
        user_taxo = "/data/projects/og_delineation/data/taxonomy/e5.taxa.sqlite"
    elif  user_taxo_file == 'Egg6':
        user_taxo = "/data/projects/og_delineation/data/taxonomy/e6.taxa.sqlite"

   
    tree, general_results,  taxo_stats,  parameters  = run_upload(ori_tree, name_file, reftree, user_counter, user_taxo, taxonomy_type, rooting)

   

    data = dict()
    data['newick'] = tree
    data['id'] = '0'
    data['name'] = 'upload_tree'
    headers = {"Authorization": 'Bearer hello'}
    res = requests.post(url, data = data, headers=headers)
    print ('response from server:',res.text)
    
    
    
    #t = Tree(newick)
    return render_template('index.html', general_results = general_results, taxonomy_result = taxo_stats, parameters = parameters, hide = 'No')
    #return render_template('index.html', success_message = success_message ) 



#RUN ANALYSIS OF DUPLICATIONS AND SPECIATION NODES
@app.route('/run_analysis', methods=['GET', 'POST'])
def run_analysis():
    
    #current_data came from run_upload()
    t = current_data["tree"]
    
    
     #Clean tree properties from previous analysis
    if 'properties' in current_data.keys():
        if isinstance(t, str):
            t = PhyloTree(t, format = 1)
            
        t = clean_tree(t)
        t, all_props = run_clean_properties(t)
        t = get_newick(t, all_props)


    outliers_node = float(request.form["out_node"])
    outliers_reftree = float(request.form["out_reft"])
    inherit_outliers = request.form["inherit_out"]
    sp_loss_perc = float(request.form["p_loss"])
    so_cell_org = float(request.form["so_cell_org"])
    so_euk = float(request.form["so_euk"])
    so_bact = float(request.form["so_bact"])
    so_arq = float(request.form["so_arq"])
    
    parameters = {
        "out_node": outliers_node,
        "out_reft": outliers_reftree,
        "inherit_outliers": inherit_outliers,
        "p_loss": sp_loss_perc,
        "so_cell_org": so_cell_org,
        "so_euk": so_euk,
        "so_bact": so_bact,
        "so_arq": so_arq
    }
   
    current_data["parameters"] = parameters

    CONTENT = current_data["tree_cache"]
    taxonomy_counter = current_data["taxo_counter"]
    taxonomy_type = current_data["taxo_type"]
    user_taxo = current_data["taxo_user"]
    SPTOTAL = current_data["total_species"]
    reftree = current_data["reftree"]
    total_mems_in_tree = current_data["total_mems"]
    
    
    #Clean properties from previous analysis
    # if 'properties' in current_data.keys():
    
        # for node in t.traverse():
            # for prop in current_data['properties']:
                # if prop not in  USER_PROPS:
                    # node.del_prop(prop)


    taxonomy_db = run_load_taxonomy(taxonomy = taxonomy_type, user_taxonomy= user_taxo)
    
    #Outliers and Dups score functions
    t =  run_outliers_dup_score(t, outliers_node, outliers_reftree, sp_loss_perc, so_cell_org, so_arq, so_bact, so_euk, CONTENT, taxonomy_counter, taxonomy_db, SPTOTAL, reftree, inherit_outliers)
         
    #Detect duplications and OGs
    total_mems_in_ogs = set()
    t, total_mems_in_ogs, ogs_info, taxid_dups_og = run_dups_and_ogs(t, outliers_node, outliers_reftree, sp_loss_perc, so_cell_org, so_arq, so_bact, so_euk, taxonomy_db, total_mems_in_tree)
   

    t, ogs_info = add_ogs_up_down(t, ogs_info)
  
    #Extended OGs at each taxid level
    taxlev2ogs =  get_taxlevel2ogs(t, taxid_dups_og, total_mems_in_ogs,taxonomy_db)
    
    
    #Flag seqs out OGs
    t = flag_seqs_out_og(t, total_mems_in_ogs, total_mems_in_tree)
    best_match = defaultdict()
    recovery_seqs = set()

    #Run recovery pipeline
    if 'aln_path' in current_data.keys():

        aln_path = current_data['aln_path']
        aln_name = current_data['aln_name']
        mode = 'fast'
        pathout = UPLOAD_FOLDER+'/user_data'
        
        clean_folder(pathout, aln_path)
        run_write_fastas(aln_path, aln_name, pathout, ogs_info, total_mems_in_ogs, total_mems_in_tree, mode)
       
        #Build HMM
        run_create_hmm_og(pathout)

        #Run hmmscan
        tblfile = run_hmmscan(pathout)

        #Get best match: for each seqs, best og
        best_match = get_best_match(tblfile)
        
        #og_info_recovery = og_name : recover_seqs
        og_info_recovery, recovery_seqs = expand_hmm(best_match, ogs_info)
        recovery_seqs = set(best_match.keys())
        
        og_info_updated = update_og_info(ogs_info, og_info_recovery)

        total_mems_in_ogs.update(recovery_seqs)

        #Update in taxlev2ogs
        taxlev2ogs_updated = update_taxlevel2ogs(taxlev2ogs, og_info_recovery, og_info_updated) 

        t = update_tree(t, og_info_recovery)
        t, ogs_info = add_ogs_up_down(t, og_info_updated)


    
    seq2ogs = get_seq2og(t, best_match)
    #write_seq2ogs(seq2ogs, pathout)

    #Clean properties
    t, all_props = run_clean_properties(t)
    
    
    #Prune tree, First I need to copy the tree
    dup_tree = t.copy("deepcopy")
    prune_t = prune_tree(dup_tree, total_mems_in_ogs)


    all_props = set()
    for n in t.traverse():
        for p in n.props:
            all_props.add(p)
     
    
    # Get tree in newick format, necessary for ete4 smartview
    t = get_newick(t, all_props)
    t_prune = get_newick(prune_t, all_props)
    
    #Save results
    current_data["tree"] = t
    
    current_data["prune_tree"] = t_prune
    current_data["properties"] = all_props
    current_data["mems_in_ogs"] = total_mems_in_ogs


    glob_og_info.clear()
    for og_name, og_info in ogs_info.items():
        glob_og_info[og_name] = og_info

    glob_taxlev2ogs.clear()
    for taxid, info in taxlev2ogs.items():
        glob_taxlev2ogs[taxid] = info
  

    stats_taxo = defaultdict(dict)
    for tax, info in taxlev2ogs.items():
        sci_name = info["sci_name"]
        name = str(sci_name)+'_'+str(tax)
        stats_taxo[name]["num_ogs"] = len(info["ogs_names"])
        stats_taxo[name]["num_mems"]= len(info["mems"])
    
    glob_stats_taxo.clear()
    for k,v in stats_taxo.items():
        glob_stats_taxo[k] = v
     
    name_file = current_data["tree_name"]


    general_results["Tree_name"] = name_file
    general_results["Total_seqs"] = len(total_mems_in_tree)
    general_results["Total_species"] = SPTOTAL
    general_results["Seqs_in_OGs"] = len(set(total_mems_in_ogs))
    general_results["Recovery_seqs"] = len(recovery_seqs)
    general_results["Seqs_out_OGs"] = len(total_mems_in_tree)-len(total_mems_in_ogs)
    general_results["Num_OGs"] = len(ogs_info)
 
   
    #Send data to ete4 smartview server
    data = dict()
    data['newick'] = t
    data['id'] = '0'
    data['name'] = 'upload_tree'
    headers = {"Authorization": 'Bearer hello'}
    res = requests.post(url, data = data, headers=headers)

    
    return render_template('index.html', general_results = general_results, taxonomy_result = stats_taxo, parameters = parameters, hide = 'no')



#RUN RECOVERY SECQUENCES OUT FROM CORE-OGS
# @app.route('/run_recover', methods = ['GET', 'POST'])
# def run_recover():

   
    # fasta = request.files['fasta']          
    # name_fasta = secure_filename(fasta.filename)
    
    # #Get fastas: Core-OGs, seqs_out
    # #Create dir to compute recovery pipelin
    # try:
        # path = UPLOAD_FOLDER+'/user_data'
        # subprocess.run("rm -r %s" % path, shell = True)
        # subprocess.run("mkdir -p %s" % path, shell = True)
        
        
    # except OSError:
        # print ("Creation of the directory user data folder failed")
    # else:
        # print ("Successfully created the directory user data folder ")
    
    # #Write fastas
    # fasta.save(path+'/'+name_fasta)
    # total_mems_in_tree = current_data["total_mems"]
    # total_mems_in_ogs = current_data["mems_in_ogs"]

    # path2fasta_server = path+'/'+name_fasta

    # run_write_fastas(path2fasta_server, name_fasta, path, glob_og_info, total_mems_in_ogs, total_mems_in_tree)
    
    
    # #Build HMM
    # run_create_hmm_og(path)

    # #Run hmmscan
    # tblfile = run_hmmscan(path)
    # #print(tblfile)

    # #Get best match
    # best_match = get_best_match(tblfile)
    
    # #Create table for og after recovery
    # og_info_recovery,  total_recovery_seqs = expand_hmm(best_match, glob_og_info)
    # print(type(og_info_recovery))
    # recovery_seqs = set(best_match.keys())

    
    # glob_og_info_updated = update_og_info(glob_og_info, og_info_recovery)
    
    # total_mems_in_ogs.update(recovery_seqs)

#    #Update in taxlev2ogs
    # glob_taxlev2ogs_updated = update_taxlevel2ogs(glob_taxlev2ogs, og_info_recovery, glob_og_info_updated) 
    

    # stats_taxo_updated = defaultdict(dict)
    # for tax, info in glob_taxlev2ogs_updated.items():
        # sci_name = info["sci_name"]
        # name = str(sci_name)+'_'+str(tax)
        # stats_taxo_updated[name]["num_ogs"] = len(info["ogs_names"])
        # stats_taxo_updated[name]["num_mems"]= len(info["mems"])
    

    # name_file = current_data["tree_name"]
    # SPTOTAL = current_data["total_species"]

    # general_results["Tree name"] = name_file
    # general_results["Total seqs"] = len(total_mems_in_tree)
    # general_results["Total_species"] = SPTOTAL
    # general_results["Seqs in OGs"] = len(total_mems_in_ogs)
    # general_results["Recovery seqs"] = len(recovery_seqs) 
    # general_results["Seqs out OGS"] = len(total_mems_in_tree)-len(total_mems_in_ogs)
    # general_results["Num Ogs"] = len(glob_og_info_updated)    

#    # print(total_mems_in_tree.difference(total_mems_in_ogs))

    # # stats_taxo = defaultdict()
    # # stats_taxo = taxo_stats
    # parameters = current_data["parameters"]
    
    
    
    # t = current_data["tree"]
    # data = dict()
    # data['newick'] = t
    # data['id'] = '0'
    # data['name'] = 'upload_tree'
    # headers = {"Authorization": 'Bearer hello'}
    # res = requests.post(url, data = data, headers=headers)
    
    # return render_template('index.html', general_results = general_results, taxonomy_result = stats_taxo_updated, parameters = parameters)
               
              

#ACTIONS ON NODES AND TAXIDS 

 
@app.route('/show_prune', methods = ['GET', 'POST'])
def show_prune():

    # results = defaultdict()
    # results = general_results

    # stats_taxo = defaultdict()
    # stats_taxo = taxo_stats
    parameters = current_data["parameters"]
    
    
    t = current_data["prune_tree"]
    data = dict()
    data['newick'] = t
    data['id'] = '0'
    data['name'] = 'upload_tree'
    headers = {"Authorization": 'Bearer hello'}
    res = requests.post(url, data = data, headers=headers)
    print('AA',res)
    return render_template('index.html', general_results = general_results, taxonomy_result = stats_taxo, parameters = parameters, check="yes")


@app.route('/show_original', methods = ['GET', 'POST'])
def show_orginal(): 

    results = defaultdict()
    results = general_results

    stats_taxo = defaultdict()
    stats_taxo = taxo_stats
    parameters = current_data["parameters"]
    print(request.form)  

    t = current_data["tree"]
    data = dict()
    data['newick'] = t
    data['id'] = '0'
    data['name'] = 'upload_tree'
    headers = {"Authorization": 'Bearer hello'}
    res = requests.post(url, data = data, headers=headers)
    print('BB', res)
    return render_template('index.html', general_results = general_results, taxonomy_result = stats_taxo, parameters = parameters)
        

@app.route('/search', methods = ['GET', 'POST'])
def search():
    input_value = request.form["search"]
    stats_taxo = defaultdict()
    parameters = current_data["parameters"]
    
    for taxid_sci_name, info in taxo_stats.items():
        if input_value == taxid_sci_name.split('_',1)[0]:
            stats_taxo[taxid_sci_name] = info
        elif input_value == taxid_sci_name.split('_',1)[1]:
            stats_taxo[taxid_sci_name] = info

    results = defaultdict()
    results = general_results

    if len(stats_taxo.keys()) == 0:
        stats_taxo = taxo_stats

    return render_template('index.html', general_results = results, taxonomy_result = stats_taxo, parameters = parameters)


@app.route('/collapse/<sci_name_taxid>', methods=['GET', 'POST'])
def collapse(sci_name_taxid):

    taxid = int(sci_name_taxid.split('_',1)[1])
    t = Tree(current_data["tree"], format = 1)

    parameters = current_data["parameters"]

    #Remove collapsed prop from previous visualizations
    for node in t.traverse():
        if node.props.get('collapsed'):
            node.props['collapsed'] = 'true'

    nodes2collapse = glob_taxlev2ogs[taxid]['ogs_names']
    for node_name in nodes2collapse:
        node = t.search_nodes(name=node_name)[0]
        node.props['collapsed'] = "true"

    t, all_props = run_clean_properties(t)

    t = get_newick(t, all_props)
    current_data["tree"] = t
    
    headers = {"Authorization": 'Bearer hello'}
    data = dict()
    data['newick'] = t
    data['id'] = '0'
    data['name'] = 'upload_tree'
    res = requests.post(url, data = data, headers=headers)
    
    results = defaultdict()
    stats_taxo = defaultdict()
    results = general_results
    stats_taxo = taxo_stats


    return render_template('index.html', general_results = results, taxonomy_result = stats_taxo, parameters = parameters)


@app.route('/uncollapse/<sci_name_taxid>', methods=['GET', 'POST'])
def uncollapse(sci_name_taxid):

    taxid = int(sci_name_taxid.split('_',1)[1])
    t = Tree(current_data["tree"], format = 1)

    nodes2collapse = glob_taxlev2ogs[taxid]['ogs_names']
    for node_name in nodes2collapse:
        node = t.search_nodes(name=node_name)[0]
        node.props['collapsed'] = "false"

    t, all_props = run_clean_properties(t)

    t = get_newick(t, all_props)
    current_data["tree"] = t
    
    headers = {"Authorization": 'Bearer hello'}
    data = dict()
    data['newick'] = t
    data['id'] = '0'
    data['name'] = 'upload_tree'
    res = requests.post(url, data = data, headers=headers)
    
    results = defaultdict()
    stats_taxo = defaultdict()
    results = general_results
    stats_taxo = taxo_stats
    return render_template('index.html', general_results = results, taxonomy_result = stats_taxo)



#DOWNLOAD DATA
@app.route('/download_full_tree', methods = ['GET', 'POST'])
def download_full_tree():
    t = current_data["tree"]
    tree_name=current_data["tree_name"].split('.')[0]
    parameters = current_data["parameters"]
    values = []
    print(parameters)
    for val in parameters.values():
        print(val)
        values.append(str(val))#.split('.')[1])
    
    file_name = tree_name+'_'+'_'.join(values)+'.nw'
    head = {'Content-Disposition':  'attachment;filename='+file_name} 
    
    return Response(t, mimetype="text/css", headers=head)

    
@app.route('/download_full_og_info', methods = ['GET', 'POST'])
def download_full_og_info():
    # og_data = jsonify(OG_INFO)
    # og_data.headers['Content-Disposition'] = 'attachment;filename=og_info.json'
    # return og_data


    def generate(og_info):
       
        data = StringIO()
        w = csv.writer(data, delimiter='\t') 

        # write header
        w.writerow(('#OG_name', 'Tax_scope_OG','Dup_name','Dup_lca',  'num_OG_mems', 'Recovery_seqs','OG_up', 'OG_down', 'total_leaves', 'sp_in_OG', 'sp_out_OG', 'members'))
        yield data.getvalue()
        data.seek(0)
        data.truncate(0) 
        # write each log item
        for og_name, info in og_info.items():
            w.writerow((
                og_name,
                info['tax_scope_og'],
                info['name_dup'],
                info['lca_dup'],
                len(info['ogs_mems']),#.split('|')),
                len(info['recovery_mems']), 
                info['ogs_up'],
                info['ogs_down'],
                info['total_leaves'],
                info['num_sp_OGs'],
                info['num_sp_out'],
                '|'.join(list(info['ogs_mems']))
                
            ))
            yield data.getvalue()
            data.seek(0)
            data.truncate(0)

    # stream the response as the data is generated
    if len(glob_og_info_updated) >0:
        response = Response(generate(glob_og_info_updated), mimetype='text/csv')
    else:
        response = Response(generate(glob_og_info), mimetype='text/csv')
        
    # add a filename
    file_name = current_data["tree_name"]+'_total_OG_info.csv'
    response.headers.set("Content-Disposition", "attachment", filename=file_name)
    return response

@app.route('/download/<sci_name_taxid>', methods = ['GET', 'POST'])
def download(sci_name_taxid):

    t = PhyloTree(current_data["tree"], format = 1)
   
    taxid = int(sci_name_taxid.split('_',1)[1])
    download_data = defaultdict(dict)

    if len(glob_taxlev2ogs_updated.keys())>0:
        ogs_names_download = glob_taxlev2ogs_updated[taxid]['ogs_names']
    else:
        ogs_names_download = glob_taxlev2ogs[taxid]['ogs_names']

   
    #ogs_names_download = OG_EXTENDED[taxid]['ogs_names']
    
    
    for og_name in ogs_names_download:

        node = t.search_nodes(name=og_name)[0]
        node_name = node.props.get('name')
       
        download_data[node_name]["OG_lca"] = node.props.get('lca_node')
        
        if node.props.get('node_is_og') and taxid == int(node.props.get('lca_dup')):
            
            download_data[node_name]["DUP_name"] = node.props.get('_dup_node_name')
            download_data[node_name]["DUP_lca"] = node.props.get('lca_dup')
            download_data[node_name]["num_leaves_OGs"] = len(node.props.get('_mems_og').split('|'))
            download_data[node_name]["num_total_leaves"] = node.props.get('total_leaves')
            download_data[node_name]["OGs_up"] = node.props.get('_ogs_up', 'None')
            download_data[node_name]["OGs_down"] = node.props.get('_ogs_down', 'None')
            download_data[node_name]["num_sp_OGs"] = node.props.get('len_sp_in')
            download_data[node_name]["num_sp_out"] = node.props.get('sp_out','0')
            download_data[node_name]["Extended_og"] = 'No'
            download_data[node_name]["members_OG"] = node.props.get('_mems_og')#.split('|')
        
        else:
            download_data[node_name]["DUP_name"] = node.props.get('name')
            download_data[node_name]["DUP_lca"] = node.props.get('lca_node')
            download_data[node_name]["num_leaves_OGs"] = len(node.props.get('_leaves_in').split('|'))
            download_data[node_name]["num_total_leaves"] = node.props.get('total_leaves')
            download_data[node_name]["OGs_up"] = node.props.get('_ogs_up', 'None')
            download_data[node_name]["OGs_down"] = node.props.get('_ogs_down', 'None')
            download_data[node_name]["num_sp_OGs"] = node.props.get('len_sp_in')
            download_data[node_name]["num_sp_out"] = node.props.get('sp_out', '0')
            download_data[node_name]["Extended_og"] = 'Yes'   
            download_data[node_name]["members_OG"] = node.props.get('_mems_og')#.split('|')     




    def generate(download_data):
        data = StringIO()
        w = csv.writer(data, delimiter='\t')

        # write header
        w.writerow(('#OG_name', 'Tax_scope_OG','Dup_lca','Dup_name', 'Extended', 'OG_up', 'OG_down', 'num_leaves_OG', 'total_leaves', 'sp_in_OG', 'sp_out_OG', 'members'))
        yield data.getvalue()
        data.seek(0)
        data.truncate(0) 
        # write each log item
        for og_name, og_info in download_data.items():
            w.writerow((
                og_name,
                og_info['OG_lca'],
                og_info['DUP_lca'],
                og_info['DUP_name'],
                og_info['Extended_og'],
                og_info['OGs_up'],
                og_info['OGs_down'],
                og_info['num_leaves_OGs'],
                og_info['num_total_leaves'],
                og_info['num_sp_OGs'],
                og_info['num_sp_out'],
                og_info['members_OG']  
            ))
            yield data.getvalue()
            data.seek(0)
            data.truncate(0)

    # stream the response as the data is generated
    response = Response(generate(download_data), mimetype='text/csv')
    # add a filename
    file_name = current_data["tree_name"]+'_OG_info_'+str(taxid)+'.csv'
    response.headers.set("Content-Disposition", "attachment", filename=file_name)
    return response
    




@app.route('/')
def index():
    return render_template('index.html')



if __name__ == "__main__":
    
    app.run(port=5001, host = '138.4.138.141')