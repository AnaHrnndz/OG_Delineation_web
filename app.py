from flask import Flask, render_template, request, json, Response, jsonify, send_file
import requests
from werkzeug.utils import secure_filename
import os
import sys
from collections import defaultdict
from ete4 import Tree, PhyloTree
from ete4.smartview import TreeStyle, NodeStyle, TreeLayout
import csv
from io import StringIO

sys.path.append('/data/projects/og_delineation')
from  og_delineation import run_load_tree, run_load_taxonomy, run_load_reftree, run_load_taxonomy_counter, run_preanalysis_annot_tree, get_og_info, expand_ogs_annotated_tree
from  og_delineation import run_preanalysis, run_outliers_dup_score, run_dups_and_ogs, run_clean_properties, get_newick,  expand_ogs, add_ogs_up_down, get_analysis_parameters, flag_seqs_out_og



UPLOAD_FOLDER = '/data/projects/find_ogs_web/tmp'
url = 'http://138.4.138.141:5000/trees'

app = Flask(__name__)
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER

current_data = defaultdict()
general_results = defaultdict()
taxo_stats = defaultdict(dict)
OG_INFO = defaultdict(dict)
OG_expand = defaultdict(dict)

def run_upload(tree, name_tree, reftree, user_counter, user_taxo, taxonomy_type, midpoint):
    #Load files and DBs:
    t = run_load_tree(tree= tree)
    taxonomy_db = run_load_taxonomy(taxonomy = taxonomy_type, user_taxonomy= user_taxo)
    reftree = run_load_reftree(rtree= reftree, t = t, taxonomy_db = taxonomy_db)
    taxonomy_counter = run_load_taxonomy_counter(reftree=reftree, user_taxonomy_counter = user_counter)

    #Pre-analysis: rooting, annotate tree, etc
    t , sp_set, total_mems_in_tree, SPTOTAL, CONTENT = run_preanalysis(t, name_tree, taxonomy_db, midpoint)

    current_data["tree_name"] = name_tree
    current_data["reftree"] = reftree
    current_data["total_species"] = SPTOTAL
    current_data["tree_cache"] = CONTENT
    current_data["taxo_counter"] = taxonomy_counter
    current_data["taxo_db"] = taxonomy_db
    current_data["taxo_type"] = taxonomy_type
    current_data["taxo_user"] = user_taxo
    current_data["total_mems"] = total_mems_in_tree
    
    #Clean properties
    t, all_props = run_clean_properties(t)
    t = get_newick(t, all_props)
    current_data["tree"] = t

    return t

    

def upload_annotated_tree(tree, name_tree, reftree, user_counter, user_taxo, taxonomy_type, midpoint):
    #Load files and DBs:
    t = run_load_tree(tree= tree)
    taxonomy_db = run_load_taxonomy(taxonomy = taxonomy_type, user_taxonomy= user_taxo)
    reftree = run_load_reftree(rtree= reftree, t = t, taxonomy_db = taxonomy_db)
    taxonomy_counter = run_load_taxonomy_counter(reftree=reftree, user_taxonomy_counter = user_counter)


    t , sp_set, total_mems_in_tree, SPTOTAL, CONTENT = run_preanalysis_annot_tree(t, name_tree)
    current_data["tree_name"] = name_tree
    current_data["reftree"] = reftree
    current_data["total_species"] = SPTOTAL
    current_data["tree_cache"] = CONTENT
    current_data["taxo_counter"] = taxonomy_counter
    current_data["taxo_db"] = taxonomy_db
    current_data["taxo_type"] = taxonomy_type
    current_data["taxo_user"] = user_taxo
    current_data["total_mems"] = total_mems_in_tree
    
    parameters = get_analysis_parameters(t) 
    current_data["parameters"] = parameters
    ogs_info, taxid_dups_og, total_mems_in_ogs = get_og_info(t)
    
    ogs_expanded =  expand_ogs_annotated_tree(t, taxid_dups_og, total_mems_in_ogs,taxonomy_db)

    #Clean properties
    t, all_props = run_clean_properties(t)

    t = get_newick(t, all_props)

    current_data["tree"] = t

    for og_name, og_info in ogs_info.items():
        OG_INFO[og_name] = og_info

    for taxid, info in ogs_expanded.items():
        OG_expand[int(taxid)] = info
    
    
    name_file = current_data["tree_name"]


    stats_taxo = defaultdict(dict)
    for tax, info in ogs_expanded.items():
        sci_name = info["sci_name"]
        name = str(sci_name)+'_'+str(tax)
        stats_taxo[name]["num_ogs"] = len(info["ogs_names"])
        stats_taxo[name]["num_mems"]= len(info["mems"])
    
    for k,v in stats_taxo.items():
        taxo_stats[k] = v
        

    general_results["Tree name"] = name_file
    general_results["Total seqs"] = len(total_mems_in_tree)
    general_results["Seqs in OGs"] = len(total_mems_in_ogs)
    general_results["Seqs out OGS"] = len(total_mems_in_tree)-len(total_mems_in_ogs)
    general_results["Num Ogs"] = len(ogs_info)


    
    return  t, general_results,  stats_taxo,  parameters

#UPLOAD TREE WITH OR WITHOUT ANNOTATATIONS
@app.route('/upload_tree', methods=['GET', 'POST'])
def upload_tree():    
        
    
    if not request.files['tree'].filename == "":
        f = request.files['tree']
        name_file = secure_filename(f.filename)
        ori_tree = f.read().decode("utf-8")
        success_message = "Load tree: "+ name_file
    else:
        error_message = "Upload tree"
        return render_template('index.html', error_message = error_message)

    
    taxonomy_type = request.form["taxo_type"]
    midpoint = request.form["midpoint"]
    annot_tree = request.form["annotated_tree"]

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
    

    user_taxo_file = request.files["user_taxo"]
    if user_taxo_file.filename == '':    
        user_taxo = "/data/projects/og_delineation/egg5_pfamA/data/taxonomy/e5.taxa.sqlite"
    else:
        db_name = secure_filename(user_taxo_file.filename)
        user_taxo_file.save(os.path.join(app.config['UPLOAD_FOLDER'], db_name))
        user_taxo = os.path.join(app.config['UPLOAD_FOLDER'], db_name)
   

   

    if annot_tree == 'Yes':
        tree, general_results,  taxo_stats,  parameters = upload_annotated_tree(ori_tree, name_file, reftree, user_counter, user_taxo, taxonomy_type, midpoint)
        data = dict()
        data['newick'] = tree
        data['id'] = '0'
        data['name'] = 'upload_tree'
        headers = {"Authorization": 'Bearer hello'}
        res = requests.post(url, data = data, headers=headers)
        print ('response from server:',res.text)
        return render_template('index.html', general_results = general_results, taxonomy_result = taxo_stats, parameters = parameters)
    

    tree = run_upload(ori_tree, name_file, reftree, user_counter, user_taxo, taxonomy_type, midpoint)

    data = dict()
    data['newick'] = tree
    data['id'] = '0'
    data['name'] = 'upload_tree'
    headers = {"Authorization": 'Bearer hello'}
    res = requests.post(url, data = data, headers=headers)
    print ('response from server:',res.text)
    
    
    
    #t = Tree(newick)
    return render_template('index.html', success_message = success_message ) 



#RUN ANALYSIS OF DUPLICATIONS AND SPECIATION NODES
@app.route('/run_analysis', methods=['GET', 'POST'])
def run_analysis():
    
    #t = requests.get(url + '/0' + '/newick')
    t = PhyloTree(current_data["tree"], format = 1)

    outliers_node = float(request.form["out_node"])
    outliers_reftree = float(request.form["out_reft"])
    sp_loss_perc = float(request.form["p_loss"])
    so_euk = float(request.form["so_euk"])
    so_bact = float(request.form["so_bact"])
    so_arq = float(request.form["so_arq"])

    parameters = {
        "out_node": outliers_node,
        "out_reft": outliers_reftree,
        "p_loss": sp_loss_perc,
        "so_euk": so_euk,
        "so_bact": so_bact,
        "so_arq": so_arq,
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
    if 'properties' in current_data.keys():
    
        for node in t.traverse():
            for prop in current_data['properties']:
                if prop not in  ["sci_name", "lineage", "taxid", "support", "dist", "name"]:
                    node.del_prop(prop)


    taxonomy_db = run_load_taxonomy(taxonomy = taxonomy_type, user_taxonomy= user_taxo)
    
    #Outliers and Dups score functions
    t =  run_outliers_dup_score(t, outliers_node, outliers_reftree, sp_loss_perc, so_arq, so_bact, so_euk, CONTENT, taxonomy_counter, taxonomy_db, SPTOTAL, reftree)
         
    #Detect duplications and OGs
    t, total_mems_in_ogs, ogs_info, taxid_dups_og = run_dups_and_ogs(t, outliers_node, outliers_reftree, sp_loss_perc, so_arq, so_bact, so_euk, taxonomy_db, total_mems_in_tree)
   
    
    t, ogs_info = add_ogs_up_down(t, ogs_info)
    
    
    #Expand OGs at each taxid level
    ogs_expanded =  expand_ogs(t, taxid_dups_og, total_mems_in_ogs,taxonomy_db)

    
    #Flag seqs out OGs
    t = flag_seqs_out_og(t, total_mems_in_ogs, total_mems_in_tree)

    #Clean properties
    t, all_props = run_clean_properties(t)

    t = get_newick(t, all_props)
    

    #Save results
    current_data["tree"] = t
    current_data["properties"] = all_props

    for og_name, og_info in ogs_info.items():
        OG_INFO[og_name] = og_info

    for taxid, info in ogs_expanded.items():
        OG_expand[taxid] = info


    stats_taxo = defaultdict(dict)
    for tax, info in ogs_expanded.items():
        sci_name = info["sci_name"]
        name = str(sci_name)+'_'+str(tax)
        stats_taxo[name]["num_ogs"] = len(info["ogs_names"])
        stats_taxo[name]["num_mems"]= len(info["mems"])
    
    for k,v in stats_taxo.items():
        taxo_stats[k] = v
     
    name_file = current_data["tree_name"]

    general_results["Tree name"] = name_file
    general_results["Total seqs"] = len(total_mems_in_tree)
    general_results["Seqs in OGs"] = len(total_mems_in_ogs)
    general_results["Seqs out OGS"] = len(total_mems_in_tree)-len(total_mems_in_ogs)
    general_results["Num Ogs"] = len(ogs_info)

    
    #Send data to ete4 smartview server
    data = dict()
    data['newick'] = t
    data['id'] = '0'
    data['name'] = 'upload_tree'
    headers = {"Authorization": 'Bearer hello'}
    res = requests.post(url, data = data, headers=headers)
    
    
    return render_template('index.html', general_results = general_results, taxonomy_result = stats_taxo, parameters = parameters)


            
              

#ACTIONS ON NODES AND TAXIDS    
@app.route('/search', methods = ['GET', 'POST'])
def search():
    input_value = request.form["search"]
    stats_taxo = defaultdict()
    
    for taxid_sci_name, info in taxo_stats.items():
        if input_value == taxid_sci_name.split('_',1)[0]:
            stats_taxo[taxid_sci_name] = info
        elif input_value == taxid_sci_name.split('_',1)[1]:
            stats_taxo[taxid_sci_name] = info

    results = defaultdict()
    results = general_results

    return render_template('index.html', general_results = results, taxonomy_result = stats_taxo)


@app.route('/collapse/<sci_name_taxid>', methods=['GET', 'POST'])
def collapse(sci_name_taxid):

    taxid = int(sci_name_taxid.split('_',1)[1])
    t = Tree(current_data["tree"], format = 1)

    #Remove collapsed prop from previous visualizations
    for node in t.traverse():
        if node.props.get('collapsed'):
            node.props['collapsed'] = 'true'

    nodes2collapse = OG_expand[taxid]['ogs_names']
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


    return render_template('index.html', general_results = results, taxonomy_result = stats_taxo)


@app.route('/uncollapse/<sci_name_taxid>', methods=['GET', 'POST'])
def uncollapse(sci_name_taxid):

    taxid = int(sci_name_taxid.split('_',1)[1])
    t = Tree(current_data["tree"], format = 1)

    nodes2collapse = OG_expand[taxid]['ogs_names']
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
    for val in parameters.values():
        values.append(str(val).split('.')[1])
    
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
        w = csv.writer(data) 

        # write header
        w.writerow(('#OG_name', 'Tax_scope_OG','Dup_name','Dup_lca', 'Size_OG', 'OG_up', 'OG_down'))
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
                len(info['ogs_mems']),
                info['ogs_up'],
                info['ogs_down'],
                
            ))
            yield data.getvalue()
            data.seek(0)
            data.truncate(0)

    # stream the response as the data is generated
    response = Response(generate(OG_INFO), mimetype='text/csv')
    # add a filename
    response.headers.set("Content-Disposition", "attachment", filename="og_info.csv")
    return response

@app.route('/download/<sci_name_taxid>', methods = ['GET', 'POST'])
def download(sci_name_taxid):

    t = PhyloTree(current_data["tree"], format = 1)
   
    taxid = int(sci_name_taxid.split('_',1)[1])
    download_data = defaultdict(dict)
   
    ogs_names_download = OG_expand[taxid]['ogs_names']
    
    for og_name in ogs_names_download:

        node = t.search_nodes(name=og_name)[0]
        node_name = node.props.get('name')
       
        download_data[node_name]["OG_lca"] = node.props.get('lca_node')
        
        if node.props.get('node_is_og') and taxid == int(node.props.get('lca_dup')):
            
            download_data[node_name]["DUP_name"] = node.props.get('_dup_node_name')
            download_data[node_name]["DUP_lca"] = node.props.get('lca_dup')
            download_data[node_name]["size_OGs"] = len(node.props.get('_mems_og').split('|'))
            download_data[node_name]["OGs_up"] = node.props.get('ogs_up')
            download_data[node_name]["OGs_down"] = node.props.get('ogs_down')
            download_data[node_name]["Expanded_og"] = 'No'
        
        else:
            download_data[node_name]["DUP_name"] = node.props.get('name')
            download_data[node_name]["DUP_lca"] = node.props.get('lca_node')
            download_data[node_name]["size_OGs"] = len(node.props.get('_leaves_in').split('|'))
            download_data[node_name]["OGs_up"] = node.props.get('ogs_up')
            download_data[node_name]["OGs_down"] = node.props.get('ogs_down')
            download_data[node_name]["Expanded_og"] = 'Yes'        




    def generate(download_data):
        data = StringIO()
        w = csv.writer(data) 

        # write header
        w.writerow(('#OG_name', 'OG_lca','Dup_lca','Dup_name', 'Expanded', 'OG_down', 'OG_up', 'Size_OG'))
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
                og_info['Expanded_og'],
                og_info['OGs_down'],
                og_info['OGs_up'],
                og_info['size_OGs']  
            ))
            yield data.getvalue()
            data.seek(0)
            data.truncate(0)

    # stream the response as the data is generated
    response = Response(generate(download_data), mimetype='text/csv')
    # add a filename
    response.headers.set("Content-Disposition", "attachment", filename="log.csv")
    return response
    



@app.route('/')
def index():
    return render_template('index.html')



if __name__ == "__main__":
    app.run(port=5001, host = '138.4.138.141')