from flask import Flask, render_template, request, json, Response, jsonify, send_file
import requests
from werkzeug.utils import secure_filename
import os
import sys
from collections import defaultdict
from ete4 import Tree, PhyloTree, SeqGroup
#from ete4.smartview import TreeStyle, NodeStyle, TreeLayout
import csv
import subprocess
from io import StringIO

sys.path.append('/data/projects/og_delineation_web/bin')

sys.path.append('/data/projects/og_delineation')
import og_delineation

sys.path.append('/data/projects/og_delineation_web/ogd_web')
import load_data
import util
from ogd_analysis import run_ogd_analyis



UPLOAD_FOLDER = '/data/projects/og_delineation_web/tmp'
url = 'http://138.4.138.141:5000/trees'

app = Flask(__name__)
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER

app.config['MAX_CONTENT_LENGTH'] = 250 * 1024 * 1024    # 50 Mb limit  ############

current_data = defaultdict()
general_results = defaultdict()
# USER_PROPS = list()

glob_og_info = defaultdict(dict)
glob_og_info_updated = defaultdict(dict)
glob_stats_taxo = defaultdict(dict)
glob_taxlev2ogs = defaultdict(dict)
glob_taxlev2ogs_updated = defaultdict(dict)


@app.route('/upload_data', methods=['GET', 'POST'])
def upload_data():

    """
        Upload data:
            - Tree:
                - With annotations
                - Without annotations
            - Aln
            - Taxonomy. Choose:
                - Egg6
                - Egg5
            - Species counter
            - Reference species tree
            - Root method:
                - Midpoint
                - MinVar
    """

    # Delete previous info
    if bool(current_data):
        current_data.clear()

    # Create dir for user data
    util.create_aln_dir(UPLOAD_FOLDER)
    path = UPLOAD_FOLDER+'/user_data/'

    ### Get tree from form ###
    if not request.files['tree'].filename == "":
        f = request.files['tree']
        name_file = secure_filename(f.filename)
        ori_tree = f.read().decode("utf-8").strip()
        success_message = "Load tree: "+ name_file
    else:
        error_message = "Upload tree"
        return render_template('index.html', error_message = error_message)


    ### Get Aln from form ###
    aln = request.files["aln"]
    if aln.filename == '':
        aln_file = None
    else:
        aln_file = request.files['aln']
        name_aln = secure_filename(aln.filename)
        current_data['aln_name'] = name_aln

        #Write fasta file
        aln.save(path+'/'+name_aln)

        path2fasta_server = path+'/'+name_aln
        current_data['aln_path'] = path2fasta_server


    ### Get root method for form ###
    rooting = request.form["rooting"]

    ### Load reference species tree ###
    rt = request.files["rtree"]
    if rt.filename == '':
        reftree = None
    else:
        reftree = rt.read().decode("utf-8")


    ### Load species counter ###
    user_counter_file = request.files["count_taxo"]
    if user_counter_file.filename == '':
        user_counter = None
    else:
        user_counter = json.load(user_counter_file)


    ### Load Taxonomy ###
    taxonomy_type = request.form["taxo_type"]

    user_taxo_file = request.form["user_taxonomy_database"]
    if   user_taxo_file == 'Egg5':
        user_taxo = "/data/projects/og_delineation/data/taxonomy/e5.taxa.sqlite"
    elif  user_taxo_file == 'Egg6':
        user_taxo = "/data/projects/og_delineation/data/taxonomy/e6.taxa.sqlite"


    ### Prepare data for first visualization ###
    t_nw, general_results,  taxo_stats,  parameters, SPTOTAL, taxonomy_db, taxonomy_counter, total_mems_in_tree \
          = load_data.run_upload(ori_tree, name_file, reftree, user_counter, user_taxo, taxonomy_type, rooting, UPLOAD_FOLDER)

    ### Save data ###
    current_data["tree_nw"] = t_nw
    current_data["tree_name"] = name_file
    current_data["reftree"] = reftree
    current_data["total_species"] = SPTOTAL
    current_data["taxo_counter"] = taxonomy_counter
    current_data["taxo_db"] = taxonomy_db
    current_data["taxo_type"] = taxonomy_type
    current_data["taxo_user"] = user_taxo
    current_data["user_counter"] = user_counter
    current_data["total_mems"] = total_mems_in_tree

    # Send data to ete4 smartview server
    data = dict()
    data['newick'] = t_nw
    data['id'] = '0'
    data['name'] = 'upload_tree'
    res = requests.post(url, data = data)
    print ('response from server:',res.text)

    return render_template('index.html', general_results = general_results, taxonomy_result = taxo_stats, parameters = parameters, hide = 'No')



#RUN ANALYSIS OF DUPLICATIONS AND SPECIATION NODES
@app.route('/run_analysis', methods=['GET', 'POST'])
def run_analysis():

    # Current_data came from run_upload()
    t_nw = current_data["tree_nw"]

    #Taxonomy have to be load
    taxonomy_db = og_delineation.load_taxonomy(taxonomy = current_data["taxo_type"], user_taxonomy = current_data["taxo_user"])

    args_dict = defaultdict(dict)
    args_dict["user_taxonomy_counter"] = current_data["user_counter"]
    args_dict["outliers_node"] = float(request.form["out_node"])
    args_dict["outliers_reftree"] = float(request.form["out_reft"])
    args_dict["inherit_out "] = request.form["inherit_out"]
    args_dict["so_bact"] = float(request.form["so_bact"])
    args_dict["so_euk"] = float(request.form["so_euk"])
    args_dict["so_arq"] = float(request.form["so_arq"])
    args_dict["so_cell_org"] = float(request.form["so_cell_org"])
    args_dict["sp_loss_perc"] = float(request.form["p_loss"])
    args = util.dict2args(args_dict)
    current_data["parameters"] = args_dict

    # Run al steps from OGD analysis
    t, all_props, prune_t, taxlev2ogs, base_ogs_annot, total_mems_in_ogs, recovery_seqs \
        = run_ogd_analyis(t_nw, taxonomy_db, args, current_data)

    all_props = set()
    for n in t.traverse():
        for p in n.props:
            all_props.add(p)

    # Get tree in newick format, necessary for ete4 smartview
    t_nw = util.get_newick(t, all_props)
    t_prune_nw = util.get_newick(prune_t, all_props)

    # Save results
    current_data["tree_nw"] = t_nw
    current_data["prune_tree"] = t_prune_nw
    current_data["properties"] = all_props
    current_data["mems_in_ogs"] = total_mems_in_ogs

    glob_og_info.clear()
    for og_name, og_info in base_ogs_annot.items():
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


    general_results["Tree_name"] = current_data["tree_name"]
    general_results["Total_seqs"] = len(current_data["total_mems"])
    general_results["Total_species"] = current_data["total_species"]
    general_results["Seqs_in_OGs"] = len(set(total_mems_in_ogs))
    general_results["Recovery_seqs"] = len(recovery_seqs)
    general_results["Seqs_out_OGs"] = len(current_data["total_mems"])-len(set(total_mems_in_ogs))
    general_results["Num_OGs"] = len(base_ogs_annot)


    # Send data to ete4 smartview server
    data = dict()
    data['newick'] = t_nw
    data['id'] = '0'
    data['name'] = 'upload_tree'
    res = requests.post(url, data = data)

    return render_template('index.html', general_results = general_results, taxonomy_result = stats_taxo, parameters = args_dict, hide = 'no')


#ACTIONS ON NODES AND TAXIDS
@app.route('/show_prune', methods = ['GET', 'POST'])
def show_prune():

    results = defaultdict()
    results = general_results

    stats_taxo = defaultdict()
    stats_taxo = glob_stats_taxo
    parameters = current_data["parameters"]


    t = current_data["prune_tree"]
    data = dict()
    data['newick'] = t
    data['id'] = '0'
    data['name'] = 'upload_tree'

    res = requests.post(url, data = data )

    return render_template('index.html', general_results = results, taxonomy_result = stats_taxo, parameters = parameters, check='Yes', hide = 'no')


@app.route('/show_original', methods = ['GET', 'POST'])
def show_orginal():

    results = defaultdict()
    results = general_results

    stats_taxo = defaultdict()
    stats_taxo = glob_stats_taxo
    parameters = current_data["parameters"]


    t = current_data["tree"]
    data = dict()
    data['newick'] = t
    data['id'] = '0'
    data['name'] = 'upload_tree'

    res = requests.post(url, data = data)
    return render_template('index.html', general_results = results, taxonomy_result = stats_taxo, parameters = parameters, hide = 'no')


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
    t = Tree(current_data["tree"])

    parameters = current_data["parameters"]

    #Remove collapsed prop from previous visualizations
    for node in t.traverse():
        if node.props.get('collapsed'):
            node.props['collapsed'] = 'true'

    nodes2collapse = glob_taxlev2ogs[taxid]['ogs_names']
    for node_name in nodes2collapse:
        node = list(t.search_nodes(name=node_name))[0]
        node.props['collapsed'] = "true"

    t, all_props = og_delineation.run_clean_properties(t)

    t = og_delineation.get_newick(t, all_props)
    current_data["tree"] = t


    data = dict()
    data['newick'] = t
    data['id'] = '0'
    data['name'] = 'upload_tree'
    res = requests.post(url, data = data)

    results = defaultdict()
    stats_taxo = defaultdict()
    results = general_results
    stats_taxo = glob_stats_taxo


    return render_template('index.html', general_results = results, taxonomy_result = stats_taxo, parameters = parameters, hide = 'no')


@app.route('/uncollapse/<sci_name_taxid>', methods=['GET', 'POST'])
def uncollapse(sci_name_taxid):

    taxid = int(sci_name_taxid.split('_',1)[1])
    t = Tree(current_data["tree"])

    nodes2collapse = glob_taxlev2ogs[taxid]['ogs_names']
    for node_name in nodes2collapse:
        node = list(t.search_nodes(name=node_name))[0]
        node.props['collapsed'] = "false"

    t, all_props = og_delineation.run_clean_properties(t)

    t = og_delineation.get_newick(t, all_props)
    current_data["tree"] = t


    data = dict()
    data['newick'] = t
    data['id'] = '0'
    data['name'] = 'upload_tree'
    res = requests.post(url, data = data)

    results = defaultdict()
    stats_taxo = defaultdict()

    results = general_results
    stats_taxo = glob_stats_taxo

    return render_template('index.html', general_results = results, taxonomy_result = stats_taxo, parameters = parameters, hide = 'no')



#DOWNLOAD DATA
@app.route('/download_full_tree', methods = ['GET', 'POST'])
def download_full_tree():
    t = current_data["tree_nw"]
    tree_name=current_data["tree_name"].split('.')[0]
    parameters = current_data["parameters"]
    values = []

    for val in parameters.values():
        values.append(str(val))

    file_name = tree_name+'_'+'_'.join(values)+'.nw'
    head = {'Content-Disposition':  'attachment;filename='+file_name}

    return Response(t, mimetype="text/css", headers=head)


@app.route('/download_full_og_info', methods = ['GET', 'POST'])
def download_full_og_info():


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