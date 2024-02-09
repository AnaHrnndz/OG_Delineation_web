from collections import defaultdict
import sys

sys.path.append('/data/projects/og_delineation')
import og_delineation
import util


def run_ogd_analyis(t_nw, taxonomy_db, args, current_data):

    taxonomy_counter = current_data["taxo_counter"]
    SPTOTAL = current_data["total_species"]
    total_mems_in_tree = current_data["total_mems"]

    # 3. Outliers and Dups score functions
    t, CONTENT  =  og_delineation.run_outliers_and_scores(t_nw, taxonomy_db, SPTOTAL, taxonomy_counter, args)

    # 4. Detect HQ-Duplications
    t, total_mems_in_ogs, taxid_dups_og  = og_delineation.run_get_main_dups(t, taxonomy_db, total_mems_in_tree, args)

    # 5. Get OGs for all taxonomical levels
    base_ogs = og_delineation.get_ogs(t, taxonomy_counter, taxonomy_db)

    # 6. Add info about nodes that split OGs up and down
    t = og_delineation.add_nodes_up_down(t)

    # 7. Annotate Basal-OGs
    base_ogs_annot = og_delineation.annot_ogs(t, base_ogs, taxonomy_db)

    best_match = defaultdict()
    recovery_seqs = set()
    #Run recovery pipeline
    if 'aln_path' in current_data.keys():

        aln_path = current_data['aln_path']
        aln_name = current_data['aln_name']
        mode = 'fast'
        pathout = UPLOAD_FOLDER+'/user_data'

        util.clean_folder(pathout, aln_path)
        og_delineation.run_write_fastas(aln_path, aln_name, pathout, ogs_info, total_mems_in_ogs, total_mems_in_tree, mode)

        #Build HMM
        og_delineation.run_create_hmm_og(pathout)

        #Run hmmscan
        tblfile = og_delineation.run_hmmscan(pathout)

        #Get best match: for each seqs, best og
        best_match = og_delineation.get_best_match(tblfile)

        #og_info_recovery = og_name : recover_seqs
        og_info_recovery, recovery_seqs = og_delineation.expand_hmm(best_match, ogs_info)
        recovery_seqs = set(best_match.keys())

        og_info_updated = og_delineation.update_og_info(ogs_info, og_info_recovery)

        total_mems_in_ogs.update(recovery_seqs)

        #Update in taxlev2ogs
        taxlev2ogs_updated = og_delineation.update_taxlevel2ogs(taxlev2ogs, og_info_recovery, og_info_updated)

        t = og_delineation.update_tree(t, og_info_recovery)
        t, ogs_info = og_delineation.add_ogs_up_down(t, og_info_updated)



    #Extended OGs at each taxid level
    taxlev2ogs =  util.get_taxlevel2ogs(t, taxid_dups_og, total_mems_in_ogs,taxonomy_db)

    # 13. Flag seqs out OGs
    t = og_delineation.flag_seqs_out_og(t, total_mems_in_ogs, total_mems_in_tree)

    # 14. Write output files
    seq2ogs = og_delineation.get_seq2og(t, best_match)
    #write_seq2ogs(seq2ogs, pathout)

    #Clean properties
    t, all_props = util.run_clean_properties(t)


    #Prune tree, First I need to copy the tree
    dup_tree = t.copy("deepcopy")

    prune_t = util.prune_tree(dup_tree, total_mems_in_ogs)

    return t, all_props, prune_t, taxlev2ogs, base_ogs_annot, total_mems_in_ogs, recovery_seqs

