# -*- coding: utf-8 -*-
"""
Created on Mon Feb 15 09:52:17 2021

@author: ananthansadagopan
"""

from varcode import Variant
from pyensembl import ensembl_grch37
from mhcflurry import Class1PresentationPredictor
from datetime import datetime
import pandas as pd
import json
import os
import sys
import tensorflow as tf
import ntpath

#Convert MAF into dataframe, filtering out variants that can't generate neoantigens

def read_maf(maf_filename_input, log_file_input):
    
    if not maf_filename_input or not os.path.isfile(maf_filename_input):
        with open(log_file_input, "a+") as f:
            f.write("MAF file does not exist: " + str(maf_filename_input))
        sys.exit()
    else:
        maf_df_old = pd.read_csv(maf_filename_input, sep="\t") #Read the MAF, must have HGVSp_Short for each variant

    potentially_immunogenic_mutation_types = ["frame_shift_del", "frame_shift_ins", "in_frame_del", "in_frame_ins", "missense_mutation", "nonstop_mutation"]

    maf_df = maf_df_old.drop(['Consequence'], axis=1, errors='ignore')

    maf_df['Variant_Classification'] = maf_df['Variant_Classification'].str.lower()

    maf_df = maf_df[maf_df['Variant_Classification'].isin(potentially_immunogenic_mutation_types)]

    maf_df = maf_df[maf_df['HGVSp_Short'].notna()] #Remove NaN in HGVSp_Short, this is essential for the next step

    maf_df = maf_df[~maf_df['HGVSp_Short'].str.endswith("*")] #Subset out nonsense mutations

    return maf_df

#get amino acid context from a variant

def get_AA_context(effect_input, variant_classification):
    
    context = None
    wt_context = None
    
    var_seq = effect_input.mutant_protein_sequence
    normal_seq = effect_input.original_protein_sequence
    var_aa_start = effect_input.aa_mutation_start_offset
    var_aa_end = effect_input.aa_mutation_end_offset
    if var_seq != None and var_aa_start != None and var_aa_end != None:
        if var_aa_start-10 > 0:
            start_context = var_aa_start-10
        else:
            start_context = 0
        if var_aa_end+10 < len(var_seq):
            end_context = var_aa_end+10
        else:
            end_context = len(var_seq)-1
        context = var_seq[start_context:end_context+1]
    
        if normal_seq != None:
            if variant_classification == "missense_mutation" or variant_classification == "in_frame_del": 
                wt_context = normal_seq[start_context:end_context+1]
            elif variant_classification == "in_frame_ins":
                wt_context = normal_seq[start_context:end_context+1]
            elif variant_classification == "frame_shift_del" or variant_classification == "frame_shift_ins" or variant_classification == "nonstop_mutation":
                wt_context = normal_seq[start_context:var_aa_start+1]

    return context, wt_context

#Execution of get amino acid context

def exec_get_AA_context(transcript_id, chromosome_number, start_pos, ref_allele, alt_allele, log_file_input, var_classification_list):
    
    mutant_context_list = []
    wt_context_list = []
    transcript_list = []    

    a=0
    while a<len(transcript_id):
        current_variant = Variant(contig=chromosome_number[a], start=start_pos[a], ref=ref_allele[a], alt=alt_allele[a], ensembl=ensembl_grch37)
        if transcript_id[a] != transcript_id[a]:
            """
            with open(log_file_input, "a+") as log_file:
                log_file.write("WARNING: No transcript ID for entry number: " + str(a+1) + "; all transcripts tested instead.\n")
            all_effects = current_variant.effects()
            """

            with open(log_file_input, "a+") as log_file:
                log_file.write("WARNING: No transcript ID for entry number: " + str(a+1) + "; highest priority transcript instead.\n")

            all_effects = [current_variant.effects().top_priority_effect()]

            temp_wt_context_list = []
            temp_var_context_list = []
            temp_transcript_list = []
            current_normal_context = ""
            current_variant_context = ""
            current_transcript = ""

            for effect in all_effects:
                current_variant_context, current_normal_context = get_AA_context(effect, var_classification_list[a])
                current_transcript = effect.transcript_id
                temp_wt_context_list.append(current_normal_context)
                temp_var_context_list.append(current_variant_context)
                temp_transcript_list.append(current_transcript)

            mutant_context_list.append(temp_var_context_list)
            wt_context_list.append(temp_wt_context_list)
            transcript_list.append(temp_transcript_list)

        else:
            b=0
            while b<len(current_variant.effects()):
                if current_variant.effects()[b].transcript_id == transcript_id[a]:
                    current_variant_context, current_normal_context = get_AA_context(current_variant.effects()[b], var_classification_list[a])
                    current_transcript = transcript_id[a]
                    break
                b=b+1
            
            if current_normal_context and current_variant_context:
                mutant_context_list.append(current_variant_context)
                wt_context_list.append(current_normal_context)
                transcript_list.append(current_transcript)

        if len(mutant_context_list) < a+1:
            topPriorityEffect = current_variant.effects().top_priority_effect()
            current_variant_context, current_normal_context = get_AA_context(topPriorityEffect, var_classification_list[a])
            current_transcript = topPriorityEffect.transcript_id
            with open(log_file_input, "a+") as log_file:
                log_file.write("No varcode transcript_id matches input transcript_id: " + str(transcript_id[a]) + "; This corresponds to entry " + str(a+1) + "; Most severe consequence used instead: " + str(current_transcript) + "\n")
            mutant_context_list.append(current_variant_context)
            wt_context_list.append(current_normal_context)
            transcript_list.append(current_transcript)
        a=a+1
    
    return mutant_context_list, wt_context_list, transcript_list

def output_result(predictor, affinity_cutoff, final_json_list, patient_alleles_dict, log_file_input, df_columns, current_barcode, current_Hugo_Symbol, chromosome, start, end, ref, ta1, ta2, current_transcript_tested, current_wt_context, current_mutated_context, master_df):

    invalid_chars = ['U', 'O', 'B', 'J', 'Z', 'X'] #Non-conventional amino acids failing in MHCflurry

    if current_wt_context != None and current_mutated_context != None and not any(element in current_mutated_context for element in invalid_chars) and not any(element in current_wt_context for element in invalid_chars):

        if current_barcode not in patient_alleles_dict:
            with open(log_file_input, "a+") as f:
                f.write("Barcode not in list of patient alleles or removed due to invalidity of alleles: " + str(current_barcode) + "; Skipped predictions on the variant above.\n")
            return master_df, final_json_list
        else:
            allele_list = patient_alleles_dict[current_barcode]
        t_total_binders_list, t_strong_binders_list, t_pres_sort_peptide_list, t_pres_sort_pres_score_list, t_pres_sort_pres_percentile_list, t_pres_sort_affinity_list, t_pres_sort_affinity_percentile_list, t_affinity_sort_peptide_list, t_affinity_sort_pres_score_list, t_affinity_sort_pres_percentile_list, t_affinity_sort_affinity_list, t_affinity_sort_affinity_percentile_list = mutation_immunogenicity(current_wt_context, current_mutated_context, allele_list, predictor, affinity_cutoff)        
        
        a=0
        PHBR_denom = 0

        while a<len(t_pres_sort_pres_percentile_list):
            if t_pres_sort_pres_percentile_list[a]:
                if float(t_pres_sort_pres_percentile_list[a]) != 0:
                    PHBR_denom = PHBR_denom + 1/float(t_pres_sort_pres_percentile_list[a])
                else:
                    PHBR_denom = PHBR_denom + 10000
                    with open(log_file_input, "a+") as f:
                        f.write("WARNING: Presentation Percentile is 0 for the following variant; 0.00001 was used during HBR calculations:\t" + str(current_barcode) + "\t" + str(current_Hugo_Symbol) + "\t" + str(start) + "\t" + str(end) + "\t" + str(ref) + "\t" + str(ta1) + "\t" + str(ta2) + "\t" + str(current_transcript_tested) + "\t\n")
            else:
                PHBR_denom = PHBR_denom + 1/100
            a=a+1

        PHBR = len(t_pres_sort_pres_percentile_list)/PHBR_denom
        PHBR = str(PHBR)

        a=0
        affinity_PHBR_denom = 0

        while a<len(t_affinity_sort_affinity_percentile_list):
            if t_affinity_sort_affinity_percentile_list[a]:
                if float(t_affinity_sort_affinity_percentile_list[a]) != 0:
                    affinity_PHBR_denom = affinity_PHBR_denom + 1/float(t_affinity_sort_affinity_percentile_list[a])
                else:
                    affinity_PHBR_denom = affinity_PHBR_denom + 10000
                    f.write("WARNING: Affinity Percentile is 0 for the following variant; 0.00001 was used during HBR calculations:\t" + str(current_barcode) + "\t" + str(current_Hugo_Symbol) + "\t" + str(start) + "\t" + str(end) + "\t" + str(ref) + "\t" + str(ta1) + "\t" + str(ta2) + "\t" + str(current_transcript_tested) + "\t\n")
            else:
                affinity_PHBR_denom = affinity_PHBR_denom + 1/100
            a=a+1

        affinity_PHBR = len(t_affinity_sort_affinity_percentile_list)/affinity_PHBR_denom
        affinity_PHBR = str(affinity_PHBR)

        w=0
        dataframe_rows = []
        predictions_dict_list = []
        
        while w<len(allele_list):
            curr_values = [current_barcode, current_Hugo_Symbol, chromosome, start, end, ref, ta1, ta2, current_transcript_tested, PHBR, affinity_PHBR, current_wt_context, current_mutated_context, allele_list[w], t_total_binders_list[w], t_strong_binders_list[w], t_pres_sort_peptide_list[w], t_pres_sort_pres_score_list[w], t_pres_sort_pres_percentile_list[w], t_pres_sort_affinity_list[w], t_pres_sort_affinity_percentile_list[w], t_affinity_sort_peptide_list[w], t_affinity_sort_pres_score_list[w], t_affinity_sort_pres_percentile_list[w], t_affinity_sort_affinity_list[w], t_affinity_sort_affinity_percentile_list[w]]
            dataframe_rows.append(curr_values)
            temp_dict = dict(zip(['HLA_Allele', 'N_total_neoantigens_500nM_cutoff', 'N_strong_binders_50nM_cutoff', 'Best_Presentation_Score_Peptide', 'Best_Presentation_Score_Presentation_Score', 'Best_Presentation_Score_Presentation_Percentile', 'Best_Presentation_Score_Affinity', 'Best_Presentation_Score_Affinity_Percentile', 'Best_Affinity_Peptide', 'Best_Affinity_Presentation_Score', 'Best_Affinity_Presentation_Percentile', 'Best_Affinity_Affinity', 'Best_Affinity_Affinity_Percentile'], [allele_list[w], t_total_binders_list[w], t_strong_binders_list[w], t_pres_sort_peptide_list[w], t_pres_sort_pres_score_list[w], t_pres_sort_pres_percentile_list[w], t_pres_sort_affinity_list[w], t_pres_sort_affinity_percentile_list[w], t_affinity_sort_peptide_list[w], t_affinity_sort_pres_score_list[w], t_affinity_sort_pres_percentile_list[w], t_affinity_sort_affinity_list[w], t_affinity_sort_affinity_percentile_list[w]]))
            predictions_dict_list.append(temp_dict)
            w=w+1
        
        json_predictions = json.dumps(predictions_dict_list)

        allele_output_array = {"Tumor_Sample_Barcode":current_barcode,"Hugo_Symbol":current_Hugo_Symbol,"Chromosome":chromosome,"Start_Position":start,"End_Position":end,"Reference_Allele":ref,"Tumor_Seq_Allele1":ta1,"Tumor_Seq_Allele2":ta2,"Transcript":current_transcript_tested,"Presentation_HBR":PHBR,"Affinity_HBR":affinity_PHBR, "WT_Amino_Acid_Context":current_wt_context,"Variant_Amino_Acid_Context":current_mutated_context,"Predictions":json_predictions}

        allele_output_array = json.dumps(allele_output_array)

        final_json_list.append(allele_output_array)

        temp_df = pd.DataFrame(dataframe_rows, columns=df_columns)

        master_df = pd.concat([master_df, temp_df], ignore_index=True, sort=False)
        
    elif (current_wt_context != None and current_mutated_context != None) and (any(element in current_mutated_context for element in invalid_chars) or any(element in current_wt_context for element in invalid_chars)):
        
        f.write("WARNING: Unconventional Amino Acid Present for Transcript Variant (Skipped):\t" + str(current_barcode) + "\t" + str(current_Hugo_Symbol) + "\t" + str(start) + "\t" + str(end) + "\t" + str(ref) + "\t" + str(ta1) + "\t" + str(ta2) + "\t" + str(current_transcript_tested) + "\t\n")
    
    else:
        
        f.write("WARNING: No Amino Acid Context Obtained for Transcript Variant (Skipped):\t" + str(current_barcode) + "\t" + str(current_Hugo_Symbol) + "\t" + str(start) + "\t" + str(end) + "\t" + str(ref) + "\t" + str(ta1) + "\t" + str(ta2) + "\t" + str(current_transcript_tested) + "\t\n")
        
    return master_df, final_json_list

def mutation_immunogenicity(wt_nmer, mutated_nmer, list_of_alleles, predictor, affinity_cutoff):
            
    total_binders_list = []
    strong_binders_list = []
    pres_sort_peptide_list = []
    pres_sort_pres_score_list = []
    pres_sort_pres_percentile_list = []
    pres_sort_affinity_list = []
    pres_sort_affinity_percentile_list = []
    affinity_sort_peptide_list = []
    affinity_sort_pres_score_list = []
    affinity_sort_pres_percentile_list = []
    affinity_sort_affinity_list = []
    affinity_sort_affinity_percentile_list = []
    
    allele_number = 0
    while allele_number < len(list_of_alleles):
    #500 nM minimum affinity threshold

        single_allele = list_of_alleles[allele_number]
        #mutated_best_pres_score = ""
        
        df = predictor.predict_sequences(
            sequences={
                'wt_protein': wt_nmer,
                'mutated_protein': mutated_nmer,
                },
            alleles=[single_allele],
            result="filtered",
            comparison_quantity="affinity",
            filter_value=affinity_cutoff,
            peptide_lengths=(8, 9, 10, 11),
            use_flanks=True,
            include_affinity_percentile=True,
            verbose=0)
        
        if not df.empty:
            
            #Remove same AA sequences between wt and mutated from df
            df.sort_values("sequence_name", ascending = False, inplace = True)
            df.drop_duplicates(subset ="peptide", 
                    keep = False, inplace = True)
            df = df[df.sequence_name != "wt_protein"]
            
            if not df.empty:
            
                #Sort by presentation scores in descending order (best mhcflurry model that combines antigen processing and presentation)
                presentation_scores = df['presentation_score'].tolist()
                peptides_sorted = [y for _,y in sorted(zip(presentation_scores,df['peptide'].tolist()))][::-1]
                proteins_sorted = [y for _,y in sorted(zip(presentation_scores,df['sequence_name'].tolist()))][::-1]
                affinity_sorted = [y for _,y in sorted(zip(presentation_scores,df['affinity'].tolist()))][::-1]
                affinity_percentile_sorted = [y for _,y in sorted(zip(presentation_scores,df['affinity_percentile'].tolist()))][::-1]
                presentation_score_percentile_sorted = [y for _,y in sorted(zip(presentation_scores,df['presentation_percentile'].tolist()))][::-1]
                presentation_scores_sorted = sorted(presentation_scores, reverse=True)
                
                x=0
                total_binders = 0
                strong_binders = 0
                mutated_best_peptide = ""
                mutated_best_pres_score = ""
                mutated_best_pres_score_percentile = ""
                mutated_best_affinity = ""
                mutated_best_affinity_percentile = ""
                
                while x<len(peptides_sorted):
                    if proteins_sorted[x] == "mutated_protein" and not mutated_best_peptide: 
                        mutated_best_pres_score = presentation_scores_sorted[x]
                        mutated_best_pres_score_percentile = presentation_score_percentile_sorted[x]
                        mutated_best_peptide = peptides_sorted[x]
                        mutated_best_affinity = affinity_sorted[x]
                        mutated_best_affinity_percentile = affinity_percentile_sorted[x]
                    if proteins_sorted[x] == "mutated_protein" and affinity_sorted[x] <= 50:
                        strong_binders = strong_binders + 1
                        total_binders = total_binders + 1
                    elif proteins_sorted[x] == "mutated_protein" and affinity_sorted[x] <= 500:
                        total_binders = total_binders + 1
                    x=x+1
                
                total_binders_list.append(str(total_binders))
                strong_binders_list.append(str(strong_binders))
                pres_sort_peptide_list.append(str(mutated_best_peptide))
                pres_sort_pres_score_list.append(str(mutated_best_pres_score))
                pres_sort_pres_percentile_list.append(str(mutated_best_pres_score_percentile))
                pres_sort_affinity_list.append(str(mutated_best_affinity))
                pres_sort_affinity_percentile_list.append(str(mutated_best_affinity_percentile))
                
                #Sort by affinity scores in ascending order (best mhcflurry model that combines antigen processing and presentation)
                affinities = df['affinity'].tolist()
                peptides_affinity_sorted = [y for _,y in sorted(zip(affinities,df['peptide'].tolist()))]
                proteins_affinity_sorted = [y for _,y in sorted(zip(affinities,df['sequence_name'].tolist()))]
                presentation_scores_affinity_sorted = [y for _,y in sorted(zip(affinities,df['presentation_score'].tolist()))]
                affinity_percentile_affinity_sorted = [y for _,y in sorted(zip(affinities,df['affinity_percentile'].tolist()))]
                presentation_score_percentile_affinity_sorted = [y for _,y in sorted(zip(affinities,df['presentation_percentile'].tolist()))]
                affinity_scores_sorted = sorted(affinities, reverse=False)

                affinity_sorted_mutated_best_peptide = ""
                affinity_sorted_mutated_best_pres_score = ""
                affinity_sorted_mutated_best_pres_score_percentile = ""
                affinity_sorted_mutated_best_affinity = ""
                affinity_sorted_mutated_best_affinity_percentile = ""

                x=0

                while x<len(peptides_affinity_sorted):
                    if proteins_affinity_sorted[x] == "mutated_protein" and not affinity_sorted_mutated_best_peptide:
                        affinity_sorted_mutated_best_pres_score = presentation_scores_affinity_sorted[x]
                        affinity_sorted_mutated_best_pres_score_percentile = presentation_score_percentile_affinity_sorted[x]
                        affinity_sorted_mutated_best_peptide = peptides_affinity_sorted[x]
                        affinity_sorted_mutated_best_affinity = affinity_scores_sorted[x]
                        affinity_sorted_mutated_best_affinity_percentile = affinity_percentile_affinity_sorted[x]
                        break
                    x=x+1
    
                affinity_sort_peptide_list.append(str(affinity_sorted_mutated_best_peptide))
                affinity_sort_pres_score_list.append(str(affinity_sorted_mutated_best_pres_score))
                affinity_sort_pres_percentile_list.append(str(affinity_sorted_mutated_best_pres_score_percentile))
                affinity_sort_affinity_list.append(str(affinity_sorted_mutated_best_affinity))
                affinity_sort_affinity_percentile_list.append(str(affinity_sorted_mutated_best_affinity_percentile))

            else:
                total_binders_list.append(str(0))
                strong_binders_list.append(str(0))
                pres_sort_peptide_list.append("")
                pres_sort_pres_score_list.append("")
                pres_sort_pres_percentile_list.append("")
                pres_sort_affinity_list.append("")
                pres_sort_affinity_percentile_list.append("")
                affinity_sort_peptide_list.append("")
                affinity_sort_pres_score_list.append("")
                affinity_sort_pres_percentile_list.append("")
                affinity_sort_affinity_list.append("")
                affinity_sort_affinity_percentile_list.append("")
                        
        else:
            total_binders_list.append(str(0))
            strong_binders_list.append(str(0))
            pres_sort_peptide_list.append("")
            pres_sort_pres_score_list.append("")
            pres_sort_pres_percentile_list.append("")
            pres_sort_affinity_list.append("")
            pres_sort_affinity_percentile_list.append("")
            affinity_sort_peptide_list.append("")
            affinity_sort_pres_score_list.append("")
            affinity_sort_pres_percentile_list.append("")
            affinity_sort_affinity_list.append("")
            affinity_sort_affinity_percentile_list.append("")
            
        allele_number = allele_number + 1
    
    return total_binders_list, strong_binders_list, pres_sort_peptide_list, pres_sort_pres_score_list, pres_sort_pres_percentile_list, pres_sort_affinity_list, pres_sort_affinity_percentile_list, affinity_sort_peptide_list, affinity_sort_pres_score_list, affinity_sort_pres_percentile_list, affinity_sort_affinity_list, affinity_sort_affinity_percentile_list

def annotate_maf(maf_filename_input, OUTPUT_FOLDER, key):

    key = str(key)

    if not key:
        key = str(datetime.now().strftime('%Y_%m_%d_%H_%M_%S')) + "_AA_context_annotated_"

    if not OUTPUT_FOLDER or not os.path.isdir(OUTPUT_FOLDER):
        path = os.getcwd()
        OUTPUT_FOLDER = os.path.join(path, 'output_files')

        if not os.path.isdir(OUTPUT_FOLDER):
            os.mkdir(OUTPUT_FOLDER)

        log_file_input = os.path.join(OUTPUT_FOLDER, str(key + '_progress.log'))
        
        with open(log_file_input, "w+") as erase_file:
            erase_file.close()
        
        with open(log_file_input, "a+") as f:
            f.write("No output folder or invalid output folder specified;" + " Creating and/or using the following folder instead: " + str(OUTPUT_FOLDER))

    else:

        log_file_input = os.path.join(OUTPUT_FOLDER, str(key + '_progress.log'))

        with open(log_file_input, "w+") as erase_file:
            erase_file.close()

    maf_df = read_maf(maf_filename_input, log_file_input)
    full_maf_df = pd.read_csv(maf_filename_input, sep="\t")

    full_transcript_id = full_maf_df['Transcript_ID'].tolist()
    full_start_pos = full_maf_df['Start_Position'].tolist()
    full_chromosome_number = full_maf_df['Chromosome'].tolist()
    full_ref_allele = full_maf_df['Reference_Allele'].tolist()
    full_alt_allele = full_maf_df['Tumor_Seq_Allele2'].tolist()

    var_classification_list = maf_df['Variant_Classification'].tolist()
    transcript_id = maf_df['Transcript_ID'].tolist()
    start_pos = maf_df['Start_Position'].tolist()
    chromosome_number = maf_df['Chromosome'].tolist()
    ref_allele = maf_df['Reference_Allele'].tolist()
    alt_allele = maf_df['Tumor_Seq_Allele2'].tolist() 

    #Define and create output folder if not specified or if it doesn't exist

    output_maf = os.path.join(OUTPUT_FOLDER, str(key + ".maf"))

    with open(output_maf, "w+") as erase_file:
        erase_file.close()

    mutant_context_list, wt_context_list, transcript_list = exec_get_AA_context(transcript_id, chromosome_number, start_pos, ref_allele, alt_allele, log_file_input, var_classification_list)

    full_maf_rows = full_maf_df.values.tolist() 
    column_names = full_maf_df.columns.tolist()

    column_names.append("Amino_Acid_Context_Transcript")
    column_names.append("WT_Amino_Acid_Context")
    column_names.append("Variant_Amino_Acid_Context")

    new_maf_rows = []

    a=0
    b=0
    while a<len(transcript_id):
        if transcript_id[a] == full_transcript_id[b] and start_pos[a] == full_start_pos[b] and chromosome_number[a] == full_chromosome_number[b] and ref_allele[a] == full_ref_allele[b] and alt_allele[a] == full_alt_allele[b]:

            row = full_maf_rows[b]

            if not transcript_list[a]:
                row.append("ERROR: NO TRANSCRIPT OBTAINED")
            else:
                row.append(transcript_list[a])

            if not wt_context_list[a]:
                row.append("ERROR: NO WT CONTEXT OBTAINED")
            else:
                row.append(wt_context_list[a])

            if not mutant_context_list[a]:
                row.append("ERROR: NO VARIANT CONTEXT OBTAINED")
            else:
                row.append(mutant_context_list[a])

            new_maf_rows.append(row)
            a=a+1
            b=b+1

        else:            
            row = full_maf_rows[b]
            row.append("")
            row.append("")
            row.append("")

            new_maf_rows.append(row)
            b=b+1

    new_maf_rows = pd.DataFrame(new_maf_rows, columns=column_names)

    new_maf_rows.to_csv(output_maf, sep="\t", index=False)

    with open(log_file_input, "a+") as f:
            f.write("Success!")

    return output_maf

#main predictor function, only MHCflurry supported for now to enable pip installation

def predict(allele_annot, maf_filename_input, OUTPUT_FOLDER, key):

    tf.compat.v1.logging.set_verbosity(tf.compat.v1.logging.ERROR)

    key = str(key)

    #Define key if not specified

    if not key:
        key = datetime.now().strftime('%Y_%m_%d_%H_%M_%S')

    #Define and create output folder if not specified or if it doesn't exist

    if not OUTPUT_FOLDER or not os.path.isdir(OUTPUT_FOLDER):
        path = os.getcwd()
        OUTPUT_FOLDER = os.path.join(path, 'output_files')

        if not os.path.isdir(OUTPUT_FOLDER):
            os.mkdir(OUTPUT_FOLDER)

        log_file_input = os.path.join(OUTPUT_FOLDER, str(key + '_progress.log'))
        
        with open(log_file_input, "w+") as erase_file:
            erase_file.close()
        
        with open(log_file_input, "a+") as f:
            f.write("No output folder or invalid output folder specified; " + "Creating and/or using the following folder instead: " + str(OUTPUT_FOLDER))

    else:

        log_file_input = os.path.join(OUTPUT_FOLDER, str(key + '_progress.log'))

        with open(log_file_input, "w+") as erase_file:
            erase_file.close()

    maf_df = read_maf(maf_filename_input, log_file_input)

    var_classification_list = maf_df['Variant_Classification'].tolist()
    Hugo_Symbol_list = maf_df['Hugo_Symbol'].tolist()
    start_list = maf_df['Start_Position'].tolist()
    end_list = maf_df['End_Position'].tolist()
    chr_list = maf_df['Chromosome'].tolist()
    ta_1_list = maf_df['Tumor_Seq_Allele1'].tolist()
    ta_2_list = maf_df['Tumor_Seq_Allele2'].tolist()
    refa_list = maf_df['Reference_Allele'].tolist()
    barcode_list = maf_df['Tumor_Sample_Barcode'].tolist()

    transcript_id = maf_df['Transcript_ID'].tolist()
    start_pos = maf_df['Start_Position'].tolist()
    chromosome_number = maf_df['Chromosome'].tolist()
    ref_allele = maf_df['Reference_Allele'].tolist()
    alt_allele = maf_df['Tumor_Seq_Allele2'].tolist()

    output_csv = os.path.join(OUTPUT_FOLDER, str(key + '_output.csv'))

    with open(output_csv, "w+") as erase_file:
        erase_file.close()

    output_json = os.path.join(OUTPUT_FOLDER, str(key + '_output.json'))

    with open(output_json, "w+") as erase_file:
        erase_file.close()

    predictor = Class1PresentationPredictor.load() #Load MHCflurry predictor

    if not allele_annot or not os.path.isfile(allele_annot):

        with open(log_file_input, "a+") as f:
            f.write("No alleles inputted, or allele annotation file has an invalid path; 318 common HLA class I alleles used instead.\n")

        uq_barcodes = list(set(barcode_list))

        #Analysis will be performed on 318 common HLA alleles if they are not inputted - https://www.biorxiv.org/content/10.1101/2020.12.08.416271v1.full
        list_of_alleles=['B0801', 'C0701', 'B4402', 'A0101', 'C0704', 'A0301', 'B0702', 'C0702', 'C0501', 'B4501', 'B5601', 'C0602', 'A2601', 'A3002', 'A2301', 'B0705', 'B1801', 'C1505', 'B5701', 'A2902', 'B2702', 'C0202', 'A6801', 'A0201', 'B1501', 'C0304', 'B3901', 'C1203', 'A2402', 'B1402', 'A1101', 'C1202', 'B1302', 'B5201', 'C1601', 'B4403', 'A3001', 'B4002', 'C1502', 'B5101', 'C0803', 'B4801', 'B3801', 'B3502', 'B4901', 'B2705', 'A2501', 'A1102', 'A3303', 'B5801', 'C0302', 'B4601', 'C0102', 'B4101', 'A3101', 'C0303', 'B4001', 'C1402', 'B1550', 'C0713', 'B3501', 'C0401', 'A3201', 'C0802', 'B1401', 'C0801', 'A0205', 'B7801', 'B3503', 'A2608', 'B5001', 'B3701', 'A6802', 'A3004', 'B1535', 'B3802', 'A3401', 'A3402', 'C1509', 'B4102', 'A6601', 'B5501', 'A2403', 'B1803', 'A3301', 'A0203', 'B5502', 'B5802', 'B1510', 'B5301', 'B1525', 'B1301', 'A0302', 'B1502', 'A2901', 'C1604', 'B4006', 'A0207', 'B1517', 'A7401', 'B1512', 'B4701', 'B5105', 'B8101', 'B1516', 'B3908', 'C0717', 'B2706', 'A2407', 'B4405', 'A0206', 'A0202', 'B3905', 'C0403', 'B4201', 'B1518', 'B1503', 'B1527', 'C0210', 'B5401', 'B5703', 'B3505', 'B4803', 'B1511', 'B1825', 'B4104', 'C0718', 'A6901', 'B1507', 'A0238', 'C0381', 'B3906', 'B3530', 'B4404', 'A3102', 'B5108', 'B1504', 'A8001', 'B3508', 'B4413', 'A0211', 'A0216', 'B4202', 'B3902', 'A2431', 'B3910', 'B3543', 'B3577', 'A3305', 'B5504', 'A0102', 'A3601', 'C0804', 'C1504', 'B4406', 'A7409', 'B3913', 'B3504', 'A2305', 'A0222', 'C0722', 'B5002', 'B5704', 'A0227', 'B4409', 'A0229', 'B1509', 'B0704', 'B1403', 'B5129', 'B1524', 'B7301', 'B0809', 'B2704', 'B2707', 'B8201', 'A6808', 'B4410', 'B3924', 'C0617', 'A0103', 'C0305', 'B5005', 'B4408', 'A0325', 'B1505', 'B15125', 'B3513', 'B4805', 'C0746', 'B7802', 'C0608', 'A0230', 'C0606', 'A0312', 'C0795', 'C1602', 'A1110', 'B5806', 'A0138', 'B4407', 'A2910', 'C0813', 'C0130', 'B5302', 'B2709', 'A2423', 'B2714', 'A6824', 'A6602', 'A3603', 'B0712', 'B1531', 'B5702', 'B4012', 'B5109', 'B1805', 'B1539', 'B0709', 'A2614', 'A0214', 'B5604', 'C0317', 'A2410', 'B1521', 'B1513', 'A1104', 'B5102', 'B1508', 'A0217', 'A2417', 'B3509', 'B4802', 'B1807', 'A1103', 'A6806', 'B2703', 'B0722', 'A0106', 'A0220', 'B3931', 'B4036', 'A2602', 'B6701', 'B4003', 'B4040', 'C1403', 'B5605', 'B1553', 'B4021', 'B5901', 'A0235', 'B0710', 'C0306', 'A3106', 'A6803', 'A6805', 'B3517', 'B1804', 'C1212', 'B1811', 'B3512', 'B4005', 'C0509', 'A0224', 'A3405', 'A6603', 'C0705', 'A0204', 'A7411', 'B3911', 'B1530', 'B5114', 'C0727', 'B1534', 'A2420', 'B4103', 'C1701', 'B3522', 'B1806', 'C0307', 'B2708', 'B4004', 'B1520', 'A0267', 'B3912', 'B3909', 'B3511', 'A7403', 'B5004', 'B1537', 'B1573', 'C0117', 'B5106', 'C1204', 'B3523', 'A2615', 'B1568', 'C0406', 'A2304', 'A3213', 'B3805', 'B15135', 'A0219', 'B0812', 'B1554', 'B3524', 'B5710', 'B5602', 'A2414', 'C0308', 'B1808', 'A0221', 'C0780', 'B0713', 'A3008', 'A0260', 'C0815']

        double_list_of_alleles = []

        a=0
        while a<len(uq_barcodes):
            double_list_of_alleles.append(list_of_alleles)
            a=a+1

        patient_alleles_dict = dict(zip(uq_barcodes, double_list_of_alleles))

    else:

        #Create a dictionary from the allele annotation file
        #Note: Alleles must be supported with MHCflurry, format: HLA-A*01:01

        df_allele = pd.read_csv(allele_annot)
        df_allele_columns = df_allele.columns.tolist()
        sample_ids = df_allele[df_allele_columns[0]].tolist()
        df_allele_no_sample = df_allele.drop(df_allele_columns[0], 1)
        df_allele_rows = df_allele_no_sample.values.tolist()

        index_to_remove = []

        #Remove samples with non-supported alleles

        all_alleles = predictor.supported_alleles

        

        a=0
        while a<len(df_allele_rows):
            if any(element not in all_alleles for element in df_allele_rows[a]):
                index_to_remove.append(a)
            a=a+1

        index_to_remove = sorted(index_to_remove, reverse=True)

        removed_samples = []

        for i in index_to_remove:
            removed_samples.append(sample_ids[i])
            del sample_ids[i]
            del df_allele_rows[i]

        with open(log_file_input, "a+") as log_file:
            i=0
            log_file.write("The following samples were removed due to allele invalidity: ")
            while i<len(removed_samples):
                log_file.write(str(removed_samples[i]) + "\t")
                i=i+1
            log_file.write("\n")

        patient_alleles_dict = dict(zip(sample_ids, df_allele_rows))

    affinity_cutoff = 50000

    save_cycle = 5 #number of iterations before saving

    mutant_context_list, wt_context_list, transcript_list = exec_get_AA_context(transcript_id, chromosome_number, start_pos, ref_allele, alt_allele, log_file_input, var_classification_list)

    d=0

    df_columns=['Tumor_Sample_Barcode', 'Hugo_Symbol', 'Chromosome', 'Start_Position', 'End_Position', 'Reference_Allele', 'Tumor_Seq_Allele1', 'Tumor_Seq_Allele2', 'Transcript', 'Presentation_HBR', 'Affinity_HBR', 'WT_Amino_Acid_Context', 'Variant_Amino_Acid_Context', 'HLA_Allele', 'N_total_neoantigens_500nM_cutoff', 'N_strong_binders_50nM_cutoff', 'Best_Presentation_Score_Peptide', 'Best_Presentation_Score_Presentation_Score', 'Best_Presentation_Score_Presentation_Percentile', 'Best_Presentation_Score_Affinity', 'Best_Presentation_Score_Affinity_Percentile', 'Best_Affinity_Peptide', 'Best_Affinity_Presentation_Score', 'Best_Affinity_Presentation_Percentile', 'Best_Affinity_Affinity', 'Best_Affinity_Affinity_Percentile']
    master_df = pd.DataFrame(columns=df_columns)

    final_json_list = []

    while d<len(Hugo_Symbol_list):
        with open(log_file_input, "a+") as log_file:
            log_file.write("Starting prediction for variant: " + str(d+1) + " of " + str(len(Hugo_Symbol_list)) + "\n")

        current_mutated_context = mutant_context_list[d]
        current_wt_context = wt_context_list[d]
        current_transcript_tested = transcript_list[d]
        current_Hugo_Symbol = Hugo_Symbol_list[d]
        current_start = start_list[d]
        current_end = end_list[d]
        current_chr = chr_list[d]
        current_ta_1 = ta_1_list[d]
        current_ta_2 = ta_2_list[d]
        current_refa = refa_list[d]
        current_barcode = barcode_list[d]

        if isinstance(current_mutated_context, list) and isinstance(current_wt_context, list) and isinstance(current_transcript_tested, list):
            
            r=0
            while r<len(current_mutated_context):
                master_df, final_json_list = output_result(predictor, affinity_cutoff, final_json_list, patient_alleles_dict, log_file_input, df_columns, current_barcode, current_Hugo_Symbol, current_chr, current_start, current_end, current_refa, current_ta_1, current_ta_2, current_transcript_tested[r], current_wt_context[r], current_mutated_context[r], master_df)
                r=r+1

        else: 
            master_df, final_json_list = output_result(predictor, affinity_cutoff, final_json_list, patient_alleles_dict, log_file_input, df_columns, current_barcode, current_Hugo_Symbol, current_chr, current_start, current_end, current_refa, current_ta_1, current_ta_2, current_transcript_tested, current_wt_context, current_mutated_context, master_df)

        if d % save_cycle == 0:
            master_df.to_csv(output_csv, index=False)

        d=d+1

    json_output_file = open(output_json, "a+")
    
    for x in final_json_list:
        json_output_file.write(str(x) + "\n")
    
    json_output_file.close()

    master_df.to_csv(output_csv, index=False)

    with open(log_file_input, "a+") as log_file:
        log_file.write("Success!")

    return output_json, output_csv
