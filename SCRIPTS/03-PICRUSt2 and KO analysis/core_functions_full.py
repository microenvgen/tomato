#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Last version: 5 June 2023

@description:
    This script generates a list of core functions for a given phylogenetic
    core group (PCG). This list is considered that PCG's pangenome, defined
    as annotations present in at least 90% of available genomes for each PCG.
    The input for this tool is a path containing .tsv annotation files
    generated with eggNOG-mapper, for which a consensus annotation file will
    be generated. The core reactions can be retrieved in KEGG, EC and in BiGG
    format. A file containing descriptions corresponding to EC and KEGG entries
    is also generated.
"""

# Setup =======================================================================

## modules
import os
# from carveme.reconstruction.eggnog import load_eggnog_data
from utils_core_functions import split_and_expand, load_eggnog_data, consensus_eggnog
from carveme.reconstruction.scoring import *
import pandas as pd
import csv
import Bio.KEGG.REST


wkdir = "/home/silvia/repos/get_core_reactions"
os.chdir(wkdir)

##nodes = ["Node13286", "Node15327", "Node17416", "Node21165", "Node45983", "Node54724", "Node55822", "Node70899", "Node71253", "Node101358", "Node115604", "Node116106", "Node119364"] #> 1%
##nodes = ["Node15327", "Node17416", "Node21165", "Node45983", "Node54724", "Node55822", "Node119364", "Node137058"] #> 0.5% + un octavo pcg
nodes = ["Node15327", "Node19398", "Node21166", "Node28120", "Node45985", "Node52227", "Node53994", "Node55420", "Node55904", "Node119365", "Node126898", "Node147305", "Node17416", "Node45983", "Node54724", "Node55822", "Node119364"]

for node in nodes:
    
    ## variables ==============================================================
    # annotations_folder = "/home/silvia/annots_pangenomes_2022_12_14/" + node + "_annots" #> 1% y creo que 137058
    annotations_folder = "/home/silvia/annotated_genomes_per_pcg_COHERENT/" + node + "/annotated_genomes" #> all coherent nodes 1% and 0.5%
    ## ORIGIN:
    ## /home/silviatm/micro/TFM/annotated_genomes_per_pcg/Node*/annotated_genomes
    ## launch_consenso_nodos; /home/silviatm/micro/TFM/all_0.01_nodes_plus_one.csv
    perc = 0.90
    pangenome_folder = "pangenomes" + str(perc) # 90%
    if not os.path.exists(pangenome_folder):
        os.makedirs(pangenome_folder)
        
    gprs_bigg = "/home/silvia/anaconda3/envs/carveme/lib/python3.7/site-packages/carveme/data/generated/bigg_gprs.csv.gz"
    gprs_bigg = "/home/silvia/bigg_gprs.csv.gz"
    ko_list = "./ko"
    enzyme_dat = "./EC"

    # output ==================================================================
    output_folder = "./analisis_pans" + str(100*perc) + "/anotaciones_" + str(100*perc)
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    descriptions_output = output_folder+"/"+node+'_descriptions.txt'
    bigg_reac_output    = output_folder+"/"+node+'_bigg_reac.txt'
    kegg_reac_output    = output_folder+"/"+node+'_kegg_reac.txt'
    kegg_ko_output      = output_folder+"/"+node+'_kegg_ko.txt'
    
    # Obtain consensus annotations ============================================
    consensus_eggnog(annotations_folder, pangenome_folder, outputname="pan" + node, perc=perc)
    consensus_file = pangenome_folder+"/pan" + node + ".tsv"  # input from previous line
    
    # Initial parsing =========================================================
    # Parse consensus annotations
    annotation = load_eggnog_data(consensus_file, drop_unused_cols=False, sep="\t")
    
    # take the annotations and filter best match for each gene
    gene2gene = annotation.sort_values(by='score', ascending=False) \
    	              .groupby('BiGG_gene', as_index=False).apply(lambda x: x.iloc[0])
    
    # drop duplicates (keep first==best score!)
    bef = len(gene2gene)
    gene2gene = gene2gene.drop_duplicates("query_gene", keep="first")
    print("Removed "+str(bef - len(gene2gene))+"/"+str(bef)+" duplicated entries ('query_gene') from the annotations table")
    
    # Parse BiGG database
    gprs = pd.read_csv(gprs_bigg)
    gprs['BiGG_gene'] = gprs.apply(lambda row: f"{row['model']}.{row['gene'][2:]}", axis=1)
    
    # save available BiGG_Reactions ===========================================
        # This way we can link each annotation to a BiGG reaction
        # en gprs; aquí fusiono columnas en base a su col en común, bigg_gene
        # 'left' significa que me quedo solo con las filas presentes en gene2gene
        # (annots) y añado lo que pueda de la tabla de BiGG. "inner" == solo para las
        # que pueda aportar algo de la tabla de bigg.
        # -> the length of the table can grow: some annotations might have multiple
        # reactions
    df     = pd.merge(gene2gene, gprs, how='left')
    # df=gene_scores_in_gprs.drop_duplicates("query_gene")
    # print("Further removed "+str(bef - len(df))+"/"+str(bef)+" duplicated entries after merging with BiGG GPRs table")
        # Duplicates might actually be two reactions encoded in the same annotation.
        # Like R_MCTP1App and R_MCTP2App
    with open(bigg_reac_output, 'w') as the_file:
        for r in df["reaction"]:
            if not pd.isna(r):
                the_file.write(r+"\n")
    print("Saved BiGG reactions to "+bigg_reac_output+"\n")

    
    
    # save available ECs ======================================================
    my_IDs=list(df[df["EC"].notnull()]["EC"])
    newList=[]
    for i in my_IDs:
        newList.append(i.split(',')[0])
    my_IDs = newList
    
    print("Retrieved ECs: "+str(len(my_IDs))+" (genes with no ECs: "+str(len(gene2gene[gene2gene["EC"].isnull()]["EC"]))+")")
    
    # save available EC descriptions in human language ========================
    import Bio.ExPASy.Enzyme
    handle = open(enzyme_dat)
    descriptions = {record.get("ID"):record for record in Bio.ExPASy.Enzyme.parse(handle)}
    
    with open(descriptions_output, 'w') as the_file:
        for id_ in my_IDs:
            the_file.write(descriptions[id_].get("DE")+"\n")
    print("Saved reaction descriptions to "+descriptions_output+"\n")
    
    # save some more descriptions from KOs
    # http://rest.kegg.jp/list/ko
    with open(ko_list, mode='r') as infile:
        reader = csv.reader(infile,delimiter="\t")
        mydict = {rows[0]:rows[1] for rows in reader}
    
    with open(descriptions_output, 'a') as the_file:
        for keggkos in list(gene2gene[gene2gene["EC"].isnull()]["KEGG_ko"]):
            try:
                for keggko in keggkos.split(","):
                    the_file.write(mydict[keggko]+"\n")
            except:
                print("WARNING: keggkos.split(',') failed. This is keggkos: "+str(keggkos))
    
    print("Saved additional reaction descriptions to "+descriptions_output+"\n")
    
    
    # save all possible kegg reactions ========================================
    with open(kegg_reac_output, 'w') as the_file:
        for rs in df["KEGG_Reaction"]:
            if not pd.isna(rs):
        	    for r in rs.split(","):
        	        ## reactions
        	        the_file.write(r+"\n")
        	        print("Saved available KEGG_Reactions to "+kegg_reac_output)
        	        ## compounds
        ##                kegg_compounds = Bio.KEGG.REST.kegg_link("compound", r)
        ##                [the_file.write(c.split("\tcpd:")[-1].strip()+"\n") for c in set(kegg_compounds)]
        	    
    
    
    # save KOs ================================================================
    with open(kegg_ko_output, 'w') as the_file:
    #        for rs in list(gene2gene[gene2gene["KEGG_Reaction"].isnull()]["KEGG_ko"]):
        for rs in list(gene2gene[gene2gene["KEGG_ko"].notnull()]["KEGG_ko"]):
            for r in rs.split(","):
                the_file.write(r.split("ko:")[1]+"\n")
            
    print("Saved available KEGG_ko to "+kegg_ko_output)
