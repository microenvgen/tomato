# The entire PICRUSt2 pipeline can be run using a single script, called
# picrust2_pipeline.py. This script will run each of the 4 key steps:
#   (1) sequence placement,
#   (2) hidden-state prediction of genomes,
#   (3) metagenome prediction,
#   (4) pathway-level predictions
# The nearest-sequenced taxon index (NSTI) will be calculated for each input ASV 
# and by default any ASVs with NSTI > 2 will be excluded from the output

PROCESSES=4
BIOM=tomate_bc_0.01/Tree/0.99/All_table_0.99.biom # same for all cutoffs!
FASTAFOLDER=tomate_bc_0.01/Tree/0.99/
FASTA=rep_set.fasta
PICRUSTOUTPUT=metagenome_predictions_picrust2
FINALOUTPUT=results_picrust2
OUTTABLE=$FINALOUTPUT"/METAGENOME_TABLE_PI2.tsv.gz"

#> Fix headers
awk '{print $1}' $FASTAFOLDER"/"$FASTA > $FASTAFOLDER"/"OTUs$FASTA

#> For me, this hung. So I made some modifications
# - installed 4.5.2 for dendropy.
# - edited util.py (call())
# important to use sepp to save RAM
picrust2_pipeline.py -s $FASTAFOLDER"/"OTUs$FASTA \
-i $BIOM \
-o $PICRUSTOUTPUT \
--max_nsti 2 \
-p $PROCESSES \
-t sepp \
--skip_minpath \
--verbose

#> With Picrust1, you obtained a biom table (that you could later transform to
#> .tsv) with predict_metagenomes.py. With Picrust2, this table is also in the
#> output folder. But we get 2: EC predicted and KO predicted. We want the KOs.
mkdir -p $FINALOUTPUT/METAG_BY_CORE
cp $PICRUSTOUTPUT"/KO_metagenome_out/pred_metagenome_unstrat.tsv.gz" $OUTTABLE

# The key output files are:
# 
#     > EC_metagenome_out - Folder containing unstratified EC number metagenome predictions
# (pred_metagenome_unstrat.tsv.gz), sequence table normalized by predicted 16S copy number
# abundances (seqtab_norm.tsv.gz), and the per-sample NSTI values weighted by the abundance of
# each ASV (weighted_nsti.tsv.gz).
# 
#     > KO_metagenome_out - As EC_metagenome_out above, but for KO metagenomes.
# 
#     > pathways_out - Folder containing predicted pathway abundances and coverages per-sample,
# based on predicted EC number abundances.

# Additional output files, which can be useful for advanced users are:
#   
#     > EC_predicted.tsv.gz - Predicted EC number copy numbers per ASV.
# 
#     > intermediate - Folder containing intermediate MinPath files and files used for sequence placement pipeline (including JPLACE file: intermediate/place_seqs/epa_out/epa_result.jplace).
# 
#     > KO_predicted.tsv.gz - As EC_predicted.tsv.gz above, but for KO predictions.
# 
#     > marker_nsti_predicted.tsv.gz - Predicted 16S copy numbers and NSTI per ASV.
# 
#     > out.tre - Tree of reference and study 16S sequences.
