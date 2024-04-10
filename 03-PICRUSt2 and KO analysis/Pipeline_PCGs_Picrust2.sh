#!/bin/bash
#conda activate picrust
#unset PYTHONPATH

pcg_table=tomate_bc_0.005/Tree/results.txt
#pcg_table=results_0005_coherent.txt

PROCESSES=4
BIOM=tomate_bc_0.01/Tree/0.99/All_table_0.99.biom # same for all cutoffs!
FASTAFOLDER=tomate_bc_0.01/Tree/0.99/
FASTA=rep_set.fasta
PICRUSTOUTPUT=metagenome_predictions_picrust2_PCG_all20
FINALOUTPUT=results_picrust2_PCG_all20
OUTTABLE="METAGENOME_TABLE_PI2.tsv.gz"

#pcg_table=tomate_bc_0.005/Tree/results.txt # 20 pcgs 0.005, coherent or not
#PICRUSTOUTPUT=metagenome_predictions_picrust2_PCG_all20 # 20 pcgs 0.005, coherent or not
#FINALOUTPUT=results_picrust2_PCG_all20 # 20 pcgs 0.005, coherent or not

mkdir -p $PICRUSTOUTPUT # we create the parent folder. it will contain multiple core-named folders
mkdir -p $FINALOUTPUT/METAG_BY_CORE

#> Fix headers (already done at Pipeline_Picrust2.sh)
awk '{print $1}' $FASTAFOLDER"/"$FASTA > $FASTAFOLDER"/"OTUs$FASTA

# For each PCG, compute the minimal_metagenome
tail -n +2 "$pcg_table" | while IFS='\t' read -r core
do
  Core=$(echo "$core" | awk -F'\t' '{print $1}')
  # if it has been done already, skip this PCG
  if [ ! -f "$FINALOUTPUT"/METAG_BY_CORE/"$Core"_"$OUTTABLE" ]; then
    Leaves=$(echo "$core" | awk -F'\t' '{print $9}')
    if [ -n "$Leaves" ]; then
      # Create a temporary file for the current row
      temp_file=$(mktemp)

      # Write the leaves to the temporary file
      echo "$Leaves" | tr ';' '\n' | head -n -1 > "$temp_file"

      # Subset the full OTU abundance table and normalize
      biom subset-table -i "$BIOM" -a observation -s "$temp_file" -o "$PICRUSTOUTPUT"/"$Core"_table.biom

      # With PICRUSt2 we do not need to normalize first
      picrust2_pipeline.py -s $FASTAFOLDER"/OTUs"$FASTA \
      -i $PICRUSTOUTPUT"/"$Core"_table.biom" \
      -o $PICRUSTOUTPUT"/"$Core \
      --max_nsti 2 \
      -p $PROCESSES \
      -t sepp \
      --skip_minpath \
      --verbose \
      --stratified
      
      
      #> With Picrust1, you obtained a biom table (that you could later transform to
      #> .tsv) with predict_metagenomes.py. With Picrust2, this table is also in the
      #> output folder. But we get 2: EC predicted and KO predicted. We want the KOs.
      mkdir -p $FINALOUTPUT/METAG_BY_CORE
      cp $PICRUSTOUTPUT"/"$Core"/KO_metagenome_out/pred_metagenome_unstrat.tsv.gz" $FINALOUTPUT"/METAG_BY_CORE/"$Core"_"$OUTTABLE

      # Convert to table, convert into reactions
      # biom convert -i "$output_folder"/"$Core"_metagenome_predictions.biom -o "$output_folder"/"$Core"_METAGENOME_TABLE.txt --to-tsv

      # Clean up the temporary file
      rm "$temp_file"
    fi
  fi
done
