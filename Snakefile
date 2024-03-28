rule all:


rule tree_assignment:
    input:
        "09_ref_tree/ref_{gene}_tree.done"
    output:
        "09_ref_tree/taxa_{gene}_assignment.done"
    shell:
        """
        cp RAxML_bipartitions.ref_{wildcards.gene}_db.out 09_ref_tree/RAxML_bipartitions.ref_{wildcards.gene}_db.tre
        module purge
        module load miniconda/22.11.1-1
        conda activate raxml_env
        module load trimal/1.4.1
        sed -i "s/&//g" 09_ref_tree/ref_{wildcards.gene}_db_aligned.fa
        sed -i "s/&//g" 09_ref_tree/RAxML_bipartitions.ref_{wildcards.gene}_db.tre
        trimal -in 09_ref_tree/ref_{wildcards.gene}_db_aligned.fa -out 09_ref_tree/ref_{wildcards.gene}_db_aligned.phylip -phylip
        cd 09_ref_tree
        mkdir -p {wildcards.gene} && cd {wildcards.gene}
        module load papara_nt/2.5
        papara -t ../RAxML_bipartitions.ref_{wildcards.gene}_db.tre -s ../ref_{wildcards.gene}_db_aligned.phylip  -q ../../asv_sequences/foram_seqs.fa -j 3 -r
        trimal -in papara_alignment.default -out papara_alignment.fa -fasta
        cut -d ' ' -f 1 < ../ref_{wildcards.gene}_db_aligned.fa > ../ref_{wildcards.gene}_db_aligned_namesFixed.fa
        sed -i "s/&//g" ../ref_{wildcards.gene}_db_aligned_namesFixed.fa
        sed -i "s/&//g" ../RAxML_bipartitions.ref_{wildcards.gene}_db.tre
        epa-ng -t ../RAxML_bipartitions.ref_{wildcards.gene}_db.tre -s ../ref_{wildcards.gene}_db_aligned_namesFixed.fa \
          -q papara_alignment.fa -m GTR+G --redo --filter-max 7 --filter-min-lwr 0.9
        sed $'s/ [^\t]*//' ../../database/foram_{wildcards.gene}_taxa_tabbed.tsv > ../../database/foram_{wildcards.gene}_taxa_tabbed_clean.tsv

        # at the genus level
        mkdir -p genus && cd genus
        awk -F'\t' '{{print $1,$7}}' OFS='\t' ../../../database/foram_{wildcards.gene}_taxa_tabbed_clean.tsv > ../../../database/foram_{wildcards.gene}_taxa_tabbed_genus.tsv
        mkdir -p sequences
        # Generate some stats and analysis
        gappa examine heat-tree --jplace-path ../epa_result.jplace --mass-norm absolute --write-svg-tree --write-newick-tree --allow-file-overwriting
        gappa prepare extract --clade-list-file ../../../database/foram_{wildcards.gene}_taxa_tabbed_genus.tsv --fasta-path ../papara_alignment.fa \
          --color-tree-file extract_tree --samples-out-dir samples --sequences-out-dir sequences --jplace-path ../epa_result.jplace --allow-file-overwriting
        gappa examine assign --jplace-path ../epa_result.jplace --taxon-file ../../../database/foram_{wildcards.gene}_taxa_tabbed_genus.tsv --per-query-results --allow-file-overwriting
        gappa examine lwr --jplace-path ../epa_result.jplace --allow-file-overwriting
        #gappa examine edpl --jplace-path ../epa_result.jplace --allow-file-overwriting
        input_dir="sequences/"
        output_dir="assignment"
        mkdir -p "$output_dir"
        for input_file in "$input_dir"/*; do
          if [ -f "$input_file" ]; then
            grep_output=$(grep '>' "$input_file")        
            file_name=$(basename "$input_file")        
            echo -e "$grep_output" | awk -v file="$file_name" '{{print $0 "\t" file}}' > "$output_dir/${{file_name}}_output.tsv"        
            echo "Tab-delimited file created with grep output and filename for $file_name.";     
          fi;
        done
        cat assignment/* > combined_assignment.tsv
        cd ..

        # at the family level
        mkdir -p family && cd family
        awk -F'\t' '{{print $1,$6}}' OFS='\t' ../../../database/foram_{wildcards.gene}_taxa_tabbed_clean.tsv > ../../../database/foram_{wildcards.gene}_taxa_tabbed_family.tsv
        mkdir -p sequences
        # Generate some stats and analysis
        gappa examine heat-tree --jplace-path ../epa_result.jplace --mass-norm absolute --write-svg-tree --write-newick-tree --allow-file-overwriting
        gappa prepare extract --clade-list-file ../../../database/foram_{wildcards.gene}_taxa_tabbed_family.tsv --fasta-path ../papara_alignment.fa \
          --color-tree-file extract_tree --samples-out-dir samples --sequences-out-dir sequences --jplace-path ../epa_result.jplace --allow-file-overwriting
        gappa examine assign --jplace-path ../epa_result.jplace --taxon-file ../../../database/foram_{wildcards.gene}_taxa_tabbed_family.tsv --per-query-results --allow-file-overwriting
        gappa examine lwr --jplace-path ../epa_result.jplace --allow-file-overwriting
        #gappa examine edpl --jplace-path ../epa_result.jplace --allow-file-overwriting
        input_dir="sequences/"
        output_dir="assignment"
        mkdir -p "$output_dir"
        for input_file in "$input_dir"/*; do
          if [ -f "$input_file" ]; then
            grep_output=$(grep '>' "$input_file")
            file_name=$(basename "$input_file")
            echo -e "$grep_output" | awk -v file="$file_name" '{{print $0 "\t" file}}' > "$output_dir/${{file_name}}_output.tsv"
            echo "Tab-delimited file created with grep output and filename for $file_name.";
          fi;
        done
        cat assignment/* > combined_assignment.tsv
        cd ..

        # at the order level
        mkdir -p order && cd order
        awk -F'\t' '{{print $1,$5}}' OFS='\t' ../../../database/foram_{wildcards.gene}_taxa_tabbed_clean.tsv > ../../../database/foram_{wildcards.gene}_taxa_tabbed_order.tsv
        mkdir -p sequences
        # Generate some stats and analysis
        gappa examine heat-tree --jplace-path ../epa_result.jplace --mass-norm absolute --write-svg-tree --write-newick-tree --allow-file-overwriting
        gappa prepare extract --clade-list-file ../../../database/foram_{wildcards.gene}_taxa_tabbed_order.tsv --fasta-path ../papara_alignment.fa \
          --color-tree-file extract_tree --samples-out-dir samples --sequences-out-dir sequences --jplace-path ../epa_result.jplace --allow-file-overwriting
        gappa examine assign --jplace-path ../epa_result.jplace --taxon-file ../../../database/foram_{wildcards.gene}_taxa_tabbed_order.tsv --per-query-results --allow-file-overwriting
        gappa examine lwr --jplace-path ../epa_result.jplace --allow-file-overwriting
        #gappa examine edpl --jplace-path ../epa_result.jplace --allow-file-overwriting
        input_dir="sequences/"
        output_dir="assignment"
        mkdir -p "$output_dir"
        for input_file in "$input_dir"/*; do
          if [ -f "$input_file" ]; then
            grep_output=$(grep '>' "$input_file")
            file_name=$(basename "$input_file")
            echo -e "$grep_output" | awk -v file="$file_name" '{{print $0 "\t" file}}' > "$output_dir/${{file_name}}_output.tsv"
            echo "Tab-delimited file created with grep output and filename for $file_name.";
          fi;
        done
        cat assignment/* > combined_assignment.tsv
        cd ..

        # at the class level
        mkdir -p class && cd class
        awk -F'\t' '{{print $1,$7}}' OFS='\t' ../../../database/foram_{wildcards.gene}_taxa_tabbed_clean.tsv > ../../../database/foram_{wildcards.gene}_taxa_tabbed_class.tsv
        mkdir -p sequences
        # Generate some stats and analysis
        gappa examine heat-tree --jplace-path ../epa_result.jplace --mass-norm absolute --write-svg-tree --write-newick-tree --allow-file-overwriting
        gappa prepare extract --clade-list-file ../../../database/foram_{wildcards.gene}_taxa_tabbed_class.tsv --fasta-path ../papara_alignment.fa \
          --color-tree-file extract_tree --samples-out-dir samples --sequences-out-dir sequences --jplace-path ../epa_result.jplace --allow-file-overwriting
        gappa examine assign --jplace-path ../epa_result.jplace --taxon-file ../../../database/foram_{wildcards.gene}_taxa_tabbed_class.tsv --per-query-results --allow-file-overwriting
        gappa examine lwr --jplace-path ../epa_result.jplace --allow-file-overwriting
        #gappa examine edpl --jplace-path ../epa_result.jplace --allow-file-overwriting
        input_dir="sequences/"
        output_dir="assignment"
        mkdir -p "$output_dir"
        for input_file in "$input_dir"/*; do
          if [ -f "$input_file" ]; then
            grep_output=$(grep '>' "$input_file")
            file_name=$(basename "$input_file")
            echo -e "$grep_output" | awk -v file="$file_name" '{{print $0 "\t" file}}' > "$output_dir/${{file_name}}_output.tsv"
            echo "Tab-delimited file created with grep output and filename for $file_name.";
          fi;
        done
        cat assignment/* > combined_assignment.tsv
        cd ..
        cd ../..
        touch {output}
        """
