#!/bin/bash

# extract all tax ids for bacterial species rank (Bacteria taxid=2)
grep 'B' ../db/categories.dmp | cut -d$'\t' -f2 | sort | uniq > ../db/2.species.taxids

# obtain all hits with E-value < 1E-3 from NCBI nr_v5 database restricted to bacterial species (2.taxids)
blastp -db $HOME/blast/db/nr_v5/nr_v5 -query ../data/query.fa -taxidlist ../db/2.species.taxids -out ../data/2.species.hits.tab -evalue 1e-3 -max_target_seqs 100000 -outfmt "7 qstart qend stitle sacc staxids evalue pident bitscore sstart send sseq" -num_threads 4
