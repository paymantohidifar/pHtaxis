The taxcat dump compressed archive(s) contains a single file -

  categories.dmp

categories.dmp contains a single line for each node that is at or below the 
species level in the NCBI Taxonomy database.

The first column is the top-level category -

  A = Archaea
  B = Bacteria
  E = Eukaryota
  V = Viruses and Viroids
  U = Unclassified
  O = Other

The third column is the taxid itself, and the second column is the 
corresponding species-level taxid.

These nodes in the taxonomy -

  242703 - Acidilobus saccharovorans
  666510 - Acidilobus saccharovorans 345-15 

will appear in categories.dmp as -

A       242703  242703
A       242703  666510

