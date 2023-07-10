# ConSequences

ConSequences is an R script designed to perform semi-automatic curation of NCBI GenBank data. The data after curation with ConSequences is ready for analysis in Modeltest-NG and RAxML-NG.

**Functions available in ConSequences script:**

- change.seq.names - automatically changes sequences names based on their accession numbers, when a table of sequence metadata is provided (see an examplary sequence metadata table - Umbelopsis.csv)
- fill.with.pluses - fills the missing molecular marker data for each strain with pluses based on a table of sequence metadata
- concatenate.sequences - concatenates molecular data for each strain

Functions change.seq.names and concatenate.sequences perform additional verification of your data. They print a report if the sequences data or metadata table are incomplete. 

Additional information about functions usage and the code can be found in ConSequences.R file.

# Files

ConSequences.R - contains thoroughly described script with functions for data curation, an examplary code to use with a test data provided and author's data necessary for citation.

Test - folder containing test data:
- Umbelopsis.csv - an examplary table containing sequence metadata to use with the test data
- act1.txt, COI.txt, ITS.txt, LSU.txt, mcm7.txt, SSU.txt - an examplary accession numbers files used to download test data via BatchEntrez on the GenBank website.

# Citation

If you use this code in Your work, please cite this repository along with the script name. Author's data necessary for citation are available in ConSequences.R file.
