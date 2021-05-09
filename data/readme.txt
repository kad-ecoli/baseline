Folder for preprocessed training data.

In naive.*, the columns are
[1] GO term
[2] background probability
[3] information content, regardless of parent
[4] information content, conditioned on parents
[5] GO term name

For items [2][3][4], probabilities are calculated with a pseudocount of 1.
For items [3][4], natural log instead of log2 is used.
The use of pseudocount and natural log is for consistency with the official
CAFA implementation at
https://github.com/yuxjiang/CAFA2/blob/master/matlab/pfp_eia.m
and is different from that defined in the original paper:
Clark, Wyatt T., and Predrag Radivojac. "Information-theoretic evaluation of
predicted ontological annotations." Bioinformatics 29.13 (2013): i53-i61.
