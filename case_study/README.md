## Case study ##

As a case study of the impact of different blast scoring, we check CAFA3 NK target
T37020011939 (UniProt ID [SCRK3_ARATH](https://www.uniprot.org/uniprot/Q9LNE4)),
which is a fructokinase found in the cytosol and extracellular region of Arabidopsis thaliana. As a probably fructokinase-3, this protein transfers the phosphate group from ATP to D-fructose. 

```
>T37020011939 SCRK3_ARATH
MASSTGEKGLIVSFGEMLIDFVPTVSGVSLSESPGFLKAPGGAPANVAIAVSRLGGRAAFVGKLGDDDFGHMLAGILRKN
GVDDQGINFDEGARTALAFVTLRSDGEREFMFYRNPSADMLLRPDELNLELIRSAKVFHYGSISLITEPCRSAHMKAMEV
AKEAGALLSYDPNLREPLWPSPEEARTQIMSIWDKADIIKVSDVELEFLTENKTMDDKTAMSLWHPNLKLLLVTLGEKGC
TYFTKKFHGSVETFHVDAVDTTGAGDSFVGALLQQIVDDQSVLEDEARLRKVLRFANACGAITTTKKGAIPALPTDIEAL
SFLKDQKKRQTNLKFSKWCCTASPC
```

The CAFA3 ground truth leaf terms for this target are:
GO:0005829 "cytosol" for MF, GO:0006000 "fructose metabolic process" for BP, and
GO:0008865 "fructokinase activity" for MF.

While all blast baseline methods fail to predict the leaf MF term GO:0008865, 
frequency-based blast baseline (blastfreq) can predict its parent MF term
GO:0016772 "transferase activity, transferring phosphorus-containing groups"
with a high probability of 0.89, while the official local identity-based
blast baseline (blastlocalID) only predict the term with probability of 0.32.
This is because 24 out of all 27 MF templates found by blast are annotated
with GO:0016772, even though all these templates are only weakly similar to
the target sequence.

Moreover, blastlocalID incorrectly predict that the CC of this protein is
GO:0009570 "chloroplast stroma" with a high probability of 0.64, respectively. 
This prediction is based on a single template
[Q9C524](https://www.uniprot.org/uniprot/Q9C524), which is the only template
with this term among all 26 CC templates. blastfreq on the other hand only
assigns a low probability of 0.04 to this term.
