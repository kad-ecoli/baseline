# baseline #
Implemtation of different CAFA baseline predictors

## Install ##
```bash
git clone https://github.com/kad-ecoli/baseline
```

The following commands are not necessary, because the binaries and data files
are already included as part of this package.
```bash
./download.sh
cd src/; make install ; cd .. # Not needed if global alignment is not used.
```

To run the IEA baseline, the ``input/target.map`` file must be manually
prepared to map entry name to uniprot accession. A version for CAFA3 targets
are already pre-generated.

## Usage ##
Run naive (``bin/predict_naive.py``), blast (``bin/predict_blast.py``) and
IEA (``bin/predict_IEA.py``) baseline predictors on CAFA3 targets:
```bash
./predict.sh
```
Input files are at ``input/``. Predictions are at ``prediction/``.

Assess the performance (Fmax) of different predictor:
```bash
./assess_both.sh
```
Assessment summaries are at ``results/``.
Summary graphic is at [results/Fmax_full.png](results/Fmax_full.png)
and [results/Smin_full.png](results/Smin_full.png)

There is also a predictor based on the Needleman-Wunsch global aligner 
(``bin/predict_nw.py``). It is not run by default because of very long running
time with little to no improvement over blast. To run and assess this global
alignment based predictor, use:
```bash
./predict_nw.sh
```
Assessment summaries are also at ``results/``.
Summary graphic is at [results/Fmax_full_nw.png](results/Fmax_full_nw.png)
and [results/Smin_full_nw.png](results/Smin_full_nw.png)

## Math ##
In naive baseline (``bin/predict_naive.py``), the confidence score of
predicting a GO term ``q`` for a target protein is calculated as:
```
Cscore_naive(q) = M(q) / M   ... (1)
```
Here, ``M`` is the total number of proteins in the whole training database;
and ``M(q)`` is the number of proteins annotated with ``q`` in the training
database. Therefore, equation (1) represents the background probability
of a GO term in the training database, regardless of target protein.

In the blast baseline (``bin/predict_blast.py``), 9 different scoring functions
are implemented, as shown in equations (2) to (10) below. In the official
assessment of CAFA1-3, the "blast" baseline predictor uses confidence score
based on local sequence identity:
```
Cscore_localID(q) =      max        { localID_t(q) }   ... (2)
                    t=1, ... , N(q)
```
Here, ``t`` is the index of blast hits in a default blast run; ``N(q)`` is the
number of blast hits annotated with GO term ``q``; and ``localID_t(q)`` is the
local sequence identity between target sequence and the ``t``-th blast hit 
with term ``q``. Local identity is calculated as
``localID_t=nident_t/length_t``, where ``nident_t`` is the number of identical
residues between target and ``t``-th hit while ``length_t`` is the number of
aligned residues between the target and the hit.

Another three scores based on global seuqence identities are also implemented:
```
Cscore_globalID1(q) =     max        { globalID1_t(q) }   ... (3)
                     t=1, ... , N(q)

Cscore_globalID2(q) =     max        { globalID2_t(q) }   ... (4)
                     t=1, ... , N(q)

Cscore_globalID3(q) =     max        { globalID1_t(q) }   ... (5)
                     t=1, ... , N(q)
```
Here, ``globalID1_t=nident_t/qlen`` is the global sequence identity normalized
by the target protein length ``qlen``; ``globalID2_t=nident_t/slen`` is
normalized by the blast hit length ``slen``;
``globalID3_t=nident_t/max(qlen,slen)`` is normalized by the maximum of target
and hit length.

The fifth score is based on  ``evalue_t(q)``, the evalue for the ``t``-th hit 
with term ``q``.  Since evalue is in the range of [0,inf) rather than (0,1], a
sigmoid function is used to rescale it into (0,1]:
```
Cscore_evalue(q) =     max        { 1 - 1 /( 1 + exp( - evalue_t(q) ) }   ... (6)
                  t=1, ... , N(q)
```
Here, ``evalue_t(q)`` is the evalue for the ``t``-th hit with term ``q``.

The sixth score is based on ``rank_t(q)``, the ranking of ``t``-th hit with
term ``q`` among all hits:
```
Cscore_rank(q) =     max        { 1 - ( rank_t(q) - 1 ) / N }   ... (7)
                t=1, ... , N(q)
```
Here, ``N`` is the number of all hits (with and without ``q``.)

The seventh score is based on frequency of a GO term among hits:
```
Cscore_freq(q) = N(q) / N   ... (8)
```

The next two scores are advanced version of frequency, weighted by either
``globalID1_t`` and ``bitscore_t``, which are the global sequence identity
and bitscore of the ``t``-th blast hit:
```
                     N(q)                         N
Cscore_metago(q) = ( sum ( globalID1_t(q) )) / ( sum ( globalID_t ))   ... (9)
		     t=1                         t=1

                     N(q)                        N
Cscore_netgo(q)  = ( sum ( bitscore_t(q) )) / ( sum ( bitscore_t ))   ... (10)
	             t=1                        t=1
```
Equations (9) and (10) are implemented by the sequence-based submodule of
MetaGO and NetGO, respecitvely, and probably should not be considered
"baseline" due to mathematical complexity. All scoring schemes from equation
(2) to (10) are based on blast local alignment using **default blast search
parameters**, except for a change of output format (``-outfmt 6``) for easy
parsing.

The ``bin/predict_blastbitscore.py`` script implement an alternative scoring
based on bitscore:
```
Cscore_bitscore1(q) =     max       { bitscore_t(q) / bitscore_self(target) }  ... (11)
                    t=1, ... , N(q)

Cscore_bitscore2(q) =     max       { bitscore_t(q) / bitscore_self(t) }       ... (12)
                    t=1, ... , N(q)

Cscore_bitscore3(q) =     max       { bitscore_t(q) / max{ bitscore_self(target),
                    t=1, ... , N(q)                        bitscore_self(t)} } ... (13)
```
Here, ``bitscore_self(target)`` and ``bitscore_self(t)`` is the bitscore of
aligning the target to itself and template ``t`` to itself, respectively.

The ``bin/predict_nw.py`` program implements equations (2) to (5) and (7) to
(13), using Needleman-Wunsch global sequence alignment (NWalign) rather than
blast local alignment. The implementation of NWalign is identical to the 
[Fortran version](https://zhanglab.ccmb.med.umich.edu/NW-align/),
but with input and output format changed to support multiple sequences in a
single input file. Due to time expense of NWalign, it will only be performed
on the hits identified by a blast run. Since NWalign does not report bitscore,
alignment score is used instead of bitscore for equation (10) to (13).
For (7), the ranking is based on NWalign alignment score rather than the
ranking in the blast output.


In IEA baseline (``bin/predict_iea.py``), the GO term of a protein is copied
from its full set of uniprot-goa, which mainly includes (but is not limited) 
electronically inferred annotations with IEA evidence, hence the name.
The confidence score of the GO term is determined by the evidence code from
in the uniprot-goa annotation. 
In [our previous study](https://doi.org/10.1093/bioinformatics/btaa548),
we obtained statistics on the portion of GO terms with the same evidence
code ``e`` that are either experimentally confirmed (``N_confirm(e)``) or 
rejected with a "NOT" qualifier in a later release (``N_reject(e)``). The 
confidence score of the evidence code ``e`` is therefore:
```
Cscore_iea(e) = N_confirm(e) / ( N_confirm(e) + N_reject(e) )   ... (14)
```

## Result ##
See full results at the [results/](results/) folder.

![results/Fmax_full.png](results/Fmax_full.png?raw=true "results/Fmax_full.png")
In terms of Fmax at "full" mode, different scoring in the order of accuracy
are: evalue < localID < naive < bitscore1 < bitscore2 ≈ bitscore3 ≈ globalID1
≈ globalID2 ≈ globalID3 < rank < iea < freq < metago < netgo. In particular,
the three scoring functions that consider all hits (freq, metago, netgo)
result in consistently higher accuracy than all scorings that use only the top
hit of each term (evalue, localID, bitsore, globalID, rank), including the
current official "blast" baseline (localID) used in CAFA assessment.


![results/Fmax_full_nw.png](results/Fmax_full_nw.png?raw=true "results/Fmax_full_nw.png")
Comparison between NW and blast show that global alignment does not result in
more accurate GO prediction than local alignment, except for alnscore/bitscore,
where NW slightly outperform blast. NW and blast result in the same performance
by freq because the NW implement only perform realignment of blast hits rather
than a full database search; since the hits are the same, the frequency of GO
terms among hits are also the same.

## License ##
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 2 of the License, or
(at your option) any later version.
