# baseline #
Implemtation of different CAFA baseline predictors

## Install ##
```bash
git clone https://github.com/kad-ecoli/baseline
cd baseline
git submodule init
git submodule update
```

The following commands are not necessary, because the binaries and data files
are already included as part of this package.
```bash
./download.sh
cd src/; make install ; cd .. # Not needed if global alignment is not used.
```

For assessment, this package use the 
[CAFA_assessment_tool](https://github.com/ashleyzhou972/CAFA_assessment_tool),
which dependends on python3 to generate correct Fmax. For prediction, both
python2.7 and python3 can be used, and will generate identical result.

## Usage ##
Run naive (``bin/predict_naive.py``) and blast (``bin/predict_blast.py``)
baseline predictors on CAFA3 targets:
```bash
./predict.sh
```
Input files are at ``input/``. Predictions are at ``prediction/``.

There is also a predictor based on the Needleman-Wunsch global aligner 
(``bin/predict_nw.py``), but is not run by default due to very long running
time.  To run this predictor, uncomment the following line in ``predict.sh``:
```bash
    #$rootdir/bin/predict_nw.py $rootdir/input/target.$species.fasta _${species}_go.txt
```

Assess the performance (Fmax) of different predictor:
```bash
./assess.sh
./plot.py
```
Assessment summaries are at ``CAFA_assessment_tool/results/``.
Summary graphic is at Fmax_full.png.

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
term ``q``:
```
Cscore_rank(q) =     max        { 1 - rank_t(q) / N }   ... (7)
                t=1, ... , N(q)
```
Here, ``N`` is the total number of hits (with and without ``q``.)

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
(2) to (10) are based on blast local alignment using **default** blast search
parameters, except for a change of output format (``-outfmt 6``) for easy
parsing.

The ``bin/predict_blastbitscore.py`` script implement an alternative scoring
based on bitscore:
```
Cscore_bitscore(q) =     max       { bitscore_t(q) / bitscore_self }   ... (11)
                    t=1, ... , N(q)
```
Here, ``bitscore_self`` is the bitscore of aligning target to itself. Thus, to
compute equation (11), two blast runs are needed: once for searching the target
against the template database (i.e. training database), the second time for
searching the target against itself.

The ``bin/predict_nw.py`` program implements equations (2) to (5) and (7) to
(10), using Needleman-Wunsch global sequence alignment (NWalign) rather than
blast local alignment. The implementation of NWalign is identical to the 
[Fortran version](https://zhanglab.ccmb.med.umich.edu/NW-align/),
but with input and output format changed to support multiple sequences in a
single input file. For consistency with blast, only the top 500 hits with the
highest alignment score are considered. Since NWalign does not report bitscore,
alignment score is used instead of bitscore for equation (10).

## Result ##
![Fmax_full.png](Fmax_full.png?raw=true "Fmax_full.png")
In terms of Fmax at "full" mode, different scoring in the order of accuracy
are: naive ≈ evalue < localID < bitscore < globalID1 ≈ globalID2 ≈ globalID3
< rank < freq < metago < netgo. In particular, the three scoring functions that
consider all hits (freq, metago, netgo) result in consistently higher accuracy
than all scorings that use only the top hit (evalue, localID, bitsore,
globalID, rank), including the current official "blast" baseline (localID).

## License ##
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 2 of the License, or
(at your option) any later version.
