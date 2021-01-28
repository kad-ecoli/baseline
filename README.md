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
(``bin/predict_nw.py``), but is not run by default due to long running time.
To run this predictor, uncomment the following line in ``predict.sh``:
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

In blast baseline (``bin/predict_blast.py``), 9 different scoring functions
are implemented. The "blast" baseline predictor uses confidence score based
on local sequence identity:
```
Cscore_blastlocalID(q) =       max       { localID_t(q) }   ... (2)
                         t=1, ... , N(q)
```
Here, ``t`` is the index of blast hits in a default blast run; ``N(q)`` is the
number of blast hits annotated with GO term ``q``; and ``localID_t(q)`` is the
local sequence identity between target sequence and the ``t``-th blast hit 
with term ``q``. Local identity is calculated as
``localID_t=nident_t/length_t``, where ``nident_t`` is the number of identical
residues between target and ``t``-th hit while ``length_t`` is the number of
aligned residues between the target and the hit.

## License ##
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 2 of the License, or
(at your option) any later version.
