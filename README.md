# baseline #
Implemtation of different CAFA baseline predictors

## Install ##
```bash
git clone https://github.com/kad-ecoli/baseline
cd baseline
git submodule init
git submodule update
./download.sh
cd src/
make install
cd ..
```
For assessment, this package use the 
[CAFA_assessment_tool](https://github.com/ashleyzhou972/CAFA_assessment_tool),
which dependends on python3 to generate correct Fmax. For prediction, both
python2.7 and python3 can be used, and will generate identical result.

## Usage ##
Run naive (``bin/predict_naive.py``) and blast (``predict_blast.py``) baseline
predictors on CAFA3 targets:
``bash
./predict.sh
```
Input files are at ``input/``. Predictions are at ``prediction/``.

Assess the performance (Fmax) of different predictor
```bash
./assess.sh
./plot.py
```
Assessment summaries are at ``CAFA_assessment_tool/results/``.
Summary graphic is at Fmax_full.png.

## License ##
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 2 of the License, or
(at your option) any later version.
