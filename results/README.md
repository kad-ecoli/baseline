# Results #
Folder for combined prediction and assessment results.

## Assessment statistics ##

* ``blast*_*_all_results.txt`` - baselines based on blast local sequence search
* ``nw*_*_all_results.txt``    - baselines based on NWalign global sequence alignment
* ``iea_1_all_results.txt``    - copy all annotations, including IEA
* ``naive_*_all_results.txt``  - background probabilities of GO term in the database

## Plotting scripts ##

* ``plot.py`` - plot {Fmax,Smin}_full.png
* ``plot_nw.py`` - plot {Fmax,Smin}_full_nw.png

## Result images ##

* [Fmax_full.png](Fmax_full.png?raw=true)    - Fmax for blast, iea and naive.
* [Smin_full.png](Smin_full.png?raw=true)    - Smin for blast, iea and naive.
* [wFmax_full.png](wFmax_full.png?raw=true)    - weighted Fmax for blast, iea and naive.
* [Fmax_full_nw.png](Fmax_full_nw.png?raw=true) - Fmax for blast and NWalign.
* [Smin_full_nw.png](Smin_full_nw.png?raw=true) - Smin for blast and NWalign.
* [wFmax_full_nw.png](wFmax_full_nw.png?raw=true) - weighted Fmax for blast and NWalign.
