# Open-Clusters-in-the-solar-neighbourhood-with-TGAS

Fits file CG18_OC_members_light.fits contains membership probability for stars in 128 open clusters. The description of the method can be found in Cantat-Gaudin et al. 2018 (https://arxiv.org/abs/1801.10042).

The python script plot.py reads the fits file and produces a six-panel png figure for each open cluster (three examples are provided, for instance ASCC_113_6panels.png).

It can be run as:
    python plot.py
to print all 128 figures, or as:
    python plot.py ASCC_113 Stock_2
to produce figures for a selected number of clusters (provided they are in the list)...
