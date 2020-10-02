# Bruker

This repository contains quick-and-dirty utility scripts for Bruker scanning


##### [Matt Budde's Bruker Paravision 6.0.1 pulse sequence for arbitrary diffusion encoding waveforms](https://osf.io/t9vqn/wiki/rFOV_DWIEpi/):

Experimentally, we have discovered that the sequence only detects the files if:
* Main orientation file Dir.dat need to have EXACTLY 6 decimal digit 
* The points in the waveform file need to STRICTLY have norm < 1 with the decimal truncation.
* Waveform file cannot have extension (need to remove the .gp)

Also,
* Waveform files need to be in  /opt/PV6.0.1/prog/curdir/USER/ParaVision/exp/lists/gp  
* only "-" and "_" special characters are allowed.


**make_waveform.m** is a very very very slightly modified version of [Daniel Topgaard Gradient waveforms for axisymmetric b-tensors code](https://github.com/daniel-topgaard/md-dmri/blob/master/acq/bruker/paravision/make_waveform.m).


**make_waveform.py** is a commandline wrapper around make_waveform.m using [oct2py](https://pypi.org/project/oct2py/).    
`make_waveform.py 0.66 0.01 sph /home/user/schemefile/ sphericalRAW.gp -p`  
to generate and plot a spherical b-tensor waveform with gmax = 0.66 T/m and a duration = 0.01 s and saving it as /home/user/schemefile/sphericalRAW.gp


**fix_gp_norm.py** is a heuristic script that fixes the waveform points that have norm above one (because of decimal truncation)  
`fix_gp_norm.py /home/user/schemefile/sphericalRAW.gp /home/user/schemefile/spherical`


**guess_tau_from_b.py** is a heuristic utility script to guess the correct waveform duration (tau) from the desired maximum gradient strenght and the target b-value (careful with units)  
`guess_tau_from_b.py 1 0.66`    
`>> 14.842199577008097`   
So a 14.84 ms waveform with gmax = 0.66 T/m will have a b-value of 1 ms/um^2

**guess_G_from_b.py** is a heuristic utility script to guess the correct gradient strenght from the desired waveform duration and the target b-value (careful with units)  
`guess_G_from_b.py 1 14.84`    
`>> 0.660146742717`   
So a 14.84 ms waveform with gmax = 0.66 T/m will have a b-value of 1 ms/um^2

**fslbvec2budde.py** converts a fsl bvec file into a compatible Dir.dat file without touching any of the vectors norm (assumes all norms = 1)  
`fslbvec2budde.py bvec.txt /home/user/schemefile/ dir.dat`
