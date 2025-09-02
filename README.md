# particle_tracking
python code to track particles through ocean models

The first is a set of routines that track particles though FVCOM fields. Most of this code was written by Vitalii Sheremet pre-pandemic and modified by JiM.

At the time of this writing, it has three different python files:
1) get_fvcom_mon_nc.py to get the monthly output field from their on-line archive at SMAST.
2) pt_functions.py a set of functions including the Rungekutta steps
3) dtr_jim.py the main program

This is a reduced version of Vitalii's extensive code which was designed specifically by JiM during the Spring of 2022 to look at the potential transport of shock tubing plastic showing up on Cape Cod beaches. The results of that study are posted <a href='http://studentdrifters.org/projects/ludwig'> here.</a>

Note: I had some sample input files on a compressed zip  since the files would not fit on github

The second is a simple example of using OpenDrift using John Wilkin's ROMS on DOPPIO grid.
