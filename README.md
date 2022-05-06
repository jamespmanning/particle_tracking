# particle_tracking
python code to track particles through FVCOM fields

This is a set of routines that track particles though FVCOM fields.
Most of this code was written by Vitalii Sheremet pre-pandemic and modified by JiM.

At the time of this writing, it has three different python files:
1) get_fvcom_mon_nc.py to get the monthly output field from their on-line archive at SMAST.
2) pt_functions.py a set of functions including the Rungekutta steps
3) dtr_jim.py the main program
