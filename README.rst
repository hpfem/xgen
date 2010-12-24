===============
Welcome to XGEN
===============

XGEN is a simple triangular mesh generator for X-Windows
based on the Advancing Front method, written using the 
Motif library.


Copyright
=========

Copyright (c) 2010 hp-FEM group at the University of Nevada,
Reno (UNR). Email: hpfem@unr.edu, home page: http://hpfem.org/.


User Documentation
==================

After downloading XGEN, change to xgen2d/doc/ and type

    pdflatex xgen.tex


Prerequisites (Ubuntu Linux)
============================

    sudo apt-get install libmotif-dev libmotif3 lesstif2 x11proto-print-dev 


Compilation
===========

Type "make" to compile.


Running in interactive mode
===========================

Run XGEN by typing
::

    ./xgen xgamm cfg/xgamm.cfg

Here "xgen" is the executable, "xgamm" a project name, and "cfg/xgamm.cfg"
a configuration file with geometry parameters. Other project names 
include "xcirc", "xduese", "xhole", "xlist", "xmunich", "xnozzle",
"xsep2d", "xspir2d", "xsquare", "xstep". See the Documentation for 
instructions on how to create your own projects. 

By default, initial point distribution is random. You can force it to 
be regular pattern resembling equilateral elements via the option 
"-overlay"::

    ./xgen xgamm cfg/xgamm.cfg -overlay


Running in batch mode
=====================

Xgen can be run in batch mode via the "-nogui" option::

    ./xgen xgamm cfg/xgamm.cfg -nogui 10 -overlay

The option "-nogui" has to be followed with an empty 
character and a natural number telling Xgen how many 
smoothing iterations over all grid points it should do. 
The "-overlay" option here is not mandatory.

If "-nogui" and "-overlay" are used together, they 
must come in this order.


Contact
=======

Please send questions and bugs to the mailing list 
hermes2d@googlegroups.com of the hp-FEM group, 
University of Nevada, Reno.



