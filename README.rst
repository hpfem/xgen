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


Testing
=======

Run XGEN by typing
::

    xgen xgamm cfg/xgamm.cfg

Here "xgen" is the executable, "xgamm" a project name, and "cfg/xgamm.cfg"
a configuration file with geometry parameters. Other project names 
include "xcirc", "xduese", "xhole", "xlist", "xmunich", "xnozzle",
"xsep2d", "xspir2d", "xsquare", "xstep". See the Documentation for 
instructions on how to create your own projects. 


Contact
=======

Please send questions and bugs to the mailing list 
hermes2d@googlegroups.com of the hp-FEM group, 
University of Nevada, Reno.



