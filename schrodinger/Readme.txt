WELCOME TO SCHRODINGER!


SCHRODINGER is a set of MATLAB files which provide moderately accurate 
numerical solutions of a number of standard problems familiar from many 
introductory texts on Quantum Mechanics.   They give the student an
interactive access to the solutions with graphical display of the 
results.

This is not a substitute for a textbook, it is a device for solving 
standard types of problems with a speed which will allow an interactive 
study of the effect of changes of the parameters.   The algorithms are
chosen for simplicity and speed, not optimal accuracy or stability. 
Calculations using matrix arguments, one of the key features of MATLAB, is
used wherever possible for maximal speed. 

There is a short manual in postscript form which gives some comments on
the quantum mechanical problems and on the algorithms, and a list of
references. The documentation for each problem is rudimentary, but I 
have tried to give literature references when this is possible. I have
also tried to set defaults for most variable parameters. In other places
the user has to guess  on the basis the uncertainty relation and the
spectrum of a square well. There is also some help in the M-files.

The correctness and accuracy of the algorithms defined in the M-files
should not be taken for granted. You are welcome to improve on them,
especially if you will share the  results.  The choice of conventions (in
the special functions,  for instance) most often conforms to that used in
Messiah, Quantum Mechanics. Throughout we have set hbar = mass = 1. 
 
The files should run on any computer with the MATLAB 5 program installed
and  a graphical interface. Almost everything should be compatible with
MATLAB 4.2.  Schrodinger has been tested on a PowerMac 7600/120
and on a PC running Linux.
 
The speed of solution evidently depends on the platform, but it can be 
adapted to your computer through the choice of the number of lattice
points and basis functions, by editing the M-files if necessary. The time
needed for most of the solutions are intended  to be at most a few tens of
seconds on a PowerMac or PC.

You can start the session by executing the file <start.m>, then making
choices pressing buttons, in menus or in the command window. Short
descriptions of the M-files are given in the file <Contents.m>. 
You can  also use the help function.

You may have problems with the size of the default for font size, figure
window size.  CHOOSE YOUR DEFAULTS in your startup.m file, or in start.m,
to get a window of a suitable size and position, and text of a size which
will fit in the window, eg
rect = [70 215 565 573];
set(0,'defaultfigureposition',rect);
set(0,'defaultaxesfontsize',12);
set(0,'defaulttextfontsize',12);
You may also have to set the fontsize in <text0.m>, and in the buttons.

The address of the home page of schrodinger is 

http://www.theophys.kth.se/mathphys/schrodinger.html
---------------------------------------------------------------------------
%%% © Göran Lindblad 1999 - gli@theophys.kth.se
