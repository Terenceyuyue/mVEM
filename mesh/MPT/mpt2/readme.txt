% MPT Multi-Parametric Toolbox
% Version 2.6.2 05-Dec-2006
%
% http://control.ee.ethz.ch/~mpt/
%
% Authors: Michal Kvasnica, Pascal Grieder, Mato Baotic
% kvasnica@control.ee.ethz.ch, grieder@control.ee.ethz.ch, baotic@control.ee.ethz.ch

MPT Toolbox installation notes


Section 1
=========

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Note: Remove any previous copies of MPT before installing this version!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

The MPT toolbox consists of the following directories

mpt/           - toolbox main directory
mpt/@mptctrl   - methods of the 'mptctrl' object
mpt/@polytope  - methods of the 'polytope' object
mpt/docs       - documentation
mpt/examples   - sample dynamical systems and other demos
mpt/extras     - extra functionality
mpt/solvers    - additional solvers
mpt/yalmip     - contains source-code distribution of the YALMIP interface

In order to use MPT, set a Matlab path to the "mpt/" directory and to all it's
subdirectories. If you are using Matlab for Windows, go to the 
"File - Set Path..." menu, choose "Add with Subfolders..." and pick 
up the MPT directory. Click on the "Save" button to store the updated 
path setting. Under Unix, you can either manually edit the file 
"startup.m", or to use the same procedure described above.

You can also add the path manually by running following command on the Matlab
prompt:

  addpath(genpath('/path/to/mpt/'));
  
For more information regarding installation consult following page:

  http://control.ee.ethz.ch/~mpt/downloads/install.php

Note: don't put mpt/@polytope folder to your path manually, it is not necessary
      and may lead to wrong behavior of the toolbox.

Once you install the toolbox, please consult Section 3 on how to set default
values of certain parameters.

To explore functionality of MPT, try one of the following:

mpt_studio
help mpt
help mpt/polytope
mpt_demo1
mpt_demo2
mpt_demo3
mpt_demo4
mpt_demo5
mpt_demo6
runExample

                                                  
MPT toolbox comes with a set of pre-defined examples which the user can go
through to get familiar with basic features of the toolbox.
The toolbox has been successfully tested with Matlab 6.1(R12) and 6.5(R13)
on Windows, Solaris and Linux, along with NAG solver, CPLEX 8.0 and CDDlib 0.93.
     

If you wish to be informed about new releases of the toolbox,
subscribe to our mailing list by sending an email to:

  mpt-request@list.ee.ethz.ch

and put the word

  subscribe

to the subject field. To unsubscribe, send an email to the same mail
address and specify

  unsubscribe

on the subject field.


If you have any questions or comments, or you observe buggy behaviour of the
toolbox, please send your reports to

  mpt@control.ee.ethz.ch



Section 2 (Additional software requirements)
============================================

LP and QP solvers:
------------------

The MPT toolbox is a package primary designed to tackle multi-parametric
programming problems. It relies on external Linear programming (LP) and
Quadratic programming (QP) solvers. Since the LP and QP solvers shipped together
with Matlab (linprog and quadprog) are rather slow in terms of computational
speed, the toolbox provides a unified interface to other solvers.

One of the supported LP solvers is the free CDD package from Komei Fukuda
(http://www.cs.mcgill.ca/~fukuda/soft/cdd_home/cdd.html)

The CDD is not only a fast and reliable LP solver, it can also solve many
problems from computational geometry, e.g. computing convex hulls, extreme
points of polytopes, calculating projections, etc.

A pre-compiled version of the Matlab interface to CDD is included in this
release of the MPT toolbox. The interface is available for Windows, Solaris and
Linux. You can download the source code of the interface from:
http://control.ee.ethz.ch/~hybrid/cdd.php

and compile it on your own for other platforms. Please consult Section 3 on how
to make CDD a default LP solver for the MPT toolbox.

An another supported possibility is the commercial CPLEX solver from ILOG. The
authors provide an interface to call CPLEX directly from Matlab, you can
download source codes and pre-compiled libraries for Windows, Solaris and Linux
from
http://control.ee.ethz.ch/~hybrid/cplexint.php

Please note that you need to be in possession of a valid CPLEX license in order
to use CPLEX solvers.

The free GLPK (GNU Linear Programming Kit) solver is also supported by MPT
toolbox. In order to use GLPK in Matlab, you need to download a mex interface
written by Nicolo Giorgetti from:
http://www-dii.ing.unisi.it/~giorgetti/downloads.php

Note that we have experienced several numerical inconsistencies when using GLPK.


Semi-definite optimization packages:
------------------------------------

Some routines of the MPT toolbox rely on Linear Matrix Inequalities (LMI)
theory. Certain functions therefore require solving a semidefinite optimization
problem. The YALMIP interface by Johan Lofberg
http://control.ee.ethz.ch/~joloef/yalmip.php

is included in this release of MPT toolbox. Since the interface is a wrapper and
calls external LMI solver, we strongly recommend to install one of the solvers
supported by YALMIP. You can obtain a list of free LMI solvers here:
http://control.ee.ethz.ch/~joloef/manual/htmldata/solvers.htm

YALMIP supports a large variety of Semi-Definite Programming packages. One of
them, namely the SeDuMi solver written by Jos Sturm, comes along with MPT.
Source codes as well as binaries for Windows are included directly, you can
compile the code for other operating systems by following the instructions in 
  mpt/solvers/SeDuMi105/Install.unix
For more information, visit
http://fewcal.kub.nl/sturm/software/sedumi.html


Solvers for projections:
------------------------

MPT allows to compute orthogonal projections of polytopes. To meet
this task, several methods for projections are available in the
toolbox. Two such methods - ESP and Fourier-Motzkin Elimination
are coded in C and need to be accessible as a mex library. These
libraries are already provided in compiled form for Linux and
Windows. For other architectures you will need to compile the
corresponding library on your own. To do so, follow instructions in
  mpt/solvers/esp
and 
  mpt/solvers/fourier
respectively.


Section 3 (Setting up default parameters)
========================================

Any routine of the MPT toolbox can be called with user-specified values of
certain global parameters. To make usage of MPT toolbox as user-friendly as
possible, we provide the option to store default values of the parameters in
variable "mptOptions", which is kept in MATLAB's workspace as a global variable
(i.e. it stays there unless one types 'clear all').

The variable is created when the toolbox get's initialized through a call to
mpt_init.

Default LP solver:
------------------
In order to set the default LP solver to be used, open the file "mpt_init" in
your editor. Scroll down to the following line:

  mptOptions.lpsolver = [];

Integer value on the right-hand side specifies the default LP solver. Allowed
values are:
0 - NAG Foundation LP solver
3 - CDD Criss-Cross method (default)
2 - CPLEX
4 - GLPK
5 - CDD Dual-Simplex method
1 - linprog
6 - QSOpt

If the right-hand argument is an empty matrix, the fastest solver available
will be used. In the above table, solvers are sorted in the order of
preference.

Default QP solver:
------------------

To change the default QP solver, locate and modify this line in mpt_init.m:

  mptOptions.qpsolver = [];

Allowed values for the right-hand side argument are the following:
0 - NAG Foundation QP solver
1 - quadprog (default)
2 - CPLEX

Again, if you leave this line unchanged, MPT will choose the fastest available
solver.

Note: If you don't have any QP solver available on your machine, set
  mptOptions.qpsolver = -1;

Then you will still be able to solve multi-parametric problems with linear
cost objective. To ensure that, set the proper flag of your problem structure,
like follows:

>> probStruct.norm = 1;

or

>> probStruct.norm = Inf;


Default solver for extreme points computation:
----------------------------------------------

Some of the functions in MPT toolbox require computing of extreme points of
polytopes given by their H-representation and calculating convex hulls of given
vertices respectively. Since efficient analytical methods are limited to low
dimensions only, we provide the possibility to pass this computation to an
external software package (CDD). However, if the user for any reason does not
want to use third-party tools, the problem can still be tackled in an analytical
way (with all the limitations mentioned earlier).

To change the default method for extreme points computation, locate the
following line in mpt_init.m:

  mptOptions.extreme_solver = [];

and change the right-hand side argument to one of these values:
0 - analytical computation 
3 - use CDD (faster computation, works also for higher dimensions)


Default tolerances:
-------------------

The MPT toolbox internally works with 2 types of tolerances:
- absolute tolerance
- relative tolerance

Default values for these two constants can be set by modifying the following
lines of mpt_init.m:

  mptOptions.rel_tol = 1e-6;

  mptOptions.abs_tol = 1e-7;


Default values for Multi-parametric solvers:
-------------------------------------------

Solving a given QP/LP in a multi-parametric way involves making "steps" across
given boundaries. Length of this "step" is given by the following variable:

  mptOptions.step_size = 1e-4;

Due to numerical problems tiny regions are sometimes difficult to calculate,
i.e. are not identified at all. This may create "gaps" in the computed control
law. For the exploration, these will be jumped over and the exploration in the
state space will continue.


Level of detecting those gaps is given by the following variable:

  mptOptions.debug_level = 1;

The right-hand side argument can have three values:

0: No debug done
1: A tolerance is given to find gap in the region partition,
   small empty regions inside the region partition will be discarded.
   Note that this is generally not a problem, since the feedback law
   is continuous and can therefore be interpolated easily.
   Correction to the calculation of the outer hull is performed as well.
2: Zero tolerance to find gaps in the region partition, empty regions
   if they exist, will be detected, i.e. the user will be notified.
   Correction to the calculation of the outer hull is performed.


Default Infinity-box
--------------------

The MPT toolbox internally converts the R^n to a box with large bounds. The
following parameter specifies size of this box:

  mptOptions.infbox = 1e6;

Note that also polyhedra (unbounded polytopes) are converted to bounded
polytopes by making an intersection with the "Infinity-box".


Default values for plotting:
----------------------------

The overloaded 'plot' function can be forced to open a new figure windows every
time the user calls it. If you want to disable this feature, go to the following
line in mpt_init.m :

  mptOptions.newfigure = 0;

and change the constant to 0 (zero)

1 means "enabled", 0 stands for "disabled"


Default verbosity level:
------------------------

Text output from functions can be limited or suppressed totally by changing
the following option in mpt_init.m :

  mptOptions.display = 1;

Allowed values are:
0 - no output at all (except of critical warnings and error messages)
1 - only the important information are displayed (default)
2 - no output suppression, all information are reported to the user



Once you modify the mpt_init.m file, type:

  mpt_init

to initialize the toolbox.
