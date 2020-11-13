INSTALLATION AND COMPILATION INSTRUCTIONS:

(NOTE: GMP is not being used in this build)

This package comes together with pre-compiled version of the CDD library for
MATLAB. Just unzip the archive to some folder, set a matlab path there and
you can directly use cdd. See 'help cddmex' for more details.

If you want to recompile the library, follow the next steps:

Windows platform:
=================

- download and instal cygwin
  http://www.cygwin.com
  
  (Note: be sure that you include Mingwc compiler in the installation. It is a
         part of the "Devel" section in the cygwin installator)

- download the latest CDD source code package from
  http://www.cs.mcgill.ca/~fukuda/download/cdd/
  
  (Note: Matlab CDD interface was successfully compiled and tested with version
         0.93c of the CDD package)
         
- extract CDD source codes to some directory

- start cygwin, go to a directory containing CDD sources and type:
  ./configure
  
- if you want to turn off the internal debugging feature of CDD, go to
  lib-src/
  directory and locate the following line in the source files:
  
    dd_boolean localdebug=dd_TRUE; 
  
  change it to
  
    dd_boolean localdebug=dd_FALSE;
  
  (Note: The directive appears several times in the following files:
         cddcore.c, cddlib.c, cddlp.c, cddproj.c)

- if you skipped previous point, go to the lib-src/ directory

- open the file 
    Makefile
  in your text editor. Locate the following line:
  
  DEFS = -DPACKAGE_NAME=\"\" -DPACKAGE_TARNAME=\"\" -DPACKAGE_VERSION=\"\" 
         -DPACKAGE_STRING=\"\" -DPACKAGE_BUGREPORT=\"\" -DPACKAGE=\"cddlib\"
         -DVERSION=\"0.93\" -DSTDC_HEADERS=1
  
  and append '-mwindows -mno-cygwin' to it, i.e. it should look like this:
  
  DEFS = -DPACKAGE_NAME=\"\" -DPACKAGE_TARNAME=\"\" -DPACKAGE_VERSION=\"\" 
         -DPACKAGE_STRING=\"\" -DPACKAGE_BUGREPORT=\"\" -DPACKAGE=\"cddlib\"
         -DVERSION=\"0.93\" -DSTDC_HEADERS=1 -mwindows -mno-cygwin 
  
  You can close the editor now.
  
- type
    make
    
  this should compile the source files and link them together.
  
- download GNUMEX from
  http://www.mrc-cbu.cam.ac.uk/Imaging/gnumex20.html
  
  and put it somewhere in your matlab path
  
- start matlab, change your directory to the lib-src\ subdirectory of CDD

  Note: be sure that you also have the cdd.c file stored in that directory

- still in matlab, type:
    gnumex
  
  a dialog windows will appear. Some of the fields will be already filled out.
  in the third option ('Mingw, Cygwin or Cygwin-mingw linking?'), choose:
    Cygwin-mingw
    
  choosing a proper linking method is _VERY_ important!!!
  
  save the options by clicking on the appropriate button.
  
- type the following on matlab prompt:

    mex -v -f Mexopts.bat -I. cddmex.c cddcore.o cddio.o cddlib.o cddlp.o cddmp.o cddproj.o setoper.o
    
  it will create the file cdd.dll which can be used as a standard matlab mex
  function. You can copy the dll to an another directory and use it according
  to help described in cdd.m
  
- if the compilation was not successful, try to check one of the following:

  * did you modify the CDD Makefile correctly? Be sure that you include
      -mwindows -mnocygwin directives into the makefile!
      
  * be sure that cdd.c as well as mexopts.bat are in the lib-src\ subdirectory
    of the CDD source code library
    
  * make sure that you have selected Cygwin-mingwc for linking in gnumex!
  
  * check if you have the proper cygwin installation
  
  * you can also send me an email to: kvasnica@control.ee.ethz.ch
    and i will try to help you
    
    

Unix platform:
==============

- download the latest CDD source code package from
  http://www.cs.mcgill.ca/~fukuda/download/cdd/
  
  (Note: Matlab CDD interface was successfully compiled and tested with version
         0.93c of the CDD package)
         
- extract CDD source codes to some directory

- on the unix prompt, type:
    ./configure
    
- if you want to turn off the internal debugging feature of CDD, go to
  lib-src/
  directory and locate the following line in the source files:
  
    dd_boolean localdebug=dd_TRUE; 
  
  change it to
  
    dd_boolean localdebug=dd_FALSE;
  
  (Note: The directive appears at line 72 of the file cddlib.c)
    
    
- if you skipped previous point, go to the lib-src/ directory

- on the unix prompt type
    make
  
  It should create the pre-compiled library
  

- type
    mex -v cddmex.c libcdd.a
    
- it should create a mex file for you platform
    cddmex.mexglx for linux
    cddmex.mexsol for solaris
    
- in case of problems feel free to contact MPT help desk at
  mailto:mpt@control.ee.ethz.ch