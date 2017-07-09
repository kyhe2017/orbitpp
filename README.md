--------------------------------------------------------------------------------
# Basic introduction to orbitpp
--------------------------------------------------------------------------------
This code is developed to process ORBIT-RF output, which is object oriented 
based on FORTRAN 2003+.                               PERSONAL USE ONLY !!!
--------------------------------------------------------------------------------

--------------------------------------------------------------------------------
# Design of code
--------------------------------------------------------------------------------
The code has two classes in gerneral: Data classes & Data processing classes.

Data class: 
   1. geometry     : collection of equilibrium variables and related functions
   2. trajectory   : collection of test particle history and related functions
   3. momenta      : collection of statistical data and related functions
   4. distribution : collection of particle distribution and related functions
   5. aux_distribution : derived class of distribution to collect particle 
      across spiecies
      
Data processing class:
   1. statistics     : collection of functions to make statistics on distribution 
      object based on linear schedule
   2. aux_statistics : derived class of statistics, linear schedule replaced by
      exponential weight assignment 
--------------------------------------------------------------------------------

--------------------------------------------------------------------------------
# Usage of code
--------------------------------------------------------------------------------
Taking geometry class as an example:
   1. first declare a pointer:
      class(geometry), pointer :: gg
   2. then allocation:
      allocate(geometry :: gg)
   3. then constructor:
      call gg % init(nsp, nst)
   4. then load data:
      call gg % ld_data(fid)
   5. then wrt output optionally, write out what you need:
      call gg % wrt('psip')
   6. finally deallocation, very important and never forget it:
      deallocate(gg)
Note: classes have different ways to initialize and subroutines to write output,
you need to know these by consulting the type definition in source code.
--------------------------------------------------------------------------------

--------------------------------------------------------------------------------
# graphical ouput
--------------------------------------------------------------------------------
The default tool is Gnuplot. 

basic plotting commands:
  * 1d :
       plot 'filename' u 1:2 w l
       plot 'filename' u 1:2 w p
       plot 'filename' u 1:2 w lp
       plot 'filename' u 1:2 w boxes
  * 2d map :
       splot 'filename' u 1:2:3 w pm3d
       plot 'filename' u 1:2:3 w image
       
For detailed usage: http://www.gnuplot.info/
--------------------------------------------------------------------------------

