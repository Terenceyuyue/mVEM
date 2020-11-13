% Builds a dll for the fourier elimination algorithm

!copy fourier.cc fourier.cpp
!copy support.cc support.cpp
!copy matlab_driver.cc matlab_driver.cpp

mex fourier.cpp support.cpp matlab_driver.cpp -output fourier

fprintf('The fourier elimination code is now compiled.\nType help fourier for more information\n');
