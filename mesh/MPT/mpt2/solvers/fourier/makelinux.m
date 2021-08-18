% Builds a mexglx for the fourier elimination routine

mex -o fourier fourier.cc support.cc matlab_driver.cc

fprintf('The fourier elimination code is now compiled.\nType help fourier for more information\n');
