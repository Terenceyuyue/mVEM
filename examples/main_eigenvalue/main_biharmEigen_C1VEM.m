clc; close all;
clear variables;

%% Parameters
nameV = 1000:1000:4000;
maxIt = length(nameV);
h = zeros(maxIt,1);  N = zeros(maxIt,1);


%% PDE data
%[node,elem] = squaremesh([0 1 0 1], 0.5,0.5);
%% Virtual element method
for k = 1:maxIt
    % load mesh
    load( ['meshdata', num2str(nameV(k)), '.mat'] );
    %[node,elem] = uniformrefine(node,elem);
    % get boundary information
    bdStruct = setboundary(node,elem);
    % solve
    lam = biharmEigen_C1VEM(node,elem,bdStruct)
end





