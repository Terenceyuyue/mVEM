function main_elasticityVEM_CookMembrane
clc; clear; close all;

%% PDE data
E = 250;  nu = 0.4999;
pde = elasticitydata_CookMembrane(E, nu);

%% Mesh
Th = load('meshdata_Cook1.mat');
%Th = load('meshdata_gunelve.mat');
node = Th.node;  elem = Th.elem;

%% Virtual element method
method = 3;
switch method
    case 1 % conforming VEM
        bdStruct = setbd(node,elem);
        [u,info] = elasticityVEM(node,elem,pde,bdStruct);
    case 2 % Kouhia-Stenberg: V_nc \times V_c
        bdStruct = setbd(node,elem);
        [u,info] = elasticityVEM_KouhiaStenberg(node,elem,pde,bdStruct);
    case 3
        refineType = 1;
        node0 = node;  elem0 = elem;
        [node,elem] = PolyMeshUniformRefine(node0,elem0,refineType);
        bdStruct = setbd(node,elem);
        [u,info] = elasticityVEM_NCreducedIntegration(node0,elem0,node,elem,pde,bdStruct,refineType);
end

%% Plot
[uhI,nodeI,elemI] = EllipticProjection(node,elem,u,info);
figure,
showsolution(nodeI+uhI,elemI,uhI(:,2)); view(2); colorbar;

end

%% Boundary information
function bdStruct = setbd(node,elem)
    bdStruct = setboundaryPiecewise(node,elem, 'x==0', 'x==48');
    % Dirichlet: left
    bdStruct.bdEdgeD = bdStruct.bdEdgeType{1};
    bdStruct.bdEdgeIdxD = bdStruct.bdEdgeIdxType{1};
    bdStruct.bdNodeIdxD = bdStruct.bdNodeIdxType{1};
    % Neumann: right (the remaining sides are traction-free)
    bdStruct.bdEdgeN = bdStruct.bdEdgeType{2};
    bdStruct.bdEdgeIdxN = bdStruct.bdEdgeIdxType{2};
    bdStruct.bdNodeIdxN = bdStruct.bdNodeIdxType{2};
end