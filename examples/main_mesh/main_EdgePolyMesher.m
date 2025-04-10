clc;clear;close all

id = 1;
switch id
    case 1 
        load meshex4
        aux = auxstructure(node,elem);
        bdEdge = aux.bdEdge;
        P1 = node(bdEdge(:,1),:);
        P2 = node(bdEdge(:,2),:);
        %figure,showmesh(node,elem);

    case 2 % flower
        [P1,P2] = edgeFlower(100);

    case 3 % heart
        [P1,P2] = edgeHeart(50);
end
options.plot = true;
options.refineTimes = 0;
[node,elem] = EdgePolyMesher(P1,P2,options);

