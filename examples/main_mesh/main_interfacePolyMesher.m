clc;clear;close all

options.plot = true;
options.hmax = 1;

id = 3;
switch id

    case 1 %% circle: default        
        options.refineTimes = 1;
        [node,elem] = interfacePolyMesher([],[],options);
        % figure, showmesh(node,elem);

    case 2 %% flower        
        options.refineTimes = 0;
        % rectangle
        P = [-1 -1; 1 -1; 1 1; -1 1];
        P1 = P(1:end,:);  P2 = P([2:end,1],:);
        g1 = decompGeo(P1,P2);
        % flower
        N = 80;
        [P1,P2] = edgeFlower(N);
        g2 = decompGeo(P1,P2);
        [node,elem] = interfacePolyMesher(g1,g2,options);
        % figure, showmesh(node,elem);

    case 3 %% heart        
        options.refineTimes = 0;
        % rectangle
        P = [-2 -2; 2 -2; 2 2; -2 2];
        P1 = P(1:end,:);  P2 = P([2:end,1],:);
        g1 = decompGeo(P1,P2);
        % flower
        N = 50;
        [P1,P2] = edgeHeart(N);
        g2 = decompGeo(P1,P2);
        [node,elem] = interfacePolyMesher(g1,g2,options);
        % figure, showmesh(node,elem);
end