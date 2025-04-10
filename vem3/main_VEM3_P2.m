clc;clear;close all

%% Some preparations
% Mesh
node3 = [  
    0 0 0; % 顶点1  
    1 0 0; % 顶点2  
    1 1 0; % 顶点3  
    0 1 0; % 顶点4  
    0 0 1; % 顶点5 
    1,0,1; % 顶点6 
    1,1,1; % 顶点7 
    0,1,1; % 顶点8 
];
elemf = cell(6,1);
elemf{1} = [5,8,7,6]; % 上
elemf{2} = [1,2,3,4]; % 下
elemf{3} = [4,8,5,1]; % 左
elemf{4} = [3,2,6,7]; % 右
elemf{5} = [1,5,6,2]; % 前
elemf{6} = [4,3,7,8]; % 后
elem3 = {elemf}; 
option.FaceAlpha = 0.5; 
showmesh(node3,elem3,option); 
findnode(node3)

% Geometry
centroid3 = [0.5,0.5,0.5]; diameter3 = sqrt(3); volume = 1;
xK = centroid3(1,1); yK = centroid3(1,2); zK = centroid3(1,3); 
hK = diameter3(1);

% auxiliary data structure
edge = [1 2; % 边 1
        1 4; % 边 2
        1 5; % 边 3
        2 3; % 边 4
        2 6; % 边 5
        3 4; % 边 6
        3 7; % 边 7
        4 8; % 边 8
        5 6; % 边 9
        5 8; % 边 10
        6 7; % 边 11
        7 8]; % 边 12
range = 1:12;
midEdge = (node3(edge(range,1),:)+node3(edge(range,2),:))/2;
hold 
plot3(midEdge(:,1),midEdge(:,2),midEdge(:,3),'s','LineWidth',1,'MarkerEdgeColor','k',...
    'MarkerFaceColor',[0.6 0.5 0.8],'MarkerSize',20);
text(midEdge(:,1)-0.025,midEdge(:,2),midEdge(:,3),int2str(range(:)), ...
    'FontSize',12,'FontWeight','bold','Color','k');

elemface2edge = cell(6,1);
elemface2edge{1} = [10 12 11 9];
elemface2edge{2} = [1  4  6  2];
elemface2edge{3} = [8  10 3  2];
elemface2edge{4} = [4  5  11 7];         
elemface2edge{5} = [3  9  5  1];         
elemface2edge{6} = [6  7  12 8];

% numbers
Nv = size(node3,1); 
Ne = size(edge,1); 
Nf = length(elemf); 

% Scaled monomials
Nm = 10;
m1 = @(x,y,z)  1+0*x;                      gradm1 = @(x,y,z)  [0+0*x,         0+0*x,              0+0*x];
m2 = @(x,y,z)  (x-xK)/hK;                  gradm2 = @(x,y,z)  [1/hK+0*x,      0+0*x,              0+0*x];
m3 = @(x,y,z)  (y-yK)/hK;                  gradm3 = @(x,y,z)  [0+0*x,         1/hK+0*x,           0+0*x];
m4 = @(x,y,z)  (z-zK)/hK;                  gradm4 = @(x,y,z)  [0+0*x,         0+0*x,              1/hK+0*x];    
m5 = @(x,y,z)  (x-xK).^2/hK^2;             gradm5 = @(x,y,z)  [2*(x-xK)/hK^2, 0+0*x,              0+0*x];
m6 = @(x,y,z)  (x-xK).*(y-yK)/hK^2;        gradm6 = @(x,y,z)  [(y-yK)/hK^2,   (x-xK)/hK^2,        0+0*x];
m7 = @(x,y,z)  (y-yK).^2/hK^2;             gradm7 = @(x,y,z)  [0+0*x,         2*(y-yK)/hK^2,      0+0*x];
m8 = @(x,y,z)  (x-xK).*(z-zK)/hK^2;        gradm8 = @(x,y,z)  [(z-zK)/hK^2,   0+0*x,              (x-xK)/hK^2];
m9 = @(x,y,z)  (y-yK).*(z-zK)/hK^2;        gradm9 = @(x,y,z)  [0+0*x,         (z-zK)/hK^2,        (y-yK)/hK^2];
m10 = @(x,y,z) (z-zK).^2/hK^2;             gradm10 = @(x,y,z) [0+0*x,         0+0*x,              2*(z-zK)/hK^2];
m = @(x,y,z) [m1(x,y,z), m2(x,y,z), m3(x,y,z), m4(x,y,z), m5(x,y,z),...
              m6(x,y,z), m7(x,y,z), m8(x,y,z), m9(x,y,z), m10(x,y,z)];
mc = {m1,m2,m3,m4,m5,m6,m7,m8,m9,m10};
gradmc = {gradm1,gradm2,gradm3,gradm4,gradm5,gradm6,gradm7,gradm8,gradm9,gradm10};   

%% G = BD
G = zeros(Nm);
for i = 1:Nm
    for j = 1:Nm
        fun = @(x,y,z) sum(gradmc{i}(x,y,z).*gradmc{j}(x,y,z),2);
        G(i,j) = integralPolyhedron(fun,5,node3,elemf);
    end
end

%% D
Ndof = Nv + Ne + Nf + 1;
index3 = 1:8;
x = node3(index3,1);  y = node3(index3,2);  z = node3(index3,3);
xe = (node3(edge(:,1),1) + node3(edge(:,2),1))/2;
ye = (node3(edge(:,1),2) + node3(edge(:,2),2))/2;
ze = (node3(edge(:,1),3) + node3(edge(:,2),3))/2;
D = zeros(Ndof,Nm);
% 顶点值、边中点值
D(1:Nv+Ne,:) = [m(x,y,z); m(xe,ye,ze)];  
% 面的矩量
for s = 1:size(elemf,1)
    % vertices of face
    face = elemf{s};  P = node3(face,:);  
    [Intm,areaP] = integralPolygon3(m,5,P); 
    D(Nv+Ne+s,:) = Intm/areaP;  
end
% 体的矩量
D(end,:) = 1/volume*integralPolyhedron(m,5,node3,elemf);

%% B
% 体积分
Lapm = zeros(Nm,1); Lapm([5,7,10]) = 2/hK^2;
Dof = [zeros(1,Ndof-1), volume];
I1 = -Lapm*Dof;
% 面积分
I2 = zeros(Nm,Ndof);
for s = 1:size(elemf,1)
    % face info
    face = elemf{s};  P = node3(face,:);  
    face2edge = elemface2edge{s};
    % elliptic projection on face (=L2 projection)
    Intf = faceEllipticProjectionk2(P,gradmc);
    id = [face, face2edge+Nv, s+Nv+Ne];
    I2(:,id) = I2(:,id) + Intf;
end
% B
B = I1 + I2;
% G
G1 = B*D;


aa = G-G1;










