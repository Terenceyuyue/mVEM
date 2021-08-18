function [cp,ce,pv,ev] = makedual2(pp,tt,varargin)

% pv: 对偶网格的顶点坐标
% ev: 对偶网格边的顶点序号
% 对偶网格的特点是: 每一个顶点相当于重心，对应一个单元 (可能存在融合的)

% ce: 按单元拉直模式记录边的编号, 即 elem2edge 拉直
% cp: - NT * 3
%     - 2,3 列给出每个单元边的自然序号排列，而 ce 是整体编号 （有了整体编号就可以获得边的顶点编号，即通过 ev）
%     - 1 列给出单元对应的重心编号（原节点）

% 综上，ce ----- elem2edge 拉直
%       
%   Darren Engwirda : 2014 --
%   Email           : engwirda@mit.edu
%   Last updated    : 17/12/2014

%----------------------------------- extract optional inputs
    ec = []; op = [];

%---------- x-y 平面数据嵌入到 R^3 -------------
    if (size(pp,2) == +2)
        pp = [pp, zeros(size(pp,1),1)] ;
    end
    op.etol = +1e-8; % bound merge tolerance
    op.atol = 1/2. ; % bound angle tolerance
 
    np = size(pp,1);
    nt = size(tt,1);
%------------ 三角剖分相邻单元 --------------
% ee: edge
% ep: e1, e2,e2, e3, .... 记录顺序
% et: 按照上面顺序记录单元序号，即 neighbor
    [ee,~,ep,et] = triaconn2(tt);
%------------------------------------ mark nonmanifold edges      
    [ke] = markedge(ee,ep,et,ec,pp,tt,op);
%------------------------------------ form tria circum-balls
    [cc] = miniball2(pp,tt); % xr,yr,zr,r
%-------------------------------------------- edge midpoints    
    [pm] = (pp(ee(:,1),:) + pp(ee(:,2),:))/2;
%-------------------------------------------- node positions
    [pv] = [pp; cc(:,1:3); pm];
%-------------------------------- make edges in dual complex
    [ev] = makeedge(ee,ep,et,ke,pv,np,nt);
%-------------------------------------- remove "short" edges
    [pv,ev] = removeEdges(pv,ev,op);
%-------------------------------------- remove un-used nodes
    [pv,ev] = removeNodes(pv,ev,np);
%-------------------------------- make cells in dual complex
    [cp,ce] = clipcell(ev)  ;
%--------------------------------------- dump extra indexing
    [ev] = ev(:,1:2);
   
end


