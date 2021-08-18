function [ev] = makeedge(ee,ep,et,ke,pv,np,nt)
%  assemble edge segments in cells of dual complex

nn = ep(:,2)-ep(:,1)+1 ; % 边的度
ev = zeros(sum(nn)*3,4);
% ev 有 四列，前两列对应对偶边，后两列对应对偶边的三角形垂直边 (含重复边)

ne = +0;
for i = 1 : size(ee,1)
    % adj. nodes
    ni = ee(i,1); nj = ee(i,2); 
    mi = i+nt+np; % mi 表示接着顶点和重心的顺序后编号
    % assemble dual cells
    if (ke(i) == +0) % 内部边
        % adj. trias
        ti = et(ep(i,1)+0)+np;   tj = et(ep(i,1)+1)+np; % 重心接着顶点后编号
        % add dual tria-tria segment
        ne = ne+1;
        ev(ne,1) = ti ; ev(ne,2) = tj ; % 两个三角形垂直平分线，三角形重心是顶点
        ev(ne,3) = ni ; ev(ne,4) = nj ; % 垂直边的起点终点
        
    else % 边界边
        % add dual node-edge segment
        ne = ne+1;
        ev(ne,1) = mi ; ev(ne,2) = ni ;
        ev(ne,3) = ni ; ev(ne,4) = +0 ;
        ne = ne+1;
        ev(ne,1) = nj ; ev(ne,2) = mi ;
        ev(ne,3) = nj ; ev(ne,4) = +0 ;
        for ti = ep(i,1) : ep(i,2)
            % add dual tria-edge segment
            tc = et(ti); tc = tc+np ;
            ne = ne+1;
            ev(ne,1) = tc ; ev(ne,2) = mi ;
            ev(ne,3) = ni ; ev(ne,4) = nj ;
        end
    end
end
%  trim allocation per demand
ev = ev(1:ne,:);