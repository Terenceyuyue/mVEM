function Pb = mpt_rectangle_domain(BdBox)

Nbox = length(BdBox); dim = Nbox/2;
if dim==2
    [x,y] = meshgrid(BdBox(1:2),BdBox(3:4));
    Vb = [x(:), y(:)];
end
if dim==3
    [x,y,z] = meshgrid(BdBox(1:2),BdBox(3:4),BdBox(5:6));
    Vb = [x(:), y(:), z(:)];
end
Pb = polytope(Vb);