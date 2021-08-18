% x==-1, y == 3
sdpvar x y
F = set(abs(abs(x+1)+3) < y)
sol = solvesdp(F,y)

% x==3, y==7
sdpvar x y
bounds(x,0,10)
F = set(abs(abs(x+1)+3) > y)+set(0<x<3);
sol = solvesdp(F,-y)
