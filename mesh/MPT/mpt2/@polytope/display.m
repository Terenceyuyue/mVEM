function sys = display(P)
%DISPLAY Displays details about the given polytope
%
% sys = display(P)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Displays details of the given POLYTOPE object.
% 
% This method is executed automatically if one types:
%   >> P
% without the semicolumn (;) at the end
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% P    - Polytope
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% sys  - String containing the information
%
% see also POLYTOPE, DIMENSION, NCONSTR, LENGTH
%

% Copyright is with the following author(s):
%
% (C) 2003 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
%          kvasnica@control.ee.ethz.ch
% (C) 2003 Mato Baotic, Automatic Control Laboratory, ETH Zurich,
%          baotic@control.ee.ethz.ch

% ---------------------------------------------------------------------------
% Legal note:
%          This program is free software; you can redistribute it and/or
%          modify it under the terms of the GNU General Public
%          License as published by the Free Software Foundation; either
%          version 2.1 of the License, or (at your option) any later version.
%
%          This program is distributed in the hope that it will be useful,
%          but WITHOUT ANY WARRANTY; without even the implied warranty of
%          MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%          General Public License for more details.
% 
%          You should have received a copy of the GNU General Public
%          License along with this library; if not, write to the 
%          Free Software Foundation, Inc., 
%          59 Temple Place, Suite 330, 
%          Boston, MA  02111-1307  USA
%
% ---------------------------------------------------------------------------

sys = [];

disp([inputname(1) '=']);
fprintf('\n');

if ~isempty(P.Array),
    % in case of P being a polyarray:
    fprintf('Polytope array: %d polytopes in %dD\n\n',length(P.Array),size(P.Array{1}.H,2));
        
else
    % in case of P being a polytope:
    str=sprintf('%dx%d polytope',size(P.H,1),size(P.H,2));
    sub_display(P,'');
end
return


%---------------------------------------------------------
function sub_display(P,str)

if P.H==1 & P.K==-Inf & P.RCheb==-Inf
    disp([str 'Empty polytope in R^1']);
    return
end

[nc,nx]=size(P.H);

switch P.normal
case 0
    str=[str 'Not normalized, '];
case 1
    str=[str 'Normalized, '];
otherwise
    str=[str '(?) normalized, '];
end
switch P.minrep
case 0
    str=[str 'non minimal representation '];
case 1
    str=[str 'minimal representation '];
otherwise
    str=[str '(?) minimal representation '];
end

disp([str 'polytope in R^' num2str(nx)]);

% make a copy of the polytope for displaying purposes
Q.H = P.H;
Q.K = P.K;
Q.normal = P.normal;
Q.minrep = P.minrep;
Q.xCheb = P.xCheb;
Q.RCheb = P.RCheb;
% display the structure
disp(struct(Q));

if size(P.H,1)<40 & size(P.H,2)<=7,
    % prepare display of Hx<=K
    % only if # of constraints is less than 100 and dimension is smaller than 7
    lb=repmat('[',size(P.H,1),1);
    rb=repmat(']',size(P.H,1),1);
    sp=repmat(' ',size(P.H,1),1);
    xs=sp;
    le=repmat('  ',size(P.H,1),1);
    xs(floor(size(P.H,1)/2)+1)='x';
    le(floor(size(P.H,1)/2)+1,:)='<=';
    disp([lb num2str(P.H) rb sp xs sp le sp lb num2str(P.K) rb])
end
fprintf('\n');
