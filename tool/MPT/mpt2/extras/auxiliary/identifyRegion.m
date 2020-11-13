function identifyRegion(P,Idx)
%IDENTIFYREGION plots the number of the polytopes into the current figure
%
% identifyRegion(P)
% identifyRegion(P,Idx)
%
% -------------------------------------------------------------------------
% INPUT
% -------------------------------------------------------------------------
% P                 - polytope array
% Idx               - vector of index numbers to be printed (optional)
%
% see also POLYTOPE/PLOT
%

% Copyright is with the following author(s):
%
% (C) 2004-2005 Frank J. Christophersen, Automatic Control Lab., ETH Zurich,
%               fjc@control.ee.ethz.ch

% -------------------------------------------------------------------------
% Legal note:
%          This program is free software; you can redistribute it and/or
%          modify it under the terms of the GNU General Public
%          License as published by the Free Software Foundation; either
%          version 2.1 of the License, or (at your option) any later 
%          version.
%
%          This program is distributed in the hope that it will be useful,
%          but WITHOUT ANY WARRANTY; without even the implied warranty of
%          MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
%          GNU General Public License for more details.
% 
%          You should have received a copy of the GNU General Public
%          License along with this library; if not, write to the 
%          Free Software Foundation, Inc., 
%          59 Temple Place, Suite 330, 
%          Boston, MA  02111-1307  USA
%
% -------------------------------------------------------------------------

error(nargchk(1,2,nargin));

global mptOptions

if ~isstruct(mptOptions),
    mpt_error;
end

if ~isa(P, 'polytope'),
    error('IDENTIFYREGION: input argument MUST be a polytope');
end
if nargin<2
    Idx = 1:length(P);
elseif nargin==2
    if length(P)<length(Idx) | min(Idx)<=0 | max(Idx) > length(P)
        error('IDENTIFYREGION: Idx must represent a possible list of indices of P.');
    end
end

% we use the Idx variable in a for-cycle, thus we have to make sure that it is a
% row vector
Idx = Idx(:)';

hold on;
h = gcf;

n = dimension(P);

for ii= Idx,
  [xc, Rc] = chebyball(P(ii));
  
  if n==2
    h2 = text(xc(1), xc(2), num2str(ii));
  elseif n==3
    h2 = text(xc(1), xc(2), xc(3), num2str(ii));
  else
    error('IDENTIFYREGION: only for 2D and 3D polytopes');
  end
  
  set(h2,'Color','k');
  set(h2,'FontSize',12);

end%ii

hold off;
