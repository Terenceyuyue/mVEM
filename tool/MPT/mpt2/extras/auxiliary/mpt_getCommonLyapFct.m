function [lP,feasible]=mpt_getCommonLyapFct(Fi,Ain,Bin,Options)
%MPT_GETCOMMONLYAPFCT Computes common Lyapunov function for PWA system
%
% [LP,feasible]=mpt_getCommonLyapFct(Fi,Ain,Bin,Options)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% This function attempts to compute a quadratic Lyapunov function V(x)=x'lPx
% which guarantees exponential stability.
% PWQ(x) = x'LPx
% PWQ(x(k+1)) - PWQ(x(k)) <= rho * x^2 
%
% (i.e. rho must be negative to guarantee exponential stability)
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% Fi        - optimal feedback law given by u=Fi{i}*x
% Ain, Bin  - system dynamics are given by x(k+1)=Ain{l}x(k)+Bin{l}u(k)
%             In case of polytopic uncertainty, you may pass a whole cell array
%             with different dynamics. This function will then attempt to
%             identify a Lyapunov function with negative decay for all dynamics
%             in the cell array. 
%
% Options.abs_tol     - Absolute tolerance
% Options.lpsolver    - Which LP solver to use (see help mpt_solveLP)
% Options.debug_level - If this is set to 1, the solution provided by the LMI 
%			            solver will be double-checked manually. We strongly 
%		                advise to set this to 1, since we've experienced 
%			            numerous numerical issues with certain LMI solvers.
%
% Note: If Options is missing or some of the fields are not defined, the default
%       values from mptOptions will be used
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% LP         - Quadratic Lyapunov function: Q(x)=x'LP{r}x%    
% feasible   - 1: stable 0: no statement about stability possible
%
% ---------------------------------------------------------------------------
% LITERATURE
% ---------------------------------------------------------------------------
% Automatica, Volume 38, issue 12, pp. 2139 - 2146 
% "Analysis of discrete-time piecewise affine and hybrid systems",
% Ferrari-Trecate G., F.A. Cuzzola, D. Mignone and M. Morari, 
%
% see also MPT_GETPWQLYAPFCT, MPT_GETSTABFEEDBACK

% Copyright is with the following author(s):
%
% (C) 2003 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
%          kvasnica@control.ee.ethz.ch
% (C) 2003 Pascal Grieder, Automatic Control Laboratory, ETH Zurich,
%          grieder@control.ee.ethz.ch

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

error(nargchk(3,4,nargin));

global mptOptions;

if ~isstruct(mptOptions),
    mpt_error;
    return
end
if nargin<4,
    Options=mptOptions;
end
if ~isfield(Options,'abs_tol'),
    Options.abs_tol=mptOptions.abs_tol;
end
if ~isfield(Options,'debug_level'),
    Options.debug_level=mptOptions.debug_level;
end
if ~isfield(Options,'lpsolver'),
    Options.lpsolver=mptOptions.lpsolver;
end
if ~isfield(Options,'verbose'),
    Options.verbose=mptOptions.verbose;
end

if(~iscell(Ain) | ~iscell(Bin))
    tmpA=Ain;
    tmpB=Bin;
    clear Ain Bin
    Ain{1}=tmpA;
    Bin{1}=tmpB;
    clear tmpA tmpB
end
if(length(Fi)>length(Ain))
    disp('mpt_getCommonLyapFct: The number of feedback laws is not identical to the number of dynamics !')  
    if(length(Ain)==1)
        disp('mpt_getCommonLyapFct: Associating Fi{j} to dynamics A{1}, B{1}.')
        disp('mpt_getCommonLyapFct: Augmenting dynamics vector...')
        for i=1:length(Fi)
            tmpA{i}=Ain{1};
            tmpB{i}=Bin{1};
        end
        clear Ain Bin
        Ain=tmpA;
        Bin=tmpB;
        clear tmpA tmpB
    else
        error('mpt_getCommonLyapFct: Dynamics not clearly associated to feedback laws !')
    end
end
if(length(Fi)<length(Ain))
    error('mpt_getCommonLyapFct: Dynamics not clearly associated to feedback laws !')
end
% if(~exist('sedumi'))
%     disp('mpt_getCommonLyapFct:')
%     disp('You need to download and install SeDuMi for this script to work...') %...and set the path as well
%     error('The package can be downloaded from http://fewcal.kub.nl/sturm/software/sedumi.html')
% end
try
    yalmip('clear') ;       %initialize yalmip
catch
    error('mpt_getCommonLyapFct: You need to download and install Yalmip for this function to work');
end
Acell=Ain;
Bcell=Bin;

%load parameters
n = size(Acell{1},2);   %no of states
noU= size(Bcell{1},2);  %no of inputs
binaryOne=dec2bin(1);   %binary one: initialize for combinatorial search later
lookahead=1;            %only examine reachability for 1 step ahead


%---------------------------------------------------
%Initialize LMI Variables
%Origin must be contained in region 1
%PWQ Lyapunov function: x'Q{i}x+x'L{i}+C{i}      if x \in Pn(i)
P = sdpvar(n,n,'symmetric'); %initialize PWQ Lyapunov variables

rho=sdpvar(1,1);   %optimization variable for exponential stability
%---------------------------------------------------
epsilon=Options.abs_tol;
myprog=lmi;
for j=1:length(Fi)
    ABK=Ain{j}+Bin{j}*Fi{j}(1:noU,:);
    myprog = myprog + set('ABK''*P*ABK-P<rho*eye(n)');
end
myprog = myprog + set('rho>-10'); %bound from below

%options=sdpsettings('Silent',(Options.verbose==0),'Solver','SEDUMI');
if isempty(mptOptions.sdpsettings),
    options=sdpsettings('Verbose',Options.verbose);
else
    options=mptOptions.sdpsettings;
end
%%switch off output and use SEDUMI
options.sedumi.stepdif=0;
options.savesolverinput=1;
%solution = solvesdp(myprog,[],rho,options);                   %find solution using LMI solver
solution = solvesdp(myprog,rho,options);                   %find solution using LMI solver

lP=double(P);
feasible=1;
if(solution.problem~=0 | double(rho)>=0)
    disp('Trying again...')
    myprog = myprog + set('rho<-epsilon');%objective to be minimized
    solution = solvesdp(myprog,[],[],options);  
    if(solution.problem~=0)
        fprintf('\n\n');
        disp('mpt_getCommonLyapFct: NO SOLUTION FOUND')
        feasible=0;
    end
end