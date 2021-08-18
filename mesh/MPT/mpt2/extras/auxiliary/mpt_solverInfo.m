function [out,err]=mpt_solverInfo(sclass, stype)
%MPT_SOLVERINFO returns information about a given solver
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
%
% internal function
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
%

% Copyright is with the following author(s):
%
% (C) 2004 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
%          kvasnica@control.ee.ethz.ch

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

error(nargchk(2,2,nargin));


if ischar(stype),
    % convert from string to integer
    % e.g. NAG -> 3
    
    if strcmpi(sclass, 'lp')
        convfunc = @char2int_lpsolver;
    elseif strcmpi(sclass, 'qp')
        convfunc = @char2int_qpsolver;
    elseif strcmpi(sclass, 'milp')
        convfunc = @char2int_milpsolver;
    elseif strcmpi(sclass, 'miqp')
        convfunc = @char2int_miqpsolver;
    elseif strcmpi(sclass, 'extreme')
        convfunc = @char2int_extreme;
    else
        error(sprintf('unknwon solver class ''%s'' !', sclass));
    end
else
    % convert from integer to string
    % e.g. 0 -> NAG
    if strcmpi(sclass, 'lp')
        convfunc = @int2char_lpsolver;
    elseif strcmpi(sclass, 'qp')
        convfunc = @int2char_qpsolver;
    elseif strcmpi(sclass, 'milp')
        convfunc = @int2char_milpsolver;
    elseif strcmpi(sclass, 'miqp')
        convfunc = @int2char_miqpsolver;
    elseif strcmpi(sclass, 'extreme')
        convfunc = @int2char_extreme;
    else
        error(sprintf('unknwon solver class ''%s'' !', sclass));
    end
end
    
[out, err] = feval(convfunc,stype);
return


%------------------------------------------------------------------------
function [out,err] = int2char_lpsolver(solver);

err = 0;
switch solver,
    case 0, out = 'NAG (e04naf)';
    case 1, out = 'linprog';
    case 2, out = 'CPLEX9';
    case 3, out = 'CDD Criss-Cross';
    case 4, out = 'GLPK';
    case 5, out = 'CDD Dual-Simplex';
    case 6, out = 'SeDuMi';
    case 7, out = 'QSopt';
    case 8, out = 'CPLEX8';
    case 9, out = 'NAG (e04mbf)';
    case 10, out = 'XPRESS';
    case 11, out = 'MOSEK';
    case 12, out = 'OOQP';
    case 13, out = 'CLP';
    case 14, out = 'BPMPD';
    case 15, out = 'CPLEXMEX';
    case 16, out = 'PDCO';
    case 17, out = 'GLPKCC';
    otherwise, out = 'unknown'; err=1;
end


%------------------------------------------------------------------------
function [out,err] = int2char_qpsolver(solver);

err = 0;
switch solver,
    case -1, out = 'not available';
    case 0, out = 'NAG';
    case 1, out = 'quadprog';
    case 2, out = 'CPLEX9';
    case 3, out = 'SeDuMi';
    case 4, out = 'CPLEX8';        
    case 5, out = 'XPRESS';
    case 6, out = 'MOSEK';
    case 7, out = 'OOQP';
    case 8, out = 'CLP';
    case 9, out = 'BPMPD';
    case 10, out = 'CPLEXMEX';
    otherwise, out = 'unknown'; err=1;        
end


%------------------------------------------------------------------------
function [out,err] = int2char_milpsolver(solver);

err = 0;
switch solver,
    case 0, out = 'CPLEX9';
    case 1, out = 'YALMIP';
    case 2, out = 'GLPK';
    case 3, out = 'XPRESS';
    case 4, out = 'MOSEK';
    case 5, out = 'bintprog';
    case 6, out = 'CPLEX8';
    case 7, out = 'CPLEXMEX';
    case 8, out = 'GLPKCC';
    otherwise, out = 'unknown'; err=1;        
end


%------------------------------------------------------------------------
function [out,err] = int2char_miqpsolver(solver);

err = 0;
switch solver,
    case 0, out = 'CPLEX9';
    case 1, out = 'YALMIP';
    case 2, out = 'XPRESS';
    case 3, out = 'MOSEK';        
    case 4, out = 'CPLEX8';
    case 5, out = 'CPLEXMEX';
    otherwise, out = 'unknown'; err=1;        
end


%------------------------------------------------------------------------
function [out,err] = int2char_extreme(solver);

err = 0;
switch solver,
    case 0, out = 'Analytical method';
    case 1, out = 'LRS';
    case 2, out = 'Analytical method (alternative)';
    case 3, out = 'CDD';
    case 4, out = 'CDD without reduction';
    otherwise, out = 'unknown'; err=1;        
end


%------------------------------------------------------------------------
function [out,err] = char2int_lpsolver(solver);
% converts a string description to numerical values

err = 0; out = -Inf;
if strcmpi(solver,'cdd') | strcmpi(solver, 'cdd criss-cross'),
    out=3;
elseif strcmpi(solver, 'cdd dual-simplex')
    out = 5;
elseif strcmpi(solver,'nag') | strcmpi(solver,'e04naf') | strcmpi(solver, 'NAG (e04naf)'),
    out=0;
elseif strcmpi(solver,'cplex') | strcmpi(solver,'cplex9'),
    out=2;
elseif strcmpi(solver,'glpk'),
    out=4;
elseif strcmpi(solver,'linprog'),
    out=1;
elseif strcmpi(solver,'sedumi'),
    out=6;
elseif strcmpi(solver,'qsopt'),
    out=7;
elseif strcmpi(solver,'cplex8'),
    out=8;
elseif strcmpi(solver,'e04mbf') | strcmpi(solver, 'NAG (e04mbf)'),
    out=9;
elseif strcmpi(solver,'xpress'),
    out = 10;
elseif strcmpi(solver,'mosek'),
    out = 11;
elseif strcmpi(solver,'ooqp'),
    out = 12;
elseif strcmpi(solver,'clp'),
    out = 13;
elseif strcmpi(solver,'bpmpd'),
    out = 14;
elseif strcmpi(solver,'cplexmex'),
    out = 15;
elseif strcmpi(solver,'pdco'),
    out = 16;
elseif strcmpi(solver, 'glpkcc')
    out = 17;
elseif strcmpi(solver, 'fastest available')
    out = [];
else
    err = 1;
end


%------------------------------------------------------------------------
function [out,err] = char2int_qpsolver(solver);
% converts a string description to numerical values

err = 0; out = -Inf;
if strcmpi(solver,'quadprog'),
    out=1;
elseif strcmpi(solver,'nag'),
    out=0;
elseif strcmpi(solver,'cplex') | strcmpi(solver,'cplex9'),
    out=2;
elseif strcmpi(solver,'sedumi'),
    out=3;
elseif strcmpi(solver,'cplex8'),
    out=4;
elseif strcmpi(solver,'xpress'),
    out = 5;
elseif strcmpi(solver,'mosek'),
    out = 6;
elseif strcmpi(solver,'ooqp'),
    out = 7;
elseif strcmpi(solver,'clp'),
    out = 8;
elseif strcmpi(solver,'bpmpd'),
    out = 9;
elseif strcmpi(solver,'cplexmex'),
    out = 10;
elseif strcmpi(solver, 'fastest available')
    out = [];
else
    err = 1;
end


%------------------------------------------------------------------------
function [out,err] = char2int_extreme(solver);
% converts a string description to numerical values

err = 0; out = -Inf;
if strcmpi(solver,'cdd'),
    out=3;
elseif strcmpi(solver,'lrs'),
    out = 1;
elseif strcmpi(solver, 'matlab_alt'),
    out = 2;
elseif strcmpi(solver,'matlab') | strcmpi(solver, 'analytical method'),
    out=0;
elseif strcmpi(solver, 'fastest available')
    out = [];
else
    err = 1;
end


%------------------------------------------------------------------------
function [out,err] = char2int_milpsolver(solver);
% converts a string description to numerical values

err = 0; out = -Inf;
if strcmpi(solver,'cplex') | strcmpi(solver, 'cplex9'),
    out=0;
elseif strcmpi(solver,'yalmip'),
    out=1;
elseif strcmpi(solver,'glpk'),
    out=2;
elseif strcmpi(solver,'xpress'),
    out = 3;
elseif strcmpi(solver,'mosek'),
    out = 4;
elseif strcmpi(solver,'bintprog'),
    out = 5;
elseif strcmpi(solver,'cplex8'),
    out = 6;
elseif strcmpi(solver,'cplexmex'),
    out = 7;
elseif strcmpi(solver, 'glpkcc'),
    out = 8;
elseif strcmpi(solver, 'fastest available')
    out = [];
else
    err = 1;
end


%------------------------------------------------------------------------
function [out,err] = char2int_miqpsolver(solver);
% converts a string description to numerical values

err = 0; out = -Inf;
if strcmpi(solver,'cplex') | strcmpi(solver, 'cplex9'),
    out=0;
elseif strcmpi(solver,'yalmip'),
    out=1;
elseif strcmpi(solver,'xpress'),
    out = 2;
elseif strcmpi(solver,'mosek'),
    out = 3;
elseif strcmpi(solver,'cplex8'),
    out = 4;
elseif strcmpi(solver,'cplexmex'),
    out = 5;
elseif strcmpi(solver, 'fastest available')
    out = [];
else
    err = 1;
end
