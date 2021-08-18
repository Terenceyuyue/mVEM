%===============================================================================
% function [S, W] = mpc_checkSysWeight(S, W, verbose)
%
% Title:        checkMLD.m                                             	
%                                                                       
% Project:      Control of MLD systems
%                                                             
% Purpose:      Given the MLD system in the MLD structure and the weights, 
%               checks the dimensions and adds missing 0-matrices 
%               
% Author:       Domenico Mignone, Tobias Geyer
%                                                                       
% History:      date         subject                                       
%               2003.11.19   first public release 
%               2003.03.18   initial version
% 
% Usage:        Given the MLD system in the MLDsystem structure, checks the 
%               dimensions and adds missing 0-matrices. Thus, the dimensions of 
%               all the MLD systems in the structure need to be the same.
%
% Input:        S: Structure or cell structure with MLD system descriptions
%               W: Structure or cell structure with the weight matrices
%               verbose: 0: don't display any messages
%                        1: display messages
%
% Output:       S: Updated cell structure with MLD system descriptions
%               W: Updated cell strucutre with the weights
%               Ext: additional outputs with the following fields:
%                    ... not used yet
%
% Assumptions:  A valid models contains A and B1
%               The matrices B2, B3, B5, C, D1, D2, D3, D5, E1, E2, E3, E4 and E5
%               are NOT mandatory. Therefore, the whole framework can also be used
%               for discrete-time linear systems.
%
% Contact:      Tobias Geyer
%               Automatic Control Laboratory
%               ETH Zentrum, 
%               Zurich, Switzerland
%
%               geyer@aut.ee.ethz.ch
%
%               Comments and bug reports are highly appreciated
%

%===============================================================================
%
% Legal note:   This program is free software; you can redistribute it and/or
%               modify it under the terms of the GNU General Public License as 
%               published by the Free Software Foundation; either version 2.1 of
%               the License, or (at your option) any later version. 
%
%               This library is distributed in the hope that it will be useful,
%               but WITHOUT ANY WARRANTY; without even the implied warranty of
%               MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%               Lesser General Public License for more details.
% 
%               You should have received a copy of the GNU Lesser General Public
%               License along with this library; if not, write to the 
%               Free Software Foundation, Inc., 
%               59 Temple Place, Suite 330, 
%               Boston, MA  02111-1307  USA
%
%===============================================================================

function [S, W] = mpc_checkSysWeight(S, W, verbose)

if ~iscell(S), 
    % define it as a cell structure
    S_temp = S;
    S = [];
    S{1} = S_temp;
end;

if ~iscell(W), 
    % define it as a cell structure
    W_temp = W;
    W = {};
    W{1} = W_temp;
end;


% Take the dimensions of the first MLD system:
% --------------------------------------------

% nu=dimension of u
if isfield(S{1}, 'B1')
    nu = size(S{1}.B1,2);   
else
    error('invalid MLD model: B1 is not defined')   
end

% nd=dimension of \delta
if isfield(S{1}, 'B2')
    nd = size(S{1}.B2,2);
else
    if verbose, disp('assume MLD model has no delta variables'); end;
    nd = 0;
end

% nz=dimension of z
if isfield(S{1}, 'B3')
    nz = size(S{1}.B3,2);
else
    if verbose, disp('assume MLD model has no z variables'); end;
    nz = 0;
end

% nx=dimension of x
if isfield(S{1}, 'A')
    nx = size(S{1}.A,2);    
else
    error('invalid MLD model: A is not defined')
end

% ne=dimension of E_5
if isfield(S{1}, 'E5')
    ne = size(S{1}.E5,1);   
else
    if verbose, disp('assume MLD model has no constraints'); end;
    ne = 0;
end

% ny=number of outputs
if isfield(S{1}, 'C')
    ny = size(S{1}.C,1);    
else
    if verbose, disp('assume MLD model has not outputs'); end;
    ny = 0;
end 



% * check, if all the MLD models have the above defined dimensions and
% * set missing (not mandatory) matrices to zero
% --------------------------------------------------------------------

for i = 1:length(S)

    if isfield(S{i}, 'A')
        if ~all( size(S{i}.A) == [nx nx] )
            error('invalid MLD model: A has wrong dimension');
        end;
    else
        error('invalid MLD model: A is not defined')       
    end 
    
    if isfield(S{i}, 'B1')
        if ~all( size(S{i}.B1) == [nx nu] )
            error('invalid MLD model: B1 has wrong dimension');
        end;
    else
        error('invalid MLD model: B1 is not defined')   
    end
    
    if isfield(S{i}, 'B2')
        if ~all( size(S{i}.B2) == [nx nd] )
            error('invalid MLD model: B2 has wrong dimension');
        end;
    else
        if verbose, disp('MLD model: B2 is not defined -- assume zeros'); end;
        S{i}.B2 = zeros(nx, nd);  
    end   
    
    if isfield(S{i}, 'B3')
        if ~all( size(S{i}.B3) == [nx nz] )
            error('invalid MLD model: B3 has wrong dimension');
        end;
    else
        if verbose, disp('MLD model: B3 is not defined -- assume zeros'); end;
        S{i}.B3 = zeros(nx, nz);
    end   
    
    if isfield(S{i}, 'B5')
        if ~all( size(S{i}.B5) == [nx 1] )
            error('invalid MLD model: B5 has wrong dimension');
        end;        
    else
        if verbose, disp('MLD model: B5 is not defined -- assume zeros'); end;
        S{i}.B5 = zeros(nx, 1);
    end 
    
    if isfield(S{i}, 'C')
        if ~all( size(S{i}.C) == [ny nx] )
            error('invalid MLD model: C has wrong dimension');
        end;
    else
        if verbose, disp('MLD model: C is not defined -- assume zeros'); end;
        S{i}.C = zeros(ny, nx);
    end 
    
    if isfield(S{i}, 'D1')
        if ~all( size(S{i}.D1) == [ny nu] )
            error('invalid MLD model: D1 has wrong dimension');
        end;
    else
        if verbose, disp('MLD model: D1 is not defined -- assume zeros'); end;
        S{i}.D1 = zeros(ny, nu);  
    end
    
    if isfield(S{i}, 'D2')
        if ~all( size(S{i}.D2) == [ny nd] )
            error('invalid MLD model: D2 has wrong dimension');
        end;
    else
        if verbose, disp('MLD model: D2 is not defined -- assume zeros'); end;
        S{i}.D2 = zeros(ny, nd); 
    end   
    
    if isfield(S{i}, 'D3')
        if ~all( size(S{i}.D3) == [ny nz] )
            error('invalid MLD model: D3 has wrong dimension');
        end;
    else        
        if verbose, disp('MLD model: D3 is not defined -- assume zeros'); end;
        S{i}.D3 = zeros(ny, nz);
    end   
    
    if isfield(S{i}, 'D5')
        if ~all( size(S{i}.D5) == [ny 1] )
            error('invalid MLD model: D5 has wrong dimension');
        end;
    else        
        if verbose, disp('MLD model: D5 is not defined -- assume zeros'); end;
        S{i}.D5 = zeros(ny, 1);
    end 
    
    if isfield(S{i}, 'E1')
        if ~all( size(S{i}.E1) == [ne nu] )
            error('invalid MLD model: E1 has wrong dimension');
        end;
    else        
        if verbose, disp('MLD model: E1 is not defined -- assume zeros'); end;
        S{i}.E1 = zeros(ne, nu);
    end  
    
    if isfield(S{i}, 'E2')
        if ~all( size(S{i}.E2) == [ne nd] )
            error('invalid MLD model: E2 has wrong dimension');
        end;
    else                
        if verbose, disp('MLD model: E2 is not defined -- assume zeros'); end;
        S{i}.E2 = zeros(ne, nd);
    end   
    
    if isfield(S{i}, 'E3')
        if ~all( size(S{i}.E3) == [ne nz] )
            error('invalid MLD model: E3 has wrong dimension');
        end;
    else                
        if verbose, disp('MLD model: E3 is not defined -- assume zeros'); end;
        S{i}.E3 = zeros(ne, nz);
    end
    
    if isfield(S{i}, 'E4')
        if ~all( size(S{i}.E4) == [ne nx] )
            error('invalid MLD model: E4 has wrong dimension');
        end;
    else                
        if verbose, disp('MLD model: E4 is not defined -- assume zeros'); end;
        S{i}.E4 = zeros(ne, nx);
    end  
    
    if isfield(S{i}, 'E5')
        if ~all( size(S{i}.E5) == [ne 1] )
            error('invalid MLD model: E5 has wrong dimension');
        end;
    else                
        if verbose, disp('MLD model: E5 is not defined -- assume zeros'); end;
        S{i}.E5 = zeros(ne, 1);
    end
    
end;



% * check, if all the weights have the above defined dimensions and
% * set missing weights to zero
% --------------------------------------------------------------------

for i = 1:length(W)
    
    if isfield(W{i}, 'Qu')
        if ~all( size(W{i}.Qu) == [nu nu] )
            error('invalid weight matrix: Qu has wrong dimension');
        end;
    else        
        if verbose, disp('weight matrices: Qu is not defined -- assume zeros'); end;
        W{i}.Qu = zeros(nu);
    end 
    
    if isfield(W{i}, 'Qd')
        if ~all( size(W{i}.Qd) == [nd nd] )
            error('invalid weight matrix: Qd has wrong dimension');
        end;
    else        
        if verbose, disp('weight matrices: Qd is not defined -- assume zeros'); end;
        W{i}.Qd = zeros(nd);
    end 

    if isfield(W{i}, 'Qz')
        if ~all( size(W{i}.Qz) == [nz nz] )
            error('invalid weight matrix: Qz has wrong dimension');
        end;
    else        
        if verbose, disp('weight matrices: Qz is not defined -- assume zeros'); end;
        W{i}.Qz = zeros(nz);
    end 
    
    if isfield(W{i}, 'Qx')
        if ~all( size(W{i}.Qx) == [nx nx] )
            error('invalid weight matrix: Qx has wrong dimension');
        end;
    else        
        if verbose, disp('weight matrices: Qx is not defined -- assume zeros'); end;
        W{i}.Qx = zeros(nx);
    end 
    
    if isfield(W{i}, 'Qy')
        if ~all( size(W{i}.Qy) == [ny ny] )
            error('invalid weight matrix: Qy has wrong dimension');
        end;
    else        
        if verbose, disp('weight matrices: Qy is not defined -- assume zeros'); end;
        W{i}.Qy = zeros(ny);
    end     
    
end;