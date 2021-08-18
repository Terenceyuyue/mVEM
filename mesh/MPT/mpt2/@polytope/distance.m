function [d, x1, x2] = distance(P1,P2,Options)
%DISTANCE Distance between two sets
%
% [d, x1, x2] = distance(P1,P2,Options)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Computes the distance between two sets
%
% d =  min    norm(x1-x2)
%     x1,x2
%      s.t.   x1 in P1
%             x2 in P2
%
% USAGE:
%   [d, x1, x2] = distance(P1,P2)
%   [d, x1, x2] = distance(P1,P2,Options)
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% P1             - point, polytope or polyarray
% P2             - point, polytope or polyarray
% Options.dist_norm - norm for computing the distance,  1, 2, inf
% Options.abs_tol   - absolute tolerance
% Options.lpsolver  - LP solver,   see mpt_solveLP for details
% Options.qpsolver  - QP solver,   see mpt_solveQP for details
%
% Note: At least one input set must be polytope or polyarray!
%       If Options is missing or some of the fields are not defined, the default
%       values from mptOptions will be used.
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% d         - distance between two sets
% x1, x2    - points in the polytopes for which the distance d is achieved
%
% see also POLYTOPE
%

% $Revision: 1.1 $ $Date: 2005/02/23 12:17:39 $
%
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

global mptOptions;

if ~isstruct(mptOptions),
    mpt_error;
end

if nargin<3
    Options=[];
end

if ~isfield(Options,'dist_norm')
    Options.dist_norm=2;    % Euclidean norm
end


if ~isfield(Options,'abs_tol')
    Options.abs_tol=mptOptions.abs_tol;    % absolute tolerance
end

if(~isfield(Options,'lpsolver'))
    Options.lpsolver=mptOptions.lpsolver;
end

if(~isfield(Options,'qpsolver'))
    Options.qpsolver=mptOptions.qpsolver;
end

if nargin<2
    error('DISTANCE: At least two arguments required!');
end

% decide the dimension of the input arguments
% deduce the class of inputs and put the argument that belongs to a higher class in the first place
P1_class=0; % 0 = point; 1 = polytope; 2 = polyarray;
P2_class=0; % 0 = point; 1 = polytope; 2 = polyarray;
if isa(P1,'polytope')
    P1_class=1;
    lenP1 = length(P1.Array);
    if lenP1>0
        P1_class=2;
    end
    n1=dimension(P1);
else
    lenP1 = 0;
    n1=size(P1,1);
end
if isa(P2,'polytope')
    P2_class=1;
    lenP2 = length(P2.Array);
    if lenP2>0
        P2_class=2;
    end
    n2=dimension(P2);
else
    lenP2 = 0;
    n2=size(P2,1);
end

reverse_pos = 0; % do we need to reverse the position of the arguments?
if P1_class < P2_class
    reverse_pos = 1;
    temp=P2_class;
    P2_class=P1_class;
    P1_class=temp;
    temp=P2;
    P2=P1;
    P1=temp;
    temp=lenP2;
    lenP2=lenP1;
    lenP1=temp;
    temp=n2;
    n2=n1;
    n1=temp;
end

if P1_class == 0 & P2_class == 0
    error('DISTANCE: At least one argument MUST be polytope object!');
end

if n1~=n2,
    error('DISTANCE: Arguments must have the same dimension!');
end

n=n1;

d = inf;


if P2_class==0 % second argument is point

    x2 = P2(:);

    if P1_class==1 % first argument is a polytope

        if Options.dist_norm==2
            H = 2*eye(n);
            f = zeros(n,1);
            A=P1.H;
            B=P1.K-P1.H*x2;
            [xopt,lambda,how,exitflag,objqp]=mpt_solveQP(H,f,A,B,[],[],[],Options.qpsolver);
            x1=xopt+x2;
            d=norm(xopt,Options.dist_norm);
        elseif Options.dist_norm==inf
            f = [zeros(n,1); 1];
            A=[P1.H zeros(size(P1.H,1),1); eye(n) -ones(n,1); -eye(n) -ones(n,1)];
            B=[P1.K; x2; -x2];
            [xopt,fval,lambda,exitflag,how]=mpt_solveLP(f,A,B,[],[],[],Options.lpsolver);
            x1=xopt(1:n);
            d=norm(x1-x2,Options.dist_norm);
        elseif Options.dist_norm==1
            f = [zeros(n,1); ones(n,1)];
            A=[P1.H zeros(size(P1.H,1),n); eye(n) -eye(n); -eye(n) -eye(n)];
            B=[P1.K; x2; -x2];
            [xopt,fval,lambda,exitflag,how]=mpt_solveLP(f,A,B,[],[],[],Options.lpsolver);
            x1=xopt(1:n);
            d=norm(x1-x2,Options.dist_norm);
        end
        if strcmp(how,'infeasible')
            error('DISTANCE: Infeasible problem!!!');
        elseif strcmp(how,'unbounded')
            error('DISTANCE: Unbounded problem!!!');
        end
        
    elseif P1_class==2 % first argument is a polyarray
        PA = P1.Array;
        for ii=1:lenP1
            if Options.dist_norm==2
                H = 2*eye(n);
                f = zeros(n,1);
                A=PA{ii}.H;
                B=PA{ii}.K-PA{ii}.H*x2;
                [xopt,lambda,how,exitflag,objqp]=mpt_solveQP(H,f,A,B,[],[],[],Options.qpsolver);
                aux_x1=xopt+x2;
                aux_d=norm(xopt,Options.dist_norm);
            elseif Options.dist_norm==inf
                f = [zeros(n,1); 1];
                A=[PA{ii}.H zeros(size(PA{ii}.H,1),1); eye(n) -ones(n,1); -eye(n) -ones(n,1)];
                B=[PA{ii}.K; x2; -x2];
                [xopt,fval,lambda,exitflag,how]=mpt_solveLP(f,A,B,[],[],[],Options.lpsolver);
                aux_x1=xopt(1:n);
                aux_d=norm(aux_x1-x2,Options.dist_norm);
            elseif Options.dist_norm==1
                f = [zeros(n,1); ones(n,1)];
                A=[PA{ii}.H zeros(size(PA{ii}.H,1),n); eye(n) -eye(n); -eye(n) -eye(n)];
                B=[PA{ii}.K; x2; -x2];
                [xopt,fval,lambda,exitflag,how]=mpt_solveLP(f,A,B,[],[],[],Options.lpsolver);
                aux_x1=xopt(1:n);
                aux_d=norm(aux_x1-x2,Options.dist_norm);
            end
            if strcmp(how,'infeasible')
                error('DISTANCE: Infeasible problem!!!');
            elseif strcmp(how,'unbounded')
                error('DISTANCE: Unbounded problem!!!');
            else
                if aux_d < d
                    x1=aux_x1;
                    d=aux_d;
                end
            end
            if d==0
                break;
            end
        end
    end
    
elseif P2_class == 1 % second argument is a polytope


    if P1_class==1 % first argument is a polytope

        if Options.dist_norm==2
            H=2*[eye(n) -eye(n);-eye(n) eye(n)];
            f=zeros(2*n,1);
            A=[P1.H zeros(size(P1.H,1),n); zeros(size(P2.H,1),n) P2.H];
            B=[P1.K; P2.K];
            [xopt,lambda,how,exitflag,objqp]=mpt_solveQP(H,f,A,B,[],[],[],Options.qpsolver);
        elseif Options.dist_norm==inf
            f = [zeros(n,1); zeros(n,1); 1];
            A = [P1.H zeros(size(P1.H,1),n) zeros(size(P1.H,1),1); 
                zeros(size(P2.H,1),n) P2.H zeros(size(P2.H,1),1); 
                eye(n) -eye(n) -ones(n,1); 
                -eye(n) eye(n) -ones(n,1)];
            B=[P1.K; P2.K; zeros(n,1); zeros(n,1)];
            [xopt,fval,lambda,exitflag,how]=mpt_solveLP(f,A,B,[],[],[],Options.lpsolver);
        elseif Options.dist_norm==1
            f = [zeros(n,1); zeros(n,1); ones(n,1)];
            A = [P1.H zeros(size(P1.H,1),n) zeros(size(P1.H,1),n); 
                zeros(size(P2.H,1),n) P2.H zeros(size(P2.H,1),n); 
                eye(n) -eye(n) -eye(n); 
                -eye(n) eye(n) -eye(n)];
            B=[P1.K; P2.K; zeros(n,1); zeros(n,1)];
            [xopt,fval,lambda,exitflag,how]=mpt_solveLP(f,A,B,[],[],[],Options.lpsolver);
        end
        if strcmp(how,'infeasible')
            error('DISTANCE: Infeasible problem!!!');
        elseif strcmp(how,'unbounded')
            error('DISTANCE: Unbounded problem!!!');
        else
            x1=xopt(1:n);
            x2=xopt(n+1:2*n);
            d=norm(x1-x2,Options.dist_norm);
        end
        
    elseif P1_class==2 % first argument is a polyarray
        
        PA = P1.Array;
        for ii=1:lenP1
            if Options.dist_norm==2
                H=2*[eye(n) -eye(n);-eye(n) eye(n)];
                f=zeros(2*n,1);
                A=[PA{ii}.H zeros(size(PA{ii}.H,1),n); zeros(size(P2.H,1),n) P2.H];
                B=[PA{ii}.K; P2.K];
                [xopt,lambda,how,exitflag,objqp]=mpt_solveQP(H,f,A,B,[],[],[],Options.qpsolver);
            elseif Options.dist_norm==inf
                f = [zeros(n,1); zeros(n,1); 1];
                A = [PA{ii}.H zeros(size(PA{ii}.H,1),n) zeros(size(PA{ii}.H,1),1); 
                    zeros(size(P2.H,1),n) P2.H zeros(size(P2.H,1),1); 
                    eye(n) -eye(n) -ones(n,1); 
                    -eye(n) eye(n) -ones(n,1)];
                B=[PA{ii}.K; P2.K; zeros(n,1); zeros(n,1)];
                [xopt,fval,lambda,exitflag,how]=mpt_solveLP(f,A,B,[],[],[],Options.lpsolver);
            elseif Options.dist_norm==1
                f = [zeros(n,1); zeros(n,1); ones(n,1)];
                A = [PA{ii}.H zeros(size(PA{ii}.H,1),n) zeros(size(PA{ii}.H,1),n); 
                    zeros(size(P2.H,1),n) P2.H zeros(size(P2.H,1),n); 
                    eye(n) -eye(n) -eye(n); 
                    -eye(n) eye(n) -eye(n)];
                B=[PA{ii}.K; P2.K; zeros(n,1); zeros(n,1)];
                [xopt,fval,lambda,exitflag,how]=mpt_solveLP(f,A,B,[],[],[],Options.lpsolver);
            end
            if strcmp(how,'infeasible')
                error('DISTANCE: Infeasible problem!!!');
            elseif strcmp(how,'unbounded')
                error('DISTANCE: Unbounded problem!!!');
            else
                aux_x1 = xopt(1:n);
                aux_x2 = xopt(n+1:2*n);
                aux_d = norm(aux_x1-aux_x2,Options.dist_norm);
                if aux_d < d
                    x1=aux_x1;
                    x2=aux_x2;
                    d=aux_d;
                end
            end
            if d==0
                break;
            end
        end
        
    end
    
elseif P2_class==2 % second argument is polyarray
    
    if P1_class==2 % first argument is polyarray
        
        PA1 = P1.Array;
        PA2 = P2.Array;
        for ii=1:lenP1
            for jj=1:lenP2
                if Options.dist_norm==2
                    H=2*[eye(n) -eye(n);-eye(n) eye(n)];
                    f=zeros(2*n,1);
                    A=[PA1{ii}.H zeros(size(PA1{ii}.H,1),n); zeros(size(PA2{jj}.H,1),n) PA2{jj}.H];
                    B=[PA1{ii}.K; PA2{jj}.K];
                    [xopt,lambda,how,exitflag,objqp]=mpt_solveQP(H,f,A,B,[],[],[],Options.qpsolver);
                elseif Options.dist_norm==inf
                    f = [zeros(n,1); zeros(n,1); 1];
                    A = [PA1{ii}.H zeros(size(PA1{ii}.H,1),n) zeros(size(PA1{ii}.H,1),1); 
                        zeros(size(PA2{jj}.H,1),n) PA2{jj}.H zeros(size(PA2{jj}.H,1),1); 
                        eye(n) -eye(n) -ones(n,1); 
                        -eye(n) eye(n) -ones(n,1)];
                    B=[PA1{ii}.K; PA2{jj}.K; zeros(n,1); zeros(n,1)];
                    [xopt,fval,lambda,exitflag,how]=mpt_solveLP(f,A,B,[],[],[],Options.lpsolver);
                elseif Options.dist_norm==1
                    f = [zeros(n,1); zeros(n,1); ones(n,1)];
                    A = [PA1{ii}.H zeros(size(PA1{ii}.H,1),n) zeros(size(PA1{ii}.H,1),n); 
                        zeros(size(PA2{jj}.H,1),n) PA2{jj}.H zeros(size(PA2{jj}.H,1),n); 
                        eye(n) -eye(n) -eye(n); 
                        -eye(n) eye(n) -eye(n)];
                    B=[PA1{ii}.K; PA2{jj}.K; zeros(n,1); zeros(n,1)];
                    [xopt,fval,lambda,exitflag,how]=mpt_solveLP(f,A,B,[],[],[],Options.lpsolver);
                end
                if strcmp(how,'infeasible')
                    error('DISTANCE: Infeasible problem!!!');
                elseif strcmp(how,'unbounded')
                    error('DISTANCE: Unbounded problem!!!');
                else
                    aux_x1 = xopt(1:n);
                    aux_x2 = xopt(n+1:2*n);
                    aux_d = norm(aux_x1-aux_x2,Options.dist_norm);
                    if aux_d < d
                        x1=aux_x1;
                        x2=aux_x2;
                        d=aux_d;
                    end
                end
                if d==0
                    break;
                end
            end
            
            if d==0
                break;
            end
            
        end
        
    end
    
end

if reverse_pos == 1
    temp=P2;
    P2=P1;
    P1=temp;
    temp=x2;
    x2=x1;
    x1=temp;
end