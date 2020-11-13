function blks(p,Q_in,v_in)

% Start by placing blocks in separate matrices
[v1,dummy1,r1,dummy3]=dmperm(Q_in{1}+eye(length(Q_in{1})));
for blocks = 1:length(r1)-1
    i1 = r1(blocks);
    i2 = r1(blocks+1)-1;   
    Q{blocks} = Q_in{1}(i1:i2,i1:i2);
    v{blocks} = v_in{1}(i1:i2);
    S{blocks} = eye(length(v{blocks}));
end

pdecomp = 0;
F = set([]);
obj = 0;

for blocks = 1:length(Q)
    if length(Q{blocks})> 1 
        if blocks <= 100
        [Si,index] = blkOne(Q{blocks});
        else
            Si = eye(length(Q{blocks}));
        end
    else
        Si = 1;
        index = 1;
    end    
    Qreduced{blocks} = sdpvar(size(Si,2));
    F = F + set(Qreduced{blocks} >= 0);
    S{blocks} = S{blocks}*Si;
    vi = S{blocks}'*v{blocks};
    pdecomp = pdecomp + vi'*Qreduced{blocks}*vi;
    obj = obj + trace(Qreduced{blocks});
end

solvesdp(F+set(coefficients(p-pdecomp,recover(depends(p))) == 0),obj,sdpsettings('dualize',1));

for i = 1:length(Q)
    Q{i} = double(Qreduced{i});    
end

clear ss
for i = 1:length(Q)
    ss(i) = norm(Q{i}) > 1e-4;
end
Q = {Q{ss}};
S = {S{ss}};
v = {v{ss}};

F       


function [S,index] = blkOne(T)

s = 1;
S = [];
clear v
for i = 1:size(T,1)
    for j = i+1:size(T,1)
        v(i,j) = (T(i,:)*T(j,:)'./(norm(T(i,:))*norm(T(j,:))));%var(T(i,:)./T(j,:))/abs(mean(T(i,:)./T(j,:)));       
    end
end
v = abs((abs(v)-1)) < 1e-6 ;

top = 1;
S = [];
remove = [];
for i = 1:size(T,1)    
    if any(v(:,i))
        j = find(v(:,i));
        for k = 1:length(j)
            nn = (norm(T(i,:))/norm(T(j(k),:)));
            s = (T(i,:)*T(j(k),:)'./(norm(T(i,:))*norm(T(j(k),:))));
            if j(k) <= size(S,2) 
                if S(j(k),j(k))
            S(i,j(k)) = s*(round(10*nn)/10);
            S(i,j(k)) = round(S(i,j(k))*10)/10;
            remove = [remove i];
                end
            end
        end
    else
       S(i,i) = 1;
    end
end
        
S = S(:,any(S,1));
index = setdiff(1:size(T,2),remove);

TT = T;
TT(remove,:) = [];
TT(:,remove) = [];
S*TT*S' - T
    

function [Q,v] = prune(Q,v,Qtemp)

for sosfun = 1:length(Q)   
    keep = diag(Qtemp{sosfun})>tol;
    Qtemp(:,find(~keep)) = [];
    Qtemp(find(~keep),:) = [];

    m{sosfun} = m{sosfun}(find(keep));

    Qtemp(abs(Qtemp) < tol) = 0;
    [v1,dummy1,r1,dummy3]=dmperm(Qtemp+eye(length(Qtemp)));
    lengths{sosfun} = [];
    n{sosfun} = {};
    for blocks = 1:length(r1)-1
        i1 = r1(blocks);
        i2 = r1(blocks+1)-1;
        if i2>i1
            n{sosfun}{blocks} = m{sosfun}(v1(i1:i2));
        else
            n{sosfun}{blocks} = m{sosfun}(v1(i1));
        end
        lengths{sosfun} =  [lengths{sosfun}  length(n{sosfun}{blocks})];
    end
    lengths{sosfun} = sort(lengths{sosfun});
end

