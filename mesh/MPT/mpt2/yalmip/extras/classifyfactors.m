function [constants,general,singles,pairs] = classifyfactors(L,M,R)


% Author Johan Löfberg, Rikard Falkeborn
% $Id: classifyfactors.m,v 1.8 2009-10-27 10:26:20 joloef Exp $
nfac = length(L);
used = zeros(1,nfac);

constants = 0;
general = 0;
singles = {};
pairs = {};
for i = 1:length(L)
    if isa(M{i},'double')
        constants = constants + L{i}*M{i}*R{i};
        used(i)=1;
    end
end

for i = 1:length(L)
    if ~used(i)      
        temp = struct(M{i});
        if ~(is(M{i},'sdpcone') || is(M{i},'lpcone') || isequal(temp.originalbasis,'skew') || isequal(temp.originalbasis,'diagonal'))
            general = general + L{i}*M{i}*R{i};
            used(i) = 1;
        end
    end
end

for i = 1:length(L)
    if ~used(i)
        if issymmetric(M{i}) && (isequal(L{i},R{i}')) 
            singles{end+1}.L = L{i};
            singles{end}.M = M{i};
            singles{end}.R = R{i};
            singles{end}.negated = 0;
            used(i)=1;
        elseif issymmetric(M{i}) && (isequal(L{i},-R{i}')) 
            singles{end+1}.L = L{i};
            singles{end}.M = M{i};
            singles{end}.R = R{i};
            singles{end}.negated = 1;
            used(i)=1;          
        elseif issymmetric(M{i}) && isequal(find(L{i}),find(R{i}'))           
            Rttmp = R{i}';
            tmp = L{i}(L{i}~=0)./(Rttmp(Rttmp~=0));
            
            if length(unique(tmp)) == 1                   
                  singles{end+1}.L = L{i}/sqrt(abs(tmp(1)));
                  singles{end}.M = M{i};
                  singles{end}.R = R{i}*sqrt(abs(tmp(1)));
                  if tmp(1) > 0
                      singles{end}.negated = 0;
                  else
                      singles{end}.negated = 1;
                      
                  end
                  used(i)=1;                
            end
    
        else            
            for j = i+1:length(L)
                if isa(L{i}*M{i}*R{i}-R{j}'*M{j}'*L{j}','double')
                    pairs{end+1}.L = L{i};
                    pairs{end}.M = M{i};
                    pairs{end}.R = R{i};
                    used(j)=1;
                    used(i)=1;
                    break
                end
            end
        end
    end
end


for k = find(~used)
    if isequal(size(M{k}),[1 1])
        general = general + L{k}*M{k}*R{k};
    else
        error('Couldn''t classify all factors');
    end
end


