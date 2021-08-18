function [S,BUCKETS] = table(flag,H,Data,BUCKETS)
%
% S = table(H,Data,flag)
%
% Inputs:
%   H    = Hash code
%   Data = data to store
%   BUCKETS = hash codes
%
% Action taken based on flag:
%   'init' - intialize hash table
%   'get'  - peek a random element
%   'add'  - add the data to the table
%   'num'  - return number of elements in table
%

global DATA     % data blocks
global FREEDATA % indicies of free data blocks
%global BUCKETS  % hash codes
global COLLISIONS

S = [];

if(strcmp(flag,'add'))
	S = 0;
	% Test if this element is already in the list
	I = BUCKETS(H);
	% Test for data collision
	while(I ~= 0)
		D = DATA{I};
		% Item's already here - so remove it and return
		if(D.User.lEr == Data.lEr)
      if(sum(D.User.Er - Data.Er) == 0)
        S = 1;
        % Remove element from list
        FREEDATA = [FREEDATA I];
        if(D.Next ~= 0)
          DATA{D.Next}.Prev = D.Prev;
        end;
        if(D.Prev ~= 0)
          DATA{D.Prev}.Next = D.Next;
        end;
        if(D.Prev == 0)
          BUCKETS(H) = D.Next;
        end;
        return;
      end;
    end;
    I = D.Next;
  end;
  % Item wasn't in the list - so add it to the top
	I = GetIndex;
	DATA{I}.Prev = 0;
	DATA{I}.Next = BUCKETS(H);
	if(BUCKETS(H) ~= 0)
		COLLISIONS = COLLISIONS + 1;
		DATA{BUCKETS(H)}.Prev = I;
	end;
	DATA{I}.H    = H;
	DATA{I}.User = Data;
	BUCKETS(H)   = I;
	return;
end;

if(strcmp(flag,'init'))
    DATA     = [];
    FREEDATA = [];
    BUCKETS  = [];

    N  = 1000;
    for i=[1:N]
        % H = hash code
        % Next = pointer to next item in list
        % Prev = pointer to previous item in list
        DATA{i}   = struct('H',0,'Next',0,'Prev',0,'User',[]);
    end;
    FREEDATA  = [1:N];

    COLLISIONS = 0;
    BUCKETS   = spalloc(10^9,1,N);
%    fprintf('Storage initialized\n');
    return;
end;

if(strcmp(flag,'num'))
    S = nnz(BUCKETS);
    return;
end;

if(strcmp(flag,'get'))
    s = nonzeros(BUCKETS);
    if(~isempty(s))
%        S = DATA{s(ceil(rand*length(s)))}.User;
        S = DATA{s(1)}.User;
    else
        S = [];
    end;
    return;
end;

if(strcmp(flag,'dump'))
	s = nonzeros(BUCKETS);
	for i=[1:length(s)]
		fprintf('%i ',DATA{s(i)}.User.Er);
		fprintf('\n');
	end;
	return;
end;

error('Invalid option in table');

% helper to get a new data index - increases data store if necessary
function I = GetIndex
global FREEDATA
global DATA

% If it's empty we need more space
if(isempty(FREEDATA)) 
    N = length(DATA);
    for i=[N+1:N+500]
        DATA{i}   = struct('H',0,'Next',0,'Prev',0,'User',[]);
    end;
    FREEDATA = [N+1:N+500];
end;
I = FREEDATA(1);
FREEDATA = FREEDATA(2:end);
