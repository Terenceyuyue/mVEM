function [V,R]=lrsextreme(H,K,Options)

if nargin<3,
    Options = [];
end
if ~isfield(Options,'lrs')
    % change this to 0 to get cdd input file
    Options.lrs = 1;
end

if Options.lrs,
    lrsexec = 'c:\MatlabFiles\mptdevel\lrs\lrs.exe';
else
    lrsexec = 'c:\MatlabFiles\mptdevel\lrs\lcdd.exe';
end

precision = 1e-12;

%%[H,K] = double(P);
nc = size(H,1);
nx = size(H,2);

% prepace input text file for lrs
lrsin = cell(5+nc,1);
lrsin{1} = 'mptlrs';
lrsin{2} = 'H-representation';
lrsin{3} = 'begin';
if Options.lrs,
    lrsin{4} = sprintf('%d %d rational', nc, nx+1);
else
    lrsin{4} = sprintf('%d %d real', nc, nx+1);
end
if Options.lrs,
    [nH,dH] = rat(-H, precision);   % convert doubles to rational numbers
    [nK,dK] = rat(K, precision);
end
for ii=1:nc,
    hs='';
    for jj=1:nx,
        if Options.lrs,
            hs = [hs sprintf(' %.0f/%.0f', nH(ii,jj), dH(ii,jj))];
        else
            hs = [hs sprintf(' %.17f', -H(ii,jj))];
        end
    end
    if Options.lrs,
        ks = sprintf('%.0f/%.0f', nK(ii), dK(ii));
    else
        ks = sprintf('%.17f', K(ii));
    end
    lrsin{4+ii} = sprintf('%s%s', ks, hs);
end
lrsin{4+nc+1} = 'end';

tmpfile = tempname;
tmpine = [tempname '.ine'];
tmpout = [tempname '.out'];
fid = fopen(tmpine,'w');

if fid < 0,
    error(['error opening file "' tmpine '" for writing!']);
end
for ii=1:5+nc,
    fprintf(fid,'%s\n',lrsin{ii});
end
fclose(fid);

[status, lrsout] = system([lrsexec ' < ' tmpine ' > ' tmpout]);
delete(tmpine);
if status ~= 0,

    fid = fopen('mptlrs.out','r');
    fid = fopen(tmpout,'r');
    if fid < 0,
        error(['error opening file "' tmpout '" for reading!']);
    end
    
    while 1,
        tline = fgetl(fid);
        if ~ischar(tline), break, end
        disp(tline);
    end
    fclose(fid);
    V = []; R = [];
    error(['error executing lrs ! if you do not see any listing above, try "type ' tmpout '"']);
    return
end

fid = fopen(tmpout,'r');
if fid < 0,
    error(['error opening file "' tmpout '" for reading!']);
end

lrsV = [];
if Options.lrs,
    startline = '*****';
    offset=0;
else
    startline = 'begin';
    offset=1;
end
firsttime = 1;
while 1
    tline = fgetl(fid);
    if ~ischar(tline), break, end
    if ~isempty(strfind(tline,startline)),
        while 1
            tline = fgetl(fid);
            if ~Options.lrs,
                if firsttime,
                    tline = fgetl(fid);
                    firsttime = 0;
                end
            end
            if ~(tline(2+offset)=='1' | tline(2+offset)=='0')
                break
            end
            lrsV = [lrsV; str2num(tline)];
        end
    end
end
fclose(fid);
delete(tmpout);
V = [];
R = [];
if ~isempty(lrsV),
    Vpos = find(lrsV(:,1)==1);
    Rpos = find(lrsV(:,1)==0);
    V = lrsV(Vpos,2:end);
    R = lrsV(Rpos,2:end);
end
