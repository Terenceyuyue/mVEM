function Pb = Polybnd_Rectangle(varargin)

if nargin==1
    BdBox = varargin{1};
    if BdBox==4     % [x1,x2,y1,y2]
        x1 = BdBox(1); x2 = BdBox(2);
        y1 = BdBox(3); y2 = BdBox(4);
    end
    if BbBox==2   % [x1,x2], x1=y1, x2=y2
        x1 = BdBox(1);   x2 = BdBox(2);
        y1 = x1;         y2 = x2;
    end
end

if nargin==2    % (x1,y1), x1=y1, x2=y2
    x1 = varargin{1}; x2 = varargin{2};
    y1 = x1;          y2 = x2;
end

if nargin==4   % (x1,x2,y1,y2)
    x1 = varargin{1}; x2 = varargin{2};
    y1 = varargin{3}; y2 = varargin{4};
end

Pb = [x1 y1; x2 y1; x2 y2; x1 y2];
