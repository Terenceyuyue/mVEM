function pde = elasticitydata1(para)

% ------ Lame constants ------
if nargin == 0
    lambda = 1;   mu = 1;
else
    lambda = para.lambda; mu = para.mu;
end

para = struct('lambda',lambda,'mu',mu);
    function val = f(p)
        x = p(:,1); y = p(:,2);
        val = 2*mu*pi^3*[-cos(pi*y).*sin(pi*y).*(2*cos(2*pi*x)-1), ...
            cos(pi*x).*sin(pi*x).*(2*cos(2*pi*y)-1)];
    end
    function val = uexact(p)
        x = p(:,1); y = p(:,2);
        val = [pi*cos(pi*y).*sin(pi*x).^2.*sin(pi*y), ...
            -pi*cos(pi*x).*sin(pi*x).*sin(pi*y).^2];
    end
    function val = g_D(p)
        val = uexact(p);
    end
    function val = Dw(p)
        x = p(:,1); y = p(:,2);
        val(:,1) = 2*pi^2*cos(pi*x).*cos(pi*y).*sin(pi*x).*sin(pi*y); % u1x
        val(:,2) = pi^2*cos(2*pi*y).*sin(pi*x).^2; % u1y
        val(:,3) = -pi^2*cos(2*pi*x).*sin(pi*y).^2; % u2x
        val(:,4) = -2*pi^2*cos(pi*x).*cos(pi*y).*sin(pi*x).*sin(pi*y); % u2y
    end

pde = struct('para',para, 'f', @f, 'uexact',@uexact,'g_D',@g_D,'Dw',@Dw);
end