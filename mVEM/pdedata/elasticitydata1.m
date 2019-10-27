function pde = elasticitydata1(para)

% ------ Lame constants ------
if nargin == 0
    lambda = 1;   mu = 1;
else
    lambda = para.lambda; mu = para.mu;
end

pde = struct('lambda',lambda,'mu',mu, 'f', @f, 'uexact',@uexact,'g_D',@g_D);

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
end