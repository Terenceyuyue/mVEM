function pde = elasticitydata(para)

% ------ Lame constants ------
if nargin == 0
    lambda = 1e8;   mu = 1;
else
    lambda = para.lambda; mu = para.mu;
end

pde = struct('lambda',lambda,'mu',mu, 'f', @f, 'uexact',@uexact,'g_D',@g_D);

    function val = f(p)
        val = mu*2*pi^2.*uexact(p);
    end
    function val = uexact(p)
        x = p(:,1); y = p(:,2);
        val = [cos(pi*x).*cos(pi*y), sin(pi*x).*sin(pi*y)];
    end
    function val = g_D(p)
        val = uexact(p);
    end
end