function pde = elasticitydata(para)

% ------ Lame constants ------
if nargin == 0
    lambda = 1;   mu = 1;
else
    lambda = para.lambda; mu = para.mu;
end

para = struct('lambda',lambda, 'mu',mu);
    function val = f(p)
        val = mu*2*pi^2.*uexact(p);
    end
    function val = uexact(p)
        x = p(:,1); y = p(:,2);
        u1 = cos(pi*x).*cos(pi*y); u2 = sin(pi*x).*sin(pi*y);
        val = [u1, u2];
    end
    function val = g_D(p)
        val = uexact(p);
    end
    function val = Du(p)
        x = p(:,1); y = p(:,2);
        Du1 = [-pi*sin(pi*x).*cos(pi*y), -pi*cos(pi*x).*sin(pi*y)];
        Du2 = [pi*cos(pi*x).*sin(pi*y), pi*sin(pi*x).*cos(pi*y)];
        val = [Du1, Du2];
    end
    function val = g_N(p)
        x = p(:,1); y = p(:,2);
        E11 = -pi*sin(pi*x).*cos(pi*y); E22 = pi*sin(pi*x).*cos(pi*y);
        E12 = 0.5*(-pi*cos(pi*x).*sin(pi*y) + pi*cos(pi*x).*sin(pi*y));
        sig11 = (lambda+2*mu)*E11 + lambda*E22;
        sig22 = (lambda+2*mu)*E22 + lambda*E11;
        sig12 = 2*mu*E12;
        val = [sig11,sig22,sig12];
    end
pde = struct('para',para, 'f', @f, 'uexact',@uexact,'g_D',@g_D,'g_N',@g_N,'Du',@Du);
end