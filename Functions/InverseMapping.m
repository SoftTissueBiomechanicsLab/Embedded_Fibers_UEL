function [xi,zeta,eta] = InverseMapping(r, rk, Norder)
options = optimset('Display','off');
if Norder==1
    fun = @(x)paramfun(x, r, rk); 
    x0 = [0,0,0];
    soln = fsolve(fun,x0,options);
    xi = soln(1);
    zeta = soln(2);
    eta = soln(3);
end

function F = paramfun(x, r, rk)
    a = x(1);
    b = x(2);
    c = x(3);
    
    % Shape functions
    N(1) = 1/8*(1-a)*(1-b)*(1-c);
    N(2) = 1/8*(1+a)*(1-b)*(1-c);
    N(3) = 1/8*(1+a)*(1+b)*(1-c);
    N(4) = 1/8*(1-a)*(1+b)*(1-c);
    N(5) = 1/8*(1-a)*(1-b)*(1+c);
    N(6) = 1/8*(1+a)*(1-b)*(1+c);
    N(7) = 1/8*(1+a)*(1+b)*(1+c);
    N(8) = 1/8*(1-a)*(1+b)*(1+c);
    
    F(1) = r(1) - N*rk(:,1) ;
    F(2) = r(2) - N*rk(:,2) ;
    F(3) = r(3) - N*rk(:,3) ;
end

end



