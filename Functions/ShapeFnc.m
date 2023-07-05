function [N,dN, dN2] = ShapeFnc(xi, Norder)

if Norder==1
             
    N(:,1) = (1-xi)/2;
    N(:,2) = (1+xi)/2;
    
    dN(:,1) = -1/2;
    dN(:,2) = +1/2;
   
elseif Norder==2
    N(:,1) = -(1-xi).*xi/2;
    N(:,2) = (1-xi).*(1+xi);
    N(:,3) = (1+xi).*xi/2;
    
    dN(:,1) = xi-1/2;
    dN(:,2) = -2*xi;
    dN(:,3) = xi+1/2;
    
    dN2(:,1) = 1*ones(size(xi));
    dN2(:,2) = -2*ones(size(xi));
    dN2(:,3) = 1*ones(size(xi));
    
elseif Norder==3
    % Is defined in the [0,1] bounds
    q = (1+xi)/2;
    L = 2;
    N(1) = 1-3*q^2 + 2*q^3;
    N(2) = q^2*(3 - 2*q);
    N(3) = L*q*(1 - 2*q + q^2);
    N(4) = L*q^2*(q-1);
end