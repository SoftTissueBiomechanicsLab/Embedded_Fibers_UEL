function varargout = CurvePoint(varargin)
% Returns the total arc length of the curve or the global vector [r] given
% the BeamData or/and the arc length parameter [s]

% Export data
BeamAxesLoc = varargin{1};
s =[];
if length(varargin)==2
    s = varargin{2};
elseif length(varargin)==3
    t = varargin{3};
end



% Process data
r0 = BeamAxesLoc{1}';
rf = BeamAxesLoc{2}';
axis_length = norm(rf-r0);
curve_type = BeamAxesLoc{8};

% Define local coordinate vectors
ez = (rf-r0)/axis_length;
ex_temp = cross([0,1,0]', ez);
if norm(ex_temp) > 1e-12
    ex = ex_temp/norm(ex_temp);
    ey = cross(ez,ex)/norm(cross(ez,ex));
else
    ex = [0,0,1]';
    ey = cross(ez,ex)/norm(cross(ez,ex));
end


% Consider each case
if strcmp(curve_type, 'straight')
    
    if length(varargin)==1
        % Total arc length
        stot = axis_length;
    else
        if length(varargin)==2
            % Curve parameter
            t = s;
        end

        % Local coordinates
        x = 0;
        y = 0;
        z = t;

        % Jacobian
        J = 1;
    end

elseif strcmp(curve_type, 'helix')

    %Export Data
    Ax = BeamAxesLoc{3};  
    nloops = BeamAxesLoc{7};

%     if Ax~=Ay || wx~=wy
%         error('Invalid data for helix')
%     end
    
    % Calculate c parameter
    c = axis_length/(2*nloops*pi);
    
    if length(varargin)==1
        % Total arc length
        stot = sqrt(Ax^2+c^2)*(2*nloops*pi);
    else
        if length(varargin)==2
            % Convert to curve parameters
            t = s/(sqrt(Ax^2+c^2));
        end
                
        % Local coordinates
        x = Ax*cos(t);
        y = Ax*sin(t);
        z = c*t;

        % Jacobian
        J = sqrt(Ax^2+c^2);
    end   

elseif strcmp(curve_type, 'sin_und')
    
    Ax = BeamAxesLoc{3};  
    Ay = BeamAxesLoc{4};
    k1 = BeamAxesLoc{5};
    k2 = BeamAxesLoc{6};
    
    wx = k1*pi/axis_length;
    wy = k2*pi/axis_length;
    
       
    
    if length(varargin)==1

        % Numerical Integration
        fun = @(z) ( (Ax*wx*cos(wx*z)).^2 + (Ay*wy*cos(wy*z)).^2 +1).^(1/2);
        stot = integral(fun,0,axis_length);

    else
        
        if length(varargin)==2
            % Numerical Integration to find curve parameter
            int_term = @(z_star) integral(@(z)((Ax*wx*cos(wx*z)).^2 + (Ay*wy*cos(wy*z)).^2 +1).^(1/2), 0, z_star);
            t = fzero(@(x) int_term(x) - s, 0);
        end
                
        % Local coordinates
        x = Ax*sin(wx*t);
        y = Ay*sin(wy*t);
        z = t;
        
        dxdt = Ax*wx*cos(wx*t);
        dydt = Ay*wy*cos(wy*t);
        dzdt = 1;
        drdt = [dxdt,dydt,dzdt]';
        
        dxdt2 = -Ax*wx^2*sin(wx*t);
        dydt2 = -Ay*wy^2*sin(wy*t);
        dzdt2 = 0;
        drdt2 = [dxdt2,dydt2,dzdt2];
        
        kappa = norm(cross(drdt,drdt2))/norm(drdt)^3;
        
        tangent_vec = [dxdt,dydt,dzdt]';
        tangent_vec = tangent_vec/norm(tangent_vec);
        
        nvec_x = -Ax*wx*(-Ay^2*wy^3*sin(wy*t)*cos(wy*t)*cos(wx*t)... 
                         +Ay^2*wy^2*wx*(cos(wy*t))^2*sin(wx*t)...
                         +wx*sin(wx*t))...
                         /((Ax*wx*cos(wx*t))^2+(Ay*wy*cos(wy*t))^2+1)^(3/2);
        
        nvec_y = -Ay*wy*(-Ax^2*wx^3*sin(wx*t)*cos(wx*t)*cos(wy*t)... 
                         +Ax^2*wx^2*wy*(cos(wx*t))^2*sin(wy*t)...
                         +wy*sin(wy*t))...
                        /((Ax*wx*cos(wx*t))^2+(Ay*wy*cos(wy*t))^2+1)^(3/2);
                    
        nvec_z = (Ax^2*wx^3*sin(wx*t)*cos(wx*t)+Ay^2*wy^3*sin(wy*t)*cos(wy*t))...
                 /((Ax*wx*cos(wx*t))^2+(Ay*wy*cos(wy*t))^2+1)^(3/2);
        
        normal_vec = [nvec_x, nvec_y, nvec_z]';
        normal_vec = normal_vec/norm(normal_vec);
        
        RealGeom.t = tangent_vec;
        RealGeom.n = normal_vec;
        RealGeom.kappa = kappa; % curvature
        
        % Jacobian
        J = ((Ax*wx*cos(wx*t)).^2 + (Ay*wy*cos(wy*t)).^2 +1).^(1/2);
    end

end

if length(varargin)==1
    varargout{1} = stot;
else

    % Assemble final vector
    r = r0 + x*ex + y*ey + z*ez;
    
    varargout{1} = r;
    varargout{2} = t;
    varargout{3} = J;
    
    % Provide tangent and normal vectors
    if nargout==4 
    
    end
    

end

