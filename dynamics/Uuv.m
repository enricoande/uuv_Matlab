classdef Uuv
% Uuv.m     e.anderlini@ucl.ac.uk     07/11/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This class computes the dynamics of an UUV.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% Accessible properties:
    properties
        buoyancy;   % buoyancy force
        density;    % water density
        g;          % gravitational acceleration
        mass;       % mass of the UUV
        volume;     % displaced volume of the UUV
        weight;     % weight of the UUV
        B;          % centre of buoyancy
        C_A;        % Coriolis matrix corresponding to the added mass
        C_RB;       % Coriolis matrix corresponding to the rigid body
        D_L;        % linear damping matrix
        D_q;        % quadratic damping matrix
        G;          % centre of gravity
        I_b;        % inertia matrix
        M_A;        % added mass matrix
        Mtot;       % total mass matrix
        M_RB;       % rigid body mass matrix
        S_r_g;      % skew symmetric matrix
        f_c;        % Coriolis and centripetal force vector
        f_d;        % damping force vector
        f_h;        % hydrostatic restoring force vector
        J;          % transformation matrix
        nu_r;       % relative velocity vector
    end
    
    
    %% Accessible methods:
    methods 
        %% Initialization function:
        function obj = Uuv(uuv)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Input:
            % uuv: structure with input data
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Read the desired data: - preprocessing step may be changed
            obj.D_L = uuv.D_L;
            obj.D_q = uuv.D_q;
            obj.M_A = uuv.M_A;
            obj.M_RB = uuv.M_RB;
            obj.density = uuv.density;
            obj.g = uuv.g;
            obj.B = uuv.B;
            obj.G = uuv.G;
            obj.volume = uuv.volume;
            
            % Compute the missing properties:
            obj.mass = obj.M_RB(1,1);
            obj.I_b = obj.M_RB(3:6,3:6);
            obj.buoyancy = obj.g * obj.density * obj.volume;
            obj.weight = obj.g * obj.mass;
            obj.S_r_g = obj.skew(obj.G);
            obj.Mtot = obj.M_RB+obj.M_A;
        end
        
        %%
    end
    
    %% Protected methods:
    methods (Access = protected)
        %% Return a 3D skew-symmetric matrix:
        function S = skew(obj,x)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Input:
            % x: three-dimensional vector
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Check if size is OK:
            if length(x)>3
                error('Vector should have maximum length of 3');
            end

            % Build skew symmetric matrix:
            S = [0,-x(3),x(2);x(3),0,-x(1);-x(2),x(1),0];
        end
        
        %% Compute the transformation matrix:
        function obj = transformation_matrix(obj,phi,theta,psi)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Input:
            % phi:   roll angle
            % theta: pitch angle
            % psi:   yaw angle
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Translational transformation matrix:
            R = [cos(psi)*cos(theta), cos(psi)*sin(phi)*sin(theta)...
                -cos(phi)*sin(psi), sin(phi)*sin(psi)+cos(phi)*cos(psi)...
                *sin(theta);...
                cos(theta)*sin(psi), cos(phi)*cos(psi)+sin(phi)...
                *sin(psi)*sin(theta), cos(phi)*sin(psi)*sin(theta)...
                -cos(psi)*sin(phi);
                -sin(theta), cos(theta)*sin(phi), cos(theta)*cos(phi)];
            % Rotational transformation matrix:
            T = [1, sin(phi)*tan(theta), cos(phi)*tan(theta);...
                0, cos(phi), -sin(phi);...
                0, sin(phi)/cos(theta), cos(phi)/cos(theta)];
            % 6 degrees of freedom transformation matrix:
            obj.J = [R,zeros(6,6);zeros(6,6),T];
        end
        
        %% Compute the relative velocity:
        function obj = relative_velocity(obj,nu,nu_c)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Input:
            % nu:   UUV velocity vector in body-fixed reference frame
            % nu_c: current velocity vector in inertial reference frame
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Get the current velocity in the body fixed reference frame:
            nu_c_b = obj.J\nu_c;  
            % N.B.: As an alternative, the orthonormality of J can be
            % exploited.
            % Compute the relative velocity in the body-fixed coordinates:
            obj.nu_r = nu-nu_c_b;
        end
        
        %% Compute the damping force vector:
        function obj = damping_force(obj)
            obj.f_d=obj.D_L*obj.nu_r+obj.D_q*diag(abs(obj.nu_r))*obj.nu_r;
        end
        
        %% Compute the restoring hydrostatic force vector:
        function obj = hyrostatic_force(obj,phi,theta)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Input:
            % phi:   roll angle
            % theta: pitch angle
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            obj.f_h = [(obj.weight-obj.buoyancy)*sin(theta);...
                (obj.buoyancy-obj.weight)*cos(theta)*sin(phi);...
                (obj.buoyancy-obj.weight)*cos(theta)*cos(phi);...
                (obj.B(2)*obj.buoyancy-obj.G(2)*obj.weight)*cos(theta)...
                *cos(phi)+(obj.G(3)*obj.weight-obj.B(3)*obj.buoyancy)*...
                cos(theta)*sin(phi);...
                (obj.G(3)*obj.weight-obj.B(3)*obj.buoyancy)*sin(theta)+...
                (obj.G(1)*obj.weight-obj.B(1)*obj.buoyancy)*cos(theta)*...
                cos(phi);...
                (obj.B(1)*obj.buoyancy-obj.G(1)*obj.weight)*cos(theta)*...
                sin(phi)+(obj.B(2)*obj.buoyancy-obj.G(2)*obj.weight)*...
                sin(theta)];
        end
        
        %% Compute the Coriolis and centripetal force vector:
        
        %% Compute the state derivative:
        function x_dot = derivative(obj,tau)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Input:
            % tau: 6 degrees of freedom thrust vector
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            x_dot = [obj.J*obj.nu_r;...
                obj.Mtot\(-obj.f_h-obj.f_d+tau)];
        end
    end    
end