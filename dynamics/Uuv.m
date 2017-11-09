classdef Uuv
% Uuv.m     e.anderlini@ucl.ac.uk     08/11/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This class computes the dynamics of an UUV.
% N.B.: Have a thorough look at relative and absolute velocity concepts -
% they may be inverted and this needs to be rectified.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% Accessible properties:
    properties
        % Physical properties:
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
        % Propulsion object:
        propulsion; % propulsion object
        % Simulation parameter:
        dt;         % time step
        % Output values:
        f_c;        % Coriolis and centripetal force vector
        f_d;        % damping force vector
        f_h;        % hydrostatic restoring force vector
        nu_r;       % relative velocity vector
        tau;        % thrust force vector
    end
    
    %% Protected properties:
    properties (Access = protected)
        J;          % transformation matrix
        n;          % vector of propellers' revolutions
        nu_c;       % current velocity
%         x;          % state vector
%         x_dot;      % state derivative vector
    end
    
    %% Accessible methods:
    methods 
        %% Initialization function:
        function obj = Uuv(uuv,dt)
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
            obj.dt = dt;
            %Initialize propulsion object:
            obj.propulsion = Rov_propulsion(uuv);
            % Compute the missing properties:
            obj.mass = obj.M_RB(1,1);
            obj.I_b = obj.M_RB(4:6,4:6);
            obj.buoyancy = obj.g * obj.density * obj.volume;
            obj.weight = obj.buoyancy;  %obj.g * obj.mass;
            obj.S_r_g = obj.skew(obj.G);
            obj.Mtot = obj.M_RB+obj.M_A;
        end
        
        %% Update the system dynamics:
        function [x,obj] = update_dynamics(obj,x0,n,nu_c)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Input:
            % x0:   state vector at the previous step
            % n:    propeller revolutions (rps)
            % nu_c: current velocity (m/s)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Store the current velocity & propeller revolutions in memory:
            obj.nu_c = nu_c;
            obj.n = n;
            % Compute the state vector for the next time step:
            x = obj.ode4(x0);
            % Get the force vectors for plotting:
            obj = transformation_matrix(obj,x(4),x(5),x(6));
            obj = obj.relative_velocity(x(7:12));
            obj = obj.damping_force();
            obj = obj.hyrostatic_force(x(4),x(5));
            obj = obj.coriolis_force(x(7:12));
            obj.tau = obj.propulsion.get_thrust(obj.n,obj.nu_r);
        end
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
            R = [cos(psi)*cos(theta),...
                cos(psi)*sin(phi)*sin(theta)-cos(phi)*sin(psi),...
                sin(phi)*sin(psi)+cos(phi)*cos(psi)*sin(theta);...
                cos(theta)*sin(psi),...
                cos(phi)*cos(psi)+sin(phi)*sin(psi)*sin(theta),...
                cos(phi)*sin(psi)*sin(theta)-cos(psi)*sin(phi);
                -sin(theta),cos(theta)*sin(phi),cos(theta)*cos(phi)];
            % Rotational transformation matrix:
            T = [1, sin(phi)*tan(theta), cos(phi)*tan(theta);...
                0, cos(phi), -sin(phi);...
                0, sin(phi)/cos(theta), cos(phi)/cos(theta)];
            % 6 degrees of freedom transformation matrix:
            obj.J = [R,zeros(3,3);zeros(3,3),T];
        end
        
        %% Compute the Coriolis & centripetal matrix for the rigid body:
        function obj = coriolis_matrix_rigid_body(obj,nu)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Input:
            % nu: six-dimensional velocity vector
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            mSv1 = obj.mass*obj.skew(nu(1:3));
            mSv2 = obj.mass*obj.skew(nu(4:6));
            obj.C_RB = [zeros(3,3),-mSv1-mSv2*obj.S_r_g;...
                -mSv1+obj.S_r_g*mSv2,-obj.skew(obj.I_b*nu(4:6))];
        end
        
        %% Compute the Coriolis & centripetal matrix for the added mass:
        function obj = coriolis_matrix_added_mass(obj)
            Av = obj.M_A*obj.nu_r;
            SAv1 = obj.skew(Av(1:3));
            obj.C_A = [zeros(3,3),-SAv1;...
                -SAv1,-obj.skew(Av(4:6))];
        end
        
        %% Compute the relative velocity:
        function obj = relative_velocity(obj,nu)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Input:
            % nu:   UUV velocity vector in body-fixed reference frame
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Get the current velocity in the body fixed reference frame:
            nu_c_b = obj.J\obj.nu_c;  
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
        function obj = coriolis_force(obj,nu)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Input:
            % nu: six-dimensional velocity vector
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            obj = obj.coriolis_matrix_rigid_body(nu);
            obj = obj.coriolis_matrix_added_mass();
            obj.f_c = obj.C_RB*nu + obj.C_A*obj.nu_r;
        end
        
        %% Compute the state derivative:
        function x_dot = derivative(obj,x)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Input:
            % x0: initial state vector
            % n:  vector of propeller revolutions
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Update the relative velocity vector:
            obj = transformation_matrix(obj,x(4),x(5),x(6));
            obj = obj.relative_velocity(x(7:12));
            % Update the vector forces:
            obj = obj.damping_force();
            obj = obj.hyrostatic_force(x(4),x(5));
            obj = obj.coriolis_force(x(7:12));
            % Compute the thrust vector:
            obj.tau = obj.propulsion.get_thrust(obj.n,obj.nu_r);
            % Compute the state vector derivative
            x_dot = [obj.J*x(7:12);...
                obj.Mtot\(-obj.f_h-obj.f_d-obj.f_c+obj.tau)];
        end
        
        %% Integrate the derivative with a 4th-order Runge-Kutta method:
        function x = ode4(obj,x0)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Input:
            % x0: initial state vector
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            dxdt1 = obj.derivative(x0);
            x1 = x0 + dxdt1*obj.dt/2;
            dxdt2 = obj.derivative(x1);
            x2 = x0 + dxdt2*obj.dt/2;
            dxdt3 = obj.derivative(x2);
            x3 = x0 + dxdt3*obj.dt;
            dxdt4 = obj.derivative(x3);
            
            % Estimate the state vector at the next time step:
            x = x0 + (dxdt1+2*(dxdt2+dxdt3)+dxdt4)*obj.dt/6;
        end
    end    
end