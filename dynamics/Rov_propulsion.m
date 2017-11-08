classdef Rov_propulsion
% Rov_propulsion.m     e.anderlini@ucl.ac.uk     08/11/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This class computes the thrust vector of a ROV.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% Accessible properties:
    properties
        d;        % vector of propeller diameters (size = no. propulsors)
        density;  % water density
        k_T_c;    % coefficients used to determine the thrust vector
        n;        % propeller revolutions vector (rps)
        n_prop;   % no. of propulsors
        T;        % thrust allocation matrix
        theta_c;  % coefficients used to determine the thrust loss
    end
    
    %% Protected properties:
    properties (Access = protected)
        f_th      % thrust force vector (size = no. propulsors)
        j_a;      % advance ratio vector (size = no. propulsors)
        k_T;      % thrust coefficient vector (size = no. propulsors)
        n_kT;     % no. kT items
        theta;    % thrust loss coefficient vector (size = no. propulsors)
    end
    
    %% Accessible methods:
    methods 
        %% Initialization function:
        function obj = Rov_propulsion(rov)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Input:
            % rov: structure with input data
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            obj.d       = rov.prop_diam;
            obj.density = rov.density;
            obj.k_T_c   = rov.K_T;
            obj.n_kT    = size(obj.k_T_c,2);
            obj.n_prop  = size(obj.d,1);
            obj.T       = rov.T;
            obj.theta_c = rov.theta;
            % Initialize variables:
            obj.f_th  = zeros(obj.n_prop,1);
            obj.j_a   = zeros(obj.n_prop,1);
            obj.k_T   = zeros(obj.n_prop,1);
            obj.theta = zeros(obj.n_prop,1);
        end
        
        %% Compute the thrust vector:
        function thrust = get_thrust(obj,n,nu_r)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Input:
            % n:    propeller revolutions (rps)
            % nu_r: relative velocity of the ROV in body-fixed coordinates
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            obj.n = n;
            obj = obj.advance_ratio(nu_r);
            obj = obj.thrust_coefficient();
            obj = obj.thrust_loss_coefficient();
            obj = obj.thrust_vector();
            thrust = obj.T*obj.f_th;
        end
    end 
    
    %% Protected methods:
    methods (Access = protected)
        %% Compute the advance ratio:
        function obj = advance_ratio(obj,nu_r)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Input:
            % nu_r: relative velocity of the ROV in body-fixed coordinates
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            v_a = obj.T\nu_r;
            obj.j_a = v_a ./ (obj.n .* obj.d);
        end
        
        %% Compute the thrust coefficient vector:
        function obj = thrust_coefficient(obj)
            % Iterate over all propulsors: - not very efficient
            for i=1:obj.n_prop
                obj.k_T(i) = 0.0;
                % Check the flow direction of the propulsor:
                if obj.n(i)<0
                    j = 2;
                else
                    j = 1;
                end
                % Compute the thrust coefficient for the propulsor:
                for k=1:obj.n_kT
                    obj.k_T(i) = obj.k_T(i) +...
                        obj.k_T_c(j,k)*obj.j_a^(obj.n_kT-k);
                end
            end
        end
        
        %% Compute the thrust loss coefficient vector:
        function obj = thrust_loss_coefficient(obj)
            % Iterate over all propulsors: - not very efficient
            for i=1:obj.n_prop
                % Check the flow direction of the propulsor:
                if obj.n(i)<0
                    j = 2;
                else
                    j = 1;
                end
                obj.theta(i) = obj.theta_c(j,i);
            end
        end
        
        %% Compute the thrust vector:
        function obj = thrust_vector(obj)
            obj.f_th = obj.density .* obj.d.^4 .* obj.k_T .* obj.n...
                .* abs(obj.n) .* obj.theta;
        end
    end
end