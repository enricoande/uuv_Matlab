classdef Lspi < handle
% Lspi.m     e.anderlini@ucl.ac.uk     19/10/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This class implements least-squares policy iteration (LSPI).
% This class relies on the class Policy.
% Note that this class has been inspired by the code written by
% Lagoudakis and Parr (2003). 'Least-squares policy iteration'. The Journal
% of Machine Learning Research, Vol. 4, pp 1107-1149.
%
% Copyright 2000-2002 
%
% Michail G. Lagoudakis (mgl@cs.duke.edu)
% Ronald Parr (parr@cs.duke.edu)
%
% Department of Computer Science
% Box 90129
% Duke University, NC 27708
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% Accessible properties:
    properties
        algorithm;            % LSPI algorithm type
        delta;                % accuracy setting
        discount_factor;      % discount factor
        maxiterations;        % max. no. iterations
        previous_policy;      % previous policy object
        policy;               % current policy object
    end
    
    %% Protected properties:
    properties (Access = protected)
        ns;                   % no. state variables (not states)
    end
    
    %% Accessible methods:
    methods
        %% Initialization function:
        function obj = Lspi(policy,discount_factor,maxiterations,delta,...
            algorithm)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Input:
            % algorithm - This is a number that indicates which evaluation
            %             algorithm should be used (see the paper):
            %
            %             1-lsq       : The regular LSQ (incremental)
            %             2-lsqfast   : A fast version of LSQ (more space)
            %             3-lsqbe     : LSQ with Bellman error minimization 
            %             4-lsqbefast : A fast version of LSQBE (heavier)
            %
            %             LSQ is the evaluation algorithm for regular
            %             LSPI. Use lsqfast in general, unless you have
            %             really big sample sets. LSQBE is provided for
            %             comparison purposes and completeness.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            obj.policy = policy;
            obj.ns = size(obj.policy.states,2);
            obj.discount_factor = discount_factor;
            obj.maxiterations = maxiterations;
            obj.delta = delta;
            if nargin<4
                algorithm = 2;
            end
            % Check if the algorithm index is in the recognized range:
            if algorithm==1 || algorithm==2 || algorithm==3 || algorithm==4
                obj.algorithm = algorithm;
            else
                error('Only options 1-4 supported');
            end
        end
        
        %% LSTD-Q algorithm:
        function obj = lstdq(obj,samples)
            % If no samples, return error:
            if isempty(samples)
                error('Error: Empty sample set');
            end
            
            % Initialize policy iteration:
            iteration = 0;
            distance = inf;
            obj.previous_policy = obj.policy;
            
            % Main LSPI loop:  
            while (iteration < obj.maxiterations) && (distance > obj.delta)
                % Update and print the number of iterations:
                iteration = iteration + 1;
                if (iteration==1)
                    firsttime = true;
                else
                    firsttime = false;
                end    

                % Evaluate the current policy (and implicitly improve):
                switch obj.algorithm
                    case 1
                        obj.policy = obj.policy.set_weights(...
                            obj.lsq(samples));
                    case 2
                        obj.policy = obj.policy.set_weights(...
                            obj.lsqfast(samples,firsttime));
                    case 3
                        obj.policy = obj.policy.set_weights(...
                            obj.lsqbe(samples));
                    case 4
                        obj.policy = obj.policy.set_weights(...
                            obj.lsqbefast(samples,firsttime));
                end

                % Compute the distance between the current and the previous
                % policy:
                l1 = length(obj.policy.weights);
                l2 = length(obj.previous_policy.weights);
                if (l1 == l2)
                    difference = obj.policy.weights - ...
                        obj.previous_policy.weights;
                    LMAXnorm = norm(difference,inf);
                    L2norm = norm(difference);
                else
                    LMAXnorm = abs(norm(obj.policy.weights,inf) - ...
                        norm(obj.previous_policy.weights,inf));
                    L2norm = abs(norm(obj.policy.weights) - ...
                        norm(obj.previous_policy.weights));
                end
                distance = L2norm;
                obj.previous_policy = obj.policy;
            end
        end
    end
    
    %% Protected methods:
    methods (Access = protected)
        %% Regular LSQ algorithm:
        function [w, A, b] = lsq(obj,samples)
            % Initialize variables:
            nsamples = size(samples,1);
            k = obj.policy.nbasis;
            A = zeros(k, k);
            b = zeros(k, 1);

            % Loop through the samples:
            for i=1:nsamples
                % Compute the basis for the current state and action:
                phi = obj.policy.get_features(samples(i,1:obj.ns),...
                    samples(i,obj.ns+1));

                % Compute the policy and features at the next state:
                l = obj.ns+3;
                m = 2*obj.ns+1;
                nextaction = obj.previous_policy.select_action(...
                    samples(i,l:m));
                nextphi = obj.policy.get_features(samples(i,l:m),...
                    nextaction);

                % Update the matrices A and b:
                A = A + phi * (phi - obj.discount_factor * nextphi)';
                b = b + phi * samples(i,obj.ns+2);
            end

            % Solve the system to find the weights:
            rankA = rank(A);
            if rankA==k
                w = A\b;
            else
                w = pinv(A)*b;
            end
        end
        
        %% Fast version of LSQ algorithm:
        function [w, A, b] = lsqfast(obj,samples,firsttime)
            persistent Phihat;
            persistent Rhat;

            % Initialize variables:
            nsamples = size(samples,1);
            k = feval(new_policy.basis);
            PiPhihat = zeros(nsamples,k);

            % Precompute Phihat and Rhat for all subsequent iterations:
            if firsttime
                Phihat = zeros(nsamples,k);
                Rhat = zeros(nsamples,1);

                for i=1:nsamples
                    phi = obj.policy.get_features(samples(i,1:obj.ns),...
                        samples(i,obj.ns+1));
                    Phihat(i,:) = phi';
                    Rhat(i) = samples(i,obj.ns+2);
                end
            end

            % Loop through the samples:
            for i=1:nsamples
                % Compute the policy and features at the next state:
                l = obj.ns+3;
                m = 2*obj.ns+1;
                nextaction = policy_function(obj.previous_policy,...
                    samples(i,l:m));
                nextphi = obj.policy.get_features(samples(i,l:m),...
                    nextaction);
                PiPhihat(i,:) = nextphi';
            end

            % Compute the matrices A and b:
            A = Phihat' * (Phihat - obj.discount_factor * PiPhihat);
            b = Phihat' * Rhat;

            % Solve the system to find w:
            rankA = rank(A);
            if rankA==k
                w = A\b;
            else
                w = pinv(A)*b;
            end
        end
        
        %% LSQ algorithm with Bellman error minimization:
        function [w, A, b] = lsqbe(obj,samples)
            nsamples = size(samples,1);
            k = feval(new_policy.basis);
            A = zeros(k, k);
            b = zeros(k, 1);

            % Loop through the samples:
            for i=1:nsamples
                % Compute the basis for the current state and action:
                phi = obj.policy.get_features(samples(i,1:obj.ns),...
                    samples(i,obj.ns+1));

                % Compute the policy and features at the next state:
                l = obj.ns+3;
                m = 2*obj.ns+1;
                nextaction = obj.previous_policy.select_action(...
                    samples(i,l:m));
                nextphi = obj.policy.get_features(samples(i,l:m),...
                    nextaction);

                % Update the matrices A and b:
                A = A + (phi - obj.discount_factor * nextphi) ...
                  * (phi - obj.discount_factor * nextphi)';
                b = b + (phi - obj.discount_factor * nextphi)*...
                    samples(i,obj.ns+2);
            end

            % Solve the system to find w:
            rankA = rank(A);
            if rankA==k
                w = A\b;
            else
                w = pinv(A)*b;
            end
        end
        
        %% Fast version of LSQ algorithm with Bellman error minimization:
        function [w, A, b] = lsqbefast(obj,samples,firsttime)
            persistent Phihat;
            persistent Rhat;

            % Initialize variables:
            nsamples = size(samples,1);
            k = feval(new_policy.basis);
            PiPhihat = zeros(nsamples,k);

            % Precompute Phihat and Rhat for all subsequent iterations:
            if firsttime
                Phihat = zeros(nsamples,k);
                Rhat = zeros(nsamples,1);
                for i=1:nsamples
                    phi = obj.policy.get_features(samples(i,1:obj.ns),...
                        samples(i,obj.ns+1));
                    Phihat(i,:) = phi';
                    Rhat(i) = samples(i,obj.ns+2);
                end
            end

            % Loop through the samples:
            for i=1:nsamples
                % Compute the policy and the features at the next state: 
                l = obj.ns+3;
                m = 2*obj.ns+1;
                nextaction = obj.previous_policy.select_action(...
                    samples(i,l:m));
                nextphi = obj.policy.get_features(samples(i,l:m),...
                    nextaction);
                PiPhihat(i,:) = nextphi';
            end

            % Compute the matrices A and b:
            A = (Phihat - obj.discount_factor * PiPhihat)' * ...
              (Phihat - obj.discount_factor * PiPhihat);
            b = (Phihat - obj.discount_factor * PiPhihat)' * Rhat;

            % Solve the system to find w:
            rankA = rank(A);
            if rankA==k
                w = A\b;
            else
                w = pinv(A)*b;
            end
        end
    end
end