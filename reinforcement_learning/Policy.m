classdef Policy
% Policy.m     e.anderlini@ucl.ac.uk     19/10/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This class obtains the policy for a reinforcement learning problem.
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
%
% For greater information on the use of radial basis functions with
% reinforcement learning or dynamic programming see:
% Geramifard, et al. (2013). 'A Tutorial on Linear Function Approximators 
% for Dynamic Programming and Reinforcement Learning'. Foundations and 
% Trends in Machine Learning, Vol. 6, No. 4, pp. 418-419.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% Accessible properties:
    properties
        a;                  % selected action index
        actions;            % actions list
        basis;              % type of basis function
        centres;            % list of centres of the RBFs
        exploration_rate;   % exploration rate
        mu;                 % width of each RBF
        nbasis;             % no. basis functions overall
        states;             % list of discrete states
        weights;            % weights of the linear function approximation
    end

    %% Protected properties:
    properties (Access = protected)
        na;                 % no. actions
        ns;                 % no. discrete states
        nrbf;               % no. radial basis functions per action
        state;              % current state
        Q;                  % Q value for each action-state pair
    end
    
    %% Accessible methods:
    methods %(Access = protected)
        %% Initialization function:
        function obj = Policy(actions,states,epsilon,basis,mu)  
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Input:
            % actions: actions list
            % states:  states list or centres of the RBFs
            % epsilon: initial exploration rate
            % delta:   discount factor
            % basis:   type of basis functions
            % mu:      width of each RBF
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if nargin<4
                basis = 'basis_exact';
            end
            
            obj.actions = actions;
            obj.states = states;
            obj.exploration_rate = epsilon;
            
            obj.na = size(actions,1);
            switch basis
                case 'basis_exact'
                    obj.basis = basis;
                    obj.states = states;
                    obj.ns = size(states,1);
                case 'basis_rbf'
                    obj.basis = basis;
                    obj.centres = states;
                    obj.mu = mu;
                    obj.nrbf = size(states,1);
                    obj.weights = zeros((obj.nrbf+1)*obj.na,1);
                otherwise
                    error(['Only exact and radial basis functions',...
                    'supported']);
            end
        end
        
        %% Setter function for the features' weights:
        function obj = set_weights(obj,weights)
            obj.weights = weights;
        end
        
        %% Update the exploration rate with the given rate:
        function obj = update_epsilon(obj,rate)
            obj.exploration_rate = obj.exploration_rate*rate;
        end
            
        %% Select an action:
        function obj = select_action(obj,state)            
            % Store the current state:
            obj.state = state;
            % Compute the state-action values for the current state:
            obj = obj.qvalues();
            % Get the current action with an epsilon-greedy policy:
            obj.a = obj.eGreedy();
        end
        
        %% Evaluate the features:
        function phi = get_features(obj,state,action)
            % Store the current state:
            obj.state = state;
            % Get the features:
            phi = feval(obj.basis,action);
        end
    end
    
    %% Protected methods:
    methods (Access=protected)
        %% Find the discrete state:
        function s = discretizeState(obj,x)            
            % Copy the row vector entries (continuous states) to all rows:
            x = repmat(x,obj.ns,1);
            % Select the row using the minimum Eucledian distance:
            [~,s] = min(sum((obj.states-x).^2,2).^0.5);
        end
        
        %% Get the Q-value function for current state and action:
        function q = qvalue(obj,action)
            phi = feval(obj.basis,action);
            q = phi' * obj.weights;
        end
        
        %% Get the Q-value functions for the current state:
        function obj = qvalues(obj)            
            % Initialize the Q-values for the current state:
            obj.Q = zeros(obj.na,1);
            % Calculate the state-action values for the current state:
            for a=1:obj.na
                obj.Q(a) = obj.qvalue(a);
            end
        end
        
        %% Get an action with an epsilon-greedy exploration policy:
        function a = eGreedy(obj)
            % Generate a random number:
            r = rand;

            % Select the action that maximises Q(s)
            if (r>obj.exploration_rate)             
                [~,a] = max(obj.Q); % value, action 
            % Choose a random action:
            else                       
                a = randi(obj.na);  % random integer based on a uniform
            end                     % distribution
        end
        
        %% Find the features for the exact basis functions:
        function phi = basis_exact(obj,action)
            %Initialize the features:
            phi = zeros(obj.nbasis,1);

            % Find the current discrete state:
            s = discretizeState(obj.state);

            % Find the starting position of the block:
            base = (action-1) * obj.ns;

            % Set the indicator:
            phi(base+s) = 1;   
        end
        
        %% Find the features for the radial basis functions:
        function phi = basis_rbf(obj, action)
            %Initialize the features:
            phi = zeros(obj.nbasis,1);

            % Find the starting position:
            base = (action-1) * (obj.nbasis/obj.na);
            % This is because the matrix Theta is converted into a line 
            % vector
            
            % Compute the RBFs:
            for i=1:obj.nrbf
                phi(base+i) = exp(-norm(obj.state-obj.centres(i,:))^2/...
                    (2*obj.mu));
            end
            % ... and the constant:
            phi(base+obj.nrbf+1) = 1;
        end
    end
end