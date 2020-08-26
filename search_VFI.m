
%%%%%%%%%%%%%%%%%%%%% cody by Jing Ren (01/896410) %%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%% The Bellman equation is written in the report %%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%% Discretize the state variables %%%%%%%%%%%%%%%%%%%%%%%

  x = (0 : 0.01 : 0.99)';  % state space for wage
  n = size(x , 1);         % n is the number of states for wage
  
  hhat = 0.2;
  h = (0 : hhat : 1)';     % state space for health
                           % h = 1 does not exist in reality but the limit
  
%%%%%%%%%%%%%%%%% Discretize the control variables %%%%%%%%%%%%%%%%%%%%%%%
  
  lambda = (0 : 0.01 : 1); % state space for search intensity
  m = size(lambda , 2);    % m is the number of states for search intensity
  
%%%%%%%%% Initialize the value function and the policy function %%%%%%%%%%
  
  V = ones(n , 6);          % value function for all health and wage states
  V(: , 6) = zeros(n , 1);  % value when h = 1
                            % As is mentioned above, we are including
                            % V(h = 1) = 0 because it can help us loop
                            % backwards
  
  V_temp = ones(n , 1);     % This vector temporarily stores the value 
                            % function from each loop and forms one
                            % column in the ultimate value function
                            % matrix (corresponding to one health state)
  V_temp_new = ones(n , 1); % the updated value function
  
  Pol = ones(n , 5);        % Initialize the policy function
  Pol_temp = ones(n , 1);   % This vector temporarily stores the policy 
                            % function from each loop

%%%%%%%%%%%%%%%%%%%%%%%%%%%% parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  r = 0.05;                 % discount rate
  delta = 0.2;              % separation rate

  beta0 = 2; 
  beta1 = 5;                % parameters for the Beta distribution

  c0 = 0.5; 
  c1 = 2;                   % parameters for the cost of search function
  
  d1 = 1; 
  dh = (h ./ (1 - h)).^d1;   % the death rate
  alpha = 0.075;             % health status transition probability

%%%%%%%%%%%%%%%%%%%%%%%% iteration criterion %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  MxIt = 1000;               % maximum number of iterations
  eps = 0.001;               % convergence criterion
  
%%%%%%%%%%%%%%%%%%%%%%% Set up the return matrix %%%%%%%%%%%%%%%%%%%%%%%%% 
        
  R = x * ones(1 , m) - ones(n , 1) * (c0 * lambda.^c1); 
  % utility for the current period
  % This part does not need to be updated in the iterations
                                                         
%%%%%%%%%%%%%%%%% Loop over different health states %%%%%%%%%%%%%%%%%%%%%%                                     
                                                         
  for k = 1 : 5
      
%%%%%%%%%%%%%%%%%%%%%%%%%% Start iteration  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  for l = 1 : MxIt
      
%%%%%%%%%%%%%%%%%% possible value lost from separation %%%%%%%%%%%%%%%%%%%

  S = delta * (V_temp(1) * ones(n , m) - V_temp * ones(1 , m)); 
  % Because b = 0, W(h,b) takes the value when x = b = 0, the first
  % element V_temp
  
%%%%%%%%%%%%% possible value change from health status change %%%%%%%%%%%%

  H = alpha * (V(: , 7 - k) - V_temp) * ones(1 , m); 
  % 7 - k is h + hhat, the column behind the current column, which was 
  % filled in in the previous loop
  
%%%%%%%%%%%%%%%%% possible value change from job arrivals %%%%%%%%%%%%%%%%

%%%%%% Use numerical integration to calculate the value of integrals %%%%%

  A = zeros(n , 1); % Initialize the first part of the integral

  for i = 1 : (n - 1)
    X = x(i : n); 
    % the width of the trapezoids, extending from x to 1, divided into 100 parts
    Y = V_temp(i : n) .* betapdf(X , beta0 , beta1); 
    % the heighth of the trapezoids
    % V_temp(i) is the W(h,y) term
    % betapdf is the dF(y) term
    A(i) = trapz(X , Y); 
    % the total area corresponding to each element in the state variable x
    % A(i) denotes the whole integral when the current wage is x(i)
  end

  B = zeros(n , 1); % Initialize the second part of the integral

  for j = 1 : (n - 1)
    X = x(j : n); 
    % the width of the trapezoids, extending from x to 1, divided into 100 parts
    Y = betapdf(X , beta0 , beta1); 
    % the heighth of the trapezoids
    B(j) = V_temp(j) .* trapz(X , Y); 
    % the total area corresponding to each element in the state variable x
    % V_temp(j) is the W(h,x) term, which can be taken outside of the
    % integral sign
    % As before, betapdf is the dF(y) term
    % B(j) denotes the whole integral when the current wage is x(j)
  end
  
  % For some reasons, i = n and j = n cannot be accepted by trapz, perhaps
  % because it needs at least two points to start with, which is why I
  % excluded i = n and j = n from the loops.

  W = (A - B) * lambda;

%%%%%%%%%%%%%%%%%%%%%%%%%% the discount factor %%%%%%%%%%%%%%%%%%%%%%%%%%%

  D = r + dh(6 - k);          % the discount rate
                              % The death rate is included, making the
                              % individual more impatient
  
%%%%%%%%%%%%%%% Maximize the RHS of the Bellman equation %%%%%%%%%%%%%%%%%
  
  RHS = R + S + W + H;        % terms inside the max sign
  
  [RHS , index] = max(RHS');  % Take maximum by choosing different lambda's
  
  V_temp_new = RHS' / D;      % Transpose back the RHS
  
  Pol_temp = index';          % Find the policy function

%%%%%%%%%%%%%%%%%%%%%%%% convergence criteria %%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  if max(abs(V_temp_new - V_temp)) < eps
    fprintf(1 , 'The value function has converged in %2.0f iterations.\n' , l);
    break;
  end

%%%%%%%%%%%% Update the value function if criteria not met %%%%%%%%%%%%%%%

  V_temp = V_temp_new;           % Update value if criterion not met

  end
  
%%%%%%%%% Fill in the value function and policy function matrices %%%%%%%%
  
  V(: , 6 - k) = V_temp;         % Fill in the value function matrix
  Pol(: , 6 - k) = Pol_temp;     % Fill in the policy function matrix
  
  end
  
%%%%%%%%%%%% Plot the value function and the policy function %%%%%%%%%%%%% 

  subplot(3 , 1 , 1)
  plot(x , V( : , 3 : 6))
  title('value functions');
  xlabel('current wage');
  legend('h = 0.4' , 'h = 0.6' , 'h = 0.8' , 'h = 1')
  
  subplot(3 , 1 , 2)
  plot(x , lambda(Pol(: , 3 : 5)))
  title('policy functions');
  xlabel('current wage');
  legend('h = 0.4' , 'h = 0.6' , 'h = 0.8')
 
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%% Simulate wage distribution %%%%%%%%%%%%%%%%%%%%%%%

  N = 10000;                        % number of simulated wage of job offers
  
  rn = betarnd(beta0 , beta1 , N , 1);  
  % simulated wage distribution
  % For each health condition, the distribution of wages is the same. 
  % What's more, in question 2, I come to the result that the steady state
  % population in each health state is the same.
  % Therefore I will use this data for every health condition. Thus I
  % have a data of size 10000 * 5 for the whole population.
                                        
%%%%%%%%%%%%%%%%%%%% Find out the job arriving rate %%%%%%%%%%%%%%%%%%%%%%

  index2 = round(rn * n + 1);           % Find out the index of each random 
                                        % number in the x vector

  lambda_rn = lambda(Pol(index2 , :));  % Find out the corresponding 
                                        % job arriving rate through the
                                        % policy function
                                        
  lambda_rn_mean = mean(lambda_rn)      % job-finding rate by health status
  
%%%%%%%%%%%%%%%%%%% Find out the unemployment rate %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% The equations below are given in the solution to question 2 %%%%%%

  eta = 0.015;                          % rate of new borns
  g = 0.2 * ones(n , 1);                % steady state population 
                                        % distribution by health condition 
  
  
  unrate = ones(5 , 1);
  unrate(1) = (delta * g(1) + eta) ./ (delta + lambda(Pol(1 , 1)) + dh(1) + alpha); 
  % This is the unemployment rate when h = 0
  % Here I am picking lambda(Pol(1 , 1)) because in this equation, we are
  % only considering the job arrival rate for the unemployed. The job
  % transition for the employed does not come into this equation. The
  % decision for the unemployed comes from the first row of the policy
  % function.
  
  for q = 2 : 5                         % unemployment rate when h > 0
      unrate(q) = (delta * g(q) + alpha * unrate(q - 1)) / (delta + alpha + lambda(Pol(1 , q)) + dh(q) + alpha);
      % In place of 0 it should be the job searching effort when x = 0. 
      % but since the lambda's I obtained are not reasonable, I used 0 instead.
  end
  
  % This is the unemployment rate when h = 0.2, 0.4, 0.6 and 0.8
  % As before, only lambda(Pol(1 , q)) is considered because in this 
  % equation, only the job arrival rate for the unemployed matters. The job
  % transition for the employed does not come into this equation. The
  % decision for the unemployed comes from the first row of the policy
  % function.
  
  unrate

%%%%%%%%%%%%%%%%%%%% Find out the job transition rate %%%%%%%%%%%%%%%%%%%%
  
  Fw = cumsum(rn)/sum(rn);               % CDF of the same wage data set
  
  tranrate = mean(lambda_rn .* (1 - Fw)) % job transition rate

%%%%%%%%%%%% Find out and plot the observed wage distribution %%%%%%%%%%%%

  Lw = (unrate' .* lambda_rn .* Fw) ./ ((1 - unrate)' .* (delta + lambda_rn .* (1 - Fw)));
  
  L = (0: 0.0001 : 0.9999);
  
  subplot(3 , 1 , 3)
  plot(L , Lw(: , 3 : 5))
  xlabel('current wage')
  title('CDF of observed wage distribution')
  legend('h = 0.4', 'h = 0.6' , 'h = 0.8')