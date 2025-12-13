function [Sh, Ih, Rh,Qh,Eh] = simulate_absir(Iv0, T, recovery_rate, quarantine_percent)

% Simulate agent-based transmission model. Uses a graph to represent social
% connectivity between agents.
%
% Note: You may find it helpful to summarize the state history matrices via
% summation. For instance sum(Ih, 1) will return the total number of
% infected persons at each timestep in the simulation.
%
% Inputs
%   M (square matrix): Adjacency matrix
%   Iv0 (column vector): Initial infection state
%   T (integer): Number of timesteps to simulate
%   infection_rate (float): Infection rate (probabilistic)
%   recovery_rate (float): Recovery rate (probabilistic)
%
% Returns
%   Sh (matrix): Susceptible state history
%   Ih (matrix): Infected state history
%   Rh (matrix): Recovered state history
%   Qh(matrix): Quarantined state history
%   E_h(matrix): Exposed state history

   % Setup
   dim = length(Iv0);  % Dimensions of initial state
   Ih = zeros(dim, T); % Infection history
   Ih(:, 1) = Iv0;         % Record the initial state to infection history
   Rh = zeros(dim, T); % Recovery history
   Eh = zeros(dim, T); % Exposed history
   Qh = zeros(dim, T); % Quarantine history

   % Compute infections with NB distribution
   R0 = 2.6;
   r = 0.16;
   p = (r)/(R0+r);
   pd_sse = makedist('NegativeBinomial', 'r', r, 'p', p);  % Negative Binomial distribution

   % Construct an "action" helper function
   function [In, Rn, En, Qn] = action(I, R, E, Q)

       % Compute random sample of infected individuals 
       Samples = sum(I);
       % Random realizations
       Sample = random(pd_sse,Samples,1);
       new_infect_count = sum(Sample);

       %Compute which people can go from exposed to infected. 25% change rn
       v_exposed_to_infected = rand(dim,1)<=0.25;

       %Compute which people can recover
       v_recover = rand(dim, 1) <= recovery_rate;

       %Compute which people can go from infected to quarantined
       v_quarantine = rand(dim,1) <= quarantine_percent;

       %Compute which people can go susceptible to exposed
       v_exposed(1:new_infect_count,1) = 1;
       v_exposed(new_infect_count+1:dim,1) = 0;
       v_exposed = v_exposed(randperm(length(v_exposed)));
     
       % Update all the stuff (._.)
       En = E | (v_exposed & ~I & ~E & ~R & ~Q);
       In = I | (v_exposed_to_infected & E);
       Qn = Q | (v_quarantine & I);
       Rn = R | (v_recover & (I | Q));
       En = En & ~In;
       Qn = Qn & ~Rn;
       In = In & ~Rn & ~Qn;
       
       
   end
   % Run simulation
   for i = 2:T
       [Ih(:, i), Rh(:, i), Eh(:, i),Qh(:, i)] = action(Ih(:, i-1), Rh(:, i-1), Eh(:, i-1),Qh(:, i-1));
   end
  
   % Compute susceptible history
   Sh = ones(dim, T) - Ih - Rh - Eh - Qh;

   %Verification Line (leave commented)
   %any((Sh+Ih+Rh+Eh+Qh - ones(dim,T)),'all')

end
% Implement distributions for the New Infection Count
