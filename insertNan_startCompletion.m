clear; clc;

% Simple case study
P = ones(6,6);
G_init = [1 2 3 4 5
           1 3 5 6 0
           1 2 3 4 6
           1 2 4 5 0
           1 2 5 6 0
           1 2 3 5 0
           1 3 4 6 0
           1 2 6 5 3
           1 4 6 5 3
           1 2 4 5 6
           1 3 4 5 0
           1 2 3 4 6
           1 3 4 6 0];
G_j = [1 1 2 2 3 3 3 4 5 5 6 6 6]';
S0 = zeros(6,1);

% Find feasible solution and post process it to remove unnecessary idle times
sol_feas = FJSSP_feasibleSol(G_init, G_j, P, S0);
sol_postProc = FJSSP_evaluateGammaDelta(G_init,G_j,P, S0, sol_feas.gamma, sol_feas.delta); % Align first solution

% Create vector to store values or NaN if job does not pass through machine
start = sol_postProc.s;
completion = sol_postProc.c;
for i=1:size(start,1)
    for j=1:size(start,2)
        if start(i,j) == completion(i,j)
            start(i,j) = NaN;
            completion(i,j) = NaN;
        end
    end
end

% Find optimal solution (for comparison)
sol_opt = FJSSP_optimalSol(G_init, G_j, P, S0);

% Plot
figure()
graph_Gantt(sol_postProc, G_init, G_j, P, gamma, 6)
figure()
graph_Gantt(sol_opt, G_init, G_j, P, sol_opt.gamma, 6)