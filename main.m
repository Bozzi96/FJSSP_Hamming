clear; clc; close all;
seed = rng(1); % control random number generation
load("Datasets/hurink_edata_la03.mat") % Optimal 550 according to literature, 483 my algo
warning('off','all')
% Parameters of the problem
%%% SIMPLE CASE TO TEST THE CORRECTNESS
P = ones(10,5);
%%% IEEE ACCESS CASE STUDY
% P = [9 5 7 10 4 12
%      4 7 3 7 1 10
%      5 7 6 3 10 1
%      4 3 10 6 4 5
%      2 4 7 3 5 2
%      1 6 5 3 6 8];
% 
% G_init = [1 2 3 4 5
%            1 3 5 6 0
%            1 2 3 4 6
%            1 2 4 5 0
%            1 2 5 6 0
%            1 2 3 5 0
%            1 3 4 6 0
%            1 2 6 5 3
%            1 4 6 5 3
%            1 2 4 5 6
%            1 3 4 5 0
%            1 2 3 4 6
%            1 3 4 6 0];
% G_j = [1 1 2 2 3 3 3 4 5 5 6 6 6]';
% S0 = [0 2 4 7 10 12]';

%% Pre processing of data
[G, P, M_init, aux, aux_alt] = pre_processing_graph(G_init, P);
J = length(unique(G_j)); %jobs
M = max(max(G)); %machines
A = length(G_j);%alternatives
D = compute_D_from_graph(G_init,G_j); % disjunctive connections (2 constraints per each connection)
map_duplicate = map_duplicate_machines(G,G_init);

%% Solve problem
% Build problem and initial conditions
prob = buildOptimizationProblem(G,G_j,P, S0);
x0.C = 0;
x0.c = zeros(J,M);
x0.s = zeros(J,M);
x0.gamma = zeros(A,1);
x0.delta = zeros(J,J,M);
% firstSol = solve(prob,x0);
% gamma_old = firstSol.gamma;
% delta_old = firstSol.delta;
%prob_Hamming = addHammingConstraints(prob,A,J,M, gamma_old, delta_old);
% x0.gamma = gamma_old;
% x0.delta = delta_old;
%Initial conditions
x0eval = x0;
clear x0eval.gamma;
clear x0eval.delta;
x0.diff_gamma = zeros(A,1);
x0.diff_delta = zeros(J,J,M);


MAX_ITER = 25;
solutions = {};
prob_Hamming = {};
prob_postProcess = {};
solution_postprocess = {};
completions = [];
bestC = inf;
bestInd = 0;
opts=optimoptions(@fmincon, 'Display','off');
% Find optimal solution
tic
sol_opt = FJSSP_optimalSol(G,G_j,P,S0);
time_for_opt = toc;
% Find feasible solution
tic
[solutions{1}, ~, prob_Hamming{1}] = FJSSP_feasibleSol(G,G_j,P,S0);
solution_postprocess{1} = FJSSP_evaluateGammaDelta(G,G_j,P, S0, solutions{1}.gamma, solutions{1}.delta); % Align first solution
completions(1) = solution_postprocess{1}.C;
time_for_heuristic(1) = toc;
%% Loop to iterate between solutions with a Hamming distance between each other
%oldoptions = optimoptions(@fmincon,'Algorithm','sqp','MaxIterations',5);
for i=2:MAX_ITER
  tic
  gamma_old = solutions{end}.gamma;
  delta_old = solutions{end}.delta;
  %%% Build the optimization problem and then add only Hamming constraints
  %%% for faster computation time
  prob_Hamming{i} = addHammingConstraints(prob, A, J, M, solutions, solution_postprocess);

  solutions{i} = solve(prob_Hamming{i}, x0, 'Options', opts);
  % There are no solutions within the maximum Hamming distance
  if isempty(solutions{i}.C)
      solutions(i) = [];
      break
  end
  % Post-processing of solution
  %solution_postprocess{i} = FJSSP_evaluateGammaDelta(G,G_j,P, S0, solutions{i}.gamma, solutions{i}.delta);
  prob_postProcess{i} = addEvaluationGammaDelta(prob, A, J, M, solutions{i}.gamma, solutions{i}.delta);
  solution_postprocess{i} = solve(prob_postProcess{i}, x0eval, 'Options', opts);
  completions(i) = solution_postprocess{i}.C;
  if completions(i) < bestC
      bestC = completions(i);
      bestInd = i;
  end
  time_for_heuristic(i) = toc;
end

%% Plot 
load("results.mat", "sol_opt", "time_for_opt")
%% Plot optimal solution
figure()
graph_Gantt(sol_opt, G_init, G_j, P, sol_opt.gamma, M_init);
ax = gca; 
ax.FontSize = 16;
set(gca,'LooseInset',get(gca,'TightInset'));
%% Plot best solution found with Hamming heuristic
figure()
graph_Gantt(solutions{end}, G_init, G_j, P, solutions{end}.gamma, M_init);
ax = gca; 
ax.FontSize = 16; 
set(gca,'LooseInset',get(gca,'TightInset'));
%% Plot computational time - C
figure
cumul_time_heuristic = cumsum(time_for_heuristic);
plot(cumsum(time_for_heuristic(1:end)), completions(1:end), 'o-', Color="blue", MarkerFaceColor="blue", LineWidth=2)
hold on
plot(time_for_heuristic(1), completions(1), 'o', Color='green', MarkerFaceColor="green", LineWidth=3)
scatter(time_for_opt, sol_opt.C, "filled", 'MarkerFaceColor','r', MarkerEdgeColor='r',  LineWidth=3);
hold on
grid on
text(time_for_heuristic(1)+1,completions(1)+1 ,"(21.2" +  ', ' + num2str(completions(1)) + ")",  FontSize=16)
text(cumul_time_heuristic(end)-5,completions(end)+20,"(" + num2str(round(cumul_time_heuristic(end),1)) + ', ' + num2str(completions(end)) + ")",  FontSize=16)
text(time_for_opt-5,sol_opt.C-15,"(" + num2str(round(time_for_opt,1)) + ', ' + num2str(sol_opt.C) + ")", FontSize=16)

% xticks(round([cumsum(time_for_heuristic) time_for_opt], 2))
% yticks([sol_opt.C sort(unique(completions))])
ylim([sol_opt.C-30 completions(1)+30])
ylabel("Makespan")
xlabel("Computation time [s]")
% mySpline = spline([cumsum(time_for_heuristic) time_for_opt], [completions sol_opt.C]);
% 
% % Define a fine grid for evaluation
% %time_fine = linspace(min(cumsum(time_for_heuristic)), max(cumsum(time_for_heuristic)), 100);
% time_fine = linspace(time_for_heuristic(1), time_for_opt, 100);
% 
% % Evaluate the spline at the fine grid
% completion_fine = ppval(mySpline, time_fine);
% 
% % Plot the result
% plot(time_fine, completion_fine, 'k-');
legend("Heuristic", "Initial solution",  "Optimal solution")
ax = gca; 
ax.FontSize = 24; 
set(gca,'LooseInset',get(gca,'TightInset'));