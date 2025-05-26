function compareSolutions(sol1, sol2)
    % figure()
    % subplot(2,1,1)
    % scatter(1:length(sol1.gamma), sol1.gamma)
    % hold on
    % scatter(1:length(sol2.gamma), sol2.gamma, [], "red")
    diff_gamma = sum(int8(sol1.gamma) ~= int8(sol2.gamma))
    delta_sol1 = deltaToArray(sol1.delta);
    delta_sol2 = deltaToArray(sol2.delta);
    % subplot(2,1,2)
    % scatter(1:length(delta_sol1), delta_sol1)
    % hold on
    % scatter(1:length(delta_sol2), delta_sol2, [], "red")
    diff_delta = sum(int8(delta_sol1) ~= int8(delta_sol2))
end