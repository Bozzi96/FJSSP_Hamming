function prob_Hamming = addHammingConstraints(prob, A, J, M, solutions, solution_postprocess)
diff_gamma = optimvar('diff_gamma', A, 'Type', 'integer', 'LowerBound', 0, 'UpperBound', 1);
diff_delta = optimvar('diff_delta', J*J*M,'Type', 'integer', 'LowerBound', 0, 'UpperBound', 1);
prob_Hamming = prob;
BigM = 10000;
 %% Hamming distance constraints
    % Hamming distance for gamma
    %diffgamma = optimconstr(2*A,1);
    gamma_old = solutions{end}.gamma;
    delta_old = solutions{end}.delta;
    % Add constraints to model the absolute difference |gamma(i) - gamma_old(i)|
    % Without using the abs() function (not supported by optimization variables)
    % for i = 1:A
    %     if ~gamma_old(i)
    %     diffgamma(i) = diff_gamma(i) == prob_Hamming.Variables.gamma(i) - gamma_old(i);
    %     else
    %     diffgamma(A+i) = diff_gamma(i) == gamma_old(i) - prob_Hamming.Variables.gamma(i);
    %     end
    % end
    % prob_Hamming.Constraints.diffgamma = diffgamma;
    prob_Hamming.Constraints.diffgamma = diff_gamma == [prob_Hamming.Variables.gamma(find(int8(gamma_old)<=0)) - gamma_old(find(int8(gamma_old)<=0)) ;...
        gamma_old(find(int8(gamma_old)>= 1)) - prob_Hamming.Variables.gamma(find(int8(gamma_old) >= 1))];

    % Hamming distance for delta
    % Add constraints to model the absolute difference |delta(i) - delta_old(i)|
    % Without using the abs() function (not supported by optimization variables)
    %diffdelta = optimconstr(2*J,2*J,2*M);
    % for m=1:M
    %     for i=1:J-1
    %         for j=i+1:J
    %             if ~delta_old(i,j,m)
    %                 diffdelta(i,j,m) = diff_delta(i,j,m) == prob_Hamming.Variables.delta(i,j,m) - delta_old(i,j,m);
    %             else
    %                 diffdelta(J+i,J+j,M+m) = diff_delta(i,j,m) == delta_old(i,j,m) - prob_Hamming.Variables.delta(i,j,m);
    %             end
    %         end
    %     end
    % end
    % prob_Hamming.Constraints.diffdelta = diffdelta;
    prob_Hamming.Constraints.diffdelta0 = diff_delta(find(int8(delta_old) <= 0)) == prob_Hamming.Variables.delta(find(int8(delta_old) <= 0)) - delta_old(find(int8(delta_old) <= 0));
    prob_Hamming.Constraints.diffdelta0 = diff_delta(find(int8(delta_old) >= 1)) == delta_old(find(int8(delta_old)>=1)) - prob_Hamming.Variables.delta(find(int8(delta_old)>=1));
    % Calculate the Hamming distance
    hammingDist = sum(diff_gamma) + sum(sum(sum(diff_delta)));
    % Add the constraints for the minimum and maximum Hamming distance
    prob_Hamming.Constraints.minDist = hammingDist >= 1;
    prob_Hamming.Constraints.maxDist = hammingDist <= 30;
 
    %% Constraints Hamming 2.0
    % constraints_xor = optimconstr(1,size(solutions,2));
    % for h=1:size(solutions,2)
    %     gamma_old = solutions{h}.gamma;
    %     delta_old = solutions{h}.delta;
    %     somma_xor = 0;
    %     for i=1:A
    %         if int8(gamma_old(i)) == 1
    %             somma_xor = somma_xor + prob_Hamming.Variables.gamma(i);
    %         else
    %             somma_xor = somma_xor -BigM*prob_Hamming.Variables.gamma(i);
    %         end
    %     end
    %     for m=1:M
    %         for i=1:J-1
    %             for j=i+1:J
    %                 if int8(delta_old(i,j,m)) == 1
    %                     somma_xor = somma_xor + prob_Hamming.Variables.delta(i,j,m);
    %                 else
    %                     somma_xor = somma_xor -BigM * prob_Hamming.Variables.delta(i,j,m);
    %                 end
    %             end
    %         end
    %     end
    %     constraints_xor(h) = somma_xor <= 0;
    % end
    % prob_Hamming.Constraints.constraints_xor = constraints_xor;

    prob_Hamming.Constraints.improve_Completion = prob_Hamming.Variables.C <= solution_postprocess{end}.C;
end