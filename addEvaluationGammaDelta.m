function prob_Eval = addEvaluationGammaDelta(prob, A, J, M, gamma, delta)
    prob_Eval = prob;
    cons_evalGamma = optimconstr(A,1);
    % for i=1:A
    %     cons_evalGamma(i) = prob_Eval.Variables.gamma(i) == gamma(i);
    % end
    prob_Eval.Constraints.cons_evalGamma = prob_Eval.Variables.gamma == gamma;
    % end;
    
    cons_evalDelta = optimconstr(J,J,M);
    % for i=1:J
    %     for j=1:J
    %         for m=1:M
    %             cons_evalDelta(i,j,m) = prob_Eval.Variables.delta(i,j,m) == delta(i,j,m);
    %         end
    %     end
    % end
    prob_Eval.Constraints.cons_evalDelta = prob_Eval.Variables.delta == delta;

    % Cost function
    prob_Eval.Objective = prob_Eval.Variables.C;%0.001*prob_Eval.Variables.C+sum(sum(prob_Eval.Variables.s)) + sum(sum(prob_Eval.Variables.c));
end