for task={'spatial', 'conceptual'}
    m=csvread(strcat('modelResults/', task{1}, 'diffevidence.csv'));
    %has to be the negative loss
    [alpha,exp_r,xp,protectedExceedenceProb,bor] = bms(-m);
    %probability of exceedance
    probofexp=protectedExceedenceProb;
    %Bayesian omnibus risk
    bayesianomni=bor;
    csvwrite(strcat('modelResults/',task{1},'PXP.csv'), protectedExceedenceProb)
end