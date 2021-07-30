function run_HMM_simu(X,T,options)

%%

clear hmm_options
hmm_options.order = 0;
hmm_options.covtype = 'full'; %('full' for covariance, 'uniquefull' for no covariance)
hmm_options.zeromean = 1; % (0 to model mean, 1 to model only covariance)
hmm_options.standardise = 1;
hmm_options.dropstates = 0;
hmm_options.K = options.K;

for i = 1:options.n_iter
    for ii = 1:size(options.subject_inconsistency,2)
        for iii = 1:size(options.state_inconsistency,2)
            clear HMM
            [HMM.hmm, HMM.Gamma, ~, HMM.vpath, ~, ~, HMM.fehist] = hmmmar(X{i,ii,iii}, T, hmm_options);
            maxFO = getMaxFractionalOccupancy(HMM.Gamma, T, hmm_options);
            mean_maxFO{i,ii,iii} = nanmean(maxFO);
        end
    end
end

end