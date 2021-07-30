function Results_simu = run_evaluate_HMM_simu(X, T, results_dir, options)

%%

if nargin < 3
    results_dir = 'Results_mixing_test';
end
if nargin < 4 || isempty(options)
    warning('Options not provided, setting to default')
    options = struct();
    options.these_regions = 1:10;
    options.n_subj = 100;
    options.K = 6;
    options.subject_inconsistency = [0.1:0.1:1];
    options.state_inconsistency = [0.1:0.1:1];
    options.n_iter = 1;
    options.measures = ["staticFC", "mixing"];
end

%%
if ~isfield(options, 'measures')
    options.measures = ["staticFC", "mixing"];
    warning('Measures not specified, setting to default (static FC similarity and mixing)');
end

if any(~isfield(options, {'n_iter', 'n_subj', 'subject_inconsistency', 'state_inconsistency', 'K', 'these_regions'}))
    error('Cannot run and evaluate HMMs for this simulated timeseries because necessary options are missing.')
end

if any(strcmp(options.measures, 'staticFC'))
    
    n_ts = T(1,1);
    
    for i = 1:options.n_iter
        for ii = 1:size(options.subject_inconsistency,2)
            for iii = 1:size(options.state_inconsistency,2)
                data_FC = cell(1,options.n_subj);
                for n = 1:options.n_subj
                    data_FC{1,n} = X{i,ii,iii}((n-1)*n_ts+1:(n-1)*n_ts+n_ts,:)';
                end
                condNames = cell(1,1);
                condNames{1,1} = 'REST1';
                [C_concat, groupmean] = make_connectomes(data_FC, [], condNames);
                Results_simu.statFC_similarity(i,ii,iii) = corr(C_concat(:),groupmean(:));
            end
        end
    end
end


if any(strcmp(options.measures, 'mixing'))
    
    clear hmm_options
    hmm_options.order = 0;
    hmm_options.covtype = 'full'; %('full' for covariance, 'uniquefull' for no covariance)
    hmm_options.zeromean = 1; % (0 to model mean, 1 to model only covariance)
    hmm_options.standardise = 1;
    hmm_options.dropstates = 0;
    hmm_options.K = options.K;
    hmm_options.useParallel = 0;
    
    for i = 1:options.n_iter
        for ii = 1:size(options.subject_inconsistency,2)
            for iii = 1:size(options.state_inconsistency,2)
                clear HMM
                [HMM.hmm, HMM.Gamma, ~, HMM.vpath, ~, ~, HMM.fehist] = hmmmar(X{i,ii,iii}, T, hmm_options);
                maxFO = getMaxFractionalOccupancy(HMM.Gamma, T, hmm_options);
                Results_simu.mean_maxFO(i,ii,iii) = nanmean(maxFO);
            end
        end
    end 
end

str1 = ['regions' num2str(min(options.these_regions)) '-' num2str(max(options.these_regions))];
str2 = ['subjects' num2str(options.n_subj)];
if size(options.subject_inconsistency, 2)>1
    str3 = ['bs_var' num2str(min(options.subject_inconsistency)) '-' num2str(max(options.subject_inconsistency))];
else
    str3 = ['bs_var' num2str(options.subject_inconsistency)];
end
if size(options.state_inconsistency,2)>1
    str4 = ['ws_var' num2str(min(options.state_inconsistency)) '-' num2str(max(options.state_inconsistency))];
else
    str4 = ['ws_var' num2str(options.state_inconsistency)];
end

if ~isdir(results_dir); mkdir(results_dir); end
save([results_dir '/Results_simu_' str1 '_' str2 '_' str3 '_' str4 '.mat'], 'Results_simu', 'options')

end