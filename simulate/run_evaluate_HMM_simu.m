function Results_simu = run_evaluate_HMM_simu(X, T, results_dir, options)
% This script runs and evaluates HMMs and/or static FC similarity of
% simulated timeseries with different options. Requires that
% simulate_timecourses.m has been run.
%
% INPUT:
% X:                simulated time series with these_regions, n_subj, and
%                   specified amount of subject and state variability 
%                   (saved as cell of size n_iter x
%                   size(subject_inconsistency) x
%                   size(state_inconsistency))
% T:                vector containing the number of timepoints for each
%                   subject/session
% results_dir:      Directory where results should be saved
% options:          structure containing fields:
%       these_regions:  vector of regions to be used in example time course 
%                       (e.g. [1:10] to use the first ten regions from a 
%                       parcellation with >=10 regions)
%       n_subj:         number of subjects to simulate
%       K:              number of states in HMM (in this script, K should
%                       correspond to the number of states in example HMM 
%                       provided, but parameters can be changed below)
%       subject_inconsistency: amount of variability added between subjects
%                       ([0.1:1] seems to be a sensible range)
%       state_inconsistency: amount of variability added between states (over
%                       time, similar range as between subjects)
%       n_iter:         number of iterations 
%
% OUTPUT:
% Results_simu:  structure containing results to evaluate static FC similarity
%           and/or model stasis:
%           (if options.measures='staticFC')
%           statFC_similarity: static FC similarity between subjects
%           (if options.measures='mixing')
%           mean_maxFO: mean maximum fractional occupancy (measure of model
%           stasis)
%
% Christine Ahrends
% (Aarhus University 2020)
%

%% check that options are properly set up

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

if ~isfield(options, 'measures')
    options.measures = ["staticFC", "mixing"];
    warning('Measures not specified, setting to default (static FC similarity and mixing)');
end

if any(~isfield(options, {'n_iter', 'n_subj', 'subject_inconsistency', 'state_inconsistency', 'K', 'these_regions'}))
    error('Cannot run and evaluate HMMs for this simulated timeseries because necessary options are missing.')
end

%% evaluate static FC similarity
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

%% run HMMs and evaluate mixing
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

%% save Results for specified options
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