%% Main script for paper 
% "Data and model considerations for estimating time-varying functional connectivity in fMRI" 
% (Ahrends et al., 2021)
%
% This script requires the HMM-MAR toolbox (available at
% https://github.com/OHBA-analysis/HMM-MAR)
% and FSLNets (available at https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FSLNets)
% and the remaining scripts and example data from https://github.com/ahrends/mixing


%% set up paths

scriptdir = '/path/to/scripts'; % directory containing mixing scripts, HMM-MAR toolbox, and FSLNets
addpath(genpath(scriptdir))

tc_dir = '/path/to/timecourses'; % directory containing timecourses in each 
% parcellation called ['hcp1003_REST1_LR_' parcellation '.mat']
% time courses should be cells called data with each cell corresponding to
% one subject's timecourse as matrix in the format timepoints x parcels
HMM_outputdir = '/directory/where/HMMs/should/be/saved';

%% specify options: which parameters should be run/evaluated?
% choose how secondary parameters should be manipulated before running HMM
% (combinations used in manuscript as comments behind each option)

options = struct();
options.parcellation = 'groupICA50'; % 'groupICA50', 'groupICA100', 'PROFUMO50', 'Yeo100', 'DK80'
options.k = 12; 
options.nsubs = 100; % 50, 100, 200
options.nts = 'all'; % 200, 500, 1200 (all)
options.sr = 1; % 1, 2, 3
options.nregions = 'all'; % 10, 25, 50, all
% optionally, use either PCA or HMM-PCA to reduce dimensionality:
% options.pca = 0.9; % 0.7, 0.8, 0.9, 0.93, 0.96, 0.99
% options.lowrank = 10; % use same number of components as in PCA runs
% (saved in Results.nPCs after running and evaluating HMMs with
% options.pca); additionally: 4, 8, 12, 16

%% run HMMs

HMM = run_HMM_params(tc_dir, HMM_outputdir, options);

%% evaluate HMMs

results_dir = '/directory/to/where/results/should/be/saved';

options.measures = ["staticFC", "mixing"];

Results = evaluate_HMM_params(tc_dir, HMM_outputdir, results_dir, options);

%% simulate timeseries

example_dir = '/path/to/example/timecourse/and/example/HMM';
sim_outputdir = '/path/where/simulated/timeseries/should/be/saved';

sim_options = struct();
sim_options.these_regions = 1:10; %1:50
sim_options.n_subj = 100; %20
sim_options.K = 6;
sim_options.subject_inconsistency = [0.1:0.1:1];
sim_options.state_inconsistency = [0.1:0.1:1];
sim_options.n_iter = 1;

[X,T] = simulate_timecourses(example_dir, sim_outputdir, sim_options);


%% run and evaluate HMMs on simulated timeseries

sim_options.measures = ["staticFC", "mixing"];

Results_simu = run_evaluate_HMM_simu(X, T, results_dir, sim_options);

%% DO NOT RUN
% This part runs and evaluates HMMs for all combinations of variables and
% saves them in a table, which is used as input for the SEM. Requires a
% computing cluster or similar since it runs 1890 models (excl. PCA and
% HMM-PCA runs).


% all combinations real data:

aggr_parcellations = cell(0,1);
aggr_nsubs = zeros(0,1);
aggr_nts = cell(0,1);
aggr_sr = zeros(0,1);
aggr_nregions = cell(0,1);
aggr_staticFC_similarity = zeros(0,1);
aggr_mean_maxFO = zeros(0,1);

clear options
k = 12;
measures = ["staticFC", "mixing"];
for parcellation = {'groupICA50', 'groupICA100', 'PROFUMO50', 'Yeo100', 'DK80'}
    for nsubs = [50, 100, 200]
        for nts = {200, 500, "all"}
            for sr = [1, 2, 3]
                for nregions = {10, 25, 50, "all"}
                    options = struct('k', k, 'measures', measures, ...
                        'parcellation', parcellation, 'nsubs', nsubs, ...
                        'nts', nts{1,1}, 'sr', sr, 'nregions', nregions{1,1});
                    % exclude sampling subsets of 50 regions for
                    % parcellations with 50 parcels
                    if ~strcmp(options.nregions, 'all')
                        if ~((options.nregions==50) && ...
                                (strcmp(options.parcellation, 'groupICA50') || strcmp(options.parcellation, 'PROFUMO50')))
                            % run 5 iterations each when sampling subsets from a
                            % parcellation:
                            for i = 1:5
                                HMM = run_HMM_params(tc_dir, HMM_outputdir, options);
                                Results = evaluate_HMM_params(tc_dir, HMM_outputdir, results_dir, options);
                                aggr_parcellations{end+1,1} = parcellation;
                                aggr_nsubs(end+1,1) = nsubs;
                                aggr_nts{end+1,1} = nts;
                                aggr_sr(end+1,1) = sr;
                                aggr_nregions{end+1,1} = nregions;
                                aggr_staticFC_similarity(end+1,1) = Results.statFC_similarity;
                                aggr_mean_maxFO(end+1,1) = Results.mean_maxFO;
                            end
                        end
                    else
                        HMM = run_HMM_params(tc_dir, HMM_outputdir, options);
                        Results = evaluate_HMM_params(tc_dir, HMM_outputdir, results_dir, options);
                        aggr_parcellations{end+1,1} = parcellation;
                        aggr_nsubs(end+1,1) = nsubs;
                        aggr_nts{end+1,1} = nts;
                        aggr_sr(end+1,1) = sr;
                        aggr_nregions{end+1,1} = nregions;
                        aggr_staticFC_similarity(end+1,1) = Results.staticFC_similarity;
                        aggr_mean_maxFO(end+1,1) = Results.mean_maxFO;
                    end
                end
            end
        end
    end
end

vars = {'parcellation', 'nsubs', 'nts_tmp', 'sr', 'nregions_tmp', 'staticFC_similarity', 'mean_maxFO'};
aggr_table = table(aggr_parcellations, aggr_nsubs, aggr_nts, aggr_sr, aggr_nregions, aggr_staticFC_similarity, aggr_mean_maxFO);
aggr_table.Properties.VariableNames = vars;

for j = 1:size(aggr_table,1)
    if strcmp(aggr_table.nts_tmp{j}{1,1}, 'all')
        aggr_table.nts(j) = 1200;
    else
        aggr_table.nts(j) = aggr_table.nts_tmp{j}{1,1};
    end
    if strcmp(aggr_table.nregions_tmp{j}{1,1}, 'all')
        if strcmp(aggr_table.parcellation{j}, 'groupICA50') || strcmp(aggr_table.parcellation{j}, 'PROFUMO50')
            aggr_table.nregions(j) = 50;
        elseif strcmp(aggr_table.parcellation{j}, 'DK80')
            aggr_table.nregions(j) = 80;
        elseif strcmp(aggr_table.parcellation{j}, 'groupICA100') || strcmp(aggr_table.parcellation{j}, 'Yeo100')
            aggr_table.nregions(j) = 100;
        end
    else
        aggr_table.nregions(j) = aggr_table.nregions_tmp{j}{1,1};
    end
end

aggr_table = removevars(aggr_table, {'nts_tmp', 'nregions_tmp'});

% save table with all results (used as input to SEM.R)
writetable(aggr_table,'realdata_mixing.csv','Delimiter',',');


%% DO NOT RUN


% all combinations simulated data

aggr_parcellations = cell(0,1);
aggr_nsubs = zeros(0,1);
aggr_nts = zeros(0,1);
aggr_sr = zeros(0,1);
aggr_nregions = zeros(0,1);
aggr_subject_inconsistency = zeros(0,1);
aggr_state_inconsistency = zeros(0,1);
aggr_staticFC_similarity = zeros(0,1);
aggr_mean_maxFO = zeros(0,1);

parcellation = 'groupICA50';
nts = 1200;
sr = 1;
K = 6;
subject_inconsistency = [0.1:0.1:1];
state_inconsistency = [0.1:0.1:1];
n_iter = 1;
measures = ["staticFC", "mixing"];

clear sim_options

for these_regions = {[1:10],[1:50]}
    clear subs
    if max(these_regions{1,1})==10
        subs = 100;
    elseif max(these_regions{1,1})==50
        subs = [20,100];
    end
    for n_subj = subs  
        sim_options = struct('these_regions', these_regions{1,1}, 'n_subj', n_subj, ...
            'K', K, 'subject_inconsistency', subject_inconsistency, ...
            'state_inconsistency', state_inconsistency, 'n_iter', n_iter, ...
            'measures', measures);
        [X,T] = simulate_timecourses(example_dir, sim_outputdir, sim_options);
        Results_simu = run_evaluate_HMM_simu(X, T, results_dir, sim_options);
        for ii = 1:size(subject_inconsistency,2)
            for iii = 1:size(state_inconsistency,2)
                aggr_parcellations{end+1,1} = parcellation;
                aggr_nsubs(end+1,1) = n_subj;
                aggr_nts(end+1,1) = nts;
                aggr_sr(end+1,1) = sr;
                aggr_nregions(end+1,1) = size(these_regions{1,1},2);
                aggr_subject_inconsistency(end+1,1) = subject_inconsistency(ii);
                aggr_state_inconsistency(end+1,1) = state_inconsistency(iii);
                aggr_staticFC_similarity(end+1,1) = Results_simu.statFC_similarity(n_iter,ii,iii);
                aggr_mean_maxFO(end+1,1) = Results_simu.mean_maxFO(n_iter,ii,iii);
            end
        end
    end
end

vars = {'parcellation', 'nsubs', 'nts', 'sr', 'nregions', 'subject_inconsistency', 'state_inconsistency', 'staticFC_similarity', 'mean_maxFO'};
aggr_table_simu = table(aggr_parcellations, aggr_nsubs, aggr_nts, ...
    aggr_sr, aggr_nregions, aggr_subject_inconsistency, ...
    aggr_state_inconsistency,aggr_staticFC_similarity,aggr_mean_maxFO);
aggr_table_simu.Properties.VariableNames = vars;

% save table with all results from simulated data (used as input to SEM.R)
writetable(aggr_table_simu,'simudata_mixing.csv','Delimiter',',');

