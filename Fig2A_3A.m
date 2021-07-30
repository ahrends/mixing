%% Reproduce figures 2A and 3A
% Paper "Data and model considerations for estimating time-varying 
% functional connectivity in fMRI" (Ahrends et al., 2021)
%
% This scripts requires that all combinations for simulated data have been
% run. If not, replace the part that loads results with part that runs the
% simulations and evaluates them (note that these are simulations + 300
% models, so it would need to be run on a computing cluster or similar)
%
% Christine Ahrends
% (Aarhus University 2020)

%% Figure 2A

example_dir = '/path/to/example/timecourse/and/example/HMM';
sim_outputdir = '/path/where/simulated/timeseries/should/be/saved';

sim_options = struct();
sim_options.these_regions = 1:10;
sim_options.n_subj = 100; 
sim_options.K = 6;
sim_options.subject_inconsistency = [0.1:0.1:1];
sim_options.state_inconsistency = [0.1:0.1:1];
sim_options.n_iter = 1;
sim_options.measures = ["staticFC", "mixing"];

% if not already run:
% [X,T] = simulate_timecourses(example_dir, sim_outputdir, sim_options);
% Results_simu = run_evaluate_HMM_simu(X, T, results_dir, sim_options);

% if already run:
str1 = ['regions' num2str(min(sim_options.these_regions)) '-' num2str(max(sim_options.these_regions))];
str2 = ['subjects' num2str(sim_options.n_subj)];
if size(sim_options.subject_inconsistency, 2)>1
    str3 = ['bs_var' num2str(min(sim_options.subject_inconsistency)) '-' num2str(max(sim_options.subject_inconsistency))];
else
    str3 = ['bs_var' num2str(sim_options.subject_inconsistency)];
end
if size(sim_options.state_inconsistency,2)>1
    str4 = ['ws_var' num2str(min(sim_options.state_inconsistency)) '-' num2str(max(sim_options.state_inconsistency))];
else
    str4 = ['ws_var' num2str(sim_options.state_inconsistency)];
end

Results1 = load([results_dir '/Results_simu_' str1 '_' str2 '_' str3 '_' str4 '.mat'], 'Results_simu');


% plot results from the first iteration
figure; subplot(1,2,1); 
surf(sim_options.state_inconsistency, sim_options.subject_inconsistency, squeeze(Results1.statFC_similarity));
xlabel('within-session variability'); ylabel('between-subject variability'); zlabel('static FC similarity');
title('static FC similarity');
subplot(1,2,2); 
surf(sim_options.state_inconsistency, sim_options.subject_inconsistency, squeeze(Results1.mean_maxFO));
xlabel('within-session variability'); ylabel('between-subject variability'); zlabel('mean maxFO');
title('mean maxFO');


%% Figure 3A

% simulate timeseries, run and evaluate HMMs for 50 regions instead of 10

sim_options = struct();
sim_options.these_regions = 1:50;
sim_options.n_subj = 100; 
sim_options.K = 6;
sim_options.subject_inconsistency = [0.1:0.1:1];
sim_options.state_inconsistency = [0.1:0.1:1];
sim_options.n_iter = 1;
sim_options.measures = {'staticFC', 'mixing'};

% if not already run:
% [X,T] = simulate_timecourses(example_dir, sim_outputdir, sim_options);
% Results_simu = run_evaluate_HMM_simu(X, T, results_dir, sim_options);

% if already run:
str1 = ['regions' num2str(min(sim_options.these_regions)) '-' num2str(max(sim_options.these_regions))];
str2 = ['subjects' num2str(sim_options.n_subj)];
if size(sim_options.subject_inconsistency, 2)>1
    str3 = ['bs_var' num2str(min(sim_options.subject_inconsistency)) '-' num2str(max(sim_options.subject_inconsistency))];
else
    str3 = ['bs_var' num2str(sim_options.subject_inconsistency)];
end
if size(sim_options.state_inconsistency,2)>1
    str4 = ['ws_var' num2str(min(sim_options.state_inconsistency)) '-' num2str(max(sim_options.state_inconsistency))];
else
    str4 = ['ws_var' num2str(sim_options.state_inconsistency)];
end

Results2 = load([results_dir '/Results_simu_' str1 '_' str2 '_' str3 '_' str4 '.mat'], 'Results_simu');

% simulate timeseries, run and evaluate HMMs for 50 regions instead of 10
% and 20 subjects instead of 100

sim_options = struct();
sim_options.these_regions = 1:50;
sim_options.n_subj = 20;
sim_options.K = 6;
sim_options.subject_inconsistency = [0.1:0.1:1];
sim_options.state_inconsistency = [0.1:0.1:1];
sim_options.n_iter = 1;
sim_options.measures = {'staticFC', 'mixing'};

% if not already run:
% [X,T] = simulate_timecourses(example_dir, sim_outputdir, sim_options);
% Results_simu = run_evaluate_HMM_simu(X, T, results_dir, sim_options);

str1 = ['regions' num2str(min(sim_options.these_regions)) '-' num2str(max(sim_options.these_regions))];
str2 = ['subjects' num2str(sim_options.n_subj)];
if size(sim_options.subject_inconsistency, 2)>1
    str3 = ['bs_var' num2str(min(sim_options.subject_inconsistency)) '-' num2str(max(sim_options.subject_inconsistency))];
else
    str3 = ['bs_var' num2str(sim_options.subject_inconsistency)];
end
if size(sim_options.state_inconsistency,2)>1
    str4 = ['ws_var' num2str(min(sim_options.state_inconsistency)) '-' num2str(max(sim_options.state_inconsistency))];
else
    str4 = ['ws_var' num2str(sim_options.state_inconsistency)];
end

Results3 = load([results_dir '/Results_simu_' str1 '_' str2 '_' str3 '_' str4 '.mat'], 'Results_simu');

% plot mean maxFO for 1. 10 regions, 100 subjects, 2. 50 regions, 100
% subjects, and 3. 50 regions, 20 subjects

figure; subplot(1,3,1); 
surf(sim_options.state_inconsistency, sim_options.subject_inconsistency, squeeze(Results1.mean_maxFO));
xlabel('within-session variability'); ylabel('between-subject variability'); zlabel('mean maxFO');
title('10 regions, 100 subjects');
subplot(1,3,2);
surf(sim_options.state_inconsistency, sim_options.subject_inconsistency, squeeze(Results2.mean_maxFO));
xlabel('within-session variability'); ylabel('between-subject variability'); zlabel('mean maxFO');
title('50 regions, 100 subjects');
subplot(1,3,3); 
surf(sim_options.state_inconsistency, sim_options.subject_inconsistency, squeeze(Results3.mean_maxFO));
xlabel('within-session variability'); ylabel('between-subject variability'); zlabel('mean maxFO');
title('50 regions, 20 subjects');

