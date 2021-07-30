function [X,T] = simulate_timecourses(example_dir, outputdir, options)
% simulate time series with specified between- and within-subject variability
% for mixing paper: these_regions=[1:10] and [1:50]; n_subj=20 and 100; 
% K=6, subject_inconsistency=[0.1:0.1:1]; state_inconsistency=[0.1:0.1:1]; 
% on HCP REST1 in groupICA50
%
% INPUT:
% example_dir:      directory where some example timecourse
%                   (example_tc.mat) and example HMM (example_HMM.mat) can
%                   be found
% outputdir:        directory where simulated timeseries should be saved
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
% (NOTE: change code if X does not need to be saved across the range of 
% iterations and noise vectors to retain only the last X)
% (NOTE: subject and state noise can be weighted by changing the weight
% parameter w in the script)
%
% OUTPUT:
% X:                simulated time series with these_regions, n_subj, and
%                   specified amount of subject and state variability 
%                   (saved as cell of size n_iter x
%                   size(subject_inconsistency) x
%                   size(state_inconsistency))
%
% Christine Ahrends, Diego Vidaurre 
% (Aarhus University 2020)
%

%%
rng('default') %to replicate - change to 'shuffle' to randomise

if nargin < 2
    outputdir = 'simulations_mixing_test';
end
if nargin < 3 || isempty(options)
    warning('Options not provided, setting to default')
    options = struct();
    options.these_regions = 1:10;
    options.n_subj = 100;
    options.K = 6;
    options.subject_inconsistency = [0.1:0.1:1];
    options.state_inconsistency = [0.1:0.1:1];
    options.n_iter = 1;
end

%% load example time course and example HMM

load([example_dir '/example_tc.mat']) % some actual time course in format time x regions, e.g. HCP rest in a parcellation (no. regions >= max(these_regions))
load([example_dir '/example_HMM.mat']) % some HMM, either run on that time course with k = K (here), or replace HMM-variables depending on K later

%% set up variables

if ~isfield(options, 'these_regions')
    options.these_regions = 1:10;
    warning('Regions not specified, setting to default (1:10)');
end
nregions = size(options.these_regions,2);
n_ts = size(example_tc,1); % change variable to use fewer time points
if ~isfield(options, 'n_subj')
    options.n_subj = 100;
    warning('Number of subjects not specified, setting to default (50)');
end
clear T
T = zeros(1, options.n_subj);
for i=1:options.n_subj
    T(1,i) = n_ts;
end

%% add noise and simulate time courses

% get covariance matrix of first subject
C1 = cov(example_tc(1:n_ts,options.these_regions));
[U,S,V] = svd(C1); % C1 = U*S*V'
U1 = U;

if ~isfield(options, 'n_iter')
    options.n_iter = 1;
    warning('Number of iterations not specified, setting to default (10)');
end
if ~isfield(options, 'subject_inconsistency')
    options.subject_inconsistency = [0.1:0.1:1];
    warning('Between-subject variability not specified, setting to default');
end
if ~isfield(options, 'state_inconsistency')
    options.state_inconsistency = [0.1:0.1:1];
    warning('Within-session variability not specified, setting to default');
end
if ~isfield(options, 'K')
    options.K = 6;
    warning('Number of states not specified, setting to default (K=6)');
end

if options.n_iter > 1
    rng('shuffle')
end

for r = 1:options.n_iter
    for i=1:size(options.subject_inconsistency,2)
        for n = 1:options.n_subj
            for j = 1:nregions
                U1(:,j) = U(:,j) + V(j,j) * options.subject_inconsistency(i) * randn(nregions,1); % add random noise (subject inconsistency amount) to eigenvectors
            end
            C12{n} = U1 * S * U1';
            X1((n-1)*n_ts+1:(n-1)*n_ts+n_ts,:) = mvnrnd(zeros(1,nregions),C12{n},n_ts); % generate random time course with covariance = new, noisy covariance matrices
        end
        
        for jj=1:size(options.state_inconsistency,2)
            for k = 1:options.K
                for j = 1:nregions
                    U1(:,j) = U(:,j) + V(j,j) * options.state_inconsistency(jj) * randn(nregions,1); % add random noise (state inconsistency amount) to eigenvectors
                end
                C22{k} = U1 * S * U1';
                % reshape this if K is different from k used in HMM:
                example_HMM.hmm.state(k).Omega.Gam_rate = C22{k} * (n_ts*options.n_subj); % create new Gamma with noisy covariance matrices
                example_HMM.hmm.state(k).Omega.Gam_shape = (n_ts*options.n_subj); 
            end
            w = 0.5; % this variable can be used to weight X1 (subject noise) vs. X2 (state noise) 
            X2 = simhmmmar(T, example_HMM.hmm); % simulate time course from HMM with Gamma made with noisy covariance matrices
            X{r,i,jj} = X1 * w + X2 * (1-w); % X is the simulated time course
            % get correlatedness:
            % [~,~,e] = pca(X);
            % e = cumsum(e)/sum(e);
            % correlatedness = mean(e);
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
save([outputdir '/sim_tc_' str1 '_' str2 '_' str3 '_' str4 '.mat'], 'X', 'T', 'options')

end


