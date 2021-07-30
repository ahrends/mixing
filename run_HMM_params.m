function HMM = run_HMM_params(tc_dir, outputdir, options)
% This script runs HMMs on the HCP resting state fMRI dataset with
% varying secondary parameters
%
% INPUT:
% tc_dir:   Directory where timecourses of data can be found. This script
%           assumes that there is a file called [hcp1003_REST1_LR_' ...
%           parcellation '.mat']. In this file, there should be a separate
%           cell for each subject containing their timecourses in the
%           format timepoints x parcels.
% outputdir:    Directory where HMM should be saved
% options:  Structure containing fields:
%           parcellation: the parcellation in which timecourse were
%           extracted (default groupICA50)
%           k: the number of states k for the HMM (default 12)
%           nsubs: the number of subjects to be used (default 100)
%           nts: the number of timepoints (either less than the total
%           number of available timepoints or set to 'all' to use full
%           length timecourse) (default 'all')
%           sr: sampling rate (1 to use every timepoints, 2 to use every
%           2nd, 3 to use every 3rd,...) (default 1)
%           nregions: number of regions (parcels) used from the
%           parcellation. This will randomly sample the specified the
%           number of regions from the parcellation. Only works if nregions
%           is smaller than the original number of parcels in that
%           parcellation (use subset of regions). Set nregions to 'all' to
%           use all regions from the parcellation (default 'all')
%           (OPTIONAL fields)
%           pca: use PCA for dimensionality reduction. Specify either the 
%           number of components or amount of variance explained by PCA
%           lowrank: use HMM-PCA for dimensionality reduction. Specify the
%           number of components
% 
% OUTPUT:
% HMM:  structure containing the HMM and some useful variables to evaluate
%       mixing:
%       hmm: the HMM itself
%       Gamma: the Gamma timecourses
%       vpath: the Viterbi path
%       fehist: the free energy from the different iterations
%       maxFO: maximum fractional occupancy of each subject/session
%       FO: fractional occupancy of all states in each subject/session
%       these_regions: the regions (parcels) that this HMM was fit to
%       (either the randomly sampled subset or just a vector of all regions 
%       if options.nregions='all')
% 
%
% Christine Ahrends
% (Aarhus University 2020)
%
%

%% set up necessary variables if not specified

if nargin < 2
    outputdir = 'HMM_mixing_test';
end
if nargin < 3 || isempty(options)
    warning('Options not provided, setting to default')
    options = struct();
    options.parcellation = 'groupICA50';
    options.k = 12;
    options.nsubs = 100;
    options.nts = 'all';
    options.sr = 1;
    options.nregions = 'all';
end

%% determine parcels used in this run

% set maximum number of parcels depending on the parcellation (change this
% if other parcellations used)
if ~isfield(options, 'parcellation')
    options.parcellation = 'groupICA50';
    warning('Setting parcellation to default (groupICA50)');
end

if strcmp(options.parcellation, 'groupICA50') || strcmp(options.parcellation, 'PROFUMO50')
    max_regions = 50;
elseif strcmp(options.parcellation, 'DK80')
    max_regions = 80;
elseif strcmp(options.parcellation, 'groupICA100') || strcmp(options.parcellation, 'Yeo100')
    max_regions = 100;
end

% use all regions from the parcellation or randomly sample a subset
if ~isfield(options, 'nregions')
    options.nregions = 'all';
    warning('Setting number of regions to default (all)');
end

if strcmp(options.nregions, 'all')
    these_regions = 1:max_regions;
else
    if options.nregions > max_regions
        warning('Trying to sample more regions than this parcellation originally contains. Setting to default (all)')
        options.nregions = 'all';
    else
        rng('default')
        rand_regions = randperm(max_regions);
        these_regions = rand_regions(1:options.nregions);
%         rand_regions = randi(max_regions,[(max_regions/options.nregions),options.nregions]);
%         these_regions = rand_regions(:);
    end
end

%% load data with specified parameters 
% (number of subjects, number of timepoints, sampling rate, and (if applicable) random subset of parcels

data_temp = load([tc_dir '/hcp1003_REST1_LR_' options.parcellation '.mat']);

if ~isfield(options, 'nsubs')
    options.nsubs = 100;
    warning('Setting number of subjects to default (100)');
end
if ~isfield(options, 'nts')
    options.nts = 'all';
    warning('Setting number of timepoints to default (all)');
end
if ~isfield(options, 'sr')
    options.sr = 1;
    warning('Setting sampling rate to default (1)');
end

data = cell(options.nsubs,1);

for s = 1:options.nsubs
    if strcmp(options.nts, 'all')
        data{s,1} = data_temp.data{s,1}(1:options.sr:end,these_regions);
    else
        data{s,1} = data_temp.data{s,1}(1:options.sr:options.nts,these_regions);
    end
end
clear data_temp

%% set up other variables necessary for HMM estimation

ts = size(data{1},1);

T = cell(1, options.nsubs);
for i=1:options.nsubs
    T{1,i} = ts; % change T depending on length of data
end 

if ~isfield(options, 'k')
    options.k = 12;
    warning('Setting k to default (12)');
end

clear hmm_options
hmm_options.order = 0;
hmm_options.covtype = 'full'; %('full' for covariance, 'uniquefull' for no covariance)
hmm_options.zeromean = 1; % (0 to model mean, 1 to model only covariance)
hmm_options.standardise = 1;
hmm_options.dropstates = 0;
hmm_options.K = options.k;
hmm_options.useParallel = 0;
% if applicable, reduce dimensionality using either PCA or HMM-PCA
% approach
if isfield(options, 'pca')
    hmm_options.pca = options.pca;
end
if isfield(options, 'lowrank')
    hmm_options.lowrank = options.lowrank;
end

%% run HMM and compute max. FOs

[HMM.hmm, HMM.Gamma, ~, HMM.vpath, ~, ~, HMM.fehist] = hmmmar(data, T, hmm_options);
HMM.maxFO = getMaxFractionalOccupancy(HMM.Gamma, T, hmm_options);
HMM.FO = getFractionalOccupancy(HMM.Gamma, T, hmm_options);

HMM.Gamma = single(HMM.Gamma);
HMM.vpath = single(HMM.vpath);
HMM.these_regions = these_regions;

%% save output

if ~isdir(outputdir); mkdir(outputdir); end
if strcmp(options.nts, 'all')
    str1 = options.nts;
else
    str1 = num2str(options.nts);
end
if strcmp(options.nregions, 'all')
    str2 = options.nregions;
else
    str2 = num2str(options.nregions);
end
if isfield(options, 'pca')
    str3 = ['_pcadim' num2str(options.pca) '.mat'];
elseif isfield(options, 'lowrank')
    str3 = ['_lowrankdim' num2str(options.lowrank) '.mat'];
else
    str3 = '.mat';
end
    
save([outputdir '/HMM_k' num2str(options.k) '_' options.parcellation '_nsubs' ...
    num2str(options.nsubs) '_nts' str1 '_sampling' num2str(options.sr) ...
    '_' str2 str3], 'HMM', 'options');


end

