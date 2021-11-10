function Results = evaluate_HMM_params(tc_dir, HMM_outputdir, results_dir, options)
% This script evaluates HMMs run on the HCP resting state fMRI dataset with
% varying secondary parameters. Requires that the script run_HMM_params.m has
% been run.
%
% INPUT:
% tc_dir:   Directory where timecourses of data can be found. This script
%           assumes that there is a file called [hcp1003_REST1_LR_' ...
%           parcellation '.mat']. In this file, there should be a separate
%           cell for each subject containing their timecourses in the
%           format timepoints x parcels.
% HMM_outputdir:    Directory where HMM was saved (outputdir of script
%                   run_HMM_params.m)
% results_dir:      Directory where results should be saved
% options:  Structure containing fields:
%           measures: measures to evaluate ('staticFC' to evaluate static
%           FC similarity, 'mixing' to evaluate model stasis, or both)
%           (default both)
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
% Results:  structure containing results to evaluate static FC similarity
%           and/or model stasis:
%           (if options.measures='staticFC')
%           statFC_similarity: static FC similarity between subjects
%           (if options.measures='mixing')
%           mean_maxFO: mean maximum fractional occupancy (measure of model
%           stasis)
%           maxFO: maximum fractional occupancy of each subject/session
%           FO: fractional occupancy of all states in each subject/session
%           (if options.pca was specified)
%           nPCs: number of principal components used in PCA
%
%
% Christine Ahrends
% (Aarhus University 2020)
%
%

%% set up necessary variables if not specified

if nargin < 2
    warning('HMM output directory not specified, setting to default')
    HMM_outputdir = 'HMM_mixing_test';
end
if nargin < 3
    results_dir = 'Results_mixing_test';
end
if nargin < 4 || isempty(options)
    warning('Options not provided, setting to default')
    options = struct();
    options.parcellation = 'groupICA50';
    options.k = 12;
    options.nsubs = 100;
    options.nts = 'all';
    options.sr = 1;
    options.nregions = 'all';
    options.measures = ["staticFC", "mixing"];
end

%% specify where results are/will be saved

if ~isdir(results_dir); mkdir(results_dir); end

if ~isfield(options, 'k'); warning('Setting k to default (12)'); options.k = 12; end
if ~isfield(options, 'parcellation'); warning('Setting parcellation to default (groupICA50)'); options.parcellation = 'groupICA50'; end
if ~isfield(options, 'nsubs'); warning('Setting number of subjects to default (100)'); options.nsubs = 100; end
if ~isfield(options, 'nts'); warning('Setting number of timepoints to default (all)'); options.nts = 'all'; end
if ~isfield(options, 'sr'); warning('Setting sampling rate to default (1)'); options.sr = 1; end
if ~isfield(options, 'nregions'); warning('Setting number of regions to default (all)'); options.nregions = 'all'; end
    
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
    
    results_file = [results_dir '/Results_k' num2str(options.k) '_' options.parcellation '_nsubs' ...
        num2str(options.nsubs) '_nts' str1 '_sampling' num2str(options.sr) ...
        '_' str2 str3];
    
    %% check if results have already been computed for this combination of parameters
    if exist(results_file, 'file')
        d = dir(results_file);
        filesize = d.bytes;
        if filesize > 0
            try
                load(d.name);
            catch ME
                if strcmp(ME.identifier, 'MATLAB:load:unableToReadMatFile')
                    warning('File is corrupt')
                    new_run = true;
                else
                    disp('File is fine')
                    new_run = false;
                end
            end
        elseif filesize == 0
            warning('File is size 0')
            new_run = true;
        end
    else
        warning('File does not exist')
        new_run = true;
    end
    
    %% load required data and evaluate measures
    if new_run
        if ~isfield(options, 'measures')
            options.measures = ["staticFC", "mixing"];
            warning('Setting evaluative measures to default (static FC similarity and mixing)');
        end
        %% evaluate static FC similarity
        if any(strcmp(options.measures, 'staticFC'))
            % load time courses
            data_temp = load([tc_dir '/hcp1003_REST1_LR_' options.parcellation '.mat']);
            % load HMM
            load([HMM_outputdir '/HMM_k' num2str(options.k) '_' options.parcellation '_nsubs' ...
                num2str(options.nsubs) '_nts' str1 '_sampling' num2str(options.sr) ...
                '_' str2 str3], 'HMM');
            these_regions = HMM.these_regions;
            % load only time courses from regions that were run (important for the runs with
            % randomly sampled subsets of regions)
            
            data = cell(options.nsubs,1);
            
            for s = 1:options.nsubs
                if strcmp(options.nts, 'all')
                    data{s,1} = data_temp.data{s,1}(1:options.sr:end,these_regions);
                else
                    data{s,1} = data_temp.data{s,1}(1:options.sr:options.nts,these_regions);
                end
            end
            clear data_temp
            
            %% compute static FC with FSLnets and between-subject similarity
            
            % reshape data to fit script to compute static FC matrices
            data_FC = cell(1,options.nsubs);
            
            for i = 1:options.nsubs
                data_FC{1,i} = data{i}(:,:)';
            end
            
            condNames = cell(1,1);
            condNames{1,1} = 'REST1';
            
            [C_concat, groupmean] = make_connectomes(data_FC, [], condNames);
            Results.staticFC_similarity = corr(C_concat(:),groupmean(:));
            
        end
        
        %% evaluate mixing
        if any(strcmp(options.measures, 'mixing'))
            if ~exist('HMM','var') % if HMM has not already been loaded, do it now
                HMM = load([HMM_outputdir '/HMM_k' num2str(options.k) '_' options.parcellation '_nsubs' ...
                    num2str(options.nsubs) '_nts' str1 '_sampling' num2str(options.sr) ...
                    '_' str2 str3], 'HMM');
            end
            
            if ~exist('ttrial', 'var')
                ttrial = size(HMM.Gamma,1)/options.nsubs; % time points per session
            end
            
            if ~exist('T','var')
                T = cell(1, options.nsubs);
                for i=1:options.nsubs
                    T{1,i} = ttrial;
                end
            end
            
            if isfield(options, 'pca')
                Results.nPCs = HMM.hmm.train.pca; % save the number of PCs that explained the desired amount of variance (pcadim)
            end
            
            %% get maximum fractional occupancy (max. FO)
            
            Results.mean_maxFO = nanmean(HMM.maxFO);
            Results.maxFO = HMM.maxFO;
            Results.FO = HMM.FO;
            
        end 
        save(results_file, 'Results');
    end
    
end


