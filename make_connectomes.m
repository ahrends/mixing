function [C_concat, groupmean] = make_connectomes(TC, settings_netmat, condNames)
% Function to compute FC networks using FSLnets
%   INPUTS:
%           TC: cell array with fMRI timecourses {conditions x subjects};
%               TC{x,y}: [nROIs x nTimepoints]
%           settings_netmat (optional): struct with fields setting options for netmats
%               .partialCorrMethod: if empty (default) or non-existing no partial correlation
%                                   will be computed
%           condNames (optional): cell including string names of conditions
%                      corresponding to the order in TC
%   OUTPUTS:
%           Cnetmats: cell array with static FC networks {nConds}.FCmeasure(nROI x nROI x nSubs)
%
% Angus Stevner (edits Christine Ahrends)

%% get static FC matrices

if nargin < 2 || isempty(settings_netmat)
    settings_netmat = struct;
end

if nargin > 2
    if ~isempty(condNames)
        condLabels = condNames;
    else
        for ii = 1:size(TC,1)
            condLabels{ii} = ['condition ' num2str(ii)];
        end
    end
else
    for ii = 1:size(TC,1)
        condLabels{ii} = ['condition ' num2str(ii)];
    end
end
nConds = size(TC,1);

Cnetmats = cell(1,nConds);

for iConds = 1:nConds
    % count subs in condition
    condSubs = [];
    countSubs = 0;
    for iCount = 1:size(TC,2)
        if ~isempty(TC{iConds,iCount})
            countSubs = countSubs + 1;
            condSubs(countSubs) = iCount;
        end
    end
    nSubs = countSubs;
    for iSubs = 1:nSubs
        nNodes = size(TC{iConds,iSubs},1);
        Cnetmats{iConds}.fullCorr(:,:,iSubs) = reshape(nets_netmats(squeeze(TC{iConds,iSubs})',0,'corr'),nNodes,nNodes);
        Cnetmats{iConds}.fullCorr_z(:,:,iSubs) = reshape(nets_netmats(squeeze(TC{iConds,iSubs})',1,'corr'),nNodes,nNodes);
        Cnetmats{iConds}.Cov(:,:,iSubs) = reshape(nets_netmats(squeeze(TC{iConds,iSubs})',0,'cov'),nNodes,nNodes);
        if isfield(settings_netmat,'partialCorrMethod')
            Cnetmats{iConds}.(genvarname(settings_netmat.partialCorrMethod))(:,:,iSubs) = ...
                reshape(nets_netmats(squeeze(TC{iConds,iSubs})',1,settings_netmat.partialCorrMethod),nNodes,nNodes);
        end
    end
    
end

%% investigate static FC

clear C_tmp C_tmptmp C_tmptmptmp C_tmptmptmptmp

FCmeasures = {'fullCorr_z'};
FClabels = {'z-transformed corr'};
condLabels = condNames;
ZnetKeep = [];

for iMeasure = 1:length(FCmeasures)
    for iConds = 1:length(Cnetmats)
        Znet = [];
        C_tmp = []; C_tmp = permute(Cnetmats{iConds}.(genvarname(FCmeasures{iMeasure})),[3,1,2]);
        C_tmptmp = reshape(C_tmp,size(C_tmp,1),size(C_tmp,2)^2);
        Znet = nets_groupmean(C_tmptmp,0);
        nNodes = size(C_tmp,2);
        Znet = reshape(Znet,[nNodes nNodes]);
        C_concat = zeros(size(C_tmp,1),length(find(~tril(ones(size(C_tmp,2))))));
        for ii = 1:size(C_tmp,1)
            C_tmptmptmptmp = []; C_tmptmptmptmp = squeeze(C_tmp(ii,:,:));
            C_concat(ii,:) = C_tmptmptmptmp(~tril(ones(size(C_tmp,2))));
        end
        groupmean = [];
        groupmean = repmat(mean(C_concat,1),size(C_concat,1),1);
    end
end

end