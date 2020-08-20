function [icc_summary,var,stats,sigmask] = run_reliability(varargin)
% does sig masking and computes 1, 2- or 3-factor ICC (detects from nfactors in factors) 
% Multiple factor (2 or 3) ICC calculated via G-Theory (Webb and Shavelson, 2005)
%
% Recommended usage: [icc,var,stats,sigmask]=run_reliability('data',data,'factors',ftbl 
%
% Required input: data and ftbl
% Optional input: datadir, procedure, correctiontype
%   data : cell array (length = number of scans) of 1-, 2-, or 3-D matrices (any type, e.g., raw data, edges, masked data &c)
%   datadir and procedure : used together to load data if not already given directly
%   correctiontype : NOT RECOMMENDED FOR NORMAL USE - used to mask variables that are not significantly different than 0 across the group


%% Parse input
p = inputParser;

defaultdata{1}=0; % no default
defaultfactors=0; % no default
defaultprocedure='NA'; % default=NA, can use e.g., {'PCC', 'RMC',...} to load from directory
defaultdatadir='~/'; % automatically adds results_%procedure%
defaultcorrectiontype='none'; % multiple comparison corr for sigmask options: {'none','Bonf','FDR'}

addParameter(p,'data',defaultdata,@iscell);
addParameter(p,'factors',defaultfactors,@isnumeric);
addParameter(p,'procedure',defaultprocedure,@isnumeric);
addParameter(p,'datadir',defaultdatadir,@ischar);
addParameter(p,'correctiontype',defaultcorrectiontype,@ischar);

parse(p,varargin{:});

data = p.Results.data;
ftbl = p.Results.factors;
procedure = p.Results.procedure;
datadir=p.Results.datadir;
correctiontype=p.Results.correctiontype;

clearvars p varargin

%% setup data

scripts_path=fileparts(which(mfilename)); 
addpath(genpath(scripts_path));

if ~strcmp(procedure,'NA')
    warning('Loading in via file here might have to be updated.')
    % TODO: either manually replace or use pwd (after cd to datadir)
    datadir=strcat(datadir,sprintf('results_%s/',procedure));
    [data,ftbl]=load_reliability_data(datadir,'*nii*');     % TODO: update input to load_reliability_data based on new version

else
    if sum(sum(ftbl))==0 | sum(sum(data{1}))==0;
        error('Data unspecified. Either pass data/ftbl or specify procedure.')
    end
    procedure='tmp';
    
    % check for missing data (per sub)
    for i=1:size(ftbl,1);
        t(i)=sum(ftbl(:,1)==i);
    end
    if length(unique(t))>2  % 2 bc unique includes "0"
        disp('Runs per subject:')
        disp(sprintf('%d ',unique(t)))
        error('Missing data. Please remove partial data.')
    end
end

%% main

% determine whether doing voxelwise image or matrix

if(ndims(data{1})==3)
    doing_voxelwise=1;
elseif (ndims(data{1})==2)
    doing_voxelwise=0;
else 
    error('weird dimensions')
end

if(doing_voxelwise) % mask voxels that are zero across subjects (non-brain)
%    gmmask = load_untouch_nii('5thpass_template_GM_WM_CSF_resl_crop_resampled.nii.gz');
%     (if mask is in uint8, do: im2single(gmmask) )
%    gmmask=gmmask.img;
%    gmmask(gmmask<3)=0;
%    gmmask(gmmask==3)=1;

	gmmask=ones(size(data{1})); % no a priori mask
    [masked_data,gmmask_composite]=get_nonzero_voxels_across_group(data,gmmask);  % get GM; remove zero entries
else
    s=size(data{1},1);
    % get lower triangle if matrix and not just 1 entry per sub
    if size(data{1},2)==s && size(data{1},2)>2 
        for i=1:length(data)
            masked_data{i}=data{i}(logical(tril(ones(s,s),-1)));
        end
    else % single entry per subject
        initmask=ones(size(data{1}));
        masked_data=data;
        %[masked_data,zeromask]=get_masked_data_special_trav(data,initmask); % may contain redundant edges - dof
%       masked_data=data; % SMN
    end
    
end

sigmask = create_sigedge_mask(masked_data,0.05,correctiontype); % "none"->all ones
data_sig = get_masked_data(masked_data,sigmask);
[icc_summary,var,stats]=calc_roi_iccs(data_sig,ftbl,'all');

% if using no sigmask, still return a real sigmask in case it's needed for future analyses
if strcmp(correctiontype,'none')
    sigmask = create_sigedge_mask(masked_data,0.05,'Bonf');
end


%% save data
clockinfo{1}=datestr(now, 'HH'); clockinfo{2}=datestr(now, 'MM'); clockinfo{3}=datestr(now, 'SS');


if(doing_voxelwise) % save mask bc voxelwise
%     save(sprintf('%s_%s.mat',procedure,correctiontype),'sigmask','stats','icc_summary','var','gmmask_composite')
    save(sprintf('%s_icc_%sH%sM%sS.mat',procedure,clockinfo{1},clockinfo{2},clockinfo{3}),'sigmask','stats','icc_summary','var','gmmask_composite')
else
%     save(sprintf('%s_%s.mat',procedure,correctiontype),'sigmask','stats','icc_summary','var')
    save(sprintf('%s_icc_%sH%sM%sS.mat',procedure,clockinfo{1},clockinfo{2},clockinfo{3}),'sigmask','stats','icc_summary','var')
end

