%% This is the ERP function embedded above
%function [cw, ccw, kna, knb, MyInfo] = NLSM_GetGraspERPs(partNum, nBlocks, eventalign, windsorize, normalize)
%partNum = 1;
%can turn this into a loop once we have worked out all the details with the
%behav files
nBlocks =20;
windsorize = 5
normalize = true
selectedtrials = 6 %the number of trials you wish to use to create the erps

for partNum = 1:10 %add as many as you have
%% Description:
% Get ERP from a participant, windsorize,
% & normalize, then save them into a ERP file
% Do not normalize and windsorize if using for RSA
% Naming convention is ERP_mnn_nonorm.mat
% whichCond: 
%   1: Pointy; Steel
%   2: Pointy; Wood
%   3: Round; Steel
%   4: Round; Wood
% blocks
%   1: CW;  2: CCW;  3: KnA (CW knuckling); 4: KnB (CCW knuckling);
eeglab;

partstrg = ['p',num2str(partNum)]
datapath = ['/home/kinkini/Downloads/Neuroscience/leenina1/Experiments/SurfaceMaterials/',partstrg,'/MoveOn/iclabel_Apr25/'];
savepath = ['/home/kinkini/Downloads/Neuroscience/Matthias/Leenina/Data/',partstrg,'/']
%% Paths & Variables
%if eventalign == 1
%    datapath = '/home/kinkini/Downloads/Neuroscience/leenina1/Experiments/SurfaceMaterials/',num2str(partNum),'/postICA_seg';
%    savepath = '/home/kinkini/Downloads/Neuroscience/Matthias/NimroV/Jihoon/Experiments/ClassificationTask/EEG/ERPs/';
%elseif eventalign == 2
%    datapath = '/home/kinkini/Downloads/Neuroscience/Matthias/NimroV/Jihoon/Experiments/ClassificationTask/EEG/processed/postICA/';
%    savepath = '/home/kinkini/Downloads/Neuroscience/Matthias/NimroV/Jihoon/Experiments/ClassificationTask/EEG/ERPs/';
%end

%% Get some variables
EEGpath = [datapath];
load(['/home/kinkini/Downloads/Neuroscience/leenina1/behavdatafiles/',['behavData', num2str(partNum)],'.mat']);
MyInfo.RT = theData.RT; MyInfo.MT = theData.MT; MyInfo.blocks = theData.blocks;
MyInfo.whichCond = theData.whichCond; MyInfo.badTrials = zeros(nBlocks,1);
MyInfo.badTrials(isnan(MyInfo.badTrials))=0;


%% Loop through files
count_cw = 0; count_ccw = 0; count_kna = 0; count_knb = 0;

for iBlock = 1:nBlocks
    if iBlock < 13
        filename = ['s1b',num2str(iBlock),'.mat'];
    else
        filename = ['s2b',num2str(iBlock-12),'.mat'];
    end
    %filename = ['b',num2str(iBlock),'.mat'];
    load([EEGpath, filename]);
    MyInfo.badTrials = theData.badTrials(:,iBlock);
    MyInfo.badTrials(isnan(MyInfo.badTrials))=0;
    bad_trial_indices = find(MyInfo.badTrials);
    % Remove entries where MyInfo.badTrials is 1
    % Logical index for keeping only the good trials (0s in MyInfo.badTrials)
    goodIndices = (MyInfo.badTrials == 0);
    
    % Filter MyInfo.whichCond to keep only the good trials
    Info.whichCond = Info.whichCond(goodIndices,:);

    % Obtain grasp types
    blockType = MyInfo.blocks(iBlock);
    
    % record timepoints info
    MyInfo.timepoints = EEG.times;
    
    % interpolate any missing electrodes if they're not already done
    EEG = pop_interp(EEG, Info.chanlocs, 'spherical');
    %nElectrodes = size(EEG.data,1); nTimes = size(EEG.data,2);
    nEpochs = size(EEG.data,3);
	% check for mismatch in number of trials in this block
    
    if nEpochs ~= length(Info.whichCond)
		disp(['Trial Number mismatch at block number ', num2str(iBlock), '. Subject Number ', partNum]);
    end
 
    % perform mnn
        clear allcov normdata alldata ndata;
        alldata = permute(EEG.data,[3 1 2]); %concatenate all observations
        
        for i = 1:size(alldata,3)%parfor possible here
            allcov(:,:,i) = covCor(alldata(:,:,i)); % calculate covariance matrix at each time point
        end
        sig = mean(allcov,3); % average covariance matrix across time
        for i = 1:size(alldata,1)
            for j = 1:size(alldata,3) %parfor possible here
                idata = squeeze(alldata(i,:,j));
                ndata(:,j) = idata * (sig^(-1/2)); % normalize ERP at each time point
            end
            normdata(i,:,:) = ndata;
        end
        EEG.data = permute(normdata,[2 3 1]);

        for condition = 1:4;

        % Logical array for condition 1
        conditionLogical = (Info.whichCond == condition);
        
        % Find the indices of trials where condition is 1
        conditionIndices = find(conditionLogical);
        
        % Keep only the first three true values
        conditionLogical(conditionIndices(selectedtrials+1:end)) = false;

    % Load into block-averaged variables (observations)
    if blockType == 1
            count_cw = count_cw+1;
            fieldname = ['ob',num2str(condition)];
            cw.fieldname(count_cw,:,:)= nanmean(EEG.data(:,:,conditionLogical),3);
    elseif blockType == 2
            count_ccw = count_ccw+1;
            fieldname = ['ob',num2str(condition)];
            ccw.fieldname(count_ccw,:,:)= nanmean(EEG.data(:,:,conditionLogical),3);
    elseif blockType == 3
            count_kna = count_kna+1;
            fieldname = ['ob',num2str(condition)];   
            kna.fieldname(count_kna,:,:)= nanmean(EEG.data(:,:,conditionLogical),3);
    elseif blockType == 4
            count_knb = count_knb+1;
            fieldname = ['ob',num2str(condition)];
            kna.fieldname(count_knb,:,:)= nanmean(EEG.data(:,:,conditionLogical),3);
    end
        end
end
%% windsorize & normalize
for iVar = 1:4
    if iVar == 1
        ERPdata = cw;
    elseif iVar == 2
        ERPdata = ccw;
    elseif iVar == 3
        ERPdata = kna;
    elseif iVar == 4
        ERPdata = knb;
    end
    fields=fieldnames(ERPdata);
    for categoryidx=1:length(fields)
        categoryname=fields{categoryidx};
            data = ERPdata.(categoryname);
            
            % Standardize, windsorize, & Normalize ERPs
            if normalize
                for j = 1:size(data,1)
                    temp = data(j,:,:);
                    ztemp = (temp - mean(temp(:))) ./ std(temp(:));
                    if windsorize
                        ztemp(ztemp>3)=3; ztemp(ztemp<-3)=-3; % Windsorize outliers
                    end
                    nztemp = (ztemp - min(ztemp(:))) ./ (max(ztemp(:)) - min(ztemp(:)));
                    data(j,:,:) = nztemp;
                end
            end
            ERPdata.(categoryname) = data;
    end
    if iVar == 1
        cw = ERPdata;
    elseif iVar == 2
        ccw = ERPdata;
    elseif iVar == 3
        kna = ERPdata;
    elseif iVar == 4
        knb = ERPdata;
    end
end

%% save ERPs
try
    cd([savepath,num2str(partNum)]);
catch
    cd(savepath);
    mkdir(num2str(partNum));
    cd(partNum);
end
save ERPs_mnn_Aug27.mat -mat cw ccw kna knb MyInfo

%end
end