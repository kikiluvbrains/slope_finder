
clc
clear
%% This is the ERP function embedded above
%function [cw, ccw, kna, knb, MyInfo] = NLSM_GetGraspERPs(partNum, nBlocks, eventalign, windsorize, normalize)
%partNum = 1;
%can turn this into a loop once we have worked out all the details with the
%behav files
addpath /home/kinkini/Downloads/Neuroscience/Matthias/Leenina/Scripts

nBlocks =20;
windsorize = 5
normalize = true

for partNum = 1:2 %add as many as you have


eeglab;

partstrg = ['p',num2str(partNum)]
datapath = ['/home/kinkini/Downloads/Neuroscience/leenina1/Experiments/SurfaceMaterials/',partstrg,'/MoveOn/iclabel_RTto1s/'];
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
    MyInfo.badTrials(iBlock) = sum(isnan(Info.whichCond));
    Info.whichCond(isnan(Info.whichCond)) = [];

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
        
        parfor i = 1:size(alldata,3)%parfor possible here
            allcov(:,:,i) = covCor(alldata(:,:,i)); % calculate covariance matrix at each time point
        end
        sig = mean(allcov,3); % average covariance matrix across time
        for i = 1:size(alldata,1)
            parfor j = 1:size(alldata,3) %parfor possible here
                idata = squeeze(alldata(i,:,j));
                ndata(:,j) = idata * (sig^(-1/2)); % normalize ERP at each time point
            end
            normdata(i,:,:) = ndata;
        end
        EEG.data = permute(normdata,[2 3 1]);
     
    % Load into block-averaged variables (observations)
    if blockType == 1
            count_cw = count_cw+1;
            cw.ob1(count_cw,:,:)= nanmean(EEG.data(:,:,(Info.whichCond==1)),3);
            cw.ob2(count_cw,:,:) = nanmean(EEG.data(:,:,(Info.whichCond==2)),3);
            cw.ob3(count_cw,:,:)= nanmean(EEG.data(:,:,(Info.whichCond==3)),3);
            cw.ob4(count_cw,:,:) = nanmean(EEG.data(:,:,(Info.whichCond==4)),3);
    elseif blockType == 2
            count_ccw = count_ccw+1;
            ccw.ob1(count_ccw,:,:)= nanmean(EEG.data(:,:,(Info.whichCond==1)),3);
            ccw.ob2(count_ccw,:,:) = nanmean(EEG.data(:,:,(Info.whichCond==2)),3);
            ccw.ob3(count_ccw,:,:)= nanmean(EEG.data(:,:,(Info.whichCond==3)),3);
            ccw.ob4(count_ccw,:,:) = nanmean(EEG.data(:,:,(Info.whichCond==4)),3);
    elseif blockType == 3
            count_kna = count_kna+1;
            kna.ob1(count_kna,:,:)= nanmean(EEG.data(:,:,(Info.whichCond==1)),3);
            kna.ob2(count_kna,:,:) = nanmean(EEG.data(:,:,(Info.whichCond==2)),3);
            kna.ob3(count_kna,:,:)= nanmean(EEG.data(:,:,(Info.whichCond==3)),3);
            kna.ob4(count_kna,:,:) = nanmean(EEG.data(:,:,(Info.whichCond==4)),3);
    elseif blockType == 4
            count_knb = count_knb+1;
            knb.ob1(count_knb,:,:)= nanmean(EEG.data(:,:,(Info.whichCond==1)),3);
            knb.ob2(count_knb,:,:) = nanmean(EEG.data(:,:,(Info.whichCond==2)),3);
            knb.ob3(count_knb,:,:)= nanmean(EEG.data(:,:,(Info.whichCond==3)),3);
            knb.ob4(count_knb,:,:) = nanmean(EEG.data(:,:,(Info.whichCond==4)),3);
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
save ERPs_mnn_Aug27_OG.mat -mat cw ccw kna knb MyInfo

end

