%function LRB_SVM(subNum,events,fname)
%% Description
% RSA_LDt(subNum,events,analysisType,fname)
% Use Linear discriminant t to conduct RDM
clc
clear
%% Initialization
addpath('/home/kinkini/eeglab/'); eeglab; close; % obtain eeglab functions
%addpath '/home/kinkini/Downloads/Neuroscience/Lawrence/2021_LRB'
fname = 'LRB_ClassifyFeatures_ERP_Preview'
subNum = 1:19
events = 1
%% Paths and Variables
% Get participants
for i = 1:subNum(end)
    subjects{i} = ['p', num2str(subNum(i))]
end

subjects(12) = [];
% ERP processing options
windsorize = 0; normalize = 0; mnn = 0; %mahal = 1;
% paths
path = '/home/kinkini/Downloads/Neuroscience/Lawrence/EEG_2019_2021_2024/study2_guo2021/Data/';
ERPpath = '/home/kinkini/Downloads/Neuroscience/Lawrence/EEG_2019_2021_2024/study2_guo2021/Analyzed/';
% Electrode selection
%Electrodes = [1:3, 33:37]; % Frontal
electrodeSel = 0;
%% Analyais Loop
% select event alignment
for ievent = 1:length(events)
    clear accurate
    eventalign = events(ievent);
    accurate = [];
    if eventalign == 1
        maxTime = 500;
    elseif eventalign == 2
        maxTime = 800;
    end
    % select subjects
    for iSub = 1:length(subjects)
        disp(iSub)
        disp(['Subject Number: ',num2str(iSub)]);
        disp('Getting ERPs...');
        tic
        % Get nBlocks
        load([path,subjects{iSub}, '/behavData.mat']);
        % Get completed Block Number
        [~, c] = find(~isnan(theData.RT)==1);
        nBlocks = max(c);
        eventalign = 1
        %Get ERPs
        %try
            if eventalign == 1
                load([ERPpath,'Preview/ERPs/',subjects{iSub},'_iTrialMinusOne.mat']);
            elseif eventalign == 2
                load([ERPpath,'Go/ERPs/',subjects{iSub},'_iTrialMinusOne.mat']);
            elseif eventalign == 3
                load([ERPpath,'MoveOn/ERPs/',subjects{iSub},'_iTrialMinusOne.mat']);
            end
        %catch
            %[L, R, B, MyInfo] = LRB_GetERPs_iTrialMinusOne_Obj(subjects{iSub}, nBlocks, eventalign,0,windsorize,normalize,mnn);
            %[L, R, B, MyInfo] = LRB_GetERPs_iTrialMinusOne_Obj(subjects{iSub}, nBlocks, eventalign,0,windsorize,normalize,mnn);
             
            %if eventalign == 1
            %    save([ERPpath,'Preview/ERPs/',subjects{iSub},'_iTrialMinusOne.mat'], 'L','R','B','MyInfo');
            %elseif eventalign == 2
            %    save([ERPpath,'Go/ERPs/',subjects{iSub},'_iTrialMinusOne.mat'], 'L','R','B','MyInfo');
            %elseif eventalign == 3
            %    save([ERPpath,'MoveOn/ERPs/',subjects{iSub},'_iTrialMinusOne.mat'], 'L','R','B','MyInfo');
            %end
        %end
        toc
        %% SVM Classification
        tic
        % Shape
        nElectrodes = size(L.cw.ob1,2); nTimes = size(L.cw.ob1,3);
        
         s1 = cat(4,L.cw.ob1, L.cw.ob2, L.ccw.ob1, L.ccw.ob2, ...
             R.cw.ob1, R.cw.ob2, R.ccw.ob1, R.ccw.ob2, ...
             B.cw.ob1, B.cw.ob2, B.ccw.ob1, B.ccw.ob2);
         s2 = cat(4,L.cw.ob3, L.cw.ob4, L.ccw.ob3, L.ccw.ob4, ...
             R.cw.ob3, R.cw.ob4, R.ccw.ob3, R.ccw.ob4, ...
             B.cw.ob3, B.cw.ob4, B.ccw.ob3, B.ccw.ob4);
         s1 = squeeze(nanmean(s1,4));
         s2 = squeeze(nanmean(s2,4));
         sData = cat(1,s1,s2);
         
         accurate.shape(iSub,:) = ClassifySVM(sData,sData);
%         
%         % Orientation
        o1 = cat(4,L.cw.ob1,L.cw.ob2,L.cw.ob3,L.cw.ob4, ...
             R.cw.ob1,R.cw.ob2,R.cw.ob3,R.cw.ob4, ...
             B.cw.ob1,B.cw.ob2,B.cw.ob3,B.cw.ob4);
         o2 = cat(4,L.ccw.ob1,L.ccw.ob2,L.ccw.ob3,L.ccw.ob4, ...
             R.ccw.ob1,R.ccw.ob2,R.ccw.ob3,R.ccw.ob4, ...
             B.ccw.ob1,B.ccw.ob2,B.ccw.ob3,B.ccw.ob4);
         o1 = squeeze(nanmean(o1,4));
         o2 = squeeze(nanmean(o2,4));
         oData = cat(1,o1,o2);
         
         accurate.orientation(iSub,:) = ClassifySVM(oData,oData);
        % Hand
        
        L = cat(4, L.cw.ob1, L.cw.ob2, L.cw.ob3, L.cw.ob4, ...
            L.ccw.ob1, L.ccw.ob2, L.ccw.ob3, L.ccw.ob4);
        R = cat(4, R.cw.ob1, R.cw.ob2, R.cw.ob3, R.cw.ob4, ...
           R.ccw.ob1, R.ccw.ob2, R.ccw.ob3, R.ccw.ob4);
        B = cat(4, B.cw.ob1, B.cw.ob2, B.cw.ob3, B.cw.ob4, ...
            B.ccw.ob1, B.ccw.ob2, B.ccw.ob3, B.ccw.ob4);
        L = squeeze(nanmean(L,4));
        R = squeeze(nanmean(R,4));
        B = squeeze(nanmean(B,4));
        
        LR = cat(1,L,R); LB = cat(1,L,B); RB = cat(1,R,B);
        
        accurate.LR(iSub,:) = ClassifySVM(LR,LR);
        accurate.LB(iSub,:) = ClassifySVM(LB,LB);
        accurate.RB(iSub,:) = ClassifySVM(RB,RB);
        accurate.hand(iSub,:) = (accurate.LR(iSub,:) + accurate.LB(iSub,:) + accurate.RB(iSub,:))./3;
        toc
        
        %% save results
        if eventalign == 1
            savepath = '/home/kinkini/Downloads/Neuroscience/Lawrence/2021_LRB/2021_LRB_Analyzed/Preview/';
        elseif eventalign == 2
            savepath = '/psyhome8/guolin1/Experiments/LRBimanual/Analyzed/Go/';
        elseif eventalign == 3
            savepath = '/psyhome8/guolin1/Experiments/LRBimanual/Analyzed/MoveOn/';
        end
        
        cd(savepath);
        save(fname, 'accurate', 'MyInfo', '-v7.3');
    end % subs
end % event
%end
