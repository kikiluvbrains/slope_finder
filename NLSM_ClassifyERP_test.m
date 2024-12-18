function NLSM_ClassifyERP_test_mani
%% Description
% Decode shape and "colour"(material) of objects across time

%% Define analysis
noerps = 6
subNum = 1:14; %should be 16 but you forgot to transfer some data for 15 and 16
events = 1;%1=align to Preview Onset; 3=align to Movement Onset
eeglab;
%% Paths and Variables
% Get participants
nSubs = length(subNum);
for i = 1:nSubs
    subjects{i} = ['p', num2str(subNum(i))];
end

% Set temporal window
nCat = 7; windsorize = 1; normalize = 1;
%% Analysis Loop
% select event alignment
for ievent = 1:length(events)
    eventalign = events(ievent);
    accurate = [];
    % select subjects
for iSub = 1:length(subjects)
        disp(['Subject Number: ',num2str(iSub)]);
        disp('Getting ERPs...');
        nBlocks = 24;
        
        %% Get ERPs
            if eventalign == 1
                load(['/home/kinkini/Downloads/Neuroscience/Matthias/Leenina/Data/',subjects{iSub},'/ERPs_mnn_Aug27.mat']);
            end
        %% Concatenate data
        disp('concatenating data...');
        tic
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

           ERPdata.ob1 = ERPdata.ob1(1:noerps,:,:)
           ERPdata.ob2 = ERPdata.ob2(1:noerps,:,:)
           ERPdata.ob3 = ERPdata.ob3(1:noerps,:,:)
           ERPdata.ob4 = ERPdata.ob4(1:noerps,:,:)
           fields=fieldnames(ERPdata);
           for categoryidx=1:length(fields)
               categoryname=fields{categoryidx};
                   data = ERPdata.(categoryname);
                    % concatenate according to nCat
                    tempdata = []; newdata = [];
                    for i = 1:size(data,3)
                        tempdata = cat(2, tempdata, data(:,:,i));
                        if mod(i,nCat) == 0
                            newdata = cat(3, newdata, tempdata);
                            tempdata = [];
                        end
                    end
                   ERPdata.(categoryname)=newdata;
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
        toc
        %% Classification
        disp('Classifying...');
        tic    

        grasp_shape_data = cat(1, cw.ob1, cw.ob2, ccw.ob1, ccw.ob2, cw.ob3, cw.ob4, ccw.ob3, ccw.ob4);
        knuckle_shape_data = cat(1, kna.ob1, kna.ob2, knb.ob1, knb.ob2, kna.ob3, kna.ob4, knb.ob3, knb.ob4);
        
        grasp_colour_data = cat(1, cw.ob1, cw.ob3, ccw.ob1, ccw.ob3, cw.ob2, cw.ob4, ccw.ob2, ccw.ob4);
        knuckle_colour_data = cat(1, kna.ob1, kna.ob3, knb.ob1, knb.ob3, kna.ob2, kna.ob4, knb.ob2, knb.ob4);
        
        % Classify and store results
        accurate.Grasp.Shape(iSub,:) = ClassifyERP(grasp_shape_data, grasp_shape_data);
        accurate.Knuckle.Shape(iSub,:) = ClassifyERP(knuckle_shape_data, knuckle_shape_data);
        
        accurate.Grasp.Colour(iSub,:) = ClassifyERP(grasp_colour_data, grasp_colour_data);
        accurate.Knuckle.Colour(iSub,:) = ClassifyERP(knuckle_colour_data, knuckle_colour_data);

        toc        
        %% save results
        if eventalign == 1
            savepath = '/home/kinkini/Downloads/Neuroscience/Matthias/Leenina/Analyzed';
        end
        
        cd(savepath);
        save NLSM_ClassifyFeaturesAction_ERP_7.mat -mat accurate MyInfo
end % subs
end % event
end

function [acc] = ClassifyERP(ERP1, ERP2)
%% Description
% Classifies binary categories in ERP. ERP should be organized in
% nObservations x nFeatures x nTimepoints
% Observations should be organized such that across nColumns,
% ERP(1:nColumns/2, :, :) should be observations of category 1, and
% ERP(nColumns/2+1, :, :) should be observations of category 2. ERP1 and
% ERP2 should have the same structure, and so cross-decoding is possible if
% ERP1 and ERP2 are different data. Otherwise, ERP1 and ERP2 should be the
% same data.
 
nFold = size(ERP1,1)/2;
grouplabels = [ones(nFold,1); ones(nFold,1).*(-1)]; % changed from 1 & 2 but makes no difference
nTimes = size(ERP1, 3);
cvlabels = [1:nFold, 1:nFold];

% Train & Classify at each timepoint
for itime = 1:nTimes
    for icv = 1:nFold % can change for to parfor for parallel processing
        trainfold = ERP1(cvlabels~=icv, :, itime);
        testfold = ERP2(cvlabels==icv, :, itime);
        trainlabel = grouplabels(cvlabels~=icv);
        testlabel = grouplabels(cvlabels==icv);        
        model = svmtrain(trainlabel, double(trainfold), sprintf('-q -t 0 -c %f', 1));
        [predictedlabels] = svmpredict(testlabel, double(testfold), model, '-q');
        iAcc(icv) = sum(predictedlabels(:) == testlabel(:)) / length(predictedlabels);
    end
    acc(itime) = mean(iAcc);
end
end
