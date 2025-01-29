function NLSM_ClassifyERP_test_mani
%% Description
% Decode shape and "colour"(material) of objects across time

%% Define analysis
subNum = 1:14; %should be 16 but you forgot to transfer some data for 15 and 16
events = 1;%1=align to Preview Onset; 3=align to Movement Onset
eeglab;
%% Paths and Variables
% Get participants
nSubs = length(subNum);
for i = 1:nSubs
    if i ~= 0
        subjects{i} = ['p', num2str(subNum(i))];
    end
end
noerps = 3
% Set temporal window
nCat = 5; windsorize = 1; normalize = 1; sloper = 0;
%% Analysis Loop
% select event alignment
for ievent = 1:length(events)
    eventalign = events(ievent);
    accurate = [];
    % select subjects
for iSub = 1:length(subjects)
    if iSub ~= 0       
        disp(['Subject Number: ',num2str(iSub)]);
        disp('Getting ERPs...');
        nBlocks = 24;
        
        %% Get ERPs
        if eventalign == 1
            load(['/Users/Owner/Downloads/Niemier - SVM analysis/Matlab_Nina/data/premade_erps/',subjects{iSub},'/ERPs_mnn_Aug27.mat']);
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
           fields=fieldnames(ERPdata);
           ERPdata.ob1 = ERPdata.ob1(1:noerps,:,:)
           ERPdata.ob2 = ERPdata.ob2(1:noerps,:,:)
           ERPdata.ob3 = ERPdata.ob3(1:noerps,:,:)
           ERPdata.ob4 = ERPdata.ob4(1:noerps,:,:)

           for categoryidx = 1:length(fields)
                categoryname = fields{categoryidx};
                data = ERPdata.(categoryname); % Original data: (trials x channels x timepoints)
                tempdata = []; newdata = [];
                if sloper ~= 1
                    for i = 1:size(data,3)
                        tempdata = cat(2, tempdata, data(:,:,i));
                        if mod(i,nCat) == 0
                            newdata = cat(3, newdata, tempdata);
                            tempdata = [];
                        end
                    end
                else
                    % Loop over each trial
                    for iTrial = 1:size(data, 1)
                        trialData = squeeze(data(iTrial, :, :)); % Extract data: (channels x timepoints)
                        
                        % Calculate the number of segments based on the timepoints
                        numSegments = floor(size(trialData, 2) / nCat); % Number of non-overlapping segments
                        
                        % Initialize transformed data (features x segments)
                        numChannels = size(trialData, 1);
                        numFeaturesPerChannel = 3; % 3 coefficients for a quadratic polynomial
                        transformedData = zeros(numChannels * numFeaturesPerChannel, numSegments);
                        
                        % Loop over each channel
                        for iChan = 1:numChannels
                            channelData = trialData(iChan, :); % Extract 1D array: (1 x timepoints)
                            featureMatrix = []; % Store polynomial coefficients for this channel
                    
                            % Loop over non-overlapping groups of 5 timepoints
                            for iTime = 1:nCat:length(channelData) - nCat + 1
                                % Get the current window of timepoints (5 points)
                                timeWindow = channelData(iTime:iTime + nCat - 1); % Window of 5 points
                    
                                % Fit a quadratic polynomial (degree = 2)
                                x = 1:nCat; % Time indices for the window
                                p = polyfit(x, timeWindow, 2); % Polynomial coefficients (3 values)
                    
                                % Store the polynomial coefficients (3 values for quadratic)
                                featureMatrix = [featureMatrix, p]; % Append coefficients (3 values)
                            end
                            
                            % Reshape the feature matrix to 3 rows (coefficients) and numSegments columns (windows)
                            featureMatrix = reshape(featureMatrix, numFeaturesPerChannel, numSegments);
                    
                            % Store coefficients for this channel: 3 coefficients per window, numSegments windows
                            transformedData((iChan - 1) * numFeaturesPerChannel + 1:iChan * numFeaturesPerChannel, :) = featureMatrix;
                        end
                    
                        % Store the transformed trial data
                        newdata(iTrial, :, :) = transformedData; % Shape: (1 trial x (numChannels * numFeaturesPerChannel) features x numSegments)
                    end
                end
                % Update ERPdata structure with the transformed data
                ERPdata.(categoryname) = newdata; % Final shape: (6 trials x 192 features x 133 segments)
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
            savepath = ['/Users/Owner/Downloads/Niemier - SVM analysis/Matlab_Nina/results/time', num2str(nCat),'/'];
        end
        
        cd(savepath);
        if sloper ~= 1
            save NLSM_ClassifyFeaturesAction_ERP_3.mat -mat accurate MyInfo
        else
            save NLSM_ClassifyFeaturesAction_ERP_3_slope.mat -mat accurate MyInfo

        end
    end % subs
end % event
end
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
% 
% Apply slope finder to ERP1 and ERP2

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

function [transformedERP] = applySlope(ERP, degree)
% Apply slope finder transformation to reduce the second dimension
% from 320 channels to 64 feature vectors by fitting a polynomial of
% degree 'degree' and using its parameters as features.

nObs = size(ERP, 1); % Number of observations (subjects)
nChans = size(ERP, 2); % Number of channels
nTime = size(ERP, 3); % Number of timepoints

% Number of groups for downsampling the 320 channels
groupSize = nChans / 64; % Assuming fixed reduction to 64
numFeatures = degree + 1; % Number of polynomial coefficients for each fit

% Initialize transformed ERP matrix
transformedERP = zeros(nObs, 64 * numFeatures, nTime); % Extracting all coefficients

for iObs = 1:nObs
    for iTime = 1:nTime
        featureVector = []; % Store features for this observation and timepoint
        for iNewChan = 1:64
            % Indices of channels in the current group
            startIdx = (iNewChan - 1) * groupSize + 1;
            endIdx = iNewChan * groupSize;
            
            % Data from the current group of channels
            channelData = ERP(iObs, startIdx:endIdx, iTime);
            
            % Use polyfit to compute polynomial coefficients
            x = 1:groupSize; % x-axis for channel indices
            p = polyfit(x, channelData, degree); % Polynomial of degree 'degree'
            
            % Store all polynomial coefficients as features
            featureVector = [featureVector, p]; % Concatenate the coefficients
        end
        transformedERP(iObs, :, iTime) = featureVector; % Store the feature vector
    end
end
end
