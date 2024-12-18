clc
clear

 % Open EEGLAB & Import plugins
cd('/home/kinkini/eeglab/'); 
eeglab; 
%addpath /home/kinkini/eeglab/plugins/ERPLAB12.00
%estudio;
%partNum = 5
for partNum = 1:1
    for objectproperty = 4:6 %change to 2:6 please append the remainder as condition from line 190 and onwards
    %set object property to 2 if you want toclassify based on object type
    %column 3: convexity
    %column 4: grasp orientation
    %column 5: colour centre
    %column 6: color grasp point
    %column 7: colour non-grasp point
    %1-6 
    for condition = 1:2
    eventalign = 1
    windsorize = 1
    normalize = 1
    %% Paths & Variables
    if eventalign == 1
        savepath = '/home/kinkini/Downloads/Neuroscience/Matthias/NimroV/Jihoon/Experiments/ClassificationTask/ERPs/';
    end

    directory = '/home/kinkini/Downloads/Neuroscience/Matthias/NimroV/Jihoon/Experiments/ClassificationTask/Behavior/';
    if condition == 1
    participant_files = dir(fullfile(directory, ['p', num2str(partNum), 's*b*_Grasp.mat'])); %if grasp change to that
    else
    participant_files = dir(fullfile(directory, ['p', num2str(partNum), 's*b*_Knuckle.mat'])); %if grasp change to that
    end

    %% Loop through files
    count_object_1 = 0;
    count_object_2 = 0;
    count_object_3 = 0;
    count_object_4 = 0;
    count_object_5 = 0;
    count_object_6 = 0;

    for file_idx = 1: length(participant_files)
        file_name = participant_files(file_idx).name;
        disp(['Processing file: ', file_name]);

        if condition == 1
        matches = regexp(file_name, 'p(\d+)s(\d+)b(\d+)_Grasp.mat', 'tokens');
        else
        matches = regexp(file_name, 'p(\d+)s(\d+)b(\d+)_Knuckle.mat', 'tokens');
        end

        if ~isempty(matches)
            session = str2double(matches{1}{2});
            block = str2double(matches{1}{3});
        end

        %% Get some variables
        EEGpath = '/home/kinkini/Downloads/Neuroscience/Matthias/NimroV/Jihoon/Experiments/ClassificationTask/EEG/processed/Stage2/';
        load([directory,file_name]);
        MyInfo.whichCond = theData.whichObj(block,:)';
        MyInfo.badTrials = theData.badTrials(:,block);
        MyInfo.badTrials(isnan(MyInfo.badTrials))=0;
        bad_trial_indices = find(MyInfo.badTrials);
        matrix = MyInfo.whichCond{1}; 
        % Update the cell array with the modified matrix
        MyInfo.whichCond{1} = matrix;
        % Remove the rows corresponding to bad_trial_indices
        matrix(bad_trial_indices, :) = [];
        % Update the cell array with the modified matrix
        MyInfo.whichCond{1} = matrix;        
        %trial deletion since now we have alignment

        EEG_file_name = [EEGpath,'p',num2str(partNum),'s',num2str(session),'b',num2str(block),'_processed_final.mat']
        try
            load(EEG_file_name);
        catch ME
            fprintf('Error loading file: %s\n', ME.message);
            continue;
        end
        MyInfo.timepoints = EEG.times;

        % perform mnn
        clear allcov normdata alldata ndata;
        alldata = permute(EEG.data,[3 1 2]); %concatenate all observations

        for i = 1:size(alldata,3)
            allcov(:,:,i) = covCor(alldata(:,:,i)); % calculate covariance matrix at each time point
        end
        sig = nanmean(allcov,3);

        for i = 1:size(alldata,1)
                for j = 1:size(alldata,3)
                    idata = squeeze(alldata(i,:,j));
                    ndata(:,j) = idata * (sig^(-1/2)); % normalize ERP at each time point
                end
                normdata(i,:,:) = ndata;
        end
        EEG.data = permute(normdata,[2 3 1]);

https://cdn.discordapp.com/attachments/1300158355859574898/1318374516291801179/NLSM_ClassifyERP_test.m?ex=67621772&is=6760c5f2&hm=0edc64c875cb8ce6783da4bd5cb74588b62213e397bee2a9aad3dc7f3c447662&
        % Get all unique conditions
        all_conditions = unique(MyInfo.whichCond{1}(:, objectproperty)); % Unique condition values
        num_conditions = numel(all_conditions); % Number of unique conditions
        
        % Define the desired total number of ERPs
        desired_erps = 6;
        
        % Initialize the output structure
        object = struct();
        
        % Loop through each unique condition
        for cond_idx = 1:num_conditions
            % Generate the field name dynamically
            field_name = sprintf('ob%d', cond_idx); % e.g., 'ob1', 'ob2', etc.
            
            % Get the current condition value
            current_cond = all_conditions(cond_idx);
            
            % Find trials that match the current condition
            matching_trials = find(MyInfo.whichCond{1}(:, objectproperty) == current_cond);
            
            % Get the EEG data corresponding to these trials
            condition_data = EEG.data(:, :, matching_trials); % Subset of data for this condition
            
            % Determine the number of trials for this condition
            num_trials = size(condition_data, 3);
            
            % Calculate the chunk size for this condition to achieve the desired ERPs
            chunk_size = ceil(num_trials / desired_erps); % Dynamic chunk size
            
            % Calculate the actual number of chunks
            num_chunks = ceil(num_trials / chunk_size);
            
            % Initialize counter for storing data
            count_object_1 = 1;
            
            % Loop through chunks within the current condition
            for chunk_idx = 1:num_chunks
                % Define the start and end indices for the current chunk
                chunk_start = (chunk_idx - 1) * chunk_size + 1;
                chunk_end = min(chunk_idx * chunk_size, num_trials);
                
                % Get the trials in the current chunk
                chunk_trials = condition_data(:, :, chunk_start:chunk_end);
                
                % Compute the nanmean for this chunk
                averaged_data = nanmean(chunk_trials, 3);
                
                % Assign the averaged data to the dynamically created field
                % Maintain the dimension order
                object.(field_name)(count_object_1, :, :) = averaged_data;
                
                % Increment the counter
                count_object_1 = count_object_1 + 1;
            end
        end



    %% windsorize & normalize
        ERPdata = object;
        fields=fieldnames(ERPdata);
        for categoryidx=1:length(fields)Info.whichCond{1,1}(:, 1)
            categoryname=fields{categoryidx};
            data = ERPdata.(categoryname);
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
            %object = ERPdata;
        if condition == 1 && objectproperty == 6
            G_Grasp_Point_Colour = ERPdata;
            %MyInfo.name = "G_Grasp_Point" %for grasp
        elseif condition == 2 && objectproperty == 6
            K_Grasp_Point_Colour = ERPdata;
            %MyInfo.name = "K_Grasp_Point" %for kuckle
        elseif condition == 1 && objectproperty == 4
            G_Grasp_Orientation = ERPdata;
            %MyInfo.name = "G_Grasp_Point" %for grasp
        elseif condition == 2 && objectproperty == 4
            K_Grasp_Orientation = ERPdata;
            %MyInfo.name = "K_Grasp_Point" %for kuckle
        elseif condition == 1 && objectproperty == 5
            G_Grasp_Colour_Centre = ERPdata;
            %MyInfo.name = "G_Grasp_Point" %for grasp
        elseif condition == 2 && objectproperty == 5
            K_Grasp_Colour_Centre = ERPdata;
            %MyInfo.name = "K_Grasp_Point" %for kuckle
        end
        try
            cd(savepath);
        catch
            cd(savepath);
            mkdir(num2str(partNum));
            cd(partNum);
        end

    
    end
    end
    end
    % Create the filename
    filename = sprintf('exp_%d.mat', partNum);
    % Save the variables into the .mat file
    save(filename,'G_Grasp_Point_Colour','K_Grasp_Point_Colour','G_Grasp_Orientation','K_Grasp_Orientation','G_Grasp_Colour_Centre','G_Grasp_Colour_Centre','MyInfo');
    %save filename -mat G_Grasp_Point K_Grasp_Point MyInfo

    
end
