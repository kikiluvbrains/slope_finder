function accuracy_permutation_tester(path, file_1, file_2, name_1, name_2, Class, Object_property, win)
    % Function to perform a two-tailed paired one sample permutation test based on
    % Lawerences original functions
    % cubes in red shows areas of signficance
    %for function calls please ensure that both mat files have the items
    %accurate and Myinfo as the same name, the function assumes that these are
    %the data variables that are in your script, alteration of said varaibles
    %will lead to the function terminating
    % Author: @ Kinkini M and Mani S, based on Nina Lee and Lawrence gao
    %
    % Inputs:
    % - data, data2: matrices of data for accuracy from svm
    % - name_1, name_2: legend Names, e.g., 'slope' 'intercept' 
    % - Class: e.g., 'Grasp'
    % - Object_property: e.g., 'Shape'
    % - win: Time concatenation window
    %
    %   example input
    %   path = '/home/kinkini/Downloads/'
    %   file_1 = 'NLSM_ClassifyFeaturesAction_ERP_6_slope_1_slope.mat'
    %   file_2 = 'NLSM_ClassifyFeaturesAction_ERP_6_slope_1_intercept.mat'
    %   name_1 = 'slope'
    %   name_2 = 'intercept'
    %   Class = 'Grasp'
    %   Object_property = 'Shape'
    %   win = 5 %change time concat window
    
    data_1 = load([path, file_1])
    accurate_1 = data_1.accurate
    My_info_1 = data_1.MyInfo
    times_1 = My_info_1.timepoints(1:win:end-2); % to accomodate concatenated 5 consecutive timepoints during classification; there's currently 61 tps
    clear accurate MyInfo
    
    data_2 = load([path, file_2])
    accurate_2 = data_2.accurate
    My_info_2 = data_2.MyInfo
    times_2 = My_info_2.timepoints(1:win:end-2); % to accomodate concatenated 5 consecutive timepoints during classification; there's currently 61 tps
    clear accurate MyInfo
    
    if times_1 == times_2
        
    %% Set Variables
    pthresh = 0.05; nperms = 10000; tailtest = 0;
    
    
    %% Select which Participants to Include
    subs = str2num('1:14'); %Define subs here
    
    %% Object_property Significance Test
    %[hgraspsObject_property_1, ~, ~] = perm_1dcluster_onesample((accurate_1.(Class).(Object_property)(:,:))-.5,pthresh,nperms, tailtest); %hObject_property(hObject_property==0) = NaN;
    %timesSig_GraspsObject_property_1 = (times(logical(hgraspsObject_property_1)));
    %[hgraspsObject_property_2, ~, ~] = perm_1dcluster_onesample((accurate_2.Grasp.Object_property(:,:))-.5,pthresh,nperms, tailtest); %hObject_property(hObject_property==0) = NaN;
    %timesSig_GraspsObject_property_2 = (times(logical(hgraspsObject_property_2)));
    
    %[hknucklesObject_property_1, ~, ~] = perm_1dcluster_onesample((accurate_1.Knuckle.Object_property(:,:))-.5,pthresh,nperms, tailtest); %hObject_property(hObject_property==0) = NaN;
    %timesSig_KnuckleObject_property_1 = (times(logical(hknucklesObject_property_1)));
    %[hknucklesObject_property_2, ~, ~] = perm_1dcluster_onesample((accurate_2.Knuckle.Object_property(:,:))-.5,pthresh,nperms, tailtest); %hObject_property(hObject_property==0) = NaN;
    %timesSig_KnuckleObject_property_2 = (times(logical(hknucklesObject_property_2)));
    
    [h_acc1_vs_acc2_property, ~, ~] = perm_1dcluster_onesample_NH_paired_new((accurate_1.(Class).(Object_property)),(accurate_2.(Class).(Object_property)),pthresh,nperms,tailtest); %hObject_property(hObject_property==0) = NaN;
    timesSig_h_acc1_vs_acc2_property = (times_1(logical(h_acc1_vs_acc2_property)));
    %[hgraspsvsknObject_property, ~, ~] = perm_1dcluster_onesample_NH_paired_new(accurate_2.Knuckle.(Object_property),accurate_1.Knuckle.(Object_property),pthresh,nperms,tailtest); %hObject_property(hObject_property==0) = NaN;
    %timesSig_GraspsvsKnObject_property = (times(logical(hgraspsvsknObject_property)));
    
    
    %% Plot Figure: Shape
    figure(); hold on; set(gcf, 'color', 'w');
    title( Class + " Classification for Object " + Object_property);
    numbins = size(times_1,2);
    data = accurate_1.(Class).(Object_property)(subs,:);%data 1 classification
    data2 = accurate_2.(Class).(Object_property)(subs,:);%data 2 classification
    clr = [0 0.2 0.6]; %dark blue
    clr2 = [0.4 0.6 0.8]; %light blue
    clr3 = [0.5 0.1 0.1]; %light blue
    hbar_height = 0.45;
    set(findobj(gcf,'type','axes'),'FontSize',12,'FontWeight','Bold','LineWidth',0.5);
    %axlimit = [-100 500 0.4 .7];
    axlimit = [-100 500 0.1 .9];
    axis(axlimit);
    meandata = mean(data);
    semdata = std(data)/sqrt(size(data,1));
    h1 = shadedErrorBar(times_1, meandata, semdata);
    h1.mainLine.LineWidth = 1.5;
    h1.mainLine.Color = clr;
    h1.patch.FaceColor = clr;
    h1.patch.FaceAlpha = 0.3; hold on;
    meandata2 = mean(data2);
    semdata2 = std(data2)/sqrt(size(data2,1));
    h2 = shadedErrorBar(times_1, meandata2, semdata2, 'k');
    h2.mainLine.LineWidth = 1.5;
    h2.mainLine.Color = clr2;
    h2.patch.FaceColor = clr2;
    h2.patch.FaceAlpha = 0.3;
    %scatter2=scatter(times,hgraspsshape*0.72,40,'s','filled','markerfacecolor',clr,'markeredgecolor', 'none' );
    %hknucklesshape(hknucklesshape==0) = NaN;
    %scatter3=scatter(times,hknucklesshape*0.7,40,'s','filled','markerfacecolor',clr2,'markeredgecolor', 'none' );
    scatter4=scatter(times_1,h_acc1_vs_acc2_property*0.7,40,'s','filled','markerfacecolor',clr3,'markeredgecolor', 'none' );
    %scatter1.MarkerFaceAlpha = .8;
    %scatter2.MarkerFaceAlpha = .8;
    %scatter3.MarkerFaceAlpha = .8;
    scatter4.MarkerFaceAlpha = .8;
    plot(times_1, ones(length(times_1))*.5,'k--','LineWidth',0.3,'MarkerSize',5); % baseline
    startline = line([0 0],[axlimit(3:4)],'color',[0.7 0.7 0.7],'LineWidth', 1);
    hold off;
    %legend(cl2, {name_1, name_2}, 'Location', 'best');
    legend([h1.mainLine, h2.mainLine], {name_1, name_2}, 'Location', 'best');
    
    
    else
    disp("your times do not match between the two files")
    end
    
    
    % These are the cluster-based stats functions
    %Author: @Lawrence Gao
    %% Description h = perm_cluster_onesample(data,thresh,nperms) 
    % returns the cluster corrected outcomes
    % data input: nSubs x nTime
    function [h_corrected, crit_h_size, h_true] = perm_1dcluster_onesample(data, thresh, nperms, tail)
    rng('Shuffle');
    if tail == -1
        tail = 'left';
    elseif tail == 0
        tail = 'both';
    elseif tail == 1
        tail = 'right';
    end
    
    nsubs = size(data,1);
    ntimes = size(data,2);
    
    h_true = zeros(1,ntimes);
    % true t-test outcome
    [~, p, ~, ~] = ttest(data, zeros(size(data)),'tail', tail);
    h_true(p<thresh|p==thresh)=1;
    
    
    % permutations
    clustersize = zeros(1,nperms);
    for iperm = 1:nperms
        permdata = data;
        h = zeros(1,ntimes);
        
        % sign based permutation
        selSample = randsample(nsubs, floor(nsubs/2)); % take a sample of half of subs
        permdata(selSample,:) = permdata(selSample,:).*-1;%.*-1; % take the classification results of selected sub sample (already -50% to test against 0) & flip sign
        % t-test
        [~,p,~,~] = ttest(permdata, zeros(size(data)), 'tail', tail); %test classification results against 0 note replace zeros
        h(p<thresh|p==thresh) = 1;
        % find cluster
        h = [0,h,0]; % add zero to beginning and end
        findstart = find(diff(h)==1); % find position in timepoints/data where h 0 turns to 1
        findend = find(diff(h)==-1)-1; % find position in timepoints/data where h 1 turns to 0 and go 1 back to have correct size
        
        if isempty(findstart)
            continue % skip if no cluster
        else
            clustersize(iperm) = max(findend-findstart); % determine biggest cluster size by chance; one off point is not cluster (2 consecutive points is cluster size 1, etc.)
        end
        end
        
        crit_h_size = prctile(clustersize,95); % determine 95 percentile cluster size (smallest cluster size required to be meaningful)
        
        
        % find clusters of ones in true_h that meet the 95% largest cluster
        % criterion
        h = [0 h_true 0];
        findstart = find(diff(h)==1);
        findend = find(diff(h)==-1)-1;
        clusters = findend-findstart;
        
        nclusters = length(clusters); c = 0;
        
        
        for i = 1:nclusters % for each cluster found, determine if meets 95 perccentile threshold
            if clusters(i) > crit_h_size || clusters(i) == crit_h_size
                c = c+1;
                cluster_start(c) = findstart(i); % index significant cluster again
                cluster_end(c) = findend(i);
            end
        end
        h_corrected = zeros(1,ntimes);
        for i = 1:c
            h_corrected(cluster_start(i):cluster_end(i)) = 1;
        end
        end
        
        function [h_corrected, crit_h_size, h_true] = perm_1dcluster_onesample_NH_paired_new(data, data2, thresh, nperms, tail)
        rng('shuffle');
        if tail == -1
            tail = 'left';
        elseif tail == 0
            tail = 'both';
        elseif tail == 1
            tail = 'right';
        end
        
        nsubs = size(data,1);
        ntimes = size(data,2);
        
        h_true = zeros(1,ntimes);
        % true t-test outcome
        [~, p, ~, ~] = ttest(data, data2,'tail', tail);
        h_true(p<thresh|p==thresh)=1;
        
        % permutations
        clustersize = zeros(1,nperms);
        for iperm = 1:nperms
            permdata = data;
            permdata2 = data2;
            h = zeros(1,ntimes);
            
            % sign based permutation
            selSample = randsample(nsubs, floor(nsubs/2));
            permdata(selSample,:) = permdata(selSample,:)*-1;
            permdata2(selSample,:) = permdata2(selSample,:)*-1;
            % t-test
            [~,p,~,~] = ttest(permdata, permdata2, 'tail', tail);
            h(p<thresh|p==thresh) = 1;
            
            % find cluster
            h = [0,h,0]; % add zero to beginning and end
            findstart = find(diff(h)==1); % find position in timepoints/data where h 0 turns to 1
            findend = find(diff(h)==-1)-1; % find position in timepoints/data where h 1 turns to 0 and go 1 back to have correct size
            
            if isempty(findstart)
                continue % skip if no cluster
            else
                clustersize(iperm) = max(findend-findstart); % determine biggest cluster size by chance; one off point is not cluster (2 consecutive points is cluster size 1, etc.)
            end
        end
        
        crit_h_size = prctile(clustersize,95);
        
        % find clusters of ones in true_h that meet the 95% largest cluster
        % criterion
        h = [0 h_true 0];
        findstart = find(diff(h)==1);
        findend = find(diff(h)==-1)-1;
        clusters = findend-findstart;
        
        nclusters = length(clusters); c = 0;
        
        
        for i = 1:nclusters % for each cluster found, determine if meets 95 perccentile threshold
            if clusters(i) > crit_h_size || clusters(i) == crit_h_size
                c = c+1;
                cluster_start(c) = findstart(i); % index significant cluster again
                cluster_end(c) = findend(i);
            end
        end
        h_corrected = zeros(1,ntimes);
        for i = 1:c
            h_corrected(cluster_start(i):cluster_end(i)) = 1;
        end
        end
end