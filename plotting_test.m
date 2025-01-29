clc
clear
%% Set Variables
pthresh = 0.05; nperms = 10000; tailtest = 1;
% Define or filter the specific range
%% Load Data
load('/home/kinkini/Downloads/Neuroscience/Matthias/Leenina/Analyzed/NLSM_ClassifyFeaturesAction_ERP_7.mat'); %your results file here
win = 7 %change time concat window
%times = MyInfo.timepoints(1:win:end-2); % to accomodate concatenated 5 consecutive timepoints during classification; there's currently 61 tps
times = MyInfo.timepoints(1:win:end); % to accomodate concatenated 5 consecutive timepoints during classification; there's currently 61 tps
times = times(1:end)% will have to play with this a bit to get a perfect match since there looks to be a rounding erro rof sorts
xlab = 'Time Relative to Preview [ms]';

%% Select which Participants to Include
subs = str2num('1:14'); %Define subs here

%% Shape Significance Test
[hgraspsshape, ~, ~] = perm_1dcluster_onesample((accurate.Grasp.Shape(:,:))-.5,pthresh,nperms, tailtest); %hShape(hShape==0) = NaN;
timesSig_GraspsShape = (times(logical(hgraspsshape)));

[hknucklesshape, ~, ~] = perm_1dcluster_onesample((accurate.Knuckle.Shape(:,:))-.5,pthresh,nperms, tailtest); %hShape(hShape==0) = NaN;
timesSig_KnuckleShape = (times(logical(hknucklesshape)));

[hgraspsvsknshape, ~, ~] = perm_1dcluster_onesample_NH_paired_new(accurate.Grasp.Shape,accurate.Knuckle.Shape,pthresh,nperms,tailtest); %hShape(hShape==0) = NaN;
timesSig_GraspsvsKnShape = (times(logical(hgraspsvsknshape)));

%% Material ("colour") Significant Test
[hgraspscolour, ~, ~] = perm_1dcluster_onesample((accurate.Grasp.Colour)-.5,pthresh,nperms, tailtest); %hShape(hShape==0) = NaN;
timesSig_GraspsColour = (times(logical(hgraspscolour)));

[hknucklescolour, ~, ~] = perm_1dcluster_onesample((accurate.Knuckle.Colour)-.5,pthresh,nperms, tailtest); %hShape(hShape==0) = NaN;
timesSig_KnuckleColour = (times(logical(hknucklescolour)));

[hgraspsvskncolour, ~, ~] = perm_1dcluster_onesample_NH_paired_new(accurate.Grasp.Colour,accurate.Knuckle.Colour,pthresh,nperms,tailtest); %hShape(hShape==0) = NaN;
timesSig_GraspsvsKnColour = (times(logical(hgraspsvskncolour)));

%% Orientation Significant Test
[hgraspsorientation, ~, ~] = perm_1dcluster_onesample((accurate.Grasp.Orientation)-.5,pthresh,nperms, tailtest); %hShape(hShape==0) = NaN;
timesSig_GraspsColour = (times(logical(hgraspsorientation)));

[hknucklesorientation, ~, ~] = perm_1dcluster_onesample((accurate.Knuckle.Orientation)-.5,pthresh,nperms, tailtest); %hShape(hShape==0) = NaN;
timesSig_KnuckleColour = (times(logical(hknucklesorientation)));

[hgraspsvsknorientation, ~, ~] = perm_1dcluster_onesample_NH_paired_new(accurate.Grasp.Colour,accurate.Knuckle.Colour,pthresh,nperms,tailtest); %hShape(hShape==0) = NaN;
timesSig_GraspsvsKnOrientation = (times(logical(hgraspsvsknorientation)));

%% Plot Figure: Shape
figure(); hold on; set(gcf, 'color', 'w');
title('Shape Classification Grasp v Knuckle');
numbins = size(times,2);
data = accurate.Grasp.Shape(subs,:);%grasp classification
data2 = accurate.Knuckle.Shape(subs,:);%knuckle classification
clr = [0 0.2 0.6]; %dark blue
clr2 = [0.4 0.6 0.8]; %light blue
hbar_height = 0.45;
set(findobj(gcf,'type','axes'),'FontSize',12,'FontWeight','Bold','LineWidth',0.5);
%axlimit = [-100 500 0.4 .7];
axlimit = [-100 500 0.1 .9];
axis(axlimit);
meandata = mean(data);
semdata = std(data)/sqrt(size(data,1));
h = shadedErrorBar(times, meandata, semdata);
h.mainLine.LineWidth = 1.5;
h.mainLine.Color = clr;
h.patch.FaceColor = clr;
h.patch.FaceAlpha = 0.3; hold on;
meandata2 = mean(data2);
semdata2 = std(data2)/sqrt(size(data2,1));
h = shadedErrorBar(times, meandata2, semdata2, 'k');
h.mainLine.LineWidth = 1.5;
h.mainLine.Color = clr2;
h.patch.FaceColor = clr2;
h.patch.FaceAlpha = 0.3;
scatter2=scatter(times,hgraspsshape*0.72,40,'s','filled','markerfacecolor',clr,'markeredgecolor', 'none' );
hknucklesshape(hknucklesshape==0) = NaN;
scatter3=scatter(times,hknucklesshape*0.7,40,'s','filled','markerfacecolor',clr2,'markeredgecolor', 'none' );
scatter1.MarkerFaceAlpha = .8;
scatter2.MarkerFaceAlpha = .8;
scatter3.MarkerFaceAlpha = .8;
plot(times, ones(length(times))*.5,'k--','LineWidth',0.3,'MarkerSize',5); % baseline
startline = line([0 0],[axlimit(3:4)],'color',[0.7 0.7 0.7],'LineWidth', 1);

%% Plot Material/Colour
figure(); hold on; set(gcf, 'color', 'w');
title('Material Classification Grasp v Knuckle');
numbins = size(times,2);
data = accurate.Grasp.Colour(subs,:);%grasp classification
data2 = accurate.Knuckle.Colour(subs,:);%knuckle classification
clr = [0.6 0 0]; %dark red
clr2 = [0.8 0.2 0]; %light red
hbar_height = 0.45;
set(findobj(gcf,'type','axes'),'FontSize',12,'FontWeight','Bold','LineWidth',0.5);
axlimit = [-100 500 0.1 .9];
axis(axlimit);
meandata = mean(data);
semdata = std(data)/sqrt(size(data,1));
h = shadedErrorBar(times, meandata, semdata, 'k');
h.mainLine.LineWidth = 1.5;
h.mainLine.Color = clr;
h.patch.FaceColor = clr;
h.patch.FaceAlpha = 0.3; hold on;
meandata2 = mean(data2);
semdata2 = std(data2)/sqrt(size(data2,1));
h = shadedErrorBar(times, meandata2, semdata2, 'k');
h.mainLine.LineWidth = 1.5;
h.mainLine.Color = clr2;
h.patch.FaceColor = clr2;
h.patch.FaceAlpha = 0.3;
scatter2=scatter(times,hgraspscolour*0.72,40,'s','filled','markerfacecolor',clr,'markeredgecolor', 'none' );
hknucklescolour(hknucklescolour==0) = NaN;
scatter3=scatter(times,hknucklescolour*0.7,40,'s','filled','markerfacecolor',clr2,'markeredgecolor', 'none' );
scatter1.MarkerFaceAlpha = .8;
scatter2.MarkerFaceAlpha = .8;
scatter3.MarkerFaceAlpha = .8;
plot(times, ones(length(times))*.5,'k--','LineWidth',0.3,'MarkerSize',5); % baseline
startline = line([0 0],[axlimit(3:4)],'color',[0.7 0.7 0.7],'LineWidth', 1);


%% Plot Figure: Orientation
figure(); hold on; set(gcf, 'color', 'w');
title('Orientation Classification Grasp v Knuckle');
numbins = size(times,2);
data = accurate.Grasp.Orientation(subs,:);%grasp classification
data2 = accurate.Knuckle.Orientation(subs,:);%knuckle classification
clr = [0 0.4 0]; %dark green
clr2 = [0.6 1 0]; %light green
hbar_height = 0.45;
set(findobj(gcf,'type','axes'),'FontSize',12,'FontWeight','Bold','LineWidth',0.5);
%axlimit = [-100 500 0.4 .7];
axlimit = [-100 500 0.1 .9];
axis(axlimit);
meandata = mean(data);
semdata = std(data)/sqrt(size(data,1));
h = shadedErrorBar(times, meandata, semdata);
h.mainLine.LineWidth = 1.5;
h.mainLine.Color = clr;
h.patch.FaceColor = clr;
h.patch.FaceAlpha = 0.3; hold on;
meandata2 = mean(data2);
semdata2 = std(data2)/sqrt(size(data2,1));
h = shadedErrorBar(times, meandata2, semdata2, 'k');
h.mainLine.LineWidth = 1.5;
h.mainLine.Color = clr2;
h.patch.FaceColor = clr2;
h.patch.FaceAlpha = 0.3;
scatter2=scatter(times,hgraspsorientation*0.72,40,'s','filled','markerfacecolor',clr,'markeredgecolor', 'none' );
hknucklesorientation(hknucklesorientation==0) = NaN;
scatter3=scatter(times,hknucklesorientation*0.7,40,'s','filled','markerfacecolor',clr2,'markeredgecolor', 'none' );
scatter1.MarkerFaceAlpha = .8;
scatter2.MarkerFaceAlpha = .8;
scatter3.MarkerFaceAlpha = .8;
plot(times, ones(length(times))*.5,'k--','LineWidth',0.3,'MarkerSize',5); % baseline
startline = line([0 0],[axlimit(3:4)],'color',[0.7 0.7 0.7],'LineWidth', 1);

% These are the cluster-based stats functions
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
