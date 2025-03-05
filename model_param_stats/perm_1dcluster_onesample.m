
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