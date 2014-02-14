% This is a MATLAB script for the
% CLPS1291 lab on attractiveness

% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
% Author: Thomas Serre
% The k-means code is a modified version of the
% simple_kmedias function by Mauricio Martinez-Garcia, 2003,2007
% Brown University
% CLPS Department
% email: Thomas_Serre@Brown.edu
% Website: http://serre-lab.clps.brown.edu
% February 2014;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Execute the code in each individual cell by moving the  %%
% cursor to a cell and press cmd+<enter>                  %%
% Below we will be playing with faces that I extracted    %%
% from a popular computer vision database used to test    %%
% face rec algorithms and called 'Faces in the Wild'      %%
% First you should download the data at                   %%
% https://www.dropbox.com/s/35kps5eb5j6sjc0/FACES.mat     %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
clear all;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Model neuron coding for distance from average for toy data
%% Excercise for home: Try to run the code on your faces...
clear all;
x   = 2.0 + 3.0*randn(2,30);
sim = 'rbf'; % 'rbf' or 'dot-prod' or 'sig'

% Algorithm starts
ind = randperm(size(x,2)); %% randomize the presentation order
m   = mean(x,2); %% sample mean
w   = m;

for ii = 1:size(x,2)
    subplot(1,2,1); plot(x(1,:), x(2,:), 'o');
    hold on;
    
    switch sim
        case 'rbf'
            % compute distance between current stimulus x(:,ind(ii)) and synaptic weight
            D2   = sum((x(:,ind(ii))-w).^2);
            sig2 = 100;
            % exponentiate
            y = exp(-D2/sig2);
            
            plot([m(1) x(1,ind(ii))], [m(2) x(2,ind(ii))], '-r', 'MarkerSize', 10);
            plot(m(1), m(2), 'pm', 'MarkerSize', 10);
            legend('samples', 'current','prototype')
            hold off
            subplot(1,2,2); bar(1,y); axis([0 2 0 1]); axis off
            
        case 'dot-prod'
            y = w'*x(:,ind(ii));
            plot([0 x(1,ind(ii))], [0 x(2,ind(ii))], '-r', 'MarkerSize', 10);
            plot([0 m(1)], [0  m(2)], '-g', 'MarkerSize', 10);
            legend('samples', 'current','prototype')
            hold off
            subplot(1,2,2); bar(1,y); axis([0 2 -10 10]); axis off
            
        case 'sig'
            y = logsig(w'*x(:,ind(ii)));
            
            plot([0 x(1,ind(ii))], [0 x(2,ind(ii))], '-r', 'MarkerSize', 10);
            plot([0 m(1)], [0  m(2)], '-g', 'MarkerSize', 10);
            
            
            legend('samples', 'current','prototype')
            hold off
            subplot(1,2,2); bar(1,y); axis([0 2 0 1]); axis off
    end
    pause
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% K-means demo

close all;
clear all;

K  = 3; %% K for k-means

c  = {'k.', 'k.', 'k.', 'ko', 'k<'};
c1 = {'bs', 'rp', 'm>', 'go', 'k<'};
c2 = {'b',  'r',  'm',  'g',  'k'};
c3 = {'s',  'p',  '>',  'o',  '<'};

dim     = 2;
my_var  = 1;
my_sdev = sqrt(my_var);

%%init
dist    = zeros(1,K);
maxerr  = 0;
m       = 1;


cluster1 = my_sdev*randn(200,dim) + kron(ones(200,1),[0,0]);
cluster2 = my_sdev*randn(200,dim) + kron(ones(200,1),[0,5]);
cluster3 = my_sdev*randn(300,dim) + kron(ones(300,1),[-5,0]);

% Build data matrix with corresponding labels lab
X    = [cluster1 ; cluster2 ; cluster3];
lab  = [ones(size(cluster1,1),1); ...
    2*ones(size(cluster2,1),1); 3*ones(size(cluster3,1),1)];

[Ndata, dims] = size(X);

figure
hold on;
plot(cluster1(:,1),cluster1(:,2), c{1});
plot(cluster2(:,1),cluster2(:,2), c{2});
plot(cluster3(:,1),cluster3(:,2), c{3});
axis([-10 5 -6 10]);


% Initial prototype assignment (arbitrary)
ind = randperm(size(X,1));
for i=1:K-1
    means(i,:) = X(ind(i),:);
end
means(K,:) = mean(X(K:Ndata,:));

cmp = 1 + maxerr;
while (cmp > maxerr)
    % Sums (class) and data counters (Nclass) initialization
    class  = zeros(K,dims);
    Nclass = zeros(K,1);
    m = m+1;
    
    myind = [];
    hold on;
    
    for i = 1:K
        plot(means(i,1), means(i,2), c1{i}, 'MarkerSize', 14, 'MarkerFaceColor', c2{i}, 'MarkerEdgeColor', 'k');
    end
    pause
    
    % Groups each elements to the nearest prototype
    for i=1:Ndata
        for j=1:K
            % Euclidean distance from data to each prototype
            dist(j) = norm(X(i,:)-means(j,:))^2;
        end
        % Find indices of minimum distance
        index_min = find(~(dist-min(dist)));
        % If there are multiple min distances, decide randomly
        index_min = index_min(ceil(length(index_min)*rand));
        class(index_min,:) = class(index_min,:) + X(i,:);
        Nclass(index_min) = Nclass(index_min) + 1;
        myind = [myind index_min];
    end
    err = 0;
    for i=1:K
        class(i,:) = class(i,:) / Nclass(i);
        ind = find(myind==i);
        for j = 1:length(ind)
            plot(X(ind(j),1),X(ind(j),2), [c2{i} c3{lab(ind(j))}]);
        end
        err = err+sum(sqrt(sum((X(ind,:)-repmat(class(i,:),length(ind),1)).^2,2)));
    end
    hold off;
    title(['Objective function: ' num2str(err)])
    pause
    
    % Compare results with previous iteration
    cmp = 0;
    for i=1:K
        cmp = norm(class(i,:)-means(i,:));
    end
    
    % Prototype update
    means = class;
end

Nmeans = Nclass;




