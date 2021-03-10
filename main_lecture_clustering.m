% This is a MATLAB script for the
% CLPS1291 lecture on k-means.

% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
% Author: Thomas Serre
% Brown University
% CLPS Department
% email: Thomas_Serre@Brown.edu
% Website: http://serre-lab.clps.brown.edu
% February 2014;
% The k-means code is a modified version of the
% simple_kmedias function by Mauricio Martinez-Garcia, 2003,2007

clc;
clear all;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Model neuron coding for distance from average for toy data
%% Excercise for home: Try to run the code on your faces...
clear all;
x   = 2.0 + 3.0*randn(2,30);
sim = 'dot-prod'; % 'rbf' or 'dot-prod' or 'sig'

%% Excercise: Try to implement the normalized dot-product

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

col = [1 0 0; 0 1 0; 0 0 1; 1 1 0; 1 0 1; 0 1 1];

dim     = 2;
%%init
dist    = zeros(1,K);
maxerr  = 0;
m       = 1;

clus_size = [100 100 700];
my_var    = [2 2 2];
my_sdev   = sqrt(my_var);


cluster1 = my_sdev(1)*randn(clus_size(1),dim) + kron(ones(clus_size(1),1),[0,0]);
cluster2 = my_sdev(2)*randn(clus_size(2),dim) + kron(ones(clus_size(2),1),[0,5]);
cluster3 = my_sdev(3)*randn(clus_size(3),dim) + kron(ones(clus_size(3),1),[-5,0]);

% Build data matrix with corresponding labels lab
X    = [cluster1 ; cluster2 ; cluster3];
lab  = [ones(size(cluster1,1),1); ...
    2*ones(size(cluster2,1),1); 3*ones(size(cluster3,1),1)];

[Ndata, dims] = size(X);

figure
scatter(X(:,1),X(:,2), 30, [0 0 0]);
axis([-10 5 -6 10]);
hold on;

% Initial prototype assignment (arbitrary)
ind = randperm(size(X,1));
for i=1:K-1
    means(i,:) = X(ind(i),:);
end
means(K,:) = mean(X(K:Ndata,:));

cmp = 1 + maxerr;
while (cmp > maxerr)
    
    class  = zeros(K,dims);
    Nclass = zeros(K,1);
    myind  = [];
    
    scatter(means(:,1), means(:,2), 200, col(1:K,:), 'fill' );
    
    pause(.5);
    
    % Groups each elements to the nearest prototype
    for ii = 1:Ndata
        for jj = 1:K
            % Euclidean distance from data to each prototype
            dist(jj) = norm(X(ii,:)-means(jj,:))^2;
        end
        % Find indices of minimum distance
        index_min = find(~(dist-min(dist)));
        
        % If there are multiple min distances, decide randomly
        index_min          = index_min(ceil(length(index_min)*rand));
        class(index_min,:) = class(index_min,:) + X(ii,:);
        Nclass(index_min)  = Nclass(index_min) + 1;
        myind              = [myind index_min];
    end
    
    err = 0;
    for ii = 1:K
        class(ii,:) = class(ii,:) / Nclass(ii);
        ind         = find(myind==ii);
        
        scatter(X(ind,1),X(ind,2), 30, repmat(col(ii,:), length(ind),1));
        
        err = err+sum(sqrt(sum((X(ind,:)-repmat(class(ii,:),length(ind),1)).^2,2)));
    end
    
    title(['Objective function: ' num2str(err)])
    pause(.5);
    
    % Compare results with previous iteration
    cmp = 0;
    for ii = 1:K
        cmp = norm(class(ii,:)-means(ii,:));
    end
    
    % Prototype update
    means = class;
end

Nmeans = Nclass;




