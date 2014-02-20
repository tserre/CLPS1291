% This is a MATLAB script for the
% CLPS1291 lab on clustering

% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
% Author: Thomas Serre
% Brown University
% CLPS Department
% email: Thomas_Serre@Brown.edu
% Website: http://serre-lab.clps.brown.edu
% February 2014;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FIRST PLEASE DOWNLOAD DATA AT                           %%
%% https://www.dropbox.com/s/7xp2h2kzle24vck/Archive.zip   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Segmentation example
d = dir('../Data/images/*.tif');
k = 2; % The k of kmeans
close all

for ii = 1:length(d)
    
    img = double(imresize(imread(['../Data/images/' d(ii).name]), .1));
    siz = size(img);
    
    figure(ii)
    subplot(2,2,1)
    imagesc(img)
    axis off; axis square;
        
    % Reshape your image so that individual (R,G,B) pixel
    % intensities are treated as samples when you call kmeans and
    % run kmeans on your image
    
    % A = reshape();

    [ind, C] = kmeans(A, k);
    
    % Plot a random subset of pixels together with the cluster
    % centers and the sample mean -- use plot3
    
    %     I = randperm();
    %     A = ;
    
    subplot(2,2,3)
    %     plot3( ,'.');
    hold on;
    %     plot3( , 'p', 'MarkerEdgeColor','r', ...
    %         'MarkerFaceColor', 'r' , 'MarkerSize', 12);
    
    m = mean(A);
    %     plot3( ,'o', 'MarkerEdgeColor','g', ...
    %         'MarkerFaceColor', 'g' , 'MarkerSize',12);
    hold off;
    
    xlabel('red');
    ylabel('green');
    zlabel('blue');
    
    % Show an image corresponding to the cluster assignment
    subplot(2,2,2)
    %     imagesc( );
    colorbar;
    axis off; axis square;
end


%% Additional excercise 2: Create a funtion called segment_image
% which takes as input an image and returns prototypes and
% segments for that image

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

clus_size = [300 100 700];
my_var    = [2 2 1];
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



%% Kmeans on digits
clc;
clear all;
close all;

% Load data
load ../Data/mnist_all

% Here you can play with different digits train0 -- train9
A   = cat(1, train0, train1, train2, train3, train4, train5, ...
    train6, train7, train8, train9);
A   = double(A)/255;
lab = cat(1, 0*ones(size(train0,1),1), 1*ones(size(train1,1),1), ...
    2*ones(size(train2,1),1), 3*ones(size(train3,1),1), ...
    4*ones(size(train4,1),1), 5*ones(size(train5,1),1), ...
    6*ones(size(train6,1),1), 7*ones(size(train7,1),1), ...
    8*ones(size(train8,1),1), 9*ones(size(train9,1),1));
    
% A   = double(train8)/255;

I   = randperm(size(A,1));
A   = A(I(1:1000),:);
lab = lab(I(1:1000));
siz = sqrt(size(A,2));

% run kmeans    
k = 10; % The k of kmeans

[ind, C, sumD, D] = kmeans(A, k);  % try help on kmeans, try different distance/similarity 
                                   % measures to see how they impact your results
                                                                
% Visualize the prototypes using subplots or montage
% Vary k and comment.
figure(1);
m1 = ceil(sqrt(k));
m2 = ceil(k/m1);

for ii = 1:k
    subplot(m1,m2,ii), imshow(reshape(C(ii,:), siz, siz)');
end


%% Vector quantization -- Think of coding individual digits with
% one of your learned prototypes

figure(2)
n = 24;

for ii = 1:n
    m = ceil(sqrt(3*n));
    % Below show inidviudal digits
    subplot(m, m, 3*(ii-1)+1), imshow(reshape(A(ii,:), [siz siz])');
    axis off; axis square; colormap gray;  title('Digit');
    % Below show the closest prototype
    subplot(m, m, 3*(ii-1)+2), imshow(reshape(C(ind(ii),:), [siz siz])');
    axis off; axis square; colormap gray; title(['Proto. ' num2str(ind(ii))]);
    % Below show the reconstruction error (the image difference
    % between the sample and the closest prototype
    subplot(m, m, 3*(ii-1)+3), imshow(reshape(A(ii,:)-C(ind(ii),:), [siz siz])');
    axis off; axis square; colormap gray; title(['Error: ' num2str(D(ii,ind(ii)))]);
end

suptitle([num2str(k) ' clusters (' num2str(mean(min(D,[],2))) ')']);
