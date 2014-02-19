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


%% Kmeans on digits
clc;
clear all;
close all;

% Load data
load ../Data/mnist_all

% Here you can play with different digits train0 -- train9
% these are siz*siz images where siz = 28;
A   = double(train8)/255;
siz = sqrt(size(A,2));

% Below we sample a subset of the data to speed up the
% computations
I   = randperm(size(A,1));
A   = A(I(1:1000),:);

% run kmeans
% try help on kmeans, try different distance/similarity 
% measures to see how they impact your results

k = 1; % The k of kmeans
% [ind, C, sumD, D] = kmeans( );  
                                                                
% Visualize the prototypes using subplots or montage
% Vary k and comment.
figure(1);



%% Vector quantization -- Think of coding individual digits with
% one of your learned prototypes

figure(2)
n = 12;
m = ceil(sqrt(3*n));

for ii = 1:n
    % Below show inidviudal digits
    subplot(m, m, 3*(ii-1)+1);
    
    axis off; axis square; colormap gray;  title('Digit');
    
    % Below show the closest prototype
    subplot(m, m, 3*(ii-1)+2);
    
    axis off; axis square; colormap gray; title(['Proto. ' num2str(ind(ii))]);
    
    % Below show the reconstruction error (the image difference
    % between the sample and the closest prototype
    subplot(m, m, 3*(ii-1)+3);
    
    axis off; axis square; colormap gray; title(['Error: ' num2str(D(ii,ind(ii)))]);
end

suptitle([num2str(k) ' clusters (' num2str(mean(min(D,[],2))) ')']);

%% Additional excercise 1: Get kmeans to run on your own data

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