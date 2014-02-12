% This is a MATLAB script for the 
% CLPS1291 lab on attractiveness 

% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
% Author: Thomas Serre 
% Brown University
% CLPS Department
% email: Thomas_Serre@Brown.edu
% Website: http://serre-lab.clps.brown.edu
% February 2014; 
% Execute the code in each individual cell by moving the cursor to a cell and press cmd+<enter>
% Below we will be playing with faces that I extracted from a
% popular computer vision database used to test face recognition
% algorithms and called 'Faces in the Wild'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% First you should download the data at                  %%
%% https://www.dropbox.com/s/35kps5eb5j6sjc0/FACES.mat    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
clear all;
close all;

%% load the FACES variable 
% find out what happened? what variable got loaded? what is their size?


% number of images in the database and their size
he    = size(IMG, 1); % image height
wi    = size(IMG, 2); % image width
Nimg  = size(IMG, 3);  % num images

%% use the imshow command to visualize a few faces

%% try the 'montage' command, use the zoom to see individual faces
% you will need to reshape your nxnxm IMG array into nxnx1xm array using the reshape command


%% Scan through the first couple of faces and display them
for ii =1:30
    ii
    pause(.2) % Use the 'pause' command to make sure things get displayed
end

Face1 = IMG(:,:,5);
subplot(2,2,1); imshow(Face1)
Face2 = IMG(:,:,6);
subplot(2,2,2); imshow(Face2)

a = .5; % Morphing parameter;
%% Let's try to morph between Face1 and Face2
% Morph = ;
subplot(2,2,3); imshow(Morph)

%% Now let's assess how attractive Angelina is...
close all;
imshow(Angelina); 

ind = randperm(Nimg); %% randomize the presentation order
M   = mean(IMG,3); %% sample mean for comparison
w   = zeros(he*wi,1);

% Let's run our simple average learning algorithm on the face
% database
for ii = 1:Nimg
    x  = IMG(:,:,ind(ii));
    % dw = ; %% compute update
    % w  = ; %% update w
    
    % Show the current image
    % NOTE: We use reshape to convert a row or column vector back into an
    % NxM image
    
    subplot(2,2,1)
    imagesc(x); % visualize the first face in the DB
    colormap gray;
    axis square; axis off;
    title(['Iter: ' num2str(ii)])

    subplot(2,2,2)
    imagesc(M); % visualize the actual average
    colormap gray;
    axis square; axis off;
    title('True sample mean estimate')

    subplot(2,2,3);
    imagesc(reshape(w, he, wi)); % visualize the weights
    colormap gray;
    axis square; axis off;
    title('Current mean estimate')
    
    subplot(2,2,4)
    imagesc(reshape(dw, he, wi)); % visualize the weight update
    colormap gray;
    axis square; axis off;
    title('Current weight update')
    pause(.1)
end


%% Compute the distance between Angelina and the average face
% using the pdist2 command
% DA = pdist2();

%% Compute the average distance between a face in the db and the average face
DF = pdist2();

disp(['Mean face-to-average distance: ' num2str(round(mean(DF)))]);
disp(['Angelina-to-average distance:  ' num2str(round(DA))]);


% organize people according to distance from average face using
% the sort command
[Val, Ind] = sort(DF);
IMG        = IMG(:,:,Ind);%% reorgnize people from closest to farthest

figure(3)
montage(reshape(IMG, [size(IMG,1) size(IMG,2) 1 size(IMG,3)]));

%% find people that seem more attractive than Angelina using the find command
% Ind     = find();

figure(4)
montage(reshape(IMG(:,:,Ind), [size(IMG,1) size(IMG,2) 1 nPeople]));

%% find people that seem less attractive than Angelina
% Ind     = find();

figure(5)
montage(reshape(IMG(:,:,Ind), [size(IMG,1) size(IMG,2) 1 nPeople]));


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

