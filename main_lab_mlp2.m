% This is a MATLAB script for the
% CLPS1291 lab on mlps #2
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
% Author: Thomas Serre
% Brown University
% CLPS Department
% email: Thomas_Serre@Brown.edu
% Website: http://serre-lab.clps.brown.edu
% March 2014;

% Type 'help nndemos' to get a list of builtin demos
close all;
nTr     = 10;
nTe     = 1000;

PCA = 0;

load '../Data/mnist_all.mat'

A = [];
Y = [];

% load data
for ii = 0:9
    Atmp       = eval(['train' num2str(ii)]);
    A          = cat(1, A, Atmp);
    n          = size(Atmp,1);
    
    % create labels:
    % (1, 0, ..., 0) for class 1 (digit 0)
    % (0, 1, ..., 0) for class 2 (digit 1)
    % etc
    
    Ytmp       = zeros(1,10);
    Ytmp(ii+1) = 1;
    Y          = cat(1, Y, repmat(Ytmp, n, 1));
end
A   = double(A')/255;
siz = sqrt(size(A,2));
Y   = Y';

if PCA
    Ncomp = 100;
    AVG   = mean(A,2);
    A     = A - repmat(AVG, 1, size(A,2));
    
    [PC, score, eigenvalues, tsquared, explained] = pca(A');
    
    A = score(:,1:Ncomp)';
end

I   = randperm(size(A,2));
Xtr = A(:,I(1:nTr));
Ytr = Y(:,I(1:nTr));

Xte = A(:,I(nTr+1:nTr+nTe));
Yte = Y(:,I(nTr+1:nTr+nTe));


%% Train MLP

h = [];
net.performParam.regularization = 0.5;

net = feedforwardnet(h); % initialize the network
% % net.layers{1}.transferFcn = 'tansig';
% % net.trainParam.epochs = 50;
% % net.trainParam.goal = 1e-5;
% % net.trainParam.time = 60;
% % net.trainParam.showCommandLine = 1;
% % net.trainParam.show = 1;
net = configure(net, Xtr, Ytr);
tic
net = train(net, Xtr, Ytr); % train the network
toc

%% Compute the training and test error
%  Note: that the perform function does not return the
%  classification error
%  Note: you can also independently use the plotconfusion function


%% Visualize hidden units


%% Visualize errors


%% try parrallel toolbox

