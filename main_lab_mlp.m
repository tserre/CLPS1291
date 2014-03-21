% This is a MATLAB script for the
% CLPS1291 lecture on mlps
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

problem = 'faces';

switch problem
    case 'toy'
        % Generate synthetic / toy data
        n   = 100; % number of training samples
        Xtr = linspace(1,10,n);
        Ytr = [-ones(floor(length(Xtr)/2),1); ones(ceil(length(Xtr)/2),1)]';
        Xte = linspace(1.5,10,n);
        Yte = [-ones(floor(length(Xte)/2),1); ones(ceil(length(Xte)/2),1)]';
        
    case 'toy2'
        n   = 100; % number of training samples
        Xtr = [rand(n,1)*5; 5+rand(n,1)*5]';
        Ytr = sin(Xtr);
        Xte = [rand(n,1)*5; 5+rand(n,1)*5]';
        Yte = sin(Xte);
        
    case 'mnist'
        load '../Data/mnist_all.mat'
        A   = cat(1, train0, train1, train2, train3, train4, train5, ...
            train6, train7, train8, train9);
        Xtr  = double(A')/255;
        AVG  = mean(Xtr,2);
        Xtr  = Xtr - repmat(AVG, 1, size(Xtr,2));
        
        % run PCA
        [PC, score, eigenvalues, tsquared, explained] = pca(Xtr');
        Xtr = score(:,1:100)';
        
        Ytr = cat(1, 0*ones(size(train0,1),1), 1*ones(size(train1,1),1), ...
            2*ones(size(train2,1),1), 3*ones(size(train3,1),1), ...
            4*ones(size(train4,1),1), 5*ones(size(train5,1),1), ...
            6*ones(size(train6,1),1), 7*ones(size(train7,1),1), ...
            8*ones(size(train8,1),1), 9*ones(size(train9,1),1))';
        
        
        
        %         I   = randperm(size(A,1));
        %         Xtr = Xtr(I(1:1000),:);
        %         Ytr = Ytr(I(1:1000));
        %         siz = sqrt(size(Xtr,2));
        
        A   = cat(1, test0, test1, test2, test3, test4, test5, ...
            test6, test7, test8, test9);
        Xte  = double(A')/255;
        Xte = Xte - repmat(AVG, 1, size(Xte,2));
        Xte = PC(:,1:100)'*Xte;
        
        Yte = cat(1, 0*ones(size(test0,1),1), 1*ones(size(test1,1),1), ...
            2*ones(size(test2,1),1), 3*ones(size(test3,1),1), ...
            4*ones(size(test4,1),1), 5*ones(size(test5,1),1), ...
            6*ones(size(test6,1),1), 7*ones(size(test7,1),1), ...
            8*ones(size(test8,1),1), 9*ones(size(test9,1),1))';
        
        
        
    case 'faces'
        load '../Data/faces2.mat'
        Xtr = [train_faces; train_nonfaces]';
        Xte = [test_faces; test_nonfaces]';
        
        Ytr = [ones(size(train_faces,1),1); -ones(size(train_nonfaces,1),1)]'; 
        Yte = [ones(size(test_faces,1),1); -ones(size(test_nonfaces,1),1)]'; 
        
end


% start with h  = [1]; This is essentially implementing a simple
% perceptron
% Then try to make the network more complicated by adding more
% hidden units or layers

h = [20];
net.performParam.regularization = 0;

net = feedforwardnet(h); % initialize the network
net.layers{1}.transferFcn = 'tansig';
net = train(net,Xtr,Ytr); % train the network
view(net) % view the network

% compute the output of the network for every sample in your
% training dataset (this is your prediction)
predYtr = net(Xtr);
perfTr = perform(net,predYtr,Ytr); % compute accuracy = error rate

% compute the output of the network for every sample in your
% test dataset and corresponding accuracy
predYte = sim(net, Xte);
perfTe  = perform(net,predYte,Yte)

disp(['Acuracy Training: ' num2str(perfTr) ' %'])
disp(['Acuracy Test: ' num2str(perfTe) ' %'])

plot(Xtr, Ytr, 'pr', Xtr, predYtr,'+r', Xte, Yte, 'kp', Xte, predYte, 'k+');
legend('training', 'training prediction', 'test', 'test prediction')