% This is a MATLAB script for the
% CLPS1291 lecture on perceptrons
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
% Author: Thomas Serre
% Brown University
% CLPS Department
% email: Thomas_Serre@Brown.edu
% Website: http://serre-lab.clps.brown.edu
% March 2014;

n   = 50;
dim = 2; % stick to 2D if you want to be able to visualize

% Generate toy data (normally distributed)
% you can try to make the classification problem harder by changing the
% distributions

sig = [.2 .2]; % width of the distribution
m   = [0 .7];
% input vectors (n points for each class, dim dimensional)
X   = [m(1)+sig(1)*randn(n,dim); m(2)+sig(2)*randn(n,dim)];
Y   = [ones(n,1); -ones(n,1)]; % associated labels

% create independent test set
Xte = [m(1)+sig(1)*randn(n,dim); m(2)+sig(2)*randn(n,dim)];
Yte = [ones(n,1); -ones(n,1)];

figure(1);

% initialization
w    = randn(1,dim+1); % bias trick -- add an additional dimension
eta0 = .50; % learning rate
k    = 0;

ErrTr = []; % training error
ErrTe = []; % test error
close all;

while(1)
    ind = randperm(size(X,1)); % shuffle training data presentation order
    
    for i = ind
        k = k+1; % iteration number
        
        % try to change the learning rate and comment
        % on what happens when the learning rate is very small or
        % very large
        
        eta = eta0;
        % eta = eta0/(1+(k)); %% decreasing learning rate as a
        % function of the iteration #
        % Good idea! Why?
        
        y_i = Y(i);
        x_i = [1 X(i,:)]'; % bias trick contd!
        
        % looking at the current output of the perceptron
        z = w*x_i; % linear part of the response
        y = sign(w*x_i); % Heavisde function - convert to -1/1 values
        
        
        figure(1)
        
        % showing all datapoints (red for pos class and blue
        % for neg
        plot(X(1:n,1), X(1:n,2),'or', X(n+1:2*n,1), X(n+1:2*n,2),'pb');
        hold on;
        
        % show current datapoint in black
        plot(X(i,1), X(i,2),'sk', 'MarkerSize', 12);
        
        % show decision boundary before update
        % remember equation is: w(1) * 1 + w(2) a + w(3) b = 0
        % where a is abscissa and b ordinate
        a = [-2 2]; %% abscissa
        b = (-w(1)-a*w(2))/w(3);
        plot(a, b, 'm-');
        axis([-2 2 -2 2]); title(num2str(k))
        pause(.1)
        
        % stop if there is an error
        if (y_i~=y)
            disp('error');
            
            % make weight update (nothing will actually happen
            % unlesss the perceptron makes an error (ie, y_i-y ~=0)
            dw = eta*(y_i-y)*x_i';
            w  = w+dw;
            
            % show the update
            b = (-w(1)-a*w(2))/w(3);
            
            % showing the datapoints
            plot(a, b,'m--');
            axis([-2 2 -2 2]); title(num2str(k))
            pause(.1)
            
        end
        hold off;
        
        % evaluate error on training set
        Pred  = sign(w*[ones(size(X,1), 1) X]');
        Cor   = (Pred == Y');
        ErrTr = [ErrTr 1-mean(Cor)];
        
        figure(2);
        plot([1:length(ErrTr)], ErrTr, '--m');
        legend('Training error')
    end
    % %     % end of epoch
    % %     % evaluate error on training set
    % %     Pred  = sign(w*[ones(size(X,1), 1) X]');
    % %     Cor   = (Pred == Y');
    % %     ErrTr = [ErrTr 1-mean(Cor)];
    % %
    % %     disp(['eta: ' num2str(eta) ' | Tr Err: ' num2str(1-mean(Cor))])
    % %
    % %
    % %
    % %     % evaluate error on Test set
    % %     Pred = sign(w*[ones(size(Xte,1), 1) Xte]');
    % %     Cor  = (Pred == Yte');
    % %     ErrTe  = [ErrTe 1-mean(Cor)];
    % %     figure(2);
    % %     plot([1:length(ErrTr)], ErrTr, '-b', [1:length(ErrTe)], ErrTe, '--m');
    % %     legend('Training error', 'Test error')
end
