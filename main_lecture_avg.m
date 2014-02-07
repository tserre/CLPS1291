% This is a MATLAB script for the 
% CLPS1520 lecture on color images 

% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
% Author: Thomas Serre 
% Brown University
% CLPS Department
% email: Thomas_Serre@Brown.edu
% Website: http://serre-lab.clps.brown.edu
% February 2014; 

%% Learn the mean in 1D
clear all;

% Initialize variables
x = 2+5*randn(1,10);
w   = 0;
n   = length(x);
ind = randperm(n); %% randomize the presentation order

% True mean 
m   = mean(x); %% sample mean for comparison

% Algorithm starts
for t = 1:n
    plot(x, zeros(size(x)), 'p', x(ind(t)), 0, 'pr', w, 0, 'sm', m, 0, 'or', 'MarkerSize', 15);
    legend('data', 'current', 'w_{pre}', 'true mean m')
    pause;
    
    dw = 1/t * (x(ind(t))-w); %% compute update
    w  = w + dw; %% update w
    
    plot(x, zeros(size(x)), 'p', x(ind(t)), 0, 'pr', w, 0, 'sm', m, 0, 'or', 'MarkerSize', 15);
    legend('data', 'current', 'w_{post}', 'true mean m')
    pause;
end


%% Learn the mean in 2D
clear all;

% Initialize variables
x   = 2+5*randn(2,30);
w   = zeros(2,1);
n   = size(x,2);
ind = randperm(n); %% randomize the presentation order

% True mean 
m   = mean(x,2); %% sample mean for comparison

% Algorithm starts
for t = 1:n
    plot(x(1,:), x(2,:), 'p', x(1,ind(t)), x(2,ind(t)), 'pr', w(1), w(2), 'sm', m(1), m(2), 'or', 'MarkerSize', 15);
    legend('data', 'current', 'w_{pre}', 'true mean m')
    pause;
    
    dw = 1/t * (x(:,ind(t))-w); %% compute update
    w  = w + dw; %% update w
    
    plot(x(1,:), x(2,:), 'p', x(1,ind(t)), x(2,ind(t)), 'pr', w(1), w(2), 'sm', m(1), m(2), 'or', 'MarkerSize', 15);
    legend('data', 'current', 'w_{post}', 'true mean m')
    pause;
end



