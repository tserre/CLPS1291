% This is a MATLAB script for the
% CLPS1291 lecture on PCA
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
% Author: Thomas Serre
% Brown University
% CLPS Department
% email: Thomas_Serre@Brown.edu
% Website: http://serre-lab.clps.brown.edu
% February 2014;

set(0,'DefaultAxesFontSize', 24)

% generate data
DATA = mvnrnd([0 0],[1 2;2 10], 50);

% center the data
DATA = DATA - repmat(mean(DATA), size(DATA,1), 1);
[PC,score,latent,tsquared,explained] = pca(DATA);

figure(1)

C = rand(size(DATA,1),3);

scatter(DATA(:,1), DATA(:,2), 20, C);
hold on; 
plot(PC(1,1)*[-15 15], PC(2,1)*[-15 15], '-r');
plot(PC(1,2)*[-15 15], PC(2,2)*[-15 15], '-b'); hold off
axis([-10 10 -10 10]); axis square
hold off

xlabel('x1'); ylabel('x2'); 
latent
explained

% project down to 1 dimension
pDATA = DATA * PC; 

figure(2)
scatter(pDATA(:,1), pDATA(:,2), 20, C);

% Drop one component and project back in the original space
ppDATA = pDATA(:,1)*PC(:,1)'; % or equivalently (PC'*DATA')'
figure(3)
scatter(ppDATA(:,1), ppDATA(:,2), 20, C);
