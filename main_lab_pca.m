% This is a MATLAB script for the
% CLPS1291 lab on PCA
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
% Author: Thomas Serre
% Brown University
% CLPS Department
% email: Thomas_Serre@Brown.edu
% Website: http://serre-lab.clps.brown.edu
% March 2014;

close all

set(0,'DefaultAxesFontSize', 24)

load ../Data/mit_cmu_faces.mat

% Good practice to center the data first (ie to subtract off the mean)

% run PCA
[PC, score, eigenvalues, tsquared, explained] = pca(DATA);
 
% you can see that beyond the first say 100 components
% not much is left to be explained in terms of variance
% ie most of the information is contained in the first 100 components
% let's visualize them


% try to mirror the faces to make the PCs more symetrical


% project the faces on some of the individual components and
% return the largest projections 


% plot the eigenvalues and the variance explained by the 
% components both inidvidually and cumulatively


% Project original faces onto first n components 
% then project back in the original space 
% show original vs. compressed version






