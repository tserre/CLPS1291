% This is a MATLAB script for the 
% CLPS1291 lab on neural decoding ** part 1 **

% Other m-files required: DecodingToolbox
% Subfunctions: none
% MAT-files required: none
% Author: Carl Olsson & Thomas Serre 
% Brown University
% CLPS Department
% email: Thomas_Serre@Brown.edu
% Website: http://serre-lab.clps.brown.edu
% April 2014; 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Neural data and details for the competition can be                    %%
% obtained at                                                           %%
% http://compneuro.clps.brown.edu/neural-decoding-competition-2014/     %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The goal for this first lab is to get you familiarized with the
% analysis of neural data. Most of this first last will be spent
% on visualizing ERP data. 

% Below is an example on how to load the data. We are giving you
% data from 4 human Ss 004, 006, 007 and 008
clear all, close all, clc
subjectNum = '004';
trainingDirectory = '../Data/competitionRelease/train/';
load(fullfile(trainingDirectory,strcat('subject',subjectNum,'_ECOG')));

%% EXERCISE 1
% For each subject, create a plot which shows ERPs for each electrode 
% -- make sure to properly label the axes -- use EEG.times
% Also make sure to have a suptitle for each figure and
% individual subplot titles indicating the electrode location and
% the gyrus



%% EXERCISE 2
% Building on the previous question, try to organize your
% electrodes per gyrus



%% EXERCISE 3
% Building on the previous question, plot separate ERPs for
% -- animal vs. non-animal trials
% -- fast vs. slow trials
% -- mask vs. no mask



