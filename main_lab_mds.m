% This is a MATLAB script for the 
% CLPS1291 lecture on MDS 

% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
% Author: Thomas Serre 
% Data source: Data were graciously provided
% by Michael Lee and can be downloaded at
% http://faculty.sites.uci.edu/mdlee/similarity-data/
% Brown University
% CLPS Department
% email: Thomas_Serre@Brown.edu
% Website: http://serre-lab.clps.brown.edu
% February 2014; 

% Each of these .mat files contains 4 variables: n = number of objects, labs = string array of labels for the objects, s = n x n symmetric matrix of normalised pairwise similarities between the objects, d = n x n symmetric matrix of normalized pairwise proximities between the objects:
%
% Human judgments of the numbers 0-9 [abstractnumbers.mat]. From research described in Shepard, R. N., Kilpatrick, D. W., & Cunningham, J. P. (1975). The internal representation of numbers. Cognitive Psychology, 7, 82-138 (with thanks to Josh Tenenbaum).
% Auditory confusions of 25 letters (all excluding ?o?) and the numbers 0-9 [auditory.mat]. From research reported in Kuennapas, T., & Janson, A-J. (1969). Multidimensional Similarity of Letters. Perceptual and Motor Skills, 28, 3-12.
% A sociologist?s judgment of the relationships between 14 bank wiring workers [bankwiring.mat]. From research reported in Roethlisberger, F. J., & Dickson, W. J. (1939). Management and the worker. Cambridge, MA: Harvard University Press.
% Voting patterns of 14 members of congress on environmental bills [congress.mat]. From raw data presented in Romesburg, H. C. (1984). Cluster analysis for researchers. Belmont, CA: Lifetime Learning Publications.
% Human judgments of 17 dot patterns [dotpatterns.mat]. From research reported in Glushko, R. J. (1975). Pattern goodness and redundancy revisited: Multidimensional scaling and hierarchical cluster analysis. Perception & Psychophysics, 17(2), 158-162.
% Reported adolescent use of 13 drug types [druguse.mat]. From research reported in Huba, G. L., Wingard, J. A., & Bentler, P. M. (1981). A comparison of two latent variable causal models for adolescent drug use. Journal of Personality and Social Psychology, 40(1), 180-193.
% Human judgments of 16 drawings of flowerpots [flowerpots.mat]. From research reported in Gati, I., & Tversky, A. (1982). Representations of qualitative and quantitative dimensions. Journal of Experimental Psychology: Human Perception and Performance, 8(2), 325-340.
% Human judgments of 21 fruits [fruits.mat]. From research reported in Tversky, A., & Hutchinson, J. W. (1986). Nearest Neighbor Analysis of Psychological Spaces. Psychological Review, 93(1), 3-22.
% Kindergarten children?s judgment of perceptual similarity of the 26 capital letters [letters.mat]. From research reported in Gibson, E. J., Osser, H., Schiff, W., & Smith, J. (1963). An analysis of critical features of letters, tested by a confusion matrix. Cooperative Research Project No. 639, U.S. Office of Education.
% Confusion of Morse code numerals [morsenumbers.mat] and numeral and letters [morseall.mat]. From research reported in Rothkopf, E. Z. (1957). A measure of stimulus similarity and errors in some paired-associate learning tasks. Journal of Experimental Psychology, 53, 94-101.
% Auditory confusion of 16 consonant phonemes [phonemes.mat]. From research reported in Miller, G. A., & Nicely, P. E. (1955). An analysis of perceptual confusions among some English consonants. Journal of the Acoustical Society of America, 27, 338-352.
% Human judgments of 18 risks [risks.mat]. From research reported in Johnson, E. J., & Tversky, A. (1984). Representations of Perceptions of Risks. Journal of Experimental Psychology: General, 113(1), 55-70.
% Human judgments of 16 rectangles [rectangles.mat]. From research described in Chapter 15 of Borg, I., & Lingoes, J. (1987). Multidimensional similarity structure analysis. New York: Springer Verlag.
% The following .mat files also contain a variable sigma_emp, which gives an empirical estimate of the precision of the similarity and proximity data:

% Human judgments (in 1967) of 17 countries [country_robinsonhefner.mat]. From research reported in Robinson, J. P., & Hefner, R (1967). Multidimensional Differences in Public and Academic Perceptions of Nations. Journal of Personality and Social Psychology, 7(3), 251-259.
% Human judgments of 8 rectangles with interior line segments [rectangles_kruschke.mat]. From research reported in Kruschke, J. K. (1993). Human category learning: Implications for backpropagation models. Connection Science, 5, 3-36.
% Human judgments of 15 kinship terms [kinship_rosenbergkim.mat]. From research reported in Rosenberg, S., & Kim, M. P. (1975). The Method of Sorting as a Data-Generating Procedure in Multivariate Research. Multivariate Behavioral Research, 10, 489-502.
% Human judgments of 21 bird names [birds_romney.mat], 21 clothing names [clothing_romney.mat], 21 different clothing names [clothing2_romney.mat], 21 fish names [fish_romney.mat], 21 fruit names [fruit_romney.mat], 21 different fruit names [fruit2_romney.mat], 21 furniture names [furniture_romney.mat], 21 different furniture names [furniture2_romney.mat], 21 semantically unrelated words [nonsense_romney.mat], 21 sport names [sport_romney.mat], 21 tool names [tools_romney.mat], 21 toy names [toys_romney.mat], 21 vegetable names [vegetables_romney.mat], 21 different vegetable names [vegetables2_romney.mat], 21 vehicle names [vehicles_romney.mat], 21 different vehicle names [vehicles2_romney.mat], 21 weapon names [weapons_romney.mat], 21 different weapon names [weapons2_romney.mat]. All from research reported in Romney, A. K., Brewer, D. D., & Batchelder, W. H. (1993). Predicting Clustering from Semantic Structure. Psychological Science, 4(1), 28-34, with thanks to Devon Brewer.
% Human judgments of 9 lines of different lengths [lines_cohen.mat], 60 faces [faces_busey.mat], 7 ?morphed? faces [faces_steyvers.mat], 9 shapes varying in size and angle [sizeangle_treat.mat], 24 bodies varying in ?affect and body size? [bodies_viken.mat]. Mark Steyvers kindly provided me with all of these, and I have yet to chase up references (although the filenames ought to make that pretty easy).
% Human judgments of 30 Brodatz textures [texturebrodatz_heaps.mat], and 24 MIT textures [texturemit_heaps.mat]. Both from research reported in Heaps, C., & Handel, S. (1999). Similarity and Features of Natural Textures. Journal of Experimental Psychology: Human Perception and Performance, 25(2), 299-320.
% Human judgments of 10 cartoon faces [cartoonfaces.mat], and forced-choice judgments of 16 countries in a similarity condition [countriessim.mat] and a dissimilarity condition [countriesdis.mat]. From the research described in Navarro, D.J., & Lee, M.D. (2004). Common and distinctive features in stimulus representation: A modified version of the contrast model. Psychonomic Bulletin & Review, 11(6), 961?974, and Navarro, D.J., & Lee, M.D. (2002). Commonalities and distinctions in featural stimulus representations. In W.G. Gray & C. D. Schunn, (Eds.), Proceedings of the 24th Annual Conference of the Cognitive Science Society, pp. 685-690. Mahwah, NJ: Erlbaum.
% Human judgments of 21 animals (presented as pictures on a 5 point scale) [animalpictures5.mat], of 21 animals (presented as pictures on a 5 point scale) [animalpictures5.mat], of 21 animals (presented as pictures on an 11 point scale) [animalpictures11.mat], of 21 animals (presented as words on a 5 point scale) [animalnames5.mat], of 21 animals (presented as words on an 11 point scale) [animalnames11.mat], of 25 faces (5 point scale) [faces5.mat], and of 25 faces (11 point scale) [faces11.mat], together with two bitmap files with the face stimuli [faces.bmp, faces2.bmp]. From (as yet; probably never-to-be) unreported research I did a while back.
% Human judgments of 24 sounds (with different similarity collection methodologies) [sounds_harbke.txt], kindly provided by Colin Harbke.

%% DATA: You need to download the data available at
% https://www.dropbox.com/s/afov03u64es4980/MDS_data.zip
% A link is also available on canvas

%% TOPIC DISCUSSED TODAY:
% - Organize your code vs. data
% - Shortcuts
% - Command line vs. script vs. function
% - Refresher on variables
% - Basic image visualization
% - The mdscale command
% - More on indexing ? ':'
% - Basic ploting 
% - The linkage and dendrogram functions

%% Initial cleanup
close all;
clear all;
clc;

%% Select the data you want to load and load the data
% data_file is in a format called char array The data have
% already been preformated for you in a format that matlab can
% understand ? next time we will talk about data import
%% Q: You will have to modify the line below to load different kinds of data
data_file = '../MDS_data/abstractnumbers.mat';
load(data_file)

%% Q: Use the whos command to find out what variables get loaded 

%% The data provided are a little inconsistent, the files all
% contain sim matrixes, some also contain dissim matrixes, some
% have dissimilarity measures but some don't ? just to be safe we
% will create/overwrite a dissimilarity matrix from similarity
% data
d = 1-s;

%% Q: Visualize the dissimilarity matrix using the 'imagesc' command
% You should also learn about the colormap and colorbar commands
% the axis command

%% Run non-metric MDS, we will start with 2D...
% The first argument that the function takes is a dissim (or sim
% matrix) the second one is the number of dimensions you want the
% mental space to be
Y = mdscale(d,2);

%% Q: Use the 'figure' command to start a new figure


%% Q: Use the 'plot' command to visualize the recovered mental space

%% Q: Label the axes using the 'xlabel' and 'ylabel' commands and
% add a title with the 'title' commands --  Remember to always label
% the figures axes and to always add a title

%% Here we call the 'hold on' command to tell matlab to hold on
% the current plot to add onto it
hold on;

% Q: Use the 'text' command to add labels to the datapoints on your plot at location (x,y)
% Shifted the x or y location by some small amount to improve readability

%% Here we call the 'hold off' command when we are done with the
% figure (just a good habit)
hold off;

%% Q: Open another figure

%% Q: Use the function pdist to recover the distances between
% points in mental space and plot these distances against the
% original similarity data Shepard'style


%% Q: Run hierarchical clutering on your dissim data using the
% 'linkage' function -- you will have to first format your dissim
% data using the 'squareform' command provided below.
D = squareform(d);
Z = linkage(D);
dendrogram(Z, 'labels', labs);

% save your figures using the print command
 
 
 
 
 
 
 