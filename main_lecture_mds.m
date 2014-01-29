% This is a MATLAB script for the 
% CLPS1291 lecture on MDS 

% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
% Author: Thomas Serre 
% Brown University
% CLPS Department
% email: Thomas_Serre@Brown.edu
% Website: http://serre-lab.clps.brown.edu
% January 2014; 

% source: http://www.mathworks.com/help/stats/multidimensional-scaling.html
% List of city names (strings) stored as a cell array 
cities = {'Atl','Chi','Den','Hou','LA','Mia','NYC','SF','Sea','WDC'};

% D is our dissimilarity matrix corresponding to the
% distance matrix containing the pairwise distance 
% between all pairs of cities
% The matrix needs to be symmetric
D = [    0  587 1212  701 1936  604  748 2139 2182   543;
       587    0  920  940 1745 1188  713 1858 1737   597;
      1212  920    0  879  831 1726 1631  949 1021  1494;
       701  940  879    0 1374  968 1420 1645 1891  1220;
      1936 1745  831 1374    0 2339 2451  347  959  2300;
       604 1188 1726  968 2339    0 1092 2594 2734   923;
       748  713 1631 1420 2451 1092    0 2571 2408   205;
      2139 1858  949 1645  347 2594 2571    0  678  2442;
      2182 1737 1021 1891  959 2734 2408  678    0  2329;
       543  597 1494 1220 2300  923  205 2442 2329     0];
   
% Below is the command for metric mds
% Matlab returns the coordinates of the
% datapoints in the the smallest space in which 
% the n points whose interpoint distances are
% given by D can be embedded.

Y = cmdscale(D);

% Using the plot command to plot the data..
plot(Y(:,1),Y(:,2),'o')

% The text command allows you to add some text at location (x,y)
% Here the x location is shifted by 25 to
% help for readability

text(Y(:,1)+25,Y(:,2),cities)

%% Label the two axes
xlabel('Miles')
ylabel('Miles')
