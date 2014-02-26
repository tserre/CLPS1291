% This is a MATLAB script for the
% CLPS1291 lecture on random variables

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

%% 1D
figure(1)

n = 100; % number of sample points
M  = 0;    % mean
S  = 1;    % std

% Use normpdf to compute the pdf of a Gaussian random variable
% with mean M and std S

x1 = -10:.2:10;
Y  = normpdf(x1, M, S);

% Sample random points from a normal distribution with mean M and
% std S (remember that randn samples from a zero mean and unit
% standard deviation but you can shift the mean by adding a
% constant M and change the std by multiplying the points by S

X  = M + S*randn(n,1);

% Normfit estimates the mean M1 and std S1 of the normal
% distribution given the data X

[M1, S1] = normfit(X);

% Plot the corresponding pdf
Y1 = 1/(sqrt(2*pi)*S1)*exp(-(x1-M1).^2/(2*S1^2));

plot(x1, Y, 'r', x1, Y1, 'b--', X , 0,'.b')

legend('true pdf', 'estimated pdf', 'random samples')
xlabel('x'); ylabel('pdf');

%% 2D

close all

n = 500; % number of sample points

M = [0 0];     % mean
S = [5 0;
     0 1]; % covariance matrix

x1 = -3:.2:3; 
x2 = -3:.2:3;

[X1, X2] = meshgrid(x1,x2); % creates a cartesian grid

% you can see what meshgrid does by plotting X1 vs X2
figure(1)
plot(X1(:), X2(:), '.');

F = mvnpdf([X1(:) X2(:)], M, S); % multivariate pdf
F = reshape(F,length(x2),length(x1)); % reshape to plot as 3D function
xlabel('x1'); ylabel('x2'); zlabel('See what meshgrid does');

% what the pdf should be
figure(2), surf(X1,X2,F);

xlabel('x1'); ylabel('x2'); title('Probability Density Function');

% sample points
X = mvnrnd(M, S, n); % sample points

% estimated mean and cov
M1 = mean(X);
S1 = cov(X);

disp('True mean and cov:');
M
S

disp('Estimated mean and cov:');
M1
S1


% Iso probability contours
figure(3), plot(X(:,1), X(:,2),'+')
% axis([-5 5 -5 5])

hold on;

x1    = linspace(0,2*pi, 50)';
x     = cos(x1); y=sin(x1);
ap    = [x(:) y(:)]';
[v,d] = eig(S);

bp    = (v*d*ap) + repmat(M', 1, size(ap,2)); 
plot(bp(1,:), bp(2,:), '--r');
hold off;

xlabel('x1'); ylabel('x2'); title('Datapoints sampled from normal distribution');
%% image

close all

n = 500; % number of sample points

img = double(rgb2gray(imread('../Data/images/pippin_jtalon0033.tif')))/255;

figure(1), imshow(img);

figure(2), hist(img(:), [0:0.05:1]);
title('Pixel histogram');
xlabel('Intensities')
ylabel('Count');

sh_img = circshift(img, [0 10]);
I      = randperm(prod(size(img)));
X1     = img(I(1:n));
X2     = sh_img(I(1:n));
X      = [X1; X2]';

disp('True mean and cov:')
% estimated mean and cov
M = mean(X)
S = cov(X)

figure(3), plot(X1, X2, '.');
xlabel('x1'); ylabel('x2'); title('Sample pixel pairs');

%% 3D

close all

n = 500; % number of sample points

M = [0 0 0];     % mean
S = [1 0 0;0 1 0; 0 0 1]; % covariance matrix

figure(1)

% sample points
X = mvnrnd(M, S, n);

[x,y,z] = sphere(50);
ap = [x(:) y(:) z(:)]';
[v,d]=eig(S);
if any(d(:) < 0)
    fprintf('warning: negative eigenvalues\n');
    d = max(d,0);
end

bp = (v*d*ap) + repmat(M', 1, size(ap,2));
xp = reshape(bp(1,:), size(x));
yp = reshape(bp(2,:), size(y));
zp = reshape(bp(3,:), size(z));

plot3(X(:,1), X(:,2), X(:,3), '.');
hold on;
h = surf(xp,yp,zp); alpha(.3)
axis([-2 2 -2 2 -2 2])

