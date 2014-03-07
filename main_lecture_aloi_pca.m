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
% March 2014;

%% Should point to your local copy of the COIL-100
% available at http://www.cs.columbia.edu/CAVE/databases/SLAM_coil-20_coil-100/coil-100/coil-100.zip

d = dir('../Data/coil-100/*.png');

nsiz   = 32;

DATA = zeros(length(d),nsiz^2);

for ii = 1:length(d)
    img = double(imread(fullfile('../Data/coil-100', d(ii).name)))/255;
    img2 = histeq(imresize(rgb2gray(img), [nsiz nsiz]));
    DATA(ii,:) = img2(:);
%     imagesc(img);
%     pause(.1)
end

AVG  = mean(DATA);
DATA = DATA - repmat(AVG, size(DATA,1), 1);


[PC, score, eigenvalues, tsquared, explained] = pca(DATA);

figure(1)
for ii =1:25
    subplot(5,5,ii), imagesc(reshape(PC(:,ii), [nsiz nsiz]));
    axis off; axis image; colormap gray
end



%%
load ../Data/mit_cmu_faces.mat
nsiz = 19;

M     = reshape([1:361], 19, 19);
M     = fliplr(M);
DATA2 = DATA(:,M(:));

DATA  = [DATA; DATA2];
  
AVG  = mean(DATA);
DATA = DATA - repmat(AVG, size(DATA,1), 1);

% run PCA
[PC, score, eigenvalues, tsquared, explained] = pca(DATA);
siz = sqrt(size(PC,1));
 
%
colormap gray
n = 25;
s = ceil(sqrt(n));


figure(1)
for ii = 1:n
    subplot(s,s,ii); imagesc(reshape(PC(:,ii), siz, siz)/255); axis off;
    
end
%% 
ind = randperm(size(DATA,1));
ind = ind(1:50);
my_score = score(ind,:);

figure(5)
n1 = 1;
n2 = 2;

plot(my_score(:,n1), my_score(:,n2), '.'); % here change for pc3, or 4 etc
xlabel('1st component')
ylabel('2nd component')
title('Images organized by pca')

xa = my_score(:,n1) - min(my_score(:,n1));
ya = my_score(:,n2) - min(my_score(:,n2));

scaling = 2;

S = max(xa);
xa = .925*xa/S+0.025;
S = max(ya);
ya = .925*ya/S+0.025;

figure(6); colormap gray;
for n = length(xa):-1:1
    h=axes('position', [xa(n) ya(n) .05*scaling .05*scaling]);
    imagesc(reshape(DATA(ind(n),:), [nsiz nsiz]), 'parent', h)
    axis('off'); axis('equal')
end

%%
kk = 200;
my_face = DATA(kk,:); 
figure(8); colormap gray;
imagesc(reshape(my_face, nsiz, nsiz));

my_pface = score(kk,1:100);
my_ref_pfaces = my_score(:,1:100);

D = pdist2(my_pface, my_ref_pfaces);


