%% 
% April 2019
% By: Alexander Fyrdahl (alexander.fyrdahl@ki.se)
%
% Contribution to the reproducible research study group initiative to reproduce [1]
%
% [1] Pruessmann, K. P.; Weiger, M.; Boernert, P. and Boesiger, P.
% Advances in sensitivity encoding with arbitrary k-space trajectories.
% Magn Reson Med 46: 638-651 (2001)

clear variables; clc; close all; 
addpath(fullfile(pwd,'nufft_files'));
addpath(fullfile(pwd,'siv'));

%% Load data

rawdata_real    = h5read('rawdata_heart_radial_55proj_34ch.h5','/rawdata');
trajectory      = h5read('rawdata_heart_radial_55proj_34ch.h5','/trajectory');

rawdata = double(rawdata_real.r+1i*rawdata_real.i); clear rawdata_real;
rawdata = permute(rawdata,[3,2,1]);
trajectory = double(permute(trajectory,[3,2,1]));
imSize  = round([2 2].*max(trajectory(:)));
[nFE,nSpokes,nCh] = size(rawdata);

k = squeeze(complex(trajectory(1,:,:),trajectory(2,:,:)));
k = k/(2*max(abs(k(:))));

FT = NUFFT(k,abs(k),imSize,nCh);
img_grid = FT'*rawdata;
img_igrid_sos = sqrt(sum(abs(img_grid.^2),3));
csm = img_grid./repmat(img_igrid_sos,[1 1 nCh]);

%% Subsampling

cardiac_55 = do_sense_recon(rawdata,k,csm,3);
cardiac_33 = do_sense_recon(rawdata(:,1:33,:),k(:,1:33),csm,4);
cardiac_22 = do_sense_recon(rawdata(:,1:22,:),k(:,1:22),csm,5);
cardiac_11 = do_sense_recon(rawdata(:,1:11,:),k(:,1:11),csm,8);

%% Reproduce figure 6
figure(6);

imshow(flipud(abs(cat(2,cardiac_55,cardiac_33,cardiac_22,cardiac_11))),[],'border','tight');
text(210,225,'55','Color','w','FontSize',20);
text(210+240*1,225,'33','Color','w','FontSize',20);
text(210+240*2,225,'22','Color','w','FontSize',20);
text(210+240*3,225,'11','Color','w','FontSize',20);

set(gcf,'PaperPositionMode','auto');
pos = get(gcf,'Position');
pos(3:4) = [imSize(1)*4 imSize(2)];
set(gcf,'Position',pos)

print figure6.png -dpng