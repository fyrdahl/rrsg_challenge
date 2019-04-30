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

rawdata_real    = h5read('rawdata_brain_radial_96proj_12ch.h5','/rawdata');
trajectory      = h5read('rawdata_brain_radial_96proj_12ch.h5','/trajectory');

rawdata = double(rawdata_real.r+1i*rawdata_real.i); clear rawdata_real;
rawdata = permute(rawdata,[3,2,1]);
trajectory = double(permute(trajectory,[3,2,1]));
imSize  = [2 2].*max(trajectory(:));
[nFE,nSpokes,nCh] = size(rawdata);

k = squeeze(complex(trajectory(1,:,:),trajectory(2,:,:)));
k = k/(2*max(k(:)));

FT = NUFFT(k,abs(k),imSize,nCh);
img_grid = FT'*rawdata;
img_grid_sos = sqrt(sum(abs(img_grid.^2),3));

csm = img_grid./repmat(img_grid_sos,[1 1 nCh]);

%% Subsampling
img_sense_iter = [];
img_coils = [];
delta = [];
nIter = 100;

% Reference image
[~,sense_ref] = do_sense_recon(rawdata,k,csm,8);

for R = 2:4
    k_u = k(:,1:R:nSpokes);
    rawdata_u = rawdata(:,1:R:nSpokes,:);
    [~,a,c,d] = do_sense_recon(rawdata_u,k_u,csm,nIter);
    img_coils = cat(4,img_coils,c);
    img_sense_iter = cat(4,img_sense_iter,a);
    delta = cat(1,delta,d);
end

%% Reproduce figure 4
figure(4);
subplot(211);
for j = 1:3
    for i = 1:nIter
        Delta(j,i) = abs(sqrt(sum(col(img_sense_iter(:,:,i,j))-col(sense_ref(:,:,end)))^2)/sqrt(sum(col(sense_ref(:,:,end))))^2);
    end
end
plot(0:nIter-1,log10(Delta(3,:)),'k:','LineWidth',1.5); hold on;
plot(0:nIter-1,log10(Delta(2,:)),'k--','LineWidth',1.5);
plot(0:nIter-1,log10(Delta(1,:)),'k-','LineWidth',1.5); hold off;
legend({'R = 4','R = 3','R = 2'},'Box','off','Position',[0.5 0.8 0.05 0.075]);
ylabel('Log_{10} \Delta');
yticks([-3:0]);
ytickformat('%.1f');
set(gca,'FontSize',20,'LineWidth',1,'XMinorTick','on','YMinorTick','on');
axis([0 100 -3 0]);
axis square;

subplot(212);
plot(0:nIter-1,log10(delta(3,:)),'k:','LineWidth',1.5); hold on;
plot(0:nIter-1,log10(delta(2,:)),'k--','LineWidth',1.5);
plot(0:nIter-1,log10(delta(1,:)),'k-','LineWidth',1.5); hold off;
legend({'R = 4','R = 3','R = 2'},'Box','off','Position',[0.5 0.3 0.05 0.075]);
ylabel('Log_{10} \delta');
yticks([-5:0]);
ytickformat('%.1f');
set(gca,'FontSize',20,'LineWidth',1,'XMinorTick','on','YMinorTick','on');
axis([0 100 -5 0]);
axis square;

set(gcf,'PaperPositionMode','auto');
pos = get(gcf,'Position');
pos(3:4) = [800 1600];
set(gcf,'Position',pos);

print figure4.png -dpng

%% Reproduce figure 5
figure(5);

subplot(331);
imshow(abs(img_coils(:,:,2,1)),[],'Border','Tight');
title('Single coil','FontWeight','Normal','FontSize',20);
ylabel('2','Rotation',0,'FontSize',20);

subplot(332);
imshow(sqrt(sum(abs(img_coils(:,:,:,1).^2),3)),[],'Border','Tight');
text(275,275,'1','Color','w','FontSize',20);
title('Initial','FontWeight','Normal','FontSize',20);

subplot(333);
imshow(abs(img_sense_iter(:,:,6,1)),[],'Border','Tight');
text(275,275,'6','Color','w','FontSize',20);
title('Final','FontWeight','Normal','FontSize',20);

subplot(334);
imshow(abs(img_coils(:,:,2,2)),[],'Border','Tight');
ylabel('3','Rotation',0,'FontSize',20);

subplot(335);
imshow(sqrt(sum(abs(img_coils(:,:,:,2).^2),3)),[],'Border','Tight');
text(275,275,'1','Color','w','FontSize',20);

subplot(336);
imshow(abs(img_sense_iter(:,:,12,2)),[],'Border','Tight');
text(255,275,'12','Color','w','FontSize',20);

subplot(337);
imshow(abs(img_coils(:,:,2,3)),[],'Border','Tight');
ylabel('4','Rotation',0,'FontSize',20);

subplot(338);
imshow(sqrt(sum(abs(img_coils(:,:,:,3).^2),3)),[],'Border','Tight');
text(275,275,'1','Color','w','FontSize',20);

subplot(339);
imshow(abs(img_sense_iter(:,:,24,3)),[],'Border','Tight');
text(255,275,'24','Color','w','FontSize',20);

set(gcf,'PaperPositionMode','auto');
pos = get(gcf,'Position');
pos(3:4) = [800 800];
set(gcf,'Position',pos)

print figure5.png -dpng
