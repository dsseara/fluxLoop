%% JFilament workflow
% Script to analyze JFilament data to obtain probability flux loop information

% First load the data you want. Will have to manually
% enter save folder name and table name

folder = '/media/daniel/storage1/Dropbox/LLM_Danny/FactinOrder/JFilament/ctrl';
snakes = ctrlsnakes;

%% 
% Parse table of data into a cell array
filament = loadJFilamentData(snakes);

%%
% Let's plot up the filaments so we can take a look at them.

px2um= .108; % 0.108 microns per pixel
maxFrame = 75; % 75 is absolute max frame possible
savestuff = 0; % boolean to save stuff or not
n = 3; % plot every nth frame

plotJFilamentData(filament, px2um,maxFrame, savestuff, folder,n);

%%
% Now we rotate the filaments so they lie on a horizontal axis (defined as the
% line between the first and last point in the first frame), and get the
% tangent information. Use this information to first get the elastohydrodynamic
% modes of slender rods, but also store in order to use PCA on it afterwards

px2um= .108; % 0.108 microns per pixel
bcs = 'free';   % Use free boundary conditions
nmax = 6;       % get the first 6 modes
savestuff = 1;  % Save the plots?
n = 3; % plot every nth frame
aa = getJFilamentModes(filament,nmax,bcs,maxFrame,px2um,savestuff,folder,n);

%%
% See how the modes correlate with each other, expect all the even
% modes and odd modes to correlate linearly with other even or odd modes,
% respectively, showing variation is mainly an even or oddness issue.

savestuff = 1;
plotJFilamentModeCorrelations(aa,savestuff, folder)

%% Now perform PCA on the mode coefficients
[coeff,score,latent,tsquared,explained,mu] = pca(aa);

%%
% Check out different modes to see which give a good analysis

for ii=1:5
    for jj = ii+1:6
        figure
        cutoff = [];
        axis1 = ii;
        axis2 = jj;

        dbin = (max(aa(:,axis2)) - min(aa(:,axis2)))/25;
        dt = 5;

        [probMat,fluxField,xEdges,yEdges] =...
            probabilityFlux([aa(:,axis1),aa(:,axis2)],dt,dbin,cutoff);
        pcolor(xEdges(1:end-1), yEdges(1:end-1),probMat), hold on
        quiver(xEdges(1:end-1)+dbin/2, yEdges(1:end-1)+dbin/2,...
                fluxField(:,:,1),fluxField(:,:,2),'w')
        colorbar
        xlabel(['Mode' num2str(axis1)])
        ylabel(['Mode' num2str(axis2)])
        axis equal
        saveas(gcf,[folder filesep 'phaseSpace' filesep...
            'modes_' num2str(axis1) num2str(axis2)],'tif');
    end
end

%%
% Vary the bin size and see if that affects the total curl found.
dt = 5;
cutoff = [];

axis1 = 5;
axis2 = 6;

nbinArray = [1:10].*2;
dbinArray = zeros(size(nbinArray));
stdArray = [1,2,3];

curlArray = zeros(numel(dbinArray),numel(stdArray)+1);

for ii =1:numel(nbinArray)
	dbin = (max(aa(:,axis2)) - min(aa(:,axis2)))/nbinArray(ii);
	dbinArray(ii) = dbin;
	for jj =  1:numel(stdArray)+1
		if jj == numel(stdArray)+1
			figure
			[probMat,fluxField,xEdges,yEdges] =...
                probabilityFlux([aa(:,axis1),aa(:,axis2)],dt,dbin,cutoff);
			pcolor(xEdges(1:end-1), yEdges(1:end-1),probMat), hold on
			quiver(xEdges(1:end-1)+dbin/2, yEdges(1:end-1)+dbin/2,...
					fluxField(:,:,1),fluxField(:,:,2),'w')
			colorbar
			xlabel(['Mode ' num2str(axis1)])
            ylabel(['Mode ' num2str(axis2)])
			title(['All curl, dbin = ',num2str(dbin)])
			axis equal
			plotbox = 0;
			nstd = [];
			curlArray(ii,jj) =...
                fluxLoopCurl([aa(:,axis1),aa(:,axis2)], probMat,fluxField,...
                xEdges,yEdges, nstd, plotbox, dbin);
             saveas(gcf,['varyBinSize/dbin_',strrep(num2str(dbin),'.',',')],'fig')
             saveas(gcf,['varyBinSize/dbin_',strrep(num2str(dbin),'.',',')],'tif')
             saveas(gcf,['varyBinSize/dbin_',strrep(num2str(dbin),'.',',')],'epsc')
        else
			[probMat,fluxField,xEdges,yEdges] =...
                probabilityFlux([aa(:,axis1),aa(:,axis2)],dt,dbin,cutoff);
			plotbox = 0;
			nstd = stdArray(jj);
			curlArray(ii,jj) = ...
                fluxLoopCurl([aa(:,axis1),aa(:,axis2)], probMat,fluxField,...
                xEdges,yEdges, nstd, plotbox, dbin);
			if plotbox
				title(['Inner curl, nstd = ',num2str(nstd), ' dbin = ', num2str(dbin)])
				axis equal
			end
		end
	end
end

%%
% Start the bootstrapping procedure
% bootStrappingLoop.m varies the bin size and area over which the
% curl is integrated in order to find out what parameter set gives
% us a statistically significant circulation in the curl

axis1 = 5;
axis2 = 6;

ctrlfolder =...
    '/media/daniel/storage1/Dropbox/LLM_Danny/FactinOrder/JFilament/ctrl';
activefolder =...
    '/media/daniel/storage1/Dropbox/LLM_Danny/FactinOrder/JFilament/active';
load([ctrlfolder filesep,'ctrlData.mat'])
ctrlTseries = [aa(:,axis1),aa(:,axis2)];
load([activefolder filesep,'activeData.mat'])
activeTseries = [aa(:,axis1),aa(:,axis2)];
% Vary the bin size and see if that affects the total curl found.
dt = 5;
cutoff = [];
plotstuff = 1;
m=500;

nbinArray = [1:10].*2; % Array of number of bins to change bin size
dbinArray = zeros(size(nbinArray)); % array to keep 
stdArray = [1,2,3]; % number of stds to integrate over also does the entire
					% domain in addition to these

[curlNormed, curlHistActive, curlHistCtrl] =...
    bootStrappingLoop(activeTseries, ctrlTseries, dt, cutoff, plotbox, m,...
    nbinArray, stdArray,plotstuff);