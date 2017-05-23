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

plotJFilamentData(filament, px2um,maxFrame, savestuff, folder);

%%
% Now we rotate the filaments so they lie on a horizontal axis (defined as the
% line between the first and last point in the first frame), and get the
% tangent information. Use this information to first get the elastohydrodynamic
% modes of slender rods, but also store in order to use PCA on it afterwards

px2um= .108; % 0.108 microns per pixel
bcs = 'free';   % Use free boundary conditions
nmax = 6;       % get the first 6 modes
savestuff = 0;  % Save the plots?
aa = getJFilamentModes(filament,nmax,bcs,maxFrame,px2um,savestuff,folder);

%%
% See how the modes correlate with each other, expect all the even
% modes and odd modes to correlate linearly with other even or odd modes,
% respectively, showing variation is mainly an even or oddness issue.

savestuff = 0;
plotJFilamentModeCorrelations(aa,savestuff, folder)

%% Now perform PCA on the mode coefficients
[coeff,score,latent,tsquared,explained,mu] = pca(aa);

%%
% Now use the first two score columns (>99% explained) as the time series for 
% the probability flux
cutoff = [];
dbin = (max(score(:,2)) - min(score(:,2)))/25;
dt = 5;

[probMat,fluxField,xEdges,yEdges] = probabilityFlux(score(:,1:2),dt,dbin,cutoff);
pcolor(xEdges(1:end-1), yEdges(1:end-1),probMat), hold on
quiver(xEdges(1:end-1)+dbin/2, yEdges(1:end-1)+dbin/2,...
		fluxField(:,:,1),fluxField(:,:,2),'w')
colorbar
xlabel('pca 1')
ylabel('pca 2')
axis equal

%%
% Vary the bin size and see if that affects the total curl found.
dt = 5;
cutoff = [];

nbinArray = [1:10].*5;
dbinArray = zeros(size(nbinArray));
stdArray = [1,2,3];

curlArray = zeros(numel(dbinArray),numel(stdArray)+1);

for ii =1:numel(nbinArray)
	dbin = (max(score(:,2)) - min(score(:,2)))/nbinArray(ii);
	dbinArray(ii) = dbin;
	for jj =  1:numel(stdArray)+1
		if jj == numel(stdArray)+1
			figure
			[probMat,fluxField,xEdges,yEdges] = probabilityFlux(score(:,1:2),dt,dbin,cutoff);
			pcolor(xEdges(1:end-1), yEdges(1:end-1),probMat), hold on
			quiver(xEdges(1:end-1)+dbin/2, yEdges(1:end-1)+dbin/2,...
					fluxField(:,:,1),fluxField(:,:,2),'w')
			colorbar
			xlabel('pca 1')
			ylabel('pca 2')
			title(['All curl, dbin = ',num2str(dbin)])
			axis equal
			plotbox = 0;
			nstd = [];
			curlArray(ii,jj) = JFilamentFluxLoop(score(:,1:2), probMat,fluxField,xEdges,yEdges, nstd, plotbox, dbin);
            saveas(gcf,['varyBinSize/dbin_',strrep(num2str(dbin),'.',',')],'fig')
            saveas(gcf,['varyBinSize/dbin_',strrep(num2str(dbin),'.',',')],'tif')
            saveas(gcf,['varyBinSize/dbin_',strrep(num2str(dbin),'.',',')],'epsc')
        else
			[probMat,fluxField,xEdges,yEdges] = probabilityFlux(score(:,1:2),dt,dbin,cutoff);
			plotbox = 0;
			nstd = stdArray(jj);
			curlArray(ii,jj) = JFilamentFluxLoop(score(:,1:2), probMat,fluxField,xEdges,yEdges, nstd, plotbox, dbin);
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

run bootStrappingLoop
