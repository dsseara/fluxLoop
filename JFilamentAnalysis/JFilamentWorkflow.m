%% JFilament workflow modified for Howard axoneme data
% Script to analyze JFilament data to obtain probability flux loop information

% First load the data you want. Will have to manually
% enter save folder name and table name

folder = '/media/daniel/storage1/local_LLM_Danny/fActinOrder/howardAxonemeData';
axonemeData = 'wt01';
tangents = wt01;

% Subtract off linear rate of rotation
dt = 0.001;
t = cumsum(ones(size(tangents,2),1).*dt);
p = polyfit(t, tangents(1,:)',1);

tangents = tangents - repmat((p(1)*t' + p(2)), [size(tangents,1),1]);
%%
% Let's plot up the filaments so we can take a look at them.

px2um= .4417; % 0.4417 microns per pixel
maxFrame = size(tangents,2); % 75 is absolute max frame possible
savestuff = 0; % boolean to save stuff or not
n = 10; % plot every nth frame
ds = cumsum(ones(size(tangents,1)).*px2um)-px2um;

figure, hold on
for ii=1:n:maxFrame
    plot(ds, tangents(:,ii));
    xlabel('arc length (\mum)')
    ylabel('\theta (rad)')
end

if savestuff
    saveas(gcf, [folder filesep axonemeData],'fig')
    saveas(gcf, [folder filesep axonemeData],'epsc')
    saveas(gcf, [folder filesep axonemeData],'tif')
end

%%
% Now we rotate the filaments so they lie on a horizontal axis (defined as the
% line between the first and last point in the first frame), and get the
% tangent information. Use this information to first get the elastohydrodynamic
% modes of slender rods, but also store in order to use PCA on it afterwards

bcs = 'free';   % Use free boundary conditions
nmax = 6;       % get the first 6 modes
savestuff = 0;  % Save the plots?
n = 10; % plot every nth frame
aa = getJFilamentModes(tangents, nmax,bcs,maxFrame,px2um,savestuff,folder,n);

%%
% See how the modes correlate with each other, expect all the even
% modes and odd modes to correlate linearly with other even or odd modes,
% respectively, showing variation is mainly an even or oddness issue.

savestuff = 1;
plotJFilamentModeCorrelations(aa,savestuff, folder)

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
        h = pcolor(xEdges(1:end-1), yEdges(1:end-1),probMat);
        hold on
        set(h, 'EdgeColor','none')
        quiver(xEdges(1:end-1)+dbin/2, yEdges(1:end-1)+dbin/2,...
                fluxField(:,:,1),fluxField(:,:,2),'w')
        colorbar
        xlabel(['Mode' num2str(axis1)])
        ylabel(['Mode' num2str(axis2)])
% %         axis equal
        saveas(gcf,[folder filesep 'phaseSpace' filesep...
             'modes_' num2str(axis1) num2str(axis2)],'tif');
    end
end

%%
% % Vary the bin size and see if that affects the total curl found.
% dt = 5;
% cutoff = [];
% 
% axis1 = 1;
% axis2 = 2;
% 
% nbinArray = [1:10].*5;
% dbinArray = zeros(size(nbinArray));
% stdArray = [1,2,3];
% 
% curlArray = zeros(numel(dbinArray),numel(stdArray)+1);
% 
% for ii =1:numel(nbinArray)
% 	dbin = (max(aa(:,axis2)) - min(aa(:,axis2)))/nbinArray(ii);
% 	dbinArray(ii) = dbin;
% 	for jj =  1:numel(stdArray)+1
% 		if jj == numel(stdArray)+1
% 			figure
% 			[probMat,fluxField,xEdges,yEdges] =...
%                 probabilityFlux([aa(:,axis1),aa(:,axis2)],dt,dbin,cutoff);
% 			h = pcolor(xEdges(1:end-1), yEdges(1:end-1),probMat);
%             hold on
% 			set(h, 'EdgeColor', 'none')
%             quiver(xEdges(1:end-1)+dbin/2, yEdges(1:end-1)+dbin/2,...
% 					fluxField(:,:,1),fluxField(:,:,2),'w')
% 			colorbar
% 			xlabel(['Mode ' num2str(axis1)])
%             ylabel(['Mode ' num2str(axis2)])
% 			title(['All curl, dbin = ',num2str(dbin)])
% 			axis equal
% 			plotbox = 0;
% 			nstd = [];
% 			curlArray(ii,jj) =...
%                 fluxLoopCurl([aa(:,axis1),aa(:,axis2)], probMat,fluxField,...
%                 xEdges,yEdges, nstd, plotbox, dbin);
% %              saveas(gcf,['varyBinSize/dbin_',strrep(num2str(dbin),'.',',')],'fig')
% %              saveas(gcf,['varyBinSize/dbin_',strrep(num2str(dbin),'.',',')],'tif')
% %              saveas(gcf,['varyBinSize/dbin_',strrep(num2str(dbin),'.',',')],'epsc')
%         else
% 			[probMat,fluxField,xEdges,yEdges] =...
%                 probabilityFlux([aa(:,axis1),aa(:,axis2)],dt,dbin,cutoff);
% 			plotbox = 1;
% 			nstd = stdArray(jj);
% 			curlArray(ii,jj) = ...
%                 fluxLoopCurl([aa(:,axis1),aa(:,axis2)], probMat,fluxField,...
%                 xEdges,yEdges, nstd, plotbox, dbin);
% 			if plotbox
% 				title(['Inner curl, nstd = ',num2str(nstd), ' dbin = ', num2str(dbin)])
% 				axis equal
% 			end
% 		end
% 	end
% end

%%
% Start the bootstrapping procedure
% bootStrappingLoop.m varies the bin size and area over which the
% curl is integrated in order to find out what parameter set gives
% us a statistically significant circulation in the curl

dt = 0.001;
cutoff = [];
plotstuff = 1;
m=500; % number of bootstrapped trajectories generated
nbinArray = [1:2:10].*5; % Array of number of bins to change bin size
dbinArray = zeros(size(nbinArray)); % array to keep 
stdArray = []; % number of stds to integrate over also does the entire
               % domain in addition to these
plotbox = 0;


parfor ii = 1:6
    for jj = ii+1:6
        close all
        axis1 = ii;
        axis2 = jj;
        
%         % Don't do the mode combinations that have already been done
%         if axis1==1 && axis2==2
%             continue
%         elseif axis1==2 && axis2==6
%             continue
%         elseif axis1==5 && axis2==6
%             continue
%         end
        
        modeDir = [folder '/bootStrap/modes_',...
            num2str(axis1), num2str(axis2)];
        mkdir(modeDir);
        mkdir([modeDir '/histograms']);
        mkdir([modeDir '/histograms/tif']);
        mkdir([modeDir '/histograms/eps']);
        mkdir([modeDir '/histograms/fig']);
        
        savefolder = [modeDir '/histograms'];
        
        tSeries = [aa(:,axis1), aa(:,axis2)];
        
        [curlNormed, curlHist, dbinArray] =...
            bootStrappingLoop(tSeries, dt, cutoff, plotbox, m,...
            nbinArray, stdArray, savefolder, plotstuff);
        
        saveFile = matfile([modeDir filesep axonemeData '_bootstrapData.mat'],'writable',true)
        
        saveFile.tSeries = tSeries;
        saveFile.curlHist = curlHist;
        saveFile.curlNormed = curlNormed;
        saveFile.dbinArray = dbinArray;

    end
end