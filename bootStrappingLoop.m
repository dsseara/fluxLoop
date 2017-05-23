% Vary the bin size and see if that affects the total curl found.
dt = 5;
cutoff = [];
plotbox = 0;
m=500;

nbinArray = [1:10].*5;
dbinArray = zeros(size(nbinArray));
stdArray = [1,2,3];

activeTSeries = activScore;
ctrlTSeries = ctrlScore;

curlNormed = zeros(numel(dbinArray),numel(stdArray)+1,2);
curlHistActive = zeros(numel(dbinArray), numel(stdArray)+1,m);
curlHistCtrl   = zeros(numel(dbinArray), numel(stdArray)+1,m);

rng(1);
for ii =1:numel(nbinArray)
	dbin = (max(activeTSeries(:,2)) - min(activeTSeries(:,2)))/nbinArray(ii);
	dbinArray(ii) = dbin;
	for jj =  1:numel(stdArray)+1
        counter = 0;
        figure
        if jj == numel(stdArray)+1
			nstd = [];
        else
            nstd = stdArray(jj);
        end
        
        tempCurl.active = nan(m,1);
        tempCurl.ctrl = nan(m,1);
        
        tic
        while counter<m
            counter = counter+1;
            newActiveTSeries = bootStrap(activeTSeries);
            newCtrlTSeries   = bootStrap(ctrlTSeries);
		
            [activeprobMat,activefluxField,activexEdges,activeyEdges] = ...
                    probabilityFlux(newActiveTSeries,dt,dbin,cutoff);
                
            tempCurl.active(counter) = JFilamentFluxLoop(newActiveTSeries,...
                    activeprobMat,activefluxField,activexEdges,activeyEdges,...
                    nstd, plotbox, dbin);

            [ctrlprobMat,ctrlfluxField,ctrlxEdges,ctrlyEdges] = ...
                    probabilityFlux(newCtrlTSeries,dt,dbin,cutoff);
                
            tempCurl.ctrl(counter) = JFilamentFluxLoop(newCtrlTSeries,...
                    ctrlprobMat,ctrlfluxField,ctrlxEdges,ctrlyEdges,...
                    nstd, plotbox, dbin);
        end
        toc
        
        [t{ii,jj},n{ii,jj},x{ii,jj}] = nhist(tempCurl);
        title(['dbin=',num2str(dbin),' nstd=',num2str(jj)])
        
        curlHistActive(ii,jj,:) = tempCurl.active;
        curlHistCtrl(ii,jj,:)   = tempCurl.ctrl;
        curlNormed(ii,jj,1) = mean(tempCurl.active)./std(tempCurl.active);
        curlNormed(ii,jj,2) = mean(tempCurl.ctrl)./std(tempCurl.ctrl);
    end
    
end