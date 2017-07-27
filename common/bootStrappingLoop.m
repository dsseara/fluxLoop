% function [curlNormed, curlHistActive, curlHistCtrl, dbinArray] =...
%     bootStrappingLoop(activeScore, ctrlScore, dt, cutoff, plotbox,...
%     m, nbinArray, stdArray, plotstuff)
%
% Function that performs a bootstrapping procedure on an active and control
% set of time series


function [curlNormed, curlHist, dbinArray] =...
    bootStrappingLoop(tSeries, dt, cutoff, plotbox,...
    m, nbinArray, stdArray, savefolder)


    curlNormed = zeros(numel(nbinArray), numel(stdArray)+1,m);
    curlHist   = zeros(numel(nbinArray), numel(stdArray)+1,m);

    rng(1);
    for ii =1:numel(nbinArray)
        dbin = (max(tSeries(:,2)) - min(tSeries(:,2)))/nbinArray(ii);
        dbinArray(ii) = dbin;
        for jj =  1:numel(stdArray)+1
            counter = 0;
            figure
            if jj == numel(stdArray)+1
                nstd = [];
            else
                nstd = stdArray(jj);
            end
            
            tempCurl = nan(m,1);

            tic
            while counter<m
                counter = counter+1;
                newTSeries = bootStrap(tSeries);

                [probMat,fluxField,xEdges,yEdges] = ...
                        probabilityFlux(newActiveTSeries,dt,dbin,cutoff);

                tempCurl(counter) = fluxLoopCurl(newTSeries,...
                        probMat,fluxField,xEdges,yEdges,...
                        nstd, plotbox, dbin);
            end
            toc

            curlHist(ii,jj, :) = tempCurl;
            curlNormed(ii,jj, :) = tempCurl./std(tempCurl.active);
        end

    end
end