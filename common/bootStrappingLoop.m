% function [curlNormed, curlHistActive, curlHistCtrl, dbinArray] =...
%     bootStrappingLoop(activeScore, ctrlScore, dt, cutoff, plotbox,...
%     m, nbinArray, stdArray, plotstuff)
%
% Function that performs a bootstrapping procedure on an active and control
% set of time series


function [curlNormed, curlHistActive, curlHistCtrl, dbinArray] =...
    bootStrappingLoop(activeScore, ctrlScore, dt, cutoff, plotbox,...
    m, nbinArray, stdArray, plotstuff, savefolder)


    activeTSeries = activeScore;
    ctrlTSeries = ctrlScore;

    curlNormed = zeros(numel(nbinArray),numel(stdArray)+1,2);
    curlHistActive = zeros(numel(nbinArray), numel(stdArray)+1,m);
    curlHistCtrl   = zeros(numel(nbinArray), numel(stdArray)+1,m);

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

                tempCurl.active(counter) = fluxLoopCurl(newActiveTSeries,...
                        activeprobMat,activefluxField,activexEdges,activeyEdges,...
                        nstd, plotbox, dbin);

                [ctrlprobMat,ctrlfluxField,ctrlxEdges,ctrlyEdges] = ...
                        probabilityFlux(newCtrlTSeries,dt,dbin,cutoff);

                tempCurl.ctrl(counter) = fluxLoopCurl(newCtrlTSeries,...
                        ctrlprobMat,ctrlfluxField,ctrlxEdges,ctrlyEdges,...
                        nstd, plotbox, dbin);
            end
            toc

            if plotstuff
                [t{ii,jj},n{ii,jj},x{ii,jj}] = nhist(tempCurl);
                xlabel('\Omega')
                ylabel('pdf')
                title(['dbin=',num2str(dbin),', std=', num2str(nstd)]);
                
                dbinStr = strrep(num2str(dbin),'.',','); % Can't have period in filename
                
                saveas(gcf, [savefolder '/fig/dbin_' dbinStr '_nstd_' num2str(nstd)], 'fig');
                saveas(gcf, [savefolder '/eps/dbin_' dbinStr '_nstd_' num2str(nstd)], 'epsc');
                saveas(gcf, [savefolder '/tif/dbin_' dbinStr '_nstd_' num2str(nstd)], 'tif');
            end

            curlHistActive(ii,jj,:) = tempCurl.active;
            curlHistCtrl(ii,jj,:)   = tempCurl.ctrl;
            curlNormed(ii,jj,1) = mean(tempCurl.active)./std(tempCurl.active);
            curlNormed(ii,jj,2) = mean(tempCurl.ctrl)./std(tempCurl.ctrl);
        end

    end
end