% function netcurl = JFilamentFluxLoop(tseries, probMat,fluxField,xEdges,
%                                           yEdges, nstd, plotbox)
%
% This function finds the net curl of the flux vector field given by 
% probabilityFlux.
%
% INPUTS	
%		 tseries : time series data used in probabilityFlux, used if 
%				   getting the net curl inside a certain region around
%                  the mean
%		 probMat : output probability map of probabilityFLux
%	   fluxField : output vector field of probabilityFlux
%		  xEdges : output x-edges used to bin the data in probabilityFlux
%	      yEdges : output y-edges used to bin the data in probabilityFlux
%			nstd : number of times the standard deviation in each direction
%				   to consider when calculating the total curl (empty for
%                  all)
%		 plotbox : boolean of whether to replot the the phase space with
%					  the area over which the curl is summed 
% 
% Created by Daniel Seara, 05/11/2017


function netcurl = fluxLoopCurl(tseries, probMat,fluxField,xEdges,...
    yEdges, nstd, plotbox, dbin)
	% if nstd = [], find total curl over entire region
	if isempty(nstd)
		netcurl = dbin.^2.*sum(sum(curl(xEdges(1:end-1), yEdges(1:end-1),...
					  fluxField(:,:,1), fluxField(:,:,2))));
	else % find curl inside a small region defined by std
		stdX = std(tseries(:,1));
		meanX = mean(tseries(:,1));
		stdY = std(tseries(:,2));
		meanY = mean(tseries(:,2));

		insideX = xEdges(1:end-1)>meanX-nstd*(stdX) &...
                  xEdges(1:end-1)<meanX+nstd*(stdX);
		insideY = yEdges(1:end-1)>meanY-nstd*(stdY) &...
                  yEdges(1:end-1)<meanY+nstd*(stdY);

		if sum(insideX)<2 || sum(insideY)<2
			netcurl = 0;
			return
        end
        %keyboard
		netcurl = dbin.^2.*sum(sum(curl(xEdges(insideX), yEdges(insideY),...
			 		  fluxField(find(insideY),find(insideX),1),...
			 		  fluxField(find(insideY), find(insideX),2))));

		if plotbox
			figure;
			pcolor(xEdges(1:end-1), yEdges(1:end-1),probMat), hold on
			quiver(xEdges(1:end-1)+dbin/2, yEdges(1:end-1)+dbin/2,...
					fluxField(:,:,1),fluxField(:,:,2),'w')
			xlabel('pca 1')
			ylabel('pca 2')
			insideXEdges = xEdges(insideX);
			insideYEdges = yEdges(insideY);

			plot([insideXEdges(1),insideXEdges(end)+dbin],...
				 [insideYEdges(1),insideYEdges(1)],'r:','LineWidth',1)

			plot([insideXEdges(1),insideXEdges(end)+dbin],...
				 [insideYEdges(end)+dbin,insideYEdges(end)+dbin],'r:','LineWidth',1)

			plot([insideXEdges(end)+dbin,insideXEdges(end)+dbin],...
				 [insideYEdges(1),insideYEdges(end)+dbin],'r:','LineWidth',1)

			plot([insideXEdges(1),insideXEdges(1)],...
				 [insideYEdges(1),insideYEdges(end)+dbin],'r:','LineWidth',1)

			title([num2str(nstd), ' std region, curl = ', num2str(netcurl)])
		end
	end
end