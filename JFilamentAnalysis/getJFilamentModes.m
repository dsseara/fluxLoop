% function aa = getJFilamentModes(filament,nmax,bcs,savestuff)
% This function rotates the filaments so that they lie horizontally,
% with the first frame end points of the filament determining what is 
% 'horizontal'. Get tangent angle data along the arc length of the filament.
% Plots both of these things together. Then also decomposes the tangent angle
% into elastohydrodynamic modes of thin rods using elastohydroModes
%
% INPUTS	 
%			 filament : Array of filament tangent angles. Each row is a position
%						along the arc length, each column is a time point
% 			     nmax : Maximum number of modes to calculate for the filaments
% 			      bcs : String of boundary conditions of the rods (use 'free')
%			 maxFrame : The most frames a single filament can be tracked
%			    px2um : Conversion factor of pixels to microns. Units: [um/pixel]
% 		    savestuff : Boolean of if to save the plots or not
% 			   folder : String. If savestuff, then save them in this folder
%
% OUTPUTS		   aa : (N)x(nmax) array of the mode coefficients for every
% 					    each time point contained in a single row
%
% This version is made for axoneme data from the following Dryad database:
% http://datadryad.org/resource/doi:10.5061/dryad.0529j
% 
% For the following paper:
%
% Dynamic curvature regulation accounts for the symmetric and asymmetric beats of
% Chlamydomonas flagella
% Sartori P, Geyer VF, Scholich A, Jülicher F, Howard J
% Date Published: May 12, 2016
% DOI: http://dx.doi.org/10.5061/dryad.0529j
%
% Created by Daniel Seara, 07/27/2017
function aa = getJFilamentModes(filament,nmax,bcs,maxFrame,px2um,savestuff,folder, frameRate)
	
	colors = colormap(parula(maxFrame));
	close all;

    aa = [];
	nframes = size(filament,2);
	figure;
	for jj = 1:nframes
        tangents = filament(:,jj);
        ds = ones(size(tangents))*px2um;
        L = sum(ds);
        % keyboard
        elastoCoeffs = elastohydroModes(tangents', cumsum(ds)', L, nmax, bcs);
        aa = [aa; elastoCoeffs'];
        if mod(jj,frameRate)==0
	        subplot(1,1,1), hold on;
	        plot(cumsum(ds), tangents,'Color', colors(jj,:))
		end
	end % End loop over frames
    
    subplot(1,1,1)
    axis equal
    xlabel('arc length (\mum)')
    ylabel('rad')

    if savestuff
    	if ~exist([folder filesep 'filaments'],'dir')
    		mkdir([folder filesep 'filaments']);
    	end

		saveas(gcf,[folder '/filaments/filament',num2str(ii)],'fig')
	    saveas(gcf,[folder '/filaments/filament',num2str(ii)],'tif')
	    saveas(gcf,[folder '/filaments/filament',num2str(ii)],'epsc')
	end
end