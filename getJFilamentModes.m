% function aa = getJFilamentModes(filament,nmax,bcs,savestuff)
% This function rotates the filaments so that they lie horizontally,
% with the first frame end points of the filament determining what is 
% 'horizontal'. Get tangent angle data along the arc length of the filament.
% Plots both of these things together. Then also decomposes the tangent angle
% into elastohydrodynamic modes of thin rods using elastohydroModes
%
% INPUTS	 
%			 filament : Mx5 array filament positions over time, output of
% 					    loadJFilamentData. Array structure:
%						frame | segment # | x | y | all zeros (?)
% 			     nmax : Maximum number of modes to calculate for the filaments
% 			      bcs : String of boundary conditions of the rods (use 'free')
%			 maxFrame : The most frames a single filament can be tracked
%			    px2um : Conversion factor of pixels to microns. Units: [um/pixel]
% 		    savestuff : Boolean of if to save the plots or not
% 			   folder : String. If savestuff, then save them in this folder
%
% OUTPUTS		   aa : (N)x(nmax) array of the mode coefficients for every
% 					    filament at every time point, each time point contained
% 					    in a single row
%
% Created by Daniel Seara, 05/10/2017
function aa = getJFilamentModes(filament,nmax,bcs,maxFrame,px2um,savestuff,folder)
	aa = [];
	colors = colormap(parula(maxFrame));
	close all;
	for ii = 1:numel(filament) 
		filaData = filament{ii};
		frames = unique(filaData(:,1));
		figure;
		for jj = 1:numel(frames)
			f = frames(jj);
			xy = [filaData(filaData(:,1)==f,3)'.*px2um;...
				  filaData(filaData(:,1)==f,4)'.*px2um];
			npoints = size(xy,2);
			if jj==1
				theta = atan2(xy(2,end)-xy(2,1), xy(1,end)-xy(1,1));
				% rotate clockwise by angle theta
		    	R = [cos(theta), sin(theta);-sin(theta), cos(theta)];
		    end % End if statement to get rotation angle

		    % Rotate and subtract starting position
	        xy = R*(xy - [xy(1,1).*ones(1,size(xy,2));xy(2,1).*ones(1,size(xy,2))]);
	        displacements = diff(xy, 1,2); % take difference along columns
	        ds = sqrt(displacements(1,:).^2 + displacements(2,:).^2);
	        L = sum(ds);
	        positions = -L/2:(L/(npoints-2)):L/2;
	        tangents = atan2(displacements(2,:), displacements(1,:));
	        elastoCoeffs = elastohydroModes(tangents, ds, L, nmax, bcs);
	        aa = [aa; elastoCoeffs'];

	        subplot(2,1,1), hold on;
	        plot(xy(1,:)-L/2,xy(2,:),'Color', colors(jj,:))
	        subplot(2,1,2), hold on;
	        plot(positions, tangents,'Color', colors(jj,:))
		end % End loop over frames
	    subplot(2,1,1)
	    xlabel('x (\mum)')
	    ylabel('y (\mum)')
	    subplot(2,1,2)
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
	end % End loop over filaments
end