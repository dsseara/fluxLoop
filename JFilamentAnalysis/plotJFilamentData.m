% function plotJFilamentData(px2um, maxFrame, savestuff)
%
% This function takes the output from loadJFilamentData and plots each
% filament for all frames in a single figure, color coded such that as
% the frame number goes up, the color of the filament gradually goes up
% the parula colormap
%
% INPUTS	 
%			 filament : cell array of data of filament positions over time.
%					    Output of loadJFilamentData
%			    px2um : Conversion factor of pixels to microns. Units: [um/pixel]
%			 maxFrame : The most frames a single filament can be tracked
%			savestuff : Boolean of if to save the plots or not
%			   folder : If savestuff, then save them in this folder
%			frameRate : Plot every frameRate-th frame (e.g. 1 for all frames,
%						10 for every 10th frame)
%
% OUTPUTS	
%			plots of filaments
%
% Created by Daniel Seara, 05/10/2017
function plotJFilamentData(filament, px2um, maxFrame, savestuff, folder,frameRate)
	colors = colormap(parula(maxFrame));
	close all;
	for ii = 3%1:numel(filament)
		filaData = filament{ii};
		frames = unique(filaData(:,1));
		figure;
		hold on;
		for jj = 9%1:frameRate:numel(frames)
			f = frames(jj);
			plot(filaData(filaData(:,1)==f,3).*px2um,...
				 filaData(filaData(:,1)==f,4).*px2um,...
			     'Color',colors(f,:))
		end
		set(gca,'YDir','reverse','xtick',[],'ytick',[])
		%xlabel('x (\mu m)')
		%ylabel('y (\mu m)')
		%axis equal
		%title(['Filament', num2str(ii)])
		
		llmFig('plw',1.5,'font','Helvetica')
		c = colorbar;
		ylabel(c,'t/t_c')
		set(c,'ytick',[0,0.5,1])

		if savestuff
			saveas(gcf,[folder filesep 'filament_', num2str(ii)],'fig')
			saveas(gcf,[folder filesep 'filament_', num2str(ii)],'tif')
			saveas(gcf,[folder filesep 'filament_', num2str(ii)],'epsc')
		end
	end