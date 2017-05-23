% function plotJFilamentModeCorrelations(aa,savestuff, folder)
% This function plots all the mode coefficients against one another
% systematically to get an idea of which modes are correlated or not

% INPUTS	    
%				aa : Nxn array of elastohydrodynamic modes of filaments,
% 		       	     output of getJFilamentModes
% 		 savestuff : boolean of if to save the plots or not
% 			folder : If savestuff, then save them in this folder
%
% OUTPUTS
%			plots of mode correlations

function plotJFilamentModeCorrelations(aa,savestuff, folder)
	n = size(aa,2);
	for ii = 1:n
		for jj = 1:n
			if jj<=ii % Don't double count mode pairs
				continue
			end
			figure
			scatter(aa(:,ii),aa(:,jj))
			xlabel(['Mode ', num2str(ii)])
			ylabel(['Mode ', num2str(jj)])

			if savestuff
				if ~exist([folder filesep 'modes'],'dir')
	    			mkdir([folder filesep 'modes']);
	    		end
				saveas(gcf,[folder '/modes/', num2str(ii),num2str(jj)],'fig')
				saveas(gcf,[folder '/modes/', num2str(ii),num2str(jj)],'epsc')
				saveas(gcf,[folder '/modes/', num2str(ii),num2str(jj)],'tif')
			end % end saving
		end % end loop over second mode
	end % end loop over first mode
end % end function