% function output = transitionTimeSeries(tseries,dt,xEdges, yEdges)
% This function creates an Nx3 matrix that gives the transition data
% so that we can bootstrap another trajectory out of these data
% while still maintaining the temporal relationship between the transitions
%
% INPUTS
%            tseries : The 2D time series being investigated
%
% OUTPUTS
%             newTSeries : shuffled version of tseries where pairs of
%                          rows are kept adjacent 
%
% Created by Daniel Seara, 05/12/2017
function newTSeries = bootStrap(tseries)
	
    if size(tseries,1)<size(tseries,2)
        % make sure data entered in right orientation
        tseries = tseries';
    end

    % First make sure that the number of elements in tseries is even to do
    % the bootstrapping in the first place
    if mod(size(tseries,1),2)
        tseries = tseries(1:end-1,:);
    end
    
    n = size(tseries,1); % total number of time points

    % gives random rows from 2:2:size(tseries,1)
    ordering = randperm(n/2).*2; 
    % now interleave ordering with a vector that is just the odd neighbors of ordering
    ordering = reshape([ordering-1;ordering],2*size(ordering,2),[]);
    
    newTSeries = tseries(ordering,:);

end % end function