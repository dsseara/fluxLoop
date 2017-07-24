% function newTSeries = bootStrap(tseries)
% This function creates an Nx3 matrix that gives the transition data
% so that we can bootstrap another trajectory out of these data
% while still maintaining the temporal relationship between the transitions
%
% INPUTS
%            tseries : Nx2 time series being investigated
%
% OUTPUTS
%             newTSeries : new Nx2 time series sampled WITH REPLACEMENT from
%                          tseries
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

    newTSeries = zeros(size(tseries));
    
    n = size(tseries,1); % total number of time points
    counter = 0;

    while counter<n/2-1
        counter = counter + 1;

        % Random row number between 1 and end-1 of tseries
        randRow = ceil(rand*(n-1));

        % Create newTSeries using sequential rows
        newTSeries(counter, :)     = tseries(randRow, :);
        newTSeries(counter + 1, :) = tseries(randRow + 1, :);
    end


end % end function