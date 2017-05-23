% function output = transitionTimeSeries(tseries,dt,xEdges, yEdges)
% This function creates an Nx3 matrix that gives the transition data
% so that we can bootstrap another trajectory out of these data
% while still maintaining the temporal relationship between the transitions
%
% INPUTS
%            tseries : The 2D time series being investigated
%                 dt : time between data points
%           x/yEdges : the edges used to do the binning 
%                        (output of probabilityFlux)
% OUTPUTS
%             output : (N-1)x3 matrix of the form:
%                         a_1     a_2 t_1
%                         a_2     a_3 t_2
%                         ...
%                         a_(N-1) a_N t_(N-1)
%                     where the 1st column is the state the system 
%                     is in before a transition, the 2nd column is
%                     the state the system is in after a transition,
%                     and the 3rd column is the time spent in the first
%                     state before making the transition
%
% Created by Daniel Seara, 05/12/2017
function newTSeries = bootStrap(tseries,dt,xEdges, yEdges)
	output = [];
	if size(tseries,1)<size(tseries,2)
        % make sure data entered in right orientation
        tseries = tseries';
    end

    n = size(tseries,1); % total number of time points
	nbinx = numel(xEdges)-1;
    nbiny = numel(yEdges)-1;




    for jj = 2:n
        % Get prior state
        px = tseries(jj-1,1);
        py = tseries(jj-1,2);
        pState = histcounts2(px,py, xEdges, yEdges);
        [pRow, pCol] = find(pState);
        
        % Get current state
        cx = tseries(jj,1);
        cy = tseries(jj,2);
        cState = histcounts2(cx,cy, xEdges, yEdges);
        [cRow, cCol] = find(cState);
        
        % Get all points traversed in phase space using Bresenham algorithm
        [xPath, yPath] = bresenham(pRow,pCol, cRow, cCol);
        % keep track of time spent at each state, divided equally between all 
        % the states on the path. Time between frames is 5 seconds
        
        % this part was tricky, switch order of path cuz x direction changes
        % in columns not rows, similar for y
        inds = sub2ind([nbiny, nbinx],yPath, xPath); 
        tPerState = dt/numel(inds);

        output = [output;...
                inds(1:end-1), inds(2:end),tPerState.*ones(numel(inds)-1,1)];
    end
end % end function