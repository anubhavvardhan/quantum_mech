function [ind,t0,s0,t0close,s0close] = crossing_vectors(S,t,level,imeth)
% to find the crossings of a given level of a signal
% no. of inputs
error(nargchk(1,4,nargin));

% checking the time vector's consistency
if nargin < 2 || isempty(t)
	% taking index vector as default time
    t = 1:length(S);
elseif length(t) ~= length(S)
	% if S and t are not of the same length, throw an error
    error('t and S must be of identical length!');    
end

% checking the level input
if nargin < 3
	% default = 0, if level is not given
    level = 0;
end

% check interpolation method input
if nargin < 4
    imeth = 'linear';
end

% row vectors
t = t(:)';
S = S(:)';

% If we want the crossing of a threshold value "level", 
% we subtract it from the values and search for zeros.
S   = S - level;

% first look for exact zeros
ind0 = find( S == 0 ); 

% then look for zero crossings between data points
S1 = S(1:end-1) .* S(2:end);
ind1 = find( S1 < 0 );

% bring exact zeros and "in-between" zeros together 
ind = sort([ind0 ind1]);

% and pick the associated time values
t0 = t(ind); 
s0 = S(ind);

if strcmp(imeth,'linear')
    % linear interpolation of crossing
    for ii=1:length(t0)
        if abs(S(ind(ii))) > eps(S(ind(ii)))
            % interpolate only when data point is not already zero
            NUM = (t(ind(ii)+1) - t(ind(ii)));
            DEN = (S(ind(ii)+1) - S(ind(ii)));
            DELTA =  NUM / DEN;
            t0(ii) = t0(ii) - S(ind(ii)) * DELTA;
            s0(ii) = 0;
        end
    end
end

