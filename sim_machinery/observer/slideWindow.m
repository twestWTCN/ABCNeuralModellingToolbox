function [y, index, reject] = slideWindow( x, winsize, overlap)
%SLIDWINDOW is a function that created a sliding window of 1D input signal
%
%   inputs:
%       x: 1D input array that we want to move the slide window over
%       winsize: length of the sliding window
%       overlap: length of overlapping. overlap must be less than winsize
%
%   outputs:
%       y : output, which is a matrix that has winsize rows, and each
%           column is a vector that represents the sliding window output at
%           each itteration.
%       index: the index of the overlapped data in the original input array
%              x.
%              0 index elements refer to padded zeros at the end of the y
%              matrix to maitain the matrix dimentions consistency.
%       reject: number of padded zeros at the last column of y (if needed)

%
%   example:
%      [y,index,reject] = slidWindow(x, winsize); % overlapping = 0 
%      [y,index,reject] = slidWindow(x, winsize, overlapping);
%      

%% Check the type of inputs arguments you get
ng =  nargin;
    if(ng == 1)
        y = x;
        index = 1:length(x);
        reject = 0;
        return
    elseif (ng == 2)
        overlap = 0;
    elseif (ng ~= 3)
        y = [];
        index = [];
        reject = NaN;
        msg = sprintf('Incorrect number of input argument\n check: help slideWindow');
        error(msg)
        help slideWindow
        return
end

%% Check if x is 1D array (also force x to be a row array)
% check if x is a vector
if(~isvector(x))
    y = [];
    index = [];
    reject = NaN;
    error('input x must be a vector')
    return
end

% force x to be a row vector
[row, col] = size(x);
if (row > col)
    x = x.';
end

%% Check that the overlap is positive integer and is not bigger than or equal the winsize
if(~(rem(overlap,1) == 0 ))
    y = [];
    index = [];
    reject = NaN;
    disp(isinteger(overlap))
    error('overlap must be an integer, not a %s.', class(overlap))
    return
elseif (overlap < 0)
    y = [];
    index = [];
    reject = NaN;
    error('overlap must be a positive integer')
    return
elseif (overlap >= winsize)
    y = [];
    index = [];
    reject = NaN;
    error('overlap must be less that window size')
    return
end

%% slideWindow algorithm
y = [];
index = [];
rest = winsize - overlap;
reject = 0;
len_x = length(x);

% for i = 1:winsize
%     start = (i-1)*rest+1 ;
%     stop = (i-1)*rest+winsize;
%     index_temp = start:stop;
% 
%     if stop > len_x
%         stop = len_x;
%         y_temp = x(start:stop);
%         reject = winsize-length(y_temp);
%         y_temp = cat(2,y_temp,zeros(1,reject));
%         index_temp(find(index_temp > len_x)) = 0;
%         y = cat(2,y,y_temp.');
%         index = cat(2,index,index_temp.');
%         %break;
%     else
%         y_temp = x(start:stop);
%         y = cat(2,y,y_temp.');
%         index = cat(2,index,index_temp.');
%     end
%     pause();
% end
stop = 0;
i = 1; % step counter for in window elements
while(stop < len_x)
    start = (i-1)*rest+1 ;
    stop = (i-1)*rest+winsize;
    index_temp = start:stop;

    if stop > len_x
        stop = len_x;
        y_temp = x(start:stop);
        reject = winsize-length(y_temp);
        y_temp = cat(2,y_temp,zeros(1,reject));
        index_temp(find(index_temp > len_x)) = 0;
        y = cat(2,y,y_temp.');
        index = cat(2,index,index_temp.');
        %break;
    else
        y_temp = x(start:stop);
        y = cat(2,y,y_temp.');
        index = cat(2,index,index_temp.');
    end
    i = i + 1;
end

