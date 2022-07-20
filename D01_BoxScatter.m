% [Function BoxScatter]
% Update: 2021-09-19
% Randomized the placements of datapoints (randperm)
%
% First created: 2021-06-09, Jongmin Kim (jongminkmg@gmail.com)
% Get input of variable, location, markersize, color number,
% Return a scattered plot of data points with a median line.

% Please note this is intended to my own use, and has fixed input
% requirements. Feel free to modify in any way.



function BoxScatter (val, loc, size, c1, c2, c3)

n_val = numel(val);
val = val(randperm(length(val)));
jitter = 0.6/n_val;
med_val = median(val);

x = loc-0.3:jitter:loc+0.3-jitter;
hold on
plot (x, val, 'o', 'MarkerSize', size, 'Color', [c1, c2, c3]);
plot ([loc-0.4 loc+0.4], [med_val med_val], 'Linewidth', 1, 'Color', [0 0 0]);
hold off
