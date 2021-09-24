function [se, boots] = se_median(u,N)
% calculates standard error of median using bootstrapping
% USAGE: [SE BOOTS] = se_median(U) 
% U is a matrix containing N x P samples. SE, 1 x P, is the standard error 
% of median for each column of U. BOOTS include bootstrapp samples, sample
% medians and standard error of medians. 
% 
% method: bootstrap
% SE is the standard error of median
% This function generates N bootstraped samples, calculates each sample
% median, mu, calculates standard deviation of all sample medians. 
% This is the standard error of (sample) median.
% Theory: standard error of median represents our uncertainty about median
% of population, which is only obtained through a sample of population.
% In theory, we could quantify our uncertainty by sampling many times
% (e.g. 100 times, each time 15 individuals). Therefore, we use bootstrap
% to simulate many groups and obtain the median of their standard deviation
% as the standard error of median.
% Note boots should be saved for exact reproduction of results, as this is
% a stochastic solution.
% 
% Payam Piray, 
% Donders institute, center for cognitive neuroimaging, June 2015

if nargin<2, N   = 5000; end

se = nan(1,size(u,2));
for i=1:size(u,2)
[bootstat,bootsam] = bootstrp(N,@(x)nanmedian(x,1),u(:,i));
mx = bootstat(:,1);
se(i) = std(mx);
boots = struct('SE',se,'sample',bootsam,'median',bootstat);

end

end

