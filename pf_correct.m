function [pf, idx] = pf_correct(pf, measurement, varargin)

particles = pf.particles;
MeasurementLikelihoodFcn = pf.MeasurementLikelihoodFcn;

likelihood = MeasurementLikelihoodFcn(particles, measurement, varargin{:});
likelihood = likelihood+eps; % add eps to prevent a possible division by 0 problem
weights = pf.weights;
weights = weights.*likelihood;
weights = weights/sum(weights);


Neff = 1/sum(weights.^2);
resample_percentaje = 0.50;
Nt = resample_percentaje*pf.NumParticles;
idx = 1:pf.NumParticles;
if Neff < Nt   
   [idx, weights] = resample(weights, pf.resampling_method);
   particles = particles(:,idx);
end
pf.particles = particles;
pf.weights = weights;

end

function [idx, wk] = resample(wk, resampling_strategy)
wk = wk';
Ns = length(wk);  % Ns = number of particles

switch resampling_strategy
   case 'systematic'
      % this is performing latin hypercube sampling on wk
      edges = min([0 cumsum(wk)'],1); % protect against accumulated round-off      
      edges(end) = 1;                 % get the upper edge exact
      u1 = rand/Ns;
      % this works like the inverse of the empirical distribution and returns
      % the interval where the sample is to be found
      [~, idx] = histc(u1:1/Ns:1, edges);
       wk = repmat(1/Ns, 1, Ns);  
    case 'multinomial'        
       idx = resample_multinomial(wk);
       wk = wk(idx,:)';
       wk = wk/sum(wk);
       
       
   otherwise
      error('Resampling strategy not implemented')
end


end

function idx = resample_multinomial(w)
N = length(w);
% Multinomial resampling
u = rand(N,1);
qc = cumsum(w);
qc = qc(:);
qc=qc/qc(N);
[~,ind1]=sort([u;qc]);
ind2=find(ind1<=N);
idx=ind2'-(0:N-1);
end
