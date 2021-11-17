function pf = pf_initialize(StateTransitionFcn, MeasurementLikelihoodFcn,NumParticles,initial_stateBounds)
initial_method = 'uniform';
if all( (initial_stateBounds(:,2)-initial_stateBounds(:,1))==0)
    initial_method = 'fixed';
end    

pf.StateTransitionFcn = StateTransitionFcn;
pf.MeasurementLikelihoodFcn = MeasurementLikelihoodFcn;
pf.NumParticles = NumParticles;
pf.initial_method = initial_method;
pf.initial_stateBounds = initial_stateBounds;
pf.resampling_method = 'systematic';

if strcmpi(initial_method,'uniform')
    siz = size(initial_stateBounds);
    d = siz(1);
    if siz(2)~=2        
        error('Initialization error');
    end
    l = initial_stateBounds(:,2)-initial_stateBounds(:,1);
    if any(l<0)
        error('Initialization error');
    end
    x0 = rand(d,NumParticles);
    particles = bsxfun(@times, l, x0) + initial_stateBounds(:,1);  

elseif strcmpi(initial_method,'fixed')
    siz = size(initial_stateBounds);
    d = siz(1);
    if siz(2)~=2        
        error('Initialization error');
    end
    l = initial_stateBounds(:,2)-initial_stateBounds(:,1);
    if any(l<0)
        error('Initialization error');
    end
    x0 = ones(d,NumParticles);
    particles = bsxfun(@times, initial_stateBounds(:,1), x0);      
end

weights = 1/NumParticles*ones(1, NumParticles);

pf.particles = particles;
pf.weights = weights;

end
