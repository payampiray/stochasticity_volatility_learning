function pf = pf_predict(pf,varargin)

particles = pf.particles;
StateTransitionFcn = pf.StateTransitionFcn;

particles = StateTransitionFcn(particles, varargin{:});
pf.particles = particles;
end