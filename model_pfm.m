function [vol,stc,lr,val] = model_pfm(o,specs,islesioned)
if nargin<3, islesioned = 0; end

lambda_s = specs.lambda_s;
s = specs.s0;
np = specs.nparticles;

x0_unc = 1;
if isfield(specs,'x0_unc')
    x0_unc = specs.x0_unc;
end

if ~islesioned
    lambda_v = specs.lambda_v;
    v = specs.v0;
    [val,vol,stc,lr] = pfm_core(o,x0_unc,lambda_v,lambda_s,v,s,np);
else
    v = specs.v0_lesioned;
    [val,vol,stc,lr] = pfm_core_vol_lesioned(o,x0_unc,lambda_s,v,s,np);
end

dim = size(lr,2);
if size(vol,2)<dim
    vol = repmat(vol,1,dim);
end
if size(stc,2)<dim
    stc = repmat(stc,1,dim);
end

end
