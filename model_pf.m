function [vol,stc,lr,val,unc] = model_pf(o,specs,islesioned)
if nargin<3, islesioned = 'healthy'; end

np = specs.nparticles;

x0_unc = 1;
if isfield(specs,'x0_unc')
    x0_unc = specs.x0_unc;
end    

switch lower(islesioned(1:3))
    case 'hea'
        lambda_v = specs.lambda_v;
        lambda_s = specs.lambda_s;
        v = specs.v0;
        s = specs.s0;
        [val,vol,stc,lr,unc] = pf_core(o,x0_unc,lambda_v,lambda_s,v,s,np);        
    case 'vol'
        lambda_s = specs.lambda_s;
        v_lesioned = specs.v0_lesioned;        
        s = specs.s0;
        [val,vol,stc,lr,unc] = pf_core_vol_lesioned(o,x0_unc,lambda_s,v_lesioned,s,np);        
    case 'sto'
        lambda_v = specs.lambda_v;
        s_lesioned = specs.s0_lesioned;
        v = specs.v0;
        [val,vol,stc,lr,unc] = pf_core_unp_lesioned(o,x0_unc,lambda_v,v,s_lesioned,np);        
    otherwise
        error('bad 3rd input: %d',islesioned);
end


end