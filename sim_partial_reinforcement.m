function h = sim_partial_reinforcement(nr,nc,subplots,do_supp)

fname = fullfile('simulations',sprintf('%s.mat',mfilename));
fsum = fullfile('sum',sprintf('%s.mat',mfilename));
do_sim = ~exist(fsum,'file');

if do_sim
    specs = {'Full','Partial'};

    parameters = struct('nparticles',100,'x0_unc',1,'lambda_v',.2,'lambda_s',.2,'v0',.1,'s0',.1);       
    config = struct('beta',5,'specs',{specs},'N1',100,'N2',10,'omega',10^-4,'rng_id',0,'nsim',40000,'model_parameters',parameters);
    rng(config.rng_id);
    
    N = config.N1+config.N2;
    nsim = config.nsim;
    specs = config.specs;
    
    outcomes = cell(1,2);
    vols = cell(1,2);
    stcs = cell(1,2);
    lrs = cell(1,2);
    vals = cell(1,2);    
    for j=1:2    
        outcome = nan(N,nsim);
        vol = nan(N,nsim);
        stc = nan(N,nsim);
        lr  = nan(N,nsim);
        val = nan(N,nsim);        
        for i=1:nsim
            outcome(:,i) = timeseries(config,specs{j});
            [vol(:,i),stc(:,i),lr(:,i),val(:,i)]=model_pf(outcome(:,i),config.model_parameters);            
        end
        outcomes{j} = outcome;
        vols{j} = vol;
        stcs{j} = stc;
        lrs{j} = lr;
        vals{j} = val;
    end
    sim = struct('config',config,'specs',{specs},...
                 'vols',{vols},'stcs',{stcs},'lrs',{lrs},'vals',{vals});
    save(fname,'sim');    
end

if ~exist(fsum,'file')
    sim = load(fname); sim = sim.sim;

    N1 = sim.config.N1;
    N2 = sim.config.N2;
    nsim = sim.config.nsim;

    beta = sim.config.beta;

    N = N1+N2;

    mp = zeros(7,2);
    ep = zeros(7,2);
    a = nan(nsim,2);
    u = nan(nsim,2);
    v = nan(nsim,2);
    for j=1:2    
        t = N-N2+1;
        a(:,j) = sim.lrs{j}(t,:);
        v(:,j) = sim.vols{j}(t,:);
        u(:,j) = sim.stcs{j}(t,:);

        t = (N-N2)+(0:6);
        m = sim.vals{j}(t,:);
        p = 1./(1+exp(-beta*m));
        mp(:,j) = median(p,2);        
    end

    mlr = mean(a,1);
    elr = serr(a,1);
    mvol = mean(v,1);
    evol = serr(v,1);
    mstc = mean(u,1);
    estc = serr(u,1);
    
    sim = struct('config',sim.config,'specs',{sim.specs},'mvol',mvol,'mstc',mstc,'mlr',mlr,'evol',evol,'estc',estc,'elr',elr,'mp',mp); %#ok<NASGU>    
    save(fsum,'sim');
end

[dat_full, dat_partial ] = data_partial_reinforcement;

sim = load(fsum); sim = sim.sim;
mvol = sim.mvol;
mstc = sim.mstc;
mlr = sim.mlr;
evol = sim.evol;
estc = sim.estc;
elr = sim.elr;
mp = sim.mp;

%--------------------------------------------------------------------------

if nargin<1
    close all;    
    nr = 1;
    nc = 3;        
    subplots = 1:3;
    fsiz = [0 0 .5 .2]; 
    figure; set(gcf,'units','normalized'); set(gcf,'position',fsiz);
    
    do_supp = 1;    
end
labels = sim.specs;
xstr = {def('lr'), def('vol'), def('stc')};

yll2 = [0 1];
yll3 = [0 1];

h = plot_bar(nr,nc,subplots(1:3),{mlr,mvol,mstc},{elr,evol,estc},labels,xstr);
set(h(2),'ylim',yll2);
set(h(3),'ylim',yll3);

%--------------------------------------------------------------------------
% Supp fig

mx = [dat_full dat_partial];

if ~do_supp
    return;
end

fsy = def('fsy');
abc = def('abc');

nr = 1;
nc = 2;
subplots = 1:2;
fsiz = [0 0 .4 .25];
figure; set(gcf,'units','normalized'); set(gcf,'position',fsiz);

[h, hp] = plot_signal(nr,nc,subplots(1:2),{mx,mp},{mx*0,mp*0},{'Relative score','Response probability'},'',nan,[],abc);
legend(hp(1,:),labels,'fontsize',fsy,'location','northeast','box','off','AutoUpdate','off');
title(h(1),'Data');
title(h(2),'Model');
end


function [o,x]=timeseries(config,mode)
N1 = config.N1;
N2 = config.N2;
omega = 10^-4;

if strcmpi(mode,'partial')
    b = N1/N1;
    x = [];
    for i=1:b
        xi = 2*ones(N1,1);
        j = randperm(N1);
        j = j(1:(N1/2));
        xi(j) = 0;
        x = [x; xi]; %#ok<AGROW>
    end
elseif strcmpi(mode,'full')
    x = 2*repmat([.5 ;.5], N1/2, 1);
end
x = [x; zeros(N2,1)];

o = x + sqrt(omega)*randn(size(x));

end
