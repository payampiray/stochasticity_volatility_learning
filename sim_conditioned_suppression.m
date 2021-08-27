function h = sim_conditioned_suppression(nr,nc,subplots,do_supp)

fname = fullfile('simulations',sprintf('%s.mat',mfilename));
fsum = fullfile('sum',sprintf('%s.mat',mfilename));
do_sim = ~exist(fsum,'file');

if do_sim
    mkdir('simulations');
    specs = {'Control','Omission'};    
    parameters = struct('nparticles',100,'x0_unc',1,'lambda_v',.2,'lambda_s',.2,'v0',.1,'s0',.1);   
    config = struct('beta',5,'specs',{specs},'shock_mild',.3,'shock_strong',1,'N1',100,'N2',20,'N_omission',5,...
                    'rng_id',0,'nsim',40000,'model_parameters',parameters);    
                
    rng(config.rng_id);
    nsim = config.nsim;
    specs = config.specs;
    
    vols = cell(1,2);
    stcs = cell(1,2);
    lrs = cell(1,2);
    vals = cell(1,2);            
    for j=1:2        
        for i=1:nsim
            o = timeseries(config,specs{j});
            if i==1
                N = size(o,1);
                vol = nan(N,nsim);
                stc = nan(N,nsim);
                lr  = nan(N,nsim);
                val = nan(N,nsim);                
            end
            [vol(:,i),stc(:,i),lr(:,i),val(:,i)]=model_pf(o,config.model_parameters);
        end        
        vols{j} = vol;
        stcs{j} = stc;
        lrs{j} = lr;
        vals{j} = val;
    end
    sim = struct('config',config,'specs',{specs},...
                 'vols',{vols},'stcs',{stcs},'lrs',{lrs},'vals',{vals}); %#ok<NASGU>    
    save(fname,'sim');    
end

if ~exist(fsum,'file')
    sim = load(fname); sim = sim.sim;
    
    
    nsim = sim.config.nsim;
    N2 = sim.config.N2;
    beta = sim.config.beta;

    ss = 1:nsim;
    nsim = length(ss);

    u = nan(nsim,2);
    a =  nan(nsim,2);
    v = nan(nsim,2);

    mp = nan(4,2);
    for j=1:2
        t = size(sim.vols{j},1);
        t = t-N2+1;    

        tt = t+(0:3);     
        val1 = sim.vals{j}(tt,ss);    
        p = 1./(1+exp(-beta*val1));
        mp(:,j) = median(p,2);

        t = size(sim.vols{j},1);
        t = t-N2+1;
        v(:,j) = sim.vols{j}(t,ss);
        u(:,j) = sim.stcs{j}(t,ss);
        a(:,j) = sim.lrs{j}(t,ss);


    end
    mvol = mean(v);
    mstc = mean(u);
    mlr = mean(a);
    evol = serr(v);
    estc = serr(u);
    elr = serr(a);
    
    sim = struct('config',sim.config,'specs',{sim.specs},'mvol',mvol,'mstc',mstc,'mlr',mlr,'evol',evol,'estc',estc,'elr',elr,'mp',mp); %#ok<NASGU>
    save(fsum,'sim');
end
[dat_omission, dat_control]=data_conditioned_suppression;

sim = load(fsum); sim = sim.sim;
mvol = sim.mvol;
mstc = sim.mstc;
mlr = sim.mlr;
evol = sim.evol;
estc = sim.estc;
elr = sim.elr;
mp = sim.mp;

labels = sim.specs;

%--------------------------------------------------------------------------
if nargin<1
    close all;    
    nr = 1;
    nc = 3;
    fsiz = [0 0 .4 .2];
    subplots = 1:3;
    figure; set(gcf,'units','normalized'); set(gcf,'position',fsiz);
    
    do_supp = 2;
end

xstr = {def('lr'), def('vol'), def('stc')};

h(1:3) = plot_bar(nr,nc,subplots(1:3),{mlr,mvol,mstc},{elr,evol,estc},labels,xstr);
for i=1:3
    ax = ancestor(h(i), 'axes');
    xaxes = get(ax,'XAxis');
    xaxes.TickLabelRotation = 30;
end

%--------------------------------------------------------------------------
% Supp fig

mx = [dat_control; dat_omission]';

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

[h, hp] = plot_signal(nr,nc,subplots(1:2),{mx,mp},{mx*0,mp*0},{'Suppression ratio','Response probability'},'',nan,[],abc);
set(h(1),'xtick',1:4,'xlim',[.5 4]);
set(h(2),'xtick',1:4,'xlim',[.5 4]);
legend(hp(1,:),labels,'fontsize',fsy,'location','northeast','box','off','AutoUpdate','off');
title(h(1),'Data');
title(h(2),'Model');
end

function [o,x]=timeseries(config,mode)
N1 = config.N1;
N2 = config.N2;
shock1 = config.shock_mild;
shock2 = config.shock_strong;
w = 0.01;
if strcmpi(mode,'control')
    n = 0;
elseif strcmpi(mode,'omission')
    n = config.N_omission;
end
x = -[ones(N1,1)*shock1; ones(n,1)*0; ones(N2,1)*shock2];
o = x + sqrt(w)*randn(N1+n+N2,1);
end
