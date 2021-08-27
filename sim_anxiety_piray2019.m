function h = sim_anxiety_piray2019(nr,nc,subplots)


fname = fullfile('simulations',sprintf('%s.mat',mfilename));
fsum = fullfile('sum',sprintf('%s.mat',mfilename));
do_sim = ~exist(fsum,'file');

if do_sim
    mkdir('simulations');
    [p_task]=data_anxiety_piray2019;
    parameters = struct('nparticles',100,'x0_unc',1,'lambda_v',.2,'lambda_s',.2,'v0',.1,'s0',.1,'s0_lesioned',0.05);
    config = struct('N',120,'p_task',p_task,'rng_id',0,'nsim',1000,'model_parameters',parameters);                 
                
    rng(config.rng_id);
    N = config.N;
    nsim = config.nsim;    
    
    vol = nan(N,nsim);
    stc = nan(N,nsim);
    lr = nan(N,nsim);
    val = nan(N,nsim);
    
    vols = cell(1,2);
    stcs = cell(1,2);
    lrs = cell(1,2);
    vals = cell(1,2);        
    lnames = {'Healthy', sprintf('%s lesion',def('stc'))};            
    groups = {'Control','Anxious'};   
    
    for j=1:2
        for i=1:nsim
            o = timeseries(config);
            [vol(:,i),stc(:,i),lr(:,i),val(:,i)] = model_pf(o,config.model_parameters,lnames{j});
        end                
        vols{j} = vol;
        stcs{j} = stc;
        lrs{j} = lr;
        vals{j} = val;
    end
     
    sim = struct('config',config,'specs',{groups},...
                 'vols',{vols},'stcs',{stcs},'lrs',{lrs},'vals',{vals}); %#ok<NASGU>
    save(fname,'sim');    
end

if ~exist(fsum,'file')
    sim = load(fname); sim = sim.sim;

    nsim = sim.config.nsim;
    ncond = length(sim.lrs);

    mx = nan(ncond,2);
    ex = nan(ncond,2);
    x = nan(nsim,ncond);
    for j=1:ncond
        for i=1:nsim
            val = sim.vals{j}(:,i);
            [x(i,:)] = performance(sim.config.p_task,val);
        end

        mx(j,:) = median(x);
        ex(j,:) = serr(x);    
    end
    
    sim = struct('config',sim.config,'specs',{sim.specs},'mx',mx,'ex',ex); %#ok<NASGU>   
    save(fsum,'sim');    
end

[~,mdat_low,mdat_high,edat_low,edat_high]=data_anxiety_piray2019;
mdat = [mdat_low; mdat_high];
edat = [edat_low; edat_high];


sim = load(fsum); sim = sim.sim;
mx = sim.mx;
ex = sim.ex;

%--------------------------------------------------------------------------
if nargin<1
    close all;
    nr = 1;
    nc = 2;
    subplots = 1:2;
    fsiz = [.3 .3 .35 .2];      
    figure; set(gcf,'units','normalized'); set(gcf,'position',fsiz);
end

col = def('col');
fsy = def('fsy');

labels = {'Stable','Volatile'};

h(1:2) = plot_bar(nr,nc,subplots(1:2),{mdat,mx},{edat,ex},{'Low anxiety','High anxiety'},{'Performance','Performance'},'',col);
set(h(1),'ylim',[0 1]);
set(h(2),'ylim',[0 1]);
legend(h(1),labels,'fontsize',fsy,'location','north','box','off');
title(h(1),'Data');
title(h(2),'Model');
end

function [y,x]=timeseries(config)
p_task = config.p_task;
omega = 0.01;
x = binornd(1,p_task);
y = x + sqrt(omega)*randn(size(x));

end

function [x, dp, tstable, tvolatile] = performance(ptask,val)

ii50 = ptask==0.5; 
ptask(ii50) = [];
dv = [0; val(1:end-1)]-.5;
p = 1./(1+exp(-dv));
p(ii50) = [];

N = length(p);
n = 10;

change_points = find(diff(ptask)~=0);
tvolatile = false(N,1);
for i=1:length(change_points)
    tvolatile(change_points(i) + (1:n)) = 1;    
end
tstable = ~tvolatile;
tstable(1:n) = 0;

corr_action = ptask>=.5;

choice = p>=.5;
perf = choice==corr_action;


mpvol = mean(perf(tvolatile,:));
mpstab = mean(perf(tstable,:));

mpvol = mean(mpvol);
mpstab = mean(mpstab);

dp = mpstab - mpvol;

x = [mpstab mpvol];
end
