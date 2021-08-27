function h = sim_2x2_norm(nr,nc,subplots,fig_no)

fname = fullfile('simulations',sprintf('%s.mat',mfilename));
fsum = fullfile('sum',sprintf('%s.mat',mfilename));
do_sim = ~exist(fsum,'file');

if do_sim
    factors.stc = {'Small','Small','Large','large'};
    factors.vol = {'Small','Large','Large','Small'};
    
    parameters = struct('nparticles',100,'x0_unc',100,'sigma_v',.1,'sigma_s',.1,'v0',1,'s0',2);                        
                    
    config = struct('true_vol',[.5 1.5 .5 1.5],'true_stc',[1 1 3 3],'factors',factors,'N',200,...
                    'rng_id',0,'nsim',10000,'model_parameters',parameters);

    rng(config.rng_id);
    true_vol = config.true_vol;
    true_stc = config.true_stc;
    N = config.N;
    nsim = config.nsim;    
    
    outcomes = cell(1,4);
    vols = cell(1,4);
    stcs = cell(1,4);
    lrs = cell(1,4);
    vals = cell(1,4);
    
    specs = cell(3,4);
    for j=1:4
        outcome = nan(N,nsim);        
        vol = nan(N,nsim);
        stc = nan(N,nsim);
        lr  = nan(N,nsim);
        val = nan(N,nsim);    
        for i=1:nsim
            [outcome(:,i)] = timeseries(config,true_vol(j),true_stc(j));        
            [vol(:,i),stc(:,i),lr(:,i),val(:,i)]=model(outcome(:,i),config.model_parameters);
        end
        outcomes{j} = outcome;        
        specs{1,j} = sprintf('%s %s',factors.stc{j},lower(def('stc')));
        specs{2,j} = sprintf('%s %s',factors.vol{j},lower(def('vol')));
        specs{3,j} = sprintf('true_vol=%0.2f, true_stc=%0.2f',true_vol(j),true_stc(j));        
        
        vols{j} = vol;
        stcs{j} = stc;
        lrs{j} = lr;
        vals{j} = val;
    end
    config.specs = specs;    
    sim = struct('config',config,'specs',{config.specs},'outcomes',{outcomes},...
                 'vols',{vols},'stcs',{stcs},'lrs',{lrs},'vals',{vals}); %#ok<NASGU>
    save(fname,'sim');
end

if ~exist(fsum,'file')
    sim = load(fname); sim = sim.sim;
    

    N = sim.config.N;
    nsim = sim.config.nsim;

    vol = nan(N,4);
    stc = nan(N,4);
    lr =  nan(N,4);
    e_vol = nan(N,4);
    e_stc = nan(N,4);
    e_lr =  nan(N,4);

    tend = (.9*N+1):N;
    a = nan(nsim,4);
    v = nan(nsim,4);
    u = nan(nsim,4);

    for j=1:4            
        vol(:,j) = mean(sim.vols{j},2);
        stc(:,j) = mean(sim.stcs{j},2);
        lr(:,j) = mean(sim.lrs{j},2);

        e_vol(:,j) = serr(sim.vols{j},2);
        e_stc(:,j) = serr(sim.stcs{j},2);
        e_lr(:,j) = serr(sim.lrs{j},2);

        a(:,j) = mean(sim.lrs{j}(tend,:),1);
        v(:,j) = mean(sim.vols{j}(tend,:),1);
        u(:,j) = mean(sim.stcs{j}(tend,:),1);

    end

    mv = mean(v);
    ms = mean(u);
    ma = mean(a);
    ev = serr(v);
    es = serr(u);
    ea = serr(a);    
    
    sim = struct('config',sim.config,'specs',{sim.specs},'vol',vol,'stc',stc,'lr',lr,'e_vol',e_vol,'e_stc',e_stc,'e_lr',e_lr,...
                'mv',mv,'ms',ms,'ma',ma,'ev',ev,'es',es,'ea',ea); %#ok<NASGU>
    save(fsum,'sim');    
end
sim = load(fsum); sim = sim.sim;

%--------------------------------------------------------------------------
vol = sim.vol;
stc = sim.stc;
lr = sim.lr;
e_vol = sim.e_vol;
e_stc = sim.e_stc;
e_lr = sim.e_lr;

ma = sim.ma;
mv = sim.mv;
ms = sim.ms;
ea = sim.ea;
ev = sim.ev;
es = sim.es;

ii1 = [1 3];
ii2 = [2 4];
ma = [ma(ii1)' ma(ii2)'];
mv = [mv(ii1)' mv(ii2)'];
ms = [ms(ii1)' ms(ii2)'];

ea = [ea(ii1)' ea(ii2)'];
ev = [ev(ii1)' ev(ii2)'];
es = [es(ii1)' es(ii2)'];

if nargin<1
    close all;
    nr = 3;
    nc = 2;     
    subplots = 1:6;    
    fsiz = [0 0 .6 .4];
    figure; set(gcf,'units','normalized'); set(gcf,'position',fsiz);
    fig_no = 2;    
end
colstrs = {'Small','Large'};
glbl = {sprintf('Small true %s',lower(def('vol'))),sprintf('Large true %s',lower(def('vol')))};
levels = {'Small','Large'};
lgtitle = sprintf('True %s',lower(def('stc')));
xltitle = sprintf('True %s',lower(def('vol')));

xstr = {def('lr'), def('vol'), def('stc')};

ylb = {[0 .75],[0 2],[0 3.5]};
yls = {[0.2 .81],[0 2],[0 3.5]};

fsy = def('fsy');
abc = def('abc');

if fig_no==1
    ii = ii1;
    [hx, hp] = plot_signal(nr,nc,subplots(1:2:6),{lr(:,ii),vol(:,ii),stc(:,ii)},...
                                               {e_lr(:,ii),e_vol(:,ii),e_stc(:,ii)},xstr,'',nan, yls, '', [], 0);
    title(hx(1),glbl{1},'fontsize',fsy);
    xlabel(hx(3),'Trial','fontsize',fsy);
    % set(hx(4),'xlabel','Trial');
    lg = legend(hp(2,:),levels,'fontsize',fsy,'box','off','orientation','horizontal');    
    title(lg,lgtitle);

    h(1:3) = hx;

    ii = ii2;
    [hx] = plot_signal(nr,nc,subplots(2:2:6),{lr(:,ii),vol(:,ii),stc(:,ii)},...
                                               {e_lr(:,ii),e_vol(:,ii),e_stc(:,ii)},xstr,'',nan, yls, '', [], 0);
    title(hx(1),glbl{2},'fontsize',fsy);
    xlabel(hx(3),'Trial','fontsize',fsy);

    h(4:6) = hx;
end

%--------------------------------------------------------------------------
if fig_no==2
    nr = 1;
    nc = 3;     
    subplots = 1:3;    
    fsiz = [0 0 .6 .25];
    figure; set(gcf,'units','normalized'); set(gcf,'position',fsiz);

    h(1:3) = plot_bar(nr,nc,subplots(1:3),{ma',mv',ms'},{ea',ev',es'},colstrs,xstr,abc,[],ylb);
    loc = {'northwest','northwest','north'};
    for i=1:length(h)
        xlabel(h(i),xltitle,'fontsize',fsy);
        if i==2
        lg = legend(h(i),levels,'fontsize',fsy,'location',loc{i},'box','off');
        title(lg,lgtitle);    
        end
    end
end

end


function [y,x]=timeseries(config,true_vol,true_stc)
N = config.N;
x = zeros(N,1);
for t=2:N
    x(t) = x(t-1) + sqrt(true_vol)*randn;
end
y = x + sqrt(true_stc)*randn(N,1);
end

function [vol,stc,lr,val,unc] = model(o,specs)
np = specs.nparticles;
sigma_v = specs.sigma_v;
sigma_u = specs.sigma_u;
x0_unc = specs.x0_unc;
v = specs.v0;
u = specs.u0;
[val,vol,stc,lr,unc] = pf_core_norm(o,x0_unc,sigma_v,sigma_u,v,u,np);        
end