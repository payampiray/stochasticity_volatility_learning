function h = sim_lesioned(nr,nc,subplots)

fname = fullfile('simulations',sprintf('%s.mat',mfilename));
fsum = fullfile('sum',sprintf('%s.mat',mfilename));
do_sim = ~exist(fsum,'file');

if do_sim
    mkdir('simulations');
    fsim = fullfile('simulations','sim_2x2.mat');
    sim2x2 = load(fsim); sim2x2 = sim2x2.sim;
    config = sim2x2.config;
    outcomes = sim2x2.outcomes;
    
    nsim = config.nsim;
    config.model_parameters.v0_lesioned = config.v0;
    config.model_parameters.s0_lesioned = config.s0;
    
    N = config.N;    

    lnames = {'Healthy'; sprintf('%s lesion',def('stc')); sprintf('%s lesion',def('vol'))};        
    config.lnames = lnames;
    
    vols = cell(length(lnames),4);
    stcs = cell(length(lnames),4);
    lrs  = cell(length(lnames),4);
    vals = cell(length(lnames),4);   
    
    vols(1,:) = sim2x2.vols;
    stcs(1,:) = sim2x2.stcs;
    lrs(1,:)  = sim2x2.lrs;
    vals(1,:) = sim2x2.vals;

    for l=2:3
        rng(config.rng_id);
        
        vol = nan(N,nsim);
        stc = nan(N,nsim);
        lr  = nan(N,nsim);
        val = nan(N,nsim);            
        for j=1:4
            for i=1:nsim
                [vol(:,i),stc(:,i),lr(:,i),val(:,i)]=model_pf(outcomes{j}(:,i),config.model_parameters,lnames{l});
            end
            vols{l,j} = vol;
            stcs{l,j} = stc;
            lrs{l,j} = lr;
            vals{l,j} = val;
        end
    end
    
    sim = struct('config',config,'specs',{config.specs},...                                
                 'vols',{vols},'stcs',{stcs},'lrs',{lrs},'vals',{vals}); %#ok<NASGU>
    save(fname,'sim');
end

if ~exist(fsum,'file')
    sim = load(fname); sim = sim.sim;
    
    N = sim.config.N;
    [L,J] = size(sim.lrs);
    nsim = sim.config.nsim;

    ttend = (.9*N+1):N;
    mlr = cell(1,L);
    mvol = cell(1,L);
    mstc = cell(1,L);
    elr = cell(1,L);
    evol = cell(1,L);
    estc = cell(1,L);

    for l=1:L
        a = nan(nsim,J);
        v = nan(nsim,J);
        u = nan(nsim,J);
        for j=1:J
            a(:,j) = mean(sim.lrs{l,j}(ttend,:),1);
            v(:,j) = mean(sim.vols{l,j}(ttend,:),1);
            u(:,j) = mean(sim.stcs{l,j}(ttend,:),1);
        end
        ma(l,:) = mean(a);
        mv(l,:) = mean(v);
        ms(l,:) = mean(u);    

        ea(l,:) = serr(a);
        ev(l,:) = serr(v);
        es(l,:) = serr(u);    

    end
    
%     sim = struct('config',sim.config,'specs',{sim.specs},'mlr',{mlr},'mvol',{mvol},'mstc',{mstc},'elr',{elr},'evol',{evol},'estc',{estc}); %#ok<NASGU>   
    sim = struct('config',sim.config,'specs',{sim.specs},'ma',ma,'mv',mv,'ms',ms,'ea',ea,'ev',ev,'es',es); %#ok<NASGU>   
    save(fsum,'sim');
end
sim = load(fsum); sim = sim.sim;

ma = sim.ma;
mv = sim.mv;
ms = sim.ms;

ea = sim.ea;
ev = sim.ev;
es = sim.es;

ii1 = [1 3];
ii2 = [2 4];

for l=1:3
    mlr{l} = [ma(l,ii1)' ma(l,ii2)']';
    mvol{l} = [mv(l,ii1)' mv(l,ii2)']';
    mstc{l} = [ms(l,ii1)' ms(l,ii2)']';

    elr{l} = [ea(l,ii1)' ea(l,ii2)']';
    evol{l} = [ev(l,ii1)' ev(l,ii2)']';
    estc{l} = [es(l,ii1)' es(l,ii2)']';    
end

%--------------------------------------------------------------------------
if nargin<1
    close all;    
    nr = 2;
    nc = 3;
    subplots = 1:9;
    fsiz = [0 0 .7 .55];
    figure; set(gcf,'units','normalized'); set(gcf,'position',fsiz);
end

colstrs = {'Small','Large'};
levels = {'Small','Large'};
lgtitle = sprintf('True %s',lower(def('stc')));
xltitle = sprintf('True %s',lower(def('vol')));

lnames = sim.config.lnames;
xstr = {def('lr'), def('vol'), def('stc')};

ylb = {[0 .75],[0 .75],[0 .75]};

fsy = def('fsy');

ii = [1 2 3];
lnames = lnames(ii);
hx = plot_bar(nr,nc,subplots(1:3),mlr(ii),elr(ii),colstrs,repmat(xstr(1),1,3),'',[],ylb);
for j=1:3
    title(hx(j),lnames{j},'fontsize',fsy);
    xlabel(hx(j),xltitle,'fontsize',fsy);
end
lg = legend(hx(1),levels,'fontsize',fsy,'location','northwest','box','off');
title(lg,lgtitle);
pos = lg.Position;
pos(1) = pos(1)*1.1;
set(lg,'position',pos);
h(1:3) = hx;

hx = plot_bar(nr,nc,subplots(5:6),{mvol{2},mstc{3}},{evol{2},estc{3}},colstrs,xstr([2 3]));
for j=1:2
    title(hx(j),lnames{j+1},'fontsize',fsy);
    xlabel(hx(j),xltitle,'fontsize',fsy);
end

h(5:6) = hx;
end

