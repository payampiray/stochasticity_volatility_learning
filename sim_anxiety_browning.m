function h = sim_anxiety_browning(nr,nc,subplots)


fname = fullfile('simulations',sprintf('%s.mat',mfilename));
fsum = fullfile('sum',sprintf('%s.mat',mfilename));
do_sim = ~exist(fsum,'file');

if do_sim
    mkdir('simulations');
    [~,tvolatile,tstable] = timeseries;

    trait_mean = [.5 1 3];
    trait_max = 4;

    
    parameters = struct('nparticles',100,'x0_unc',1,'lambda_v',[],'lambda_s',[],'v0',.001,'s0',.001);
    
    config = struct('n',30,'tvolatile',tvolatile,'tstable',tstable,'rng_id',0,'nsim',1000,...
                    'trait_mean',trait_mean,'trait_max',trait_max,'model_parameters',parameters);
    rng(config.rng_id);    
            
    n = config.n;
    N = length(tvolatile);
    nsim = config.nsim;
    
    lrs = cell(1,nsim);
    lambda_vs = cell(1,nsim);
    lambda_ss = cell(1,nsim);
        
    rr_min = trait_max^-1;
    rr = 2*(trait_mean.^-1 -rr_min);    
    for j=1:nsim
        lambda_v = .2*rand(1,n);
        
        K = length(rr);
        s2v = nan(1,0);
        for k=1:K
            s2v = [s2v rr_min+rr(k)*rand(1,n/length(trait_mean))]; %#ok<AGROW>
        end
                
%         lambda = trait.^-1;
%         lmin = (rr_min+rr).^-1;
%         lmean = (rr_min+.5*rr).^-1;  

        lambda_s = s2v.*lambda_v;        
        
        lr = nan(N,n);
        for i=1:n            
            [o] = timeseries;
            parameters = config.model_parameters;
            parameters.lambda_v = lambda_v(i);
            parameters.lambda_s = lambda_s(i);
            [~,~,lr(:,i)] = model_pf(o,parameters);            
        end
        lrs{j} = lr;
        lambda_vs{j} = lambda_v;
        lambda_ss{j} = lambda_s;    
    end
    sim = struct('config',config,...
                 'lambda_v',{lambda_vs},'lambda_s',{lambda_ss},'lr',{lrs}); %#ok<NASGU>
    save(fname,'sim');    
end

if ~exist(fsum,'file')
    sim = load(fname); sim = sim.sim;

    tvolatile = sim.config.tvolatile;
    tstable = sim.config.tstable;

    nsim = sim.config.nsim;

    c = nan(nsim,1);
    for i=1:nsim
        lp = sim.lambda_v{i}'./sim.lambda_s{i}';
        av = log10(sim.lr{i}(tvolatile,:));
        as = log10(sim.lr{i}(tstable,:));
        lrv = mean(av,1)';
        lrs = mean(as,1)';
        lr = lrv-lrs;

        c(i) = corr(lp,lr,'type','spearman');    

        if i==1
            model_anxiety = lp;
            dlr = lr;
        end
    end
    mc = median(c);
    ec = se_median(c);
    
    sim = struct('config',sim.config,'model_anxiety',model_anxiety,'dlr',dlr,'mc',mc,'ec',ec); %#ok<NASGU>
    save(fsum,'sim');     
end
sim = load(fsum); sim = sim.sim;

mc = sim.mc;
ec = sim.ec;
model_anxiety = sim.model_anxiety;
dlr = sim.dlr;

[data_anxiety, relative_llr] = data_anxiety_browning;

%--------------------------------------------------------------------------
if nargin<1
    close all;    
    nr = 1;
    nc = 2;
    subplots = 1:2;
    fsiz = [.3 .3 .5 .4];
    figure; set(gcf,'units','normalized'); set(gcf,'position',fsiz);
end

xstr = {def('lr'), def('vol'), def('stc')};

alf = def('alf');
col = def('col_br');
fsy = def('fsy');
fs = def('fs');


ylbl = sprintf('Relative log learning rate\n (volatile block - stable block)');
xlbls = {'Trait anxiety',sprintf('Model trait anxiety\n(%s to %s update rate)',lower(xstr{2}),lower(xstr{3}))};
ylbl_inst = sprintf('Correlation');
ttls = {'Data','Model'};

x = {data_anxiety,model_anxiety};
y = {relative_llr,dlr};
xl = {[15 75],[0 4]};


for i=1:2
    h(i) = subplot(nr,nc,subplots(i));
    scatter(x{i},y{i},[],col(1,:),'filled','marker','o');    
    xlim(xl{i});
    hl = lsline;
    set(hl,'color',col(1,:));
%     scatterplot(x{i},y{i},conf);    

    ylabel(ylbl,'fontsize',fsy);
    xlabel(xlbls{i},'fontsize',fsy);
    title(ttls{i},'fontsize',fsy);
    axis square;
end
% ylim([-.5 1]);


% create smaller axes in top right, and plot on it
% axes(h(2))
h_in = axes('Position',[.8 .68 .1 .2]);
errorbarKxN(mc,ec,{''},col(2,:),.2);    
set(h_in,'box','on','LineWidth',1);

ylim([-0.5 0]);
set(gca,'fontsize',fs);
alpha(alf);
xaxes = get(gca,'XAxis');
set(xaxes,'fontsize',fsy);    
title(ylbl_inst,'fontsize',fsy,'fontweight','normal');

end

function [y,tvolatile,tstable]=timeseries
n = 20;
m = [.75 .75*ones(1,5) .2 .8 .2 .8 .2];
omega = 0.01;

dm = [0 diff(m)];
dm = abs(dm)>0.05;
N = length(m)*n;

x = nan(N,1);
z = nan(N,1);

xi = zeros(n*6,1);
j = randperm(n*6);
n1 = m(1)*n*6;
j = j(1:n1);
xi(j) = 1;
ii = (1:n*6);
z(ii) = m(1); 
x(ii) = xi; 


tvolatile = nan(N,1);
for i=7:length(m)
    xi = zeros(n,1);
    j = randperm(n);
    n1 = m(i)*n;
    j = j(1:n1);
    xi(j) = 1;
    
    ii = (i-1)*n+ (1:n);
    z(ii) = m(i); 
    x(ii) = xi; 
    tvolatile(ii) = dm(i)*ones(n,1);
end

tvolatile = tvolatile==1;
tstable = ~tvolatile;
tstable(1:n ) = 0;

y = x + sqrt(omega)*randn(N,1);

end
