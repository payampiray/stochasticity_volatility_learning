function h = sim_amygdala_costa(nr,nc,subplots,do_supp)


fname = fullfile('simulations',sprintf('%s.mat',mfilename));
fsum = fullfile('sum',sprintf('%s.mat',mfilename));
do_sim = ~exist(fsum,'file');

if do_sim
    mkdir('simulations');
    parameters = struct('nparticles',100,'x0_unc',1,'lambda_v',.1,'lambda_s',.1,'v0',.5,'s0',.5,'v0_lesioned',0.01);                            
    config = struct('p',[.6 1],'N',80,'beta',3,'beta_lesioned',1,...
                    'rng_id',0,'nsim',1000,'model_parameters',parameters);
                
    rng(config.rng_id);
    nsim = config.nsim;
    N = config.N;
    p = config.p;
    
    outcomes = cell(1,4);
    vols = cell(1,4);
    stcs = cell(1,4);
    lrs = cell(1,4);
    vals = cell(1,4);    
    k = 0;
    lnames = {'Control','Lesioned'};
    lesion_model = {'Healthy', sprintf('%s lesion',def('vol'))};                
    jnames = {'Stochastic','Deterministic'};    
    for l=1:2
        for j=1:2         
            outcome = nan(N,nsim);
            vol = nan(N,nsim);
            stc = nan(N,nsim);
            lr  = nan(N,nsim);
            val_acq = nan(N,nsim);
            for i=1:nsim                
                outcome(:,i) = timeseries(config,p(j));
                [vol(:,i),stc(:,i),lr(:,i),val_acq(:,i)]=model_pf(outcome(:,i),config.model_parameters,lesion_model{l});
            end
            k = k+1;
            specs(:,k) = [{sprintf('%s-%s',lnames{l},jnames{j})} lnames(l) jnames(j) sprintf('%d/%d',p(j)*100,(1-p(j))*100) k]; %#ok<AGROW>
            
            outcomes{k} = outcome;        
            vols{k} = vol;
            stcs{k} = stc;
            lrs{k} = lr;
            vals{k} = val_acq;
        end
    end
    sim = struct('config',config,'specs',{specs},...
                 'vols',{vols},'stcs',{stcs},'lrs',{lrs},'vals',{vals}); %#ok<NASGU>    
    save(fname,'sim');    
end

if ~exist(fsum,'file')
    sim = load(fname); sim = sim.sim;
    

    J = 2;
    beta = [sim.config.beta*ones(1,J) sim.config.beta_lesioned*ones(1,J)];

    nsim = sim.config.nsim;
    N = sim.config.N;
    N1 = N/2;
    N2 = N/2;
    vol = nan(N,2);
    stc = nan(N,2);
    lr = nan(N,2);
    val = nan(N,2);
    e_vol = nan(N,2);
    e_stc = nan(N,2);
    e_lr = nan(N,2);
    e_val = nan(N,2);

    val_acq = nan(nsim,2*J);
    val_rev = nan(nsim,2*J);
    

    t_acq = 1:N1;
    t_rev = (N1+1):(N1+N2);    

    for j=1:length(sim.vols)

        v   = sim.vals{j};
        dv  = v;
        v   = 1./(1+exp(-beta(j)*dv));    
        
        val_acq(:,j) = mean(v(t_acq,:),1)';
        val_rev(:,j) = mean(v(t_rev,:),1)';

        val(:,j) = mean(v,2);
        e_val(:,j) = serr(v,2);
        
        vol(:,j) = mean(sim.vols{j},2);
        stc(:,j) = mean(sim.stcs{j},2);
        lr(:,j) = mean(sim.lrs{j},2);    

        e_vol(:,j) = serr(sim.vols{j},2);
        e_stc(:,j) = serr(sim.stcs{j},2);
        e_lr(:,j) = serr(sim.lrs{j},2);    
    end

    mval_acq = mean(val_acq);
    eval_acq = serr(val_acq);

    mval_rev = 1-mean(val_rev);
    eval_rev = serr(1-val_rev);
    
    sim = struct('config',sim.config,'specs',{sim.specs},'vol',vol,'stc',stc,'lr',lr,'val',val,'e_vol',e_vol,'e_stc',e_stc,'e_lr',e_lr,'e_val',e_val,...
                'mval_acq',mval_acq,'eval_acq',eval_acq,'mval_rev',mval_rev,'eval_rev',eval_rev); %#ok<NASGU>
    save(fsum,'sim');
end
sim = load(fsum); sim = sim.sim;
data = data_amygdala_costa;

%--------------------------------------------------------------------------

specs = sim.specs;
vol = sim.vol;
stc = sim.stc;
lr = sim.lr;
val = sim.val;
e_val = sim.e_val;

mval_acq = sim.mval_acq;
eval_acq = sim.eval_acq;
mval_rev = sim.mval_rev;
eval_rev = sim.eval_rev;


pnames = specs(3,1:2);
gnames = specs(2,[1 3]);
% cnames = {'Acquisition','Reversal'};
cnames = {'Acq','Rev'};

mdat_acq  = data.x_acq;
edat_acq  = data.e_acq;

mdat_rev  = data.x_rev;
edat_rev  = data.e_rev;

N = sim.config.N;
N1 = N/2;

p1 = [sim.config.p(1)*ones(N1,1); (1-sim.config.p(1))*ones(N1,1); ];
p2 = [sim.config.p(2)*ones(N1,1); (1-sim.config.p(2))*ones(N1,1); ];
p = [p1 p2];
%--------------------------------------------------------------------------

if nargin<1
    close all;    
    nr = 3;
    nc = 3;            
    fsiz = [0 0 .7 .8];
    subplots = 1:9;
    figure; set(gcf,'units','normalized'); set(gcf,'position',fsiz);
    do_supp = 1;
end

col0 = def('col');
col = def('col_bp');
fn = def('fn');
fsy = def('fsy');
abc = def('abc');

h(1) = subplot(nr,nc,subplots(1));
[h(1),hp]=plot_signal(nr,nc,subplots(1),{p},{[]},{sprintf('Probability of shape 1\n being correct')},'',N1,[0 1],[],col0);
title('Reversal learning task','fontsize',fsy);
legend(hp,pnames,'fontsize',fsy,'location','northeast','box','off');


h(2) = subplot(nr,nc,subplots(2));
ylbl = 'Fraction of correct choice';
plot_fraction(mdat_acq,mdat_rev,edat_acq,edat_rev,gnames,cnames,pnames,ylbl,1);
title('Data');

h(3) = subplot(nr,nc,subplots(3));
ylbl = 'Probability of correct choice';
plot_fraction(mval_acq,mval_rev,eval_acq,eval_rev,gnames,cnames,pnames,ylbl,0);
% % text(xsA,ysA,abc(3),'fontsize',fsA,'Unit','normalized','fontname',fn);        
title('Model');

xstr = {def('lr'), def('vol'), def('stc')};

ii = [2 4];
[h(4:6),hp] = plot_signal(nr,nc,subplots(4:6),{lr(:,ii),vol(:,ii),stc(:,ii)},{[],[],[]},xstr,'Deterministic',N1,[],[],col);
legend(hp(2,:),gnames,'fontsize',fsy,'location','northeast','box','off');

ii = [1 3];
[h(7:9)] = plot_signal(nr,nc,subplots(7:9),{lr(:,ii),vol(:,ii),stc(:,ii)},{[],[],[]},xstr,'Stochastic',N1,[],[],col);

%--------------------------------------------------------------------------
% Supp fig

if do_supp
    nr = 1;
    nc = 2;
    subplots = 1:2;
    fsiz = [0 0 .35 .25];
    figure; set(gcf,'units','normalized'); set(gcf,'position',fsiz);

    ii1 = [2 4];
    ii2 = [1 3];
    ylbl = 'Probability of correct choice';
    yll = [0 1];

    [hs,hp] = plot_signal(nr,nc,subplots(1:2),{val(:,ii1),val(:,ii2)},{e_val(:,ii1),e_val(:,ii2)},{ylbl,ylbl},'',N1,{yll,yll},abc,col);
    lg = legend(hp(2,:),gnames,'fontsize',fsy,'location','northeast','box','off');
    title(hs(1),'Deterministic','fontsize',fsy,'fontname',fn);
    title(hs(2),'Stochastic','fontsize',fsy,'fontname',fn);
    % set(hs(1),'box','off');
    % set(hs(2),'box','off');
end
    
end

function [h,he, col]=plot_fraction(mval_acq,mval_rev,eval_acq,eval_rev,gnames,cnames,pnames,ylbl,do_leg)
np = 2;
col = def('col_bp'); 

fs = def('fs');
fn = def('fn');
fsy = def('fsy');

for l=1:np
    ll = (l-1)*np+(1:np);
    x = {mval_acq(ll), mval_rev(ll)};
    e = {eval_acq(ll), eval_rev(ll)};
    linetype = {'-','--'};
    for j=1:2        
        h(l,j) = plot(x{j},linetype{j},'color',col(l,:),'linewidth',2); hold on;
        he(l,j) = errorbar(x{j},e{j},linetype{j},'color',col(l,:),'linewidth',2); hold on;
        gcnames{l,j} = sprintf('%s-%s',gnames{l},cnames{j});
    end    
    xlim([1-.75 np+.75]);    
    set(gca,'xtick',[1 np],'xticklabel',pnames);       
end
set(gca,'ylim',[.5 1],'ytick',0.5:.1:1,'fontsize',fs,'fontname',fn);
ylabel(ylbl,'fontsize',fsy);
xaxes = get(ancestor(gca, 'axes'),'XAxis');
set(xaxes,'fontsize',fsy,'ticklength', [0 0]);
% lg = legend(h(:,1),gnames,'fontsize',fsy,'location','northwest');

hh = [h(1:2) h(3:4)];
if do_leg
legend(hh,gcnames,'location','northwest','fontsize',fsy,'box','off');
end
end

function [h, hp] = plot_signal(nr,nc,sub_plts,x,e,mstr,all_title,N,yll,abc,col)
if nargin<6, mstr = []; end
if nargin<7, all_title = []; end
if nargin<8, N = nan; end
if nargin<9, yll = []; end
if nargin<10, abc = ''; end
if nargin<11, col = ''; end

fs = def('fs');
fn = def('fn');
fsy = def('fsy');
alf = def('alf');
fsA = def('fsA');
xsA = def('xsA');
ysA = def('ysA');

linestyle = {'-';'-'};
if isempty(col)
    col = def('col'); %col(1,:) = [0 0 0];
    if size(x{1},2)>2
        col = [col;col];
        linestyle = [linestyle; '--'; '--'];
    end
end
if isempty(abc)
    abc = [];
%     abc = abc(sub_plts);
end

if ~iscell(yll)
    yll = repmat({yll},size(x));
end

% if size(x{1},2)==1
%     col = [0 0 0];
% end

for j=1:length(x)
    h(j) = subplot(nr,nc,sub_plts(j));
    set(gca,'fontsize',fs,'fontname',fn);
    set_yl = 0;
    if isempty(yll{j})
        if ~isempty(e{j})
            maxy = max(max(x{j}+e{j}));
            miny = min(min(x{j}-e{j}));
        else
            maxy = max(max(x{j}));
            miny = min(min(x{j}));
        end
        
        bias = 0;
        if (maxy - miny)==0
            bias = .01;
        end
        yl = [miny maxy]+[-bias bias];
    else
        yl = yll{j};
        set_yl = 1;        
    end
    if ~isnan(N)
        yl0 = 1;
        while sum(abs(yl-yl0))>0
            yl0 = yl;
            plot([N N],yl+[10^-30 -10^-30],'color',.8*[1 1 1],'linewidth',2); hold on;
            yl = get(h(j),'ylim');
        end        
    end    
    
    for i=1:size(x{j},2)
        if isempty(e{j})
            tt = 1:N;
            z = nan(size(x{j},1),1);
            z(tt) = x{j}(tt,i);
            
            hp(j,i) = plot(z,'linewidth',2,'color',col(i,:),'linestyle',linestyle{i}); hold on;
            xlim([0 length(z)]);                        
            
            tt = (N):size(x{j}(:,i),1);
            z = nan(size(x{j},1),1);
            z(tt) = x{j}(tt,i);            
            plot(z,'linewidth',2,'color',col(i,:),'linestyle','--'); hold on;
            
        else
            [hp(j,i),hf(j,i)] = plot_ci(x{j}(:,i),e{j}(:,i),col(i,:)); hold on;
            set(hp(i),'linewidth',2,'linestyle',linestyle{i});
        end
    end
    nt = size(x{1},1);
    xlim([0 nt]);
    
    if set_yl
        set(h(j),'ylim',yl);
    end
    
    if ~isempty(mstr{j})
        ylabel(mstr{j},'fontsize',fsy);
    end
    xlabel('Trial','fontsize',fsy);
    
%     if ~isempty(glbls)
%         legend(glbls,'Interpreter','latex','fontsize',fsy);
%     end
        
    if ~isempty(abc)
        text(xsA,ysA,abc(j),'fontsize',fsA,'Unit','normalized','fontname',fn);    
    end
end

if ~isempty(all_title)
    text(h(1),-.25,.5,all_title,'fontsize',fsA,'Unit','normalized','HorizontalAlignment','Center','rotation',90);
end

end

function [y]=timeseries(config,p)
N = config.N;
omega = 10^-6;

N = N/2;
n = p*N;
x01 = [ones(1,n) zeros(1,N-n)]';
x02 = 1-x01;

j = randperm(N);
x1 = x01(j);

j = randperm(N);
x2 = x02(j);

x = [x1; x2];
x(x==0) = -1;

N = length(x);
y = x + sqrt(omega)*randn(N,1);

end
