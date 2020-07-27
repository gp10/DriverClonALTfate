function [nx_basal,nx_total,ntime] = MonteCarloSimulator_SP_total_varMu_evenStart(rtime,lambda,dens,r,gamma,micetime,mice_m_t,indiv,tlag,GamShape)
%% Non-Markovian Monte Carlo simulator of Single Progenitor (SP) model dynamics
% Total clone sizes are simulated over time starting from an even initial
% induction condition where both proliferating and differentiating and
% surpabasal cells can be induced. It allows to consider gamma-distributed
% cell cycle periods and a time-dependent shedding rate (mu).

% from Herms et al, 2020

%% Input:
% rtime: horizontal vector containing desired time points (expressed in weeks)
% lambda: average division rate (/week)
% dens: fraction of proliferating (progenitor) cells in the basal layer
% r: probability of symmetric division leading to two progenitors, PP, or two differentiating cells, DD
% gamma: stratification rate (/week)
% micetime: 1xn vector of experimental time points (weeks) when the suprabasal-to-basal cell ratio has been measured
% mice_m_t: 1xn vector of experimental suprabasal-to-basal cell ratio measurements
% indiv: number of independent simulations (number of simulated clones)
% tlag: refractory period between consecutive divisions (minimum cell-cycle period) (weeks)
% GamShape: 'Shape' parameter of the gamma-distributed cell-cycle period (=1 for the default exponential distribution)

%% Output:
% nx_basal: mxn matrix of clone sizes (No. of basal cells per clone) over time (m clones x n time points)
% nx_total: mxn matrix of clone sizes (No. of total cells per clone) over time (m clones x n time points)
% ntime: horizontal vector of the n time points collected

%% Example:
% rtime = [1.4286 4 12 24]; %(weeks)
% lambda = 2.9;
% dens = 0.52;
% r = 0.08;
% gamma = 3.1417;
% micetime = [1.4286 4 12 24]; %(weeks)
% mice_m_t = [0.7935 1.0417 0.8310 0.5574];
% indiv = 1000;
% tlag = 0.5/7;
% GamShape = 8;
% [nx_basal,nx_total,ntime] = MonteCarloSimulator_SP_total_varMu_evenStart(rtime,lambda,dens,r,gamma,micetime,mice_m_t,indiv,tlag,GamShape);

%% Default parameter values:
tic
% Initial error checks:
if (nargin < 10)
    indiv=100;
    tlag = 0;
    GamShape = 1;
end

% Initial definition of parameters:
timelim=rtime(1,end); % time limit
ntime = rtime; % all collection time points
nx=zeros(indiv,length(rtime),3); % stores No. of cells of each type / clone
nx_basal=zeros(indiv,length(rtime)); % stores No. of basal cells / clone
nx_total=zeros(indiv,length(rtime)); % stores No. of total cells / clone
% GamScale (Gamma-dist scale param. that fits the observed average division rate - considering tlag)
GamScale = (1/lambda - tlag) ./ GamShape;

%% ITERATION FOR DIFFERENT INDIVIDUAL CLONES
for it=1:indiv

    % Initial variables & cell attributes:
    p=zeros(5,1);
    Id_Cell = 1; % Cell identifier within the clone
    tini = 0;
    tend = 0;

    % Initial Suprabasal-to-basal cell ratio:
    m = interp1(micetime,mice_m_t,tini,'pchip');

    % Initial cell content of the clone:
    %           P      D     S
    x=mnrnd(1,[dens (1-dens) m]./(1+m)); % only single-cell clones with just one cell (can be proliferating / differentiated basal / suprabasal).

    % Asynchronous initial cond (defining a random value for the 'time-to-next-division' of initial cells at time 0)
    GamEstTime = [0:0.001:5]; % (wide time range for which to estimate the Gam-related integral)
    GamCumF = 1-gamcdf(GamEstTime,GamShape,GamScale);
    Prob_cutoff = tlag ./ (tlag + 0.001.*trapz(GamCumF)); % defined by the density of each of the 2 cdfs
    myr1 = rand;
    if myr1 < Prob_cutoff
        tDiv0 = rand*tlag; % random sampling from a uniform dist.
    else
        tDiv0 = tlag + GamEstTime(1,find(mnrnd(1,GamCumF./sum(GamCumF))==1)) + rand*0.001; % random (multinomial) sampling from a gamma-related cum. dist.
    end 

    % ITERATION FOR EACH SINGLE CELL WITHIN THE CLONE:
    while Id_Cell <= size(x,1)

        if tini(Id_Cell,1) <= timelim

            % Calculate Suprabasal-to-basal cell ratio at given time "tini(Id_Cell,1)":
            m = interp1(micetime,mice_m_t,tini(Id_Cell,1),'pchip');
            mu = dens*lambda/m; % week-1

            % Calculation of single-event probabilities:
            p(1) = lambda*r*x(Id_Cell,1); % P -> P + P
            p(2) = lambda*(1-2*r)*x(Id_Cell,1); % P -> P + D
            p(3) = lambda*r*x(Id_Cell,1); % P -> D + D
            p(4) = gamma*x(Id_Cell,2); % D -> S
            p(5) = mu*x(Id_Cell,3); % S -> lost
            % Calculation of total probability of event:
            pt=sum(p);

            % Calculate time to new event:
            if (x(Id_Cell,1)~=0) && (Id_Cell==1) % if founder progenitor
                tend(Id_Cell,1)=tini(Id_Cell,1)+tDiv0;
            elseif (x(Id_Cell,1)~=0) && (Id_Cell~=1) % any other progenitor
                tend(Id_Cell,1)=tini(Id_Cell,1)+gamrnd(GamShape,GamScale)+tlag;
            else % if differentiating (D) or suprabasal (S) cell
                tau=-(1./pt)*log(rand);
                tend(Id_Cell,1)=tini(Id_Cell,1)+tau;
            end

            if tend(Id_Cell,1) < timelim
                % Event selection:
                event = find(cumsum(p)>(rand*pt),1);
                if (event==1)
                    x = [x; 1 0 0; 1 0 0];
                    tini = [tini; tend(Id_Cell,1); tend(Id_Cell,1)];
                elseif (event==2)
                    x = [x; 1 0 0; 0 1 0];
                    tini = [tini; tend(Id_Cell,1); tend(Id_Cell,1)];
                elseif (event==3)
                    x = [x; 0 1 0; 0 1 0];
                    tini = [tini; tend(Id_Cell,1); tend(Id_Cell,1)];
                elseif (event==4)
                    x = [x; 0 0 1];
                    tini = [tini; tend(Id_Cell,1)];
                %elseif (event==5) | loss is implicit in the time update
                end
            end

        end
        Id_Cell = Id_Cell + 1;

    end

    % Save the populations of cells at certain time points:
    for bas = 1:size(rtime,2)
        nx(it,bas,1) = size(find( (tini <= rtime(1,bas)) & (rtime(1,bas) < tend) & x(:,1)==1 ),1);
        nx(it,bas,2) = size(find( (tini <= rtime(1,bas)) & (rtime(1,bas) < tend) & x(:,2)==1 ),1);
        nx(it,bas,3) = size(find( (tini <= rtime(1,bas)) & (rtime(1,bas) < tend) & x(:,3)==1 ),1);
    end

end

%% Sum both types of basal cells (P+D) to get basal-layer clone sizes:
nx_basal = nx(:,:,1)+nx(:,:,2);
% Sum the 3 types of cells (P+D+S) to get total clone sizes:
nx_total = nx(:,:,1)+nx(:,:,2)+nx(:,:,3);
toc
