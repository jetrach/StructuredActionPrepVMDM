%% objective function for forced RT model for Juliana's hierarchical task
%% JT + SDM; New Haven, CT; 01/19/2023

function ll = func_flatnostructure_dynamics(params,sub_responses,sub_PTs)

% free parameters
mu = params(1); % the MEAN time each level is prepared [hand, couplet, finger]
sigma = params(2); % variance in time each level is prepared [hand, couplet, finger]
rho = params(3); % weight on uniform guessing preparation function

% task variables
q = linspace(1,8,8)'; % 8 actions
respLevels = [1 0 0 0 0 0 0 0]; 

PTs = 0.001:.001:1.500; % vector of potential prep times (in seconds, ms resolution)
PTs = round(PTs,3); % matching data decimal places, hacky for now

% compute response action preparation functions
    i = 1; % set to 1 bc only one level in this case
    Phi(i,:) = normcdf(PTs,mu,sigma); % probability that H level has been planned by RT
    U(i,:) = ones(length(PTs),1)' * (1./length(q)); % lapse rate

% weight each response type by the action preparation functions
    coraction = 1; % define correct action in terms of each level value, all code references action 1 as correct
    cor = respLevels == coraction; % get set of correct level-specific actions
    tmp = (1-rho)*(Phi(:) .* cor)' + rho*((1./length(q)) .* ones(8,length(PTs))); % weighted mixture of correct selection and lapsing
    prespPAR = tmp./sum(tmp); % normalize mixture action selection prep to a probablity, put it in cell K
% choice types over time
presponsePAR = prespPAR;%(prespPAR{1}.*prespPAR{2}.*prespPAR{3}); % combine all preparation functions
presponsePAR = presponsePAR./sum(presponsePAR); % normalize to probs
X = [1;0;0;0;0;0;0;0]; % matrix of response codes, properly ordered! (hacky)
%% main trial loop (eventually can vectorize this step)
sim_p = nan(1,length(sub_PTs));
for n = 1:length(sub_PTs)

    % what row of simulated response matrix is sub's response on trial t? (hacky)
    resp_idx = X==sub_responses(n,1); % get row index
    pt_idx = PTs==sub_PTs(n); % what column (i.e., prep time)of simulated response matrix is sub's PT on trial t

    % not a weird trial?
    if ~isempty(pt_idx) && ~isempty(resp_idx) && sum(pt_idx)==1
        if sum(resp_idx)==1
        sim_p(n) = presponsePAR(resp_idx,pt_idx);% accumulate model probs of choice for each trial
        elseif sum(resp_idx) >1
            sim_p(n) = unique(presponsePAR(resp_idx,pt_idx));
        end

    end
    % n
end
epsilon=0.000001; % hack for avoiding 'underflow'
p = epsilon/2 + (1-epsilon)*sim_p;
ll = -nansum(log(p)); %  ll




