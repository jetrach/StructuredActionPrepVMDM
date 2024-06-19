%% Plotting response time simulations using fitted parameters
% This script with plot simulations for all 8 response types (correct and 7
% error types) and by error level. The code will plot simulations for the
% hierarchical model w one sigma value if typeFlag = 'hier' and it will
% plot simulations for the hierarchical model with three sigma values if
% typeFlag = 'hier_sigmaall'. Code uses mean fitted parameters for
% simulations. 

clear all;%close all;clc;
typeFlag = 'hier_sigmaall';
if strcmp(typeFlag, 'hier')
    load('hier_pt_model.mat')
elseif strcmp(typeFlag,'hier_sigmaall')
    load('hier_pt_model_sigmaall.mat')
else
    sprintf('using manually inputted values')
end

% 8 actions, separated by place in hiearchy
q = linspace(1,8,8)'; % 8 actions
hand = [1 1 1 1 2 2 2 2]'; % top level of hand
couplet = [1 1 2 2 1 1 2 2]'; % second level of couplet
finger = [1 2 1 2 1 2 1 2]'; % bottom level of specific finger
Nlevels = 3;
respLevels = [hand couplet finger];

PTs = 0:1:1500; % vector of potential prep times

%% first, we define preparation times for each level of the tree
if strcmp(typeFlag, 'hier')
    mu = mean(hier_pt_model.mu)*1000; % the MEAN time each level is prepared [hand, couplet, finger]
elseif strcmp(typeFlag,'hier_sigmaall')
    mu = mean(hier_pt_model_sigmaall.mu)*1000; % the MEAN time each level is prepared [hand, couplet, finger]
else
    mu = [615,1013,1264];
end

if strcmp(typeFlag, 'hier')
    sigma = [mean(hier_pt_model.sigma),mean(hier_pt_model.sigma),mean(hier_pt_model.sigma)]*1000;
elseif strcmp(typeFlag,'hier_sigmaall')
    sigma = mean(hier_pt_model_sigmaall.sigma)*1000;
else
    sigma = [69,255,253];
end

if strcmp(typeFlag, 'hier')
    rho = mean(hier_pt_model.rho);
elseif strcmp(typeFlag,'hier_sigmaall')
    rho = mean(hier_pt_model_sigmaall.rho);
else
    rho = [.384];
end

for i = 1:Nlevels
    Phi(i,:) = normcdf(PTs,mu(i),sigma(i)); % probability that H level has been planned by RT
    U(i,:) = ones(length(PTs),1)' * (1./length(q));
end

figure;
subplot(2,3,1);hold on;
plot(Phi'); % plot prep time functions
xlabel('RT');ylabel('p(correct selected)')
legend('hand','couplet','finger');

%% second, we weight each response type by the preparation functions
corA = 1; % let's say correct action is button 1
coraction = [1 1 1]; % define correct action in terms of each level value
for k = 1:Nlevels % 3 levels
    cor = respLevels(:,k) == coraction(k); % get set of correct level-specific actions
    tmp = (1-rho)*(Phi(k,:) .* cor) + rho*((1./length(q)) .* ones(8,length(PTs))); % weighted mixture of correct selection and lapsing
    prespPAR{k} = tmp./sum(tmp); % normalize mixture action selection prep to a probablity, put it in cell K

end

%% choice types over time
presponsePAR = (prespPAR{1}.*prespPAR{2}.*prespPAR{3}); % combine all preparation functions
presponsePAR = presponsePAR./sum(presponsePAR); % normalize to probs
subplot(2,3,2);hold on;
plot(presponsePAR');
xlabel('RT');ylabel('p(correct selected)')
legend;

%% choice types over time (Three lines)
presponsePAR = (prespPAR{1}.*prespPAR{2}.*prespPAR{3}); % combine all preparation functions
presponsePAR = presponsePAR./sum(presponsePAR); % normalize to probs
presponsePAR_type = nan(4,length(presponsePAR));
presponsePAR_type(1,:) = presponsePAR(1,:);
presponsePAR_type(2,:) = presponsePAR(2,:);
presponsePAR_type(3,:) = mean(presponsePAR(3:4,:));
presponsePAR_type(4,:) = mean(presponsePAR(5:8,:));

subplot(2,3,3);hold on;
plot(presponsePAR_type');
title('errorLevel')
xlabel('RT');ylabel('p(correct selected)')
legend;

% Save values for use in R script
if strcmp(typeFlag, 'hier')
    presponsePAR_t = array2table(presponsePAR','VariableNames',{'1/1/1','1/1/0','1/0/1','1/0/0','0/1/1','0/1/0','0/0/1','0/0/0'});
    presponsePAR_type_t = array2table(presponsePAR_type','VariableNames',{'corr','low','mid','top'});

    writetable(presponsePAR_t, 'presponseSims/simdata.csv');
    writetable(presponsePAR_type_t, 'presponseSims/errorLevel_simdata.csv');

elseif strcmp(typeFlag,'hier_sigmaall')
    presponsePAR_t = array2table(presponsePAR','VariableNames',{'1/1/1','1/1/0','1/0/1','1/0/0','0/1/1','0/1/0','0/0/1','0/0/0'});
    presponsePAR_type_t = array2table(presponsePAR_type','VariableNames',{'corr','low','mid','top'});

    writetable(presponsePAR_t, 'presponseSims/sigmaall_simdata.csv');
    writetable(presponsePAR_type_t, 'presponseSims/sigmaall_errorLevel_simdata.csv');
else
    sprintf('no csv')
end

% Show average parameters, conduct sign-rank test, calculate effect size
% based on zscore from signrank
if strcmp(typeFlag,'hier')
    mu
    sigma
    rho
    [pa,ha,statsa] = signrank(hier_pt_model.mu(:,1),hier_pt_model.mu(:,2))
    Za = statsa.zval; % Extract the Z-value
    Na = length(hier_pt_model.mu(:,1)); % Number of observations
    ra = Za / sqrt(Na); % Calculate the effect size
    [pb,hb,statsb] = signrank(hier_pt_model.mu(:,2),hier_pt_model.mu(:,3))
    Zb = statsb.zval; % Extract the Z-value
    Nb = length(hier_pt_model.mu(:,1)); % Number of observations
    rb = Zb / sqrt(Nb); % Calculate the effect size
elseif strcmp(typeFlag, 'hier_sigmaall')
    mu
    sigma
    rho
    [pa,ha,statsa] = signrank(hier_pt_model_sigmaall.mu(:,1),hier_pt_model_sigmaall.mu(:,2))
    Za = statsa.zval; % Extract the Z-value
    Na = length(hier_pt_model_sigmaall.mu(:,1)); % Number of observations
    ra = Za / sqrt(Na); % Calculate the effect size

    [pb,hb,statsb] = signrank(hier_pt_model_sigmaall.mu(:,2),hier_pt_model_sigmaall.mu(:,3))
    Zb = statsb.zval; % Extract the Z-value
    Nb = length(hier_pt_model_sigmaall.mu(:,1)); % Number of observations
    rb = Zb / sqrt(Nb); % Calculate the effect size
end


