%% play with forced RT model for Juliana's hierarchical task
%% SDM; New Haven, CT; 01/13/2023
clear all;close all;clc;

load('flatDyn_flatnostructure_pt_model.mat')
% 8 actions, separated by place in hiearchy
q = linspace(1,8,8)'; % 8 actions
respLevels = [1 0 0 0 0 0 0 0];

PTs = 0:1:1500; % vector of potential prep times

%% first, we define preparation times for each level of the tree
mu = mean(flatDyn_flatnostructure_pt_model.mu)*1000; % the MEAN time each level is prepared [hand, couplet, finger]
sigma = mean(flatDyn_flatnostructure_pt_model.sigma)*1000; % variance in time each level is prepared [hand, couplet, finger]
rho = mean(flatDyn_flatnostructure_pt_model.rho); % weight on uniform guessing preparation function
i = 1;
%for i = 1:Nlevels
Phi(i,:) = normcdf(PTs,mu,sigma); % probability that H level has been planned by RT
U(i,:) = ones(length(PTs),1)' * (1./length(q));
%end

figure;
subplot(2,2,1);hold on;
plot(Phi'); % plot prep time functions
xlabel('RT');ylabel('p(correct selected)')
legend('hand','couplet','finger');

%% second, we weight each response type by the preparation functions
corA = 1; % let's say correct action is button 1
coraction = 1; % define correct action in terms of each level value
cor = respLevels == coraction; % get set of correct level-specific actions
tmp = (1-rho)*(Phi(:) .* cor)' + rho*((1./length(q)) .* ones(8,length(PTs))); % weighted mixture of correct selection and lapsing
prespPAR = tmp./sum(tmp); % normalize mixture action selection prep to a probablity, put it in cell K

%% choice types over time
presponsePAR = prespPAR; % combine all preparation functions -- only 1 in this model
presponsePAR = presponsePAR./sum(presponsePAR); % normalize to probs
subplot(2,2,2);hold on;
plot(presponsePAR');
xlabel('RT');ylabel('p(correct selected)')
legend;

%% choice types over time (Three lines)
presponsePAR = prespPAR; % combine all preparation functions
presponsePAR = presponsePAR./sum(presponsePAR); % normalize to probs
presponsePAR_type = nan(4,length(presponsePAR));
presponsePAR_type(1,:) = presponsePAR(1,:);
presponsePAR_type(2,:) = presponsePAR(2,:);
presponsePAR_type(3,:) = mean(presponsePAR(3:4,:));
presponsePAR_type(4,:) = mean(presponsePAR(5:8,:));

subplot(2,2,3);hold on;
plot(presponsePAR_type');
title('errorLevel')
xlabel('RT');ylabel('p(correct selected)')
legend;


% set up simulation data for saving
presponsePAR_t = array2table(presponsePAR','VariableNames',{'1/1/1','1/1/0','1/0/1','1/0/0','0/1/1','0/1/0','0/0/1','0/0/0'});
presponsePAR_type_t = array2table(presponsePAR_type','VariableNames',{'corr','low','mid','top'});

writetable(presponsePAR_t, 'presponseSims/nostructure_simdata.csv');
writetable(presponsePAR_type_t, 'presponseSims/nostructure_type_simdata.csv');











