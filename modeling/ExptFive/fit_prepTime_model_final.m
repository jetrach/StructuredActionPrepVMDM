%% Fitting script for response time model for "Structured action preparation in VM DM"
% Script reads in compiled forced response data and fits model based on
% what modelType is input. hier = hierarchical model (3 mu, 1 sigma),
% sigmaall = hierarchical model (3 mu, 3 sigma), feature = feature based,
% flat = flat model

clear all;close all;clc

DT = readtable('Experiment5_forced.csv');

subs = unique(DT.subject);
nsubs = length(subs);
n_iter = 50; % fitting iterations

modelType = 'hier';

if strcmp(modelType, 'flat')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% model 1: truly flat model (no internal structure during fitting) %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this model doesnt consider the shared features in the task
% subject loop to fit model to each person

for si = 1:nsubs

    %% subset data
    subidx = DT.subject==subs(si); % subset data for sub "si"
    subdata = DT(subidx,:);
    sub_responses = subdata.corr; % should be the response that the person made (not N x 3 as previously)
    sub_PTs = subdata.PT_actual;

    disp(['now fitting subject',' ',num2str(si)]);

    for iter = 1:n_iter % iterations

        % model parameter starting values (random for each fitting iteration)
        minmu = .2; maxmu = 2; minsig = .01; maxsig = .6;
        mu = [unifrnd(minmu,maxmu)]; % mean prep time for each level, ONLY ONE MU FOR THIS MODEL
        sigma = unifrnd(minsig,maxsig); % variance of prep time for all levels
        rho = unifrnd(.01,.99); % lapse weight

        params = [mu sigma rho];
        options=optimset('display','off');
        LB = [minmu minsig 0.01]; % lower bounds, ONLY ONE MU FOR THIS MODEL
        UB = [maxmu maxsig 0.99]; % upper bounds, ONLY ONE MU FOR THIS MODEL
        % optimize!
        [params, ll] = fmincon(@func_flatnostructure_dynamics,params,[],[],[],[],LB,UB,[],options,sub_responses,sub_PTs);

        % best fit params this iteration
        model_p.ll(iter) = ll;
        model_p.p(iter,:) = params;
    end

    % subject's best fit params over all iterations
    [fit.ll(si),best] = min(model_p.ll);
    fit.mu(si,:) = model_p.p(best,1); % ONLY ONE MU FOR THIS MODEL
    fit.sigma(si,:) = model_p.p(best,2);
    fit.rho(si) = model_p.p(best,3);
    [fit.AIC(si), fit.BIC(si)] = aicbic(-fit.ll(si),length(params),length(sub_PTs));
end

flatnostructure_pt_model = fit; % name model
save flatnostructure_pt_model flatnostructure_pt_model; % save
elseif strcmp(modelType, 'hier')
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% model 2: hierarchical (three mu, one sigma) %%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% subject loop to fit model to each person
for si = 1:nsubs

    %% subset data
    subidx = DT.subject==subs(si); % subset data for sub "si"
    subdata = DT(subidx,:);
    sub_responses = [subdata.corrHand subdata.corrCouplet subdata.corrFinger]; % should be N X 3 matrix
    sub_PTs = subdata.PT_actual;

    disp(['now fitting subject',' ',num2str(si)]);

    for iter = 1:n_iter % iterations

        % model parameter starting values (random for each fitting iteration)
        minmu = .2; maxmu = 2; minsig = .01; maxsig = .6;
        mu = [unifrnd(minmu,maxmu) unifrnd(minmu,maxmu) unifrnd(minmu,maxmu)]; % mean prep time for each level
        sigma = unifrnd(minsig,maxsig); % variance of prep time for all levels
        rho = unifrnd(.01,.99); % lapse weight

        params = [mu sigma rho];
        options=optimset('display','off');
        LB = [minmu minmu minmu minsig 0.01]; % lower bounds
        UB = [maxmu maxmu maxmu maxsig 0.99]; % upper bounds
        % optimize!
        [params, ll] = fmincon(@func_hier_dynamics,params,[],[],[],[],LB,UB,[],options,sub_responses,sub_PTs);

        % best fit params this iteration
        model.ll(iter) = ll;
        model.p(iter,:) = params;
    end

    % subject's best fit params over all iterations
    [fit.ll(si),best] = min(model.ll);
    fit.mu(si,:) = model.p(best,1:3);
    fit.sigma(si,:) = model.p(best,4);
    fit.rho(si) = model.p(best,5);
    [fit.AIC(si), fit.BIC(si)] = aicbic(-fit.ll(si),length(params),length(sub_PTs));
end

hier_pt_model = fit; % name model
save hier_pt_model hier_pt_model; % save
elseif strcmp(modelType, 'sigmaall')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% model 3: hierarchical (three mu, three sigma) %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % subject loop to fit model to each person
for si = 1:nsubs

    %% subset data
    subidx = DT.subject==subs(si); % subset data for sub "si"
    subdata = DT(subidx,:);
    sub_responses = [subdata.corrHand subdata.corrCouplet subdata.corrFinger]; % should be N X 3 matrix
    sub_PTs = subdata.PT_actual;

    disp(['now fitting subject',' ',num2str(si)]);

    for iter = 1:n_iter % iterations

        % model parameter starting values (random for each fitting iteration)
        minmu = .2; maxmu = 2; minsig = .01; maxsig = .6;
        mu = [unifrnd(minmu,maxmu) unifrnd(minmu,maxmu) unifrnd(minmu,maxmu)]; % mean prep time for each level
        sigma = [unifrnd(minsig,maxsig) unifrnd(minsig,maxsig) unifrnd(minsig,maxsig)]; % variance of prep time for all levels
        rho = unifrnd(.01,.99); % lapse weight

        params = [mu sigma rho];
        options=optimset('display','off');
        LB = [minmu minmu minmu minsig minsig minsig 0.01]; % lower bounds
        UB = [maxmu maxmu maxmu maxsig maxsig maxsig 0.99]; % upper bounds
        % optimize!
        [params, ll] = fmincon(@func_hier_dynamics_sigmaall,params,[],[],[],[],LB,UB,[],options,sub_responses,sub_PTs);

        % best fit params this iteration
        model.ll(iter) = ll;
        model.p(iter,:) = params;
    end

    % subject's best fit params over all iterations
    [fit.ll(si),best] = min(model.ll);
    fit.mu(si,:) = model.p(best,1:3);
    fit.sigma(si,:) = model.p(best,4:6);
    fit.rho(si) = model.p(best,7);
    [fit.AIC(si), fit.BIC(si)] = aicbic(-fit.ll(si),length(params),length(sub_PTs));
end

hier_pt_model_sigmaall = fit; % name model
save hier_pt_model_sigmaall hier_pt_model_sigmaall; % save

elseif strcmp(modelType, 'feature')

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% model 4: feature-based model (single mu) %%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% subject loop to fit model to each person
for si = 1:nsubs

    %% subset data
    subidx = DT.subject==subs(si); % subset data for sub "si"
    subdata = DT(subidx,:);
    sub_responses = [subdata.corrHand subdata.corrCouplet subdata.corrFinger]; % should be N X 3 matrix
    sub_PTs = subdata.PT_actual;

    disp(['now fitting subject',' ',num2str(si)]);

    for iter = 1:n_iter % iterations

        % model parameter starting values (random for each fitting iteration)
        minmu = .2; maxmu = 2; minsig = .01; maxsig = .6;
        mu = [unifrnd(minmu,maxmu)]; % mean prep time for each level, ONLY ONE MU FOR THIS MODEL
        sigma = unifrnd(minsig,maxsig); % variance of prep time for all levels
        rho = unifrnd(.01,.99); % lapse weight

        params = [mu sigma rho];
        options=optimset('display','off');
        LB = [minmu minsig 0.01]; % lower bounds, ONLY ONE MU FOR THIS MODEL
        UB = [maxmu maxsig 0.99]; % upper bounds, ONLY ONE MU FOR THIS MODEL
        % optimize!
        [params, ll] = fmincon(@func_feature_dynamics,params,[],[],[],[],LB,UB,[],options,sub_responses,sub_PTs);

        % best fit params this iteration
        model.ll(iter) = ll;
        model.p(iter,:) = params;
    end

    % subject's best fit params over all iterations
    [fit.ll(si),best] = min(model.ll);
    fit.mu(si,:) = model.p(best,1); % ONLY ONE MU FOR THIS MODEL
    fit.sigma(si,:) = model.p(best,2);
    fit.rho(si) = model.p(best,3);
    [fit.AIC(si), fit.BIC(si)] = aicbic(-fit.ll(si),length(params),length(sub_PTs));
end

feature_pt_model = fit; % name model
save feature_pt_model feature_pt_model; % save
else
    disp('ALERT: model not fitting. No/wrong modelType input.')
end





