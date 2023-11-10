%% MCMC AM Algorithm Shell: READ ME --------------------------------------%
	%This script provides a basic shell of the MCMC Adaptive Metropolis
	%algorithm (Haario et al, 2001). There are several locations where
	%user-specific entries must be made (e.g., model data,
	%parameter initialization, priors, simulation model function,
	%likelihood function, etc.).

%% Iniitialization of workspace ------------------------------------------%
	clear all
	clc
	tic
	
%% Load model data -------------------------------------------------------%
    %%% set the model input and output 
        Nwarmup = 91; 
        Ncal = 1095;

        load('TSS_04087030.mat');         
        Output_obs_all = TSS(1:Ncal)';
        Output_obs_cal = TSS(Nwarmup:Ncal)';
        Input_all = obsQ(1:Ncal)';    

%% Initialize parameter values -------------------------------------------%
	%%% 1 - AM algorithm parameters
		ITER = 10000;	%number of iterations to perform
		i0   = 0.10;	%percentage of initial iterations before adaptation
		SD1  = 0.50;	%initial covariance matrix scaling factor (i0)
		SD2  = 0.15;	%adaptive covariance matrix scaling factor (1-i0)

	%%% 2 - initial model parameters
        a = 0.72;
        b = 1.41;
        kappa = 0.59;
        Smax = 1901; 

	%%% 3 - Likelihood function parameters
		varp = 0.41;	%variance parameter

	%%% 4 - Define parameter matrix
		theta(1,:) = [a,b,kappa,Smax,varp]; %first row using initial values
		nPAR = size(theta,2); %determines number of parameters

	%%% 5 - Define parameter names for output
		ParName = {'a','b','kappa','Smax','Variance'}; 

%% Define parameter priors -----------------------------------------------%
	%%% Define (log) prior to be used for each calibrated parameter.
        par_low(1)= 0; par_up(1)= 5;
        par_low(2)= 0; par_up(2)= 3;
        par_low(3)= 0; par_up(3)= 2;
        par_low(4)= 0; par_up(4)= 10000;
        par_low(5)= 0; par_up(5)= 10000;
        
        Pr1 = log(unifpdf(a,par_low(1),par_up(1)));
        Pr2 = log(unifpdf(b,par_low(2),par_up(2)));
        Pr3 = log(unifpdf(kappa,par_low(3),par_up(3)));
        Pr4 = log(unifpdf(Smax,par_low(4),par_up(4)));
        Pr5 = log(unifpdf(varp,par_low(5),par_up(5)));

        Pr(1) = Pr1+Pr2+Pr3+Pr4+Pr5;

        theta_prior(:,1)=unifrnd(par_low(1),par_up(1),ITER,1);
        theta_prior(:,2)=unifrnd(par_low(2),par_up(2),ITER,1);
        theta_prior(:,3)=unifrnd(par_low(3),par_up(3),ITER,1);
        theta_prior(:,4)=unifrnd(par_low(4),par_up(4),ITER,1);
        theta_prior(:,5)=unifrnd(par_low(5),par_up(5),ITER,1);

%% Initialize covariance matrix ------------------------------------------%
	%%% 1 - Define best guess variances for each calibrated parameter
		VarPar = [1E-6,1E-6,1E-6,1E-6,1E-6];

	%%% 2 - Populate initial covariance matrix
		for i = 1:nPAR
			CovPar(i,i) = SD1*VarPar(i);
        end

	%%% 3 - Define other covariance terms needed later
		epsilon = 1e-20; %small number preventing CovPar from becoming singular during adaptation
		Id = eye(nPAR);  %identity matrix used in CovPar adaptation

%% Initial run of simulation model ---------------------------------------%
	% Replace the following line with your function as appropriate, where
	% Output_sim_all is the output, BWmod_4 is the function name, 
    % a,b,kappa,Smax are model parameters
        [Output_sim_all] = BWmod_4(a,b,kappa,Smax,Input_all);
        Output_sim_cal(1,:) = Output_sim_all(Nwarmup:Ncal);
	
	% Replace the following line with the appropriate loglikelihood
	% function, where Output_obs_cal is the observed data, 
    % Output_sim_cal is the predicted data, varp is the variance parameter.
	% likelihoodMode is corresponding to the residual transformation and 
    % likelihoodPara is the parameters of transformation.
        likelihoodMode = 'log';
        likelihoodPara = 0.1;
        logL = loglikelihood(Output_obs_cal,Output_sim_cal(1,:),varp,...
            likelihoodMode,likelihoodPara);
        L(1) = logL; %initial value for the log likelihood
	
%% Run AM algorithm ------------------------------------------------------%
	progress = waitbar(0,'MCMC calibration progress...'); %progress dialog box
	Jumps = 0; %jump counter
	for i = 2:ITER
		%%% 1 - Adapt covariance matrix
			if i > i0*ITER;	%adaptive covariance matrix routine
				CovPar = SD2*cov(theta(1:i-1,:))+SD2*epsilon*Id; %updates the covariance
            end

		%%% 2 - Generate parameter proposals
			theta(i,:) = mvnrnd(theta(i-1,:),CovPar); %proposed values
            if  theta(i,1)< par_low(1) || theta(i,1) > par_up(1) ||...
                theta(i,2)< par_low(2) || theta(i,2) > par_up(2) ||...
                theta(i,3)< par_low(3) || theta(i,3) > par_up(3) ||...
                theta(i,4)< par_low(4) || theta(i,4) > par_up(4) ||...
                theta(i,5)< par_low(5) || theta(i,5) > par_up(5) 
    
                while theta(i,1)< par_low(1) || theta(i,1) > par_up(1) ||...
                theta(i,2)< par_low(2) || theta(i,2) > par_up(2) ||...
                theta(i,3)< par_low(3) || theta(i,3) > par_up(3) ||...
                theta(i,4)< par_low(4) || theta(i,4) > par_up(4) ||...
                theta(i,5)< par_low(5) || theta(i,5) > par_up(5) 
            
                     theta(i,:) = mvnrnd(theta(i-1,:),CovPar);
                end
            end

		%%% 3 - Run simulation model
            [Output_sim_all] = BWmod_4(theta(i,1),theta(i,2),theta(i,3),theta(i,4),Input_all);
            Output_sim_cal(i,:) = Output_sim_all(Nwarmup:Ncal);            
                           
        %%% 4 - Compute log likelihood
            logL = loglikelihood(Output_obs_cal,Output_sim_cal(i,:),theta(i,5),likelihoodMode,likelihoodPara);
            L(i) = logL;
                              
		%%% 5 - Compute log prior
            Pr1 = log(unifpdf(theta(i,1),par_low(1),par_up(1)));
            Pr2 = log(unifpdf(theta(i,2),par_low(2),par_up(2)));
            Pr3 = log(unifpdf(theta(i,3),par_low(3),par_up(3)));
            Pr4 = log(unifpdf(theta(i,4),par_low(4),par_up(4)));
            Pr5 = log(unifpdf(theta(i,5),par_low(5),par_up(5)));
            Pr(i) = Pr1+Pr2+Pr3+Pr4+Pr5;

		%%% 6 - Compute Metropolis ratio
			psi(i,1) = exp((L(i)-L(i-1))+(Pr(i)-Pr(i-1))); %Metropolis ratio

		%%% 7 - Determine accept/reject of proposal
			z = unifrnd(0,1);
			if z <= psi(i,:)
				theta(i,:) = theta(i,:); %jump to next theta
				Jumps = Jumps+1;
				Jump(i) = 1;
			else
				theta(i,:) = theta(i-1,:); %remain at current theta
                Output_sim_cal(i,:) = Output_sim_cal(i-1,:);
				L(i) = L(i-1);
				Pr(i) = Pr(i-1);
				Jump(i) = 0;
			end
		%%% 8 - Iterate
			waitbar(i/ITER);
	end

    %% Re-calculate model using best values ----------------------------------%
	    [optL,opti] = max(L);
        Parameter_opti = theta(opti,:)
        Parameter_post = theta(:,:);
        Output_sim_opti = Output_sim_cal(opti,:);
        Residual_error = Output_sim_opti - Output_obs_cal;
        save result_post Parameter_opti Parameter_post Output_sim_opti ... 
        Output_sim_cal Output_obs_cal Residual_error Jumps opti; 
    
    %% Create plots ----------------------------------------------------------%
    %%% Parameter Traces    
        figure('Name','Parameter Traces')
        for i = 1:nPAR
            subplot(3,3,i);plot(theta(:,i),'LineStyle','none','Marker','.');...
            xlabel('Iteration');ylabel(ParName(i));
        end
            subplot(3,3,nPAR+1);plot(L(:),'LineStyle','none','Marker','.');...
            xlabel('Iteration');ylabel('Likelihood (L)');
        saveas(gcf,'Parameter Traces.fig');
        saveas(gcf,'Parameter Traces.png');
        
        %%% Parameter Priors and Posteriors
        startshow=0.1*ITER;
        figure('Name','Parameter Marginal 1')
        for i = 1:nPAR
            subplot(2,3,i);
            [f1,xi1] = ksdensity(theta(startshow:ITER,i));
            [f2,xi2] = ksdensity(theta_prior(startshow:ITER,i));
            [AX,H1,H2]=plotyy(xi1,f1,xi2,f2);
            set(get(AX(1),'Ylabel'),'String','Posteriors Density');
            set(get(AX(2),'Ylabel'),'String','Priors Density');
            xlabel(ParName(i));
        end
        saveas(gcf,'Parameter Marginal.fig');
        saveas(gcf,'Parameter Marginal.png');          

        %%% Residual scatterplot of concentration
        figure('Name','Residual scatterplot of concentration')
            plot(Output_obs_cal,Residual_error,'k.');...
            xlabel('observation');ylabel('model residuals');
        saveas(gcf,'Residual scatterplot.fig');
        saveas(gcf,'Residual scatterplot.png');

        %%% Quantile-Quantile plot of concentration
       figure('Name','Q-Q plot of concentration');
            qqplot(Residual_error);
            ylabel('Quantiles of simulation residuals');
        saveas(gcf,'Quantile-Quantile plot.fig');
        saveas(gcf,'Quantile-Quantile plot.png');  

        %%% Reliability and sharpness of TSS concentration
        SAMPLE = 1000;
        if strcmp(likelihoodMode,'no')
            observed = Output_obs_cal;
            modelled = Output_sim_opti;
        elseif strcmp(likelihoodMode,'log')
            observed = log(Output_obs_cal + 0.1);
            modelled = log(Output_sim_opti + 0.1);
        end
        N = Ncal - Nwarmup + 1;
        variance = theta(opti,5);  
        for j = 1 : SAMPLE
            result_uncertainty(:,j) = modelled - normrnd(0,sqrt(variance),1,N);
            A = result_uncertainty(:,j);

            if strcmp(likelihoodMode,'no')
                A(A<0) = 0;
            end            
            result_uncertainty(:,j) = A;	
        end
       
        number_in = 0;        
        for l=1:N
            lo(l) = prctile(result_uncertainty(l,:),5);
            up(l) = prctile(result_uncertainty(l,:),95);  
            if observed(l)>=lo(l) && observed(l)<=up(l)
                number_in=number_in+1;
            end
            width_CI(l)=up(l)-lo(l);
        end
        reliability = number_in/N;
        sharpness = mean(width_CI);        
                
        %%% Result plot  
        figure('Name','Result plot of con')
            X=[(1:N) fliplr(1:N)];
            Y=[up fliplr(lo)];
            fill(X,Y,[0.5 0.5 0.5],'edgealpha',0);
            hold on
            plot(modelled,'r');
            plot(observed,'b');
            legend('5%-95% simulation','simulation','observation','Location','best');
            str = {['Reliability=', num2str(reliability)];[ 'Sharpness=',num2str(sharpness)]};
            text(600,max(observed)*0.8,str);
            hold off
        saveas(gcf,'Result plot.fig');
        saveas(gcf,'Result plot.png');

        save result_uncertainty  result_uncertainty modelled observed ...
            lo up reliability sharpness;

    toc 