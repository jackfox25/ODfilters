% =========================================================================
% Function: batch_OD
% Author:   Jack Fox
% Date:     02/05/25
%     General LKF used for OD. Implementing notation from Born et al.
% =========================================================================
% Inputs:
%          t_0   (1x1) Initial time (s)
%       Xhat_0   (nx1) Initial state guess
%          P_0   (nxn) Initial state covariance
%         meas   (kx(2+p)) True measurements with rows [t ID m1 m2 ... mp]
%       params   Struct containing constants used in sysFuncs
%     sysFuncs   Struct containing specific dynamic and measurement 
%                  functions. Any helper functions must be on path.
%       - sysFuncs.computeAugStateDot
%       - sysFuncs.computeMeasResidual
%         info   (bool) Flag for output text and plots
% Outputs:
%    batchOut   Struct containing outputs from batch
%      - Xhat_hist       (kxn) Full state history for each measurement time
%      - P_hist          (kxn^2) State covariance history for each measurement time
% =========================================================================
function batchOut = batch_OD(t_0,Xhat_0,P_0,meas,params,sysFuncs,info)
    
    tic;
    if info
        f = figure;
    end

    n = size(Xhat_0,1);
    k = size(meas,1);
    t_meas = meas(:,1);
    tspan = [t_0 t_meas'];

    xbar_0 = zeros(n,1);

    % Measurement noise covariance
    R = params.cov.R;

    % (A) Initialize
    t_im1 = t_0;
    Xref_im1 = Xhat_0;
    STM_0toim1 = eye(n);

    % Save STMs to use at the end
    STMhist = zeros(k,n^2);

    % ================== OUTER ITERATIVE LOOP ================== %
    converged = false;
    MAXITER = params.settings.batch_MAXITER;
    iter = 0;
    while ~converged && iter < MAXITER

        iter = iter + 1;

        % Normal equation accumulators
        Lambda = inv(P_0);
        N = P_0 \ xbar_0;
    
        Xref_aug_0 = [Xref_im1; reshape(STM_0toim1,n^2,1)];
        [~,Xref_hist] = params.settings.integrator(@(t,x) sysFuncs.computeAugStateDot(t,x,params),...
                                     tspan,Xref_aug_0,params.settings.OPTIONS);

        % Loop over measurements
        for i=1:k
            
            % (B) Read next observation
            t_i = meas(i,1);
            Y_i = meas(i,:);
            R_i = R;
    
            % If this is a distinct time:
            if t_i > t_im1
                % Integrate reference trajectory and STM from t_im1 to t_i
                Xref_i = Xref_hist(i+1,1:n)';
                STM_0toi = reshape(Xref_hist(i+1,n+1:end),n,n);
            else
                % At same time - do not perform time update
                Xref_i = Xref_im1;
                STM_0toi = STM_0toim1;
            end
    
            % Accumulate current observation
            [y_i,Htilde_i] = sysFuncs.computeMeasResidual(t_i,Xref_i',Y_i,params);
            H_i = Htilde_i * STM_0toi;
    
            % Add to normal equation accumulators
            Lambda = Lambda + H_i' * (R_i \ H_i);
            N      = N      + H_i' * (R_i \ y_i);
    
            % Set new im1 vars
            t_im1 = t_i;
            Xref_im1 = Xref_i;
            STM_0toim1 = STM_0toi;

            % Save STM
            STMhist(i,:) = reshape(STM_0toi,n^2,1)';
        end

        % Solve normal equations
        xhat_0 = Lambda \ N;
        Phat_0 = inv(Lambda);
        Xhat_0 = Xhat_0 + xhat_0;
        xbar_0 = xbar_0 - xhat_0;

        fprintf('Batch: Iter %d, dx = %.3e\n',iter,norm(xhat_0));
        if info
            Xhat_aug_0 = [Xhat_0; reshape(eye(n),n^2,1)];
            [~,Xhat_hist] = ode45(@(t,x) sysFuncs.computeAugStateDot(t,x,params),...
                                     tspan,Xhat_aug_0,params.settings.OPTIONS);
            [y,~] = sysFuncs.computeMeasResidual(tspan(2:end),Xhat_hist(2:end,:),meas,params);
            plotResid_batch(t_meas,y,params,iter);
        end

        % Check for convergence
        if norm(xhat_0) < params.settings.batch_dxtol
            converged = true;
            fprintf('Converged!\n');
        else
            % Not converged! Reset for next iteration.
            t_im1 = t_0;
            Xref_im1 = Xhat_0;
            STM_0toim1 = eye(n);
        end

    end

    % Propagate estimates forward nonlinearly
    Xhat_0_aug = [Xhat_0; reshape(eye(n),n^2,1)];
    sol = ode45(@(t,x) sysFuncs.computeAugStateDot(t,x,params),...
                                     [t_0 t_i],Xhat_0_aug,params.settings.OPTIONS);
    
    Xhat_aug_hist = deval(sol,meas(:,1))';
    Xhat_hist = Xhat_aug_hist(:,1:n);
    STM_0toi_hist = Xhat_aug_hist(:,n+1:end);

    P_hist = zeros(k,n^2);
    for i=1:k
        STM_0toi = reshape(STM_0toi_hist(i,:)',n,n);
        P_hist(i,:) = reshape((STM_0toi * Phat_0 * STM_0toi'),n^2,1)';
    end

    toc;

    % Package outputs
    batchOut.Xhat_hist = Xhat_hist;
    batchOut.P_hist = P_hist;

end
