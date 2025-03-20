% =========================================================================
% Function: ukf_OD
% Author:   Jack Fox
% Date:     03/15/25
%     General UKF used for OD. Implementing notation from Born et al.
% =========================================================================
% Inputs:
%          t_0   (1x1) Initial time (s)
%       Xhat_0   (nx1) Initial state guess
%          P_0   (nxn) Initial state covariance
%         meas   (kx(2+p)) True measurements with rows [t ID m1 m2 ... mp]
%       params   Struct containing problem-specific parameters
%     sysFuncs   Struct containing specific dynamic and measurement 
%                  functions. Any helper functions must be on path.
%       - sysFuncs.computeStateDot (NL)
%       - sysFuncs.computeMeasurement (NL)
%    verbosity   (bool) Flag for output text
% Outputs:
%     ukfOut   Struct containing outputs from SRIF
%      - Xhat_hist       (kxn) Full state history for each measurement time
%      - P_hist          (kxn^2) State covariance history for each measurement time
% =========================================================================
function ukfOut = ukf_OD(t_0,Xhat_0,P_0,meas,params,sysFuncs,verbosity)

    n = size(Xhat_0,1);
    k = size(meas,1);
    p = size(meas,2)-2;

    % Extract UKF parameters
    alpha = params.settings.UKFalpha;
    beta = params.settings.UKFbeta;
    kappa = params.settings.UKFkappa;
    lambda = alpha^2*(n+kappa)-n;

    % Set weights
    [wts_m, wts_c] = deal(zeros(2*n+1,1));
    
    wts_m(1,1) = lambda/(n+lambda);
    wts_m(2:end,1) = 1/(2*(n+lambda))*ones(2*n,1);
    
    wts_c(1,1) = lambda/(n+lambda)+1-alpha^2+beta;
    wts_c(2:end,1) = 1/(2*(n+lambda))*ones(2*n,1);

    dwts_c = diag(wts_c);

    % Set im1 vars
    t_im1 = t_0;
    Xhat_im1 = Xhat_0;
    P_im1 = P_0;

    % Sigma point arrays
    Chi_list_im1 = zeros(n,2*n+1);
    Chi_list_i   = zeros(n,2*n+1);

    % Create containers
    Xhat_hist = zeros(k,n);
    P_hist = zeros(k,n^2);

    % Loop over measurements
    for i=1:k
        % Read measurement
        t_i = meas(i,1);
        Y_i = meas(i,:);
        ID_i = meas(i,2);

        % =============== DYNAMICS PREDICTION STEP =============== %
        % Generate sigma points at t_im1
        S_im1 = chol(P_im1,"upper");
        Chi_list_im1(:,1) = Xhat_im1;
        for j=1:n
            Chi_list_im1(:,1+j)   = Xhat_im1 + sqrt(n+lambda) * S_im1(j,:)';
            Chi_list_im1(:,1+n+j) = Xhat_im1 - sqrt(n+lambda) * S_im1(j,:)';
        end

        % Propagate each sigma point to t_i
        if t_i > t_im1
            [~, xh] = ode45(@(t,x) sysFuncs.computeStateDot(t,x,n,2*n+1,params),...
                            [t_im1 t_i], reshape(Chi_list_im1,n*(2*n+1),1), params.settings.OPTIONS);
            Chi_list_i = reshape(xh(end,:)',n,2*n+1);
        else
            Chi_list_i = Chi_list_im1;
        end

        % Add process noise
        if strcmp(params.settings.noiseComp,"SNC")
            SNC = params.SNC.Gam(t_i-t_im1) * getQ_SNC(Xhat_im1,params) * params.SNC.Gam(t_i-t_im1)';
            % Pad with zeros for extra non-accel states
            if size(SNC,1) < n
                SNC = blkdiag(SNC,zeros(n-size(SNC,1)));
            end
        else
            SNC = zeros(n);
        end

        % Compute time update
        if t_i > t_im1
            Xbar_i = sum(wts_m'.*Chi_list_i,2);
            Pbar_i = (Chi_list_i - Xbar_i) * dwts_c * (Chi_list_i - Xbar_i)' + SNC;
        else
            Xbar_i = Xhat_im1;
            Pbar_i = P_im1;
        end

        % =============== MEASUREMENT UPDATE STEP =============== %
        % Generate new sigma points at t_i based on Xbar_i and Pbar_i
        S_i = chol(Pbar_i,"upper");
        Chi_list_i(:,1) = Xbar_i;
        for j=1:n
            Chi_list_i(:,1+j)   = Xbar_i + sqrt(n+lambda) * S_i(j,:)';
            Chi_list_i(:,1+n+j) = Xbar_i - sqrt(n+lambda) * S_i(j,:)';
        end

        % Compute measurement from each sigma point
        gamma_list_i = zeros(p,2*n+1);
        for j=1:2*n+1
            m = sysFuncs.computeMeasurement(t_i,Chi_list_i(:,j),ID_i,params);
            gamma_list_i(:,j) = m(3:end)';
        end
        
        % Measurement mean and innovation covariance
        Ybar_i = sum(wts_m'.*gamma_list_i,2);
        Sbar_i = (gamma_list_i - Ybar_i) * dwts_c * (gamma_list_i - Ybar_i)' + params.cov.R;

        % Compute cross-covariance matrix
        Cxy = (Chi_list_i - Xbar_i) * dwts_c * (gamma_list_i - Ybar_i)';

        % Compute Kalman gain
        K_i = Cxy / Sbar_i;

        % Measurement update
        Xhat_i = Xbar_i + K_i * (Y_i(3:end)' - Ybar_i);
        P_i = Pbar_i - K_i * Sbar_i * K_i';

        % =============== SAVE, SET NEXT ITERATION =============== %
        Xhat_hist(i,:) = Xhat_i';
        P_hist(i,:) = reshape(P_i,n^2,1)';

        % Update im1 vars
        t_im1 = t_i;
        Xhat_im1 = Xhat_i;
        P_im1 = P_i;

        % Display status
        if verbosity && (mod(i,10)==0 || i==k)
            fprintf('UKF: Processed measurement %d of %d\n',i,k);
        end

    end

    % Package into output struct
    ukfOut.Xhat_hist = Xhat_hist;
    ukfOut.P_hist = P_hist;

end