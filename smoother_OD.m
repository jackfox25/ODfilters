% =========================================================================
% Function: smoother_OD
% Author:   Jack Fox
% Date:     03/05/25
%     Smoother used for post-processing sequential LKF results. 
%     Implementing notation from Born et al.
% =========================================================================
% Inputs:
%     lkfOut            Struct containing outputs from LKF
%       - Xhat_hist       (kxn) Full state history for each measurement time
%       - P_hist          (kxn^2) State covariance history for each measurement time
%       - xhat_hist       (kxn) State deviation history for each measurement time
%       - Pbar_hist       (kxn^2) State covariance history prior to measurements
%       - STM_im1toi_hist (kxn^2) Incremental state transition matrices (ref)
%       - Q_hist          (kxn^2) Noise compensation history
%     params            Struct containing problem-specific parameters
%  verbosity            (bool) Flag for output text
% Outputs:
%     smootherOut         Struct containing smoothed data
%       - Xhat_hist       (kxn) Full state history (smoothed)
%       - P_hist          (kxn^2) State covariance history (smoothed)
%       - xhat_hist       (kxn) State deviation history (smoothed)
% =========================================================================
function smootherOut = smoother_OD(lkfOut,params,verbosity)

    % Extract data from filter output
    Xref_hist = lkfOut.Xhat_hist - lkfOut.xhat_hist;
    P_hist = lkfOut.P_hist;
    xhat_hist = lkfOut.xhat_hist;
    Pbar_hist = lkfOut.Pbar_hist;
    STM_im1toi_hist = lkfOut.STM_im1toi_hist;
    Q_hist = lkfOut.Q_hist;

    % Number of states/measurements
    n = size(Xref_hist,2);
    k = size(Xref_hist,1);

    % Create smoother containers
    xhat_s_hist = zeros(size(xhat_hist));
    P_s_hist = zeros(size(P_hist));

    % Initialize
    xhat_s_hist(end,:) = xhat_hist(end,:);
    P_s_hist(end,:) = P_hist(end,:);

    % Loop backwards over data to produce better estimates using all info
    for i=k-1:-1:1
        
        xhat_s_ip1 = xhat_s_hist(i+1,:)';
        P_s_ip1 = reshape(P_s_hist(i+1,:)',n,n);

        STM_itoip1 = reshape(STM_im1toi_hist(i+1,:)',n,n);

        % Process noise:
        if strcmp(params.settings.noiseComp,"SNC")
            
            P_i = reshape(P_hist(i,:)',n,n); 
            xhat_i = xhat_hist(i,:)';
            Pbar_ip1 = reshape(Pbar_hist(i+1,:)',n,n);
            Q_itoip1 = reshape(Q_hist(i+1,:)',n,n);
            
            S_i = P_i * STM_itoip1' / (STM_itoip1 * P_i * STM_itoip1' + Q_itoip1);
            xhat_s_i = xhat_i + S_i * (xhat_s_ip1 - STM_itoip1 * xhat_i);
            P_s_i = P_i + S_i * (P_s_ip1 - Pbar_ip1) * S_i';

        % No process noise:
        else

            S_i = inv(STM_itoip1);
            xhat_s_i = S_i * xhat_s_ip1;
            P_s_i = S_i * P_s_ip1 * S_i';

        end

        % Save to containers
        xhat_s_hist(i,:) = xhat_s_i';
        P_s_hist(i,:) = reshape(P_s_i,n^2,1)';

        % Status
        if verbosity && (mod(i,10)==0 || i==1)
            fprintf('Smoother: Smoothed measurement %d\n',i);
        end

    end

    smootherOut.Xhat_hist = Xref_hist + xhat_s_hist;
    smootherOut.P_hist = P_s_hist;
    smootherOut.xhat_hist = xhat_s_hist;

end
