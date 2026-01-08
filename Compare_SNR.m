K = 6; % UE Number
J = 6; % Eve Number
M = 6; % BS Antenna Number
Lk = 3; % Number of NLoS paths of legitimate users
Lj = 3; % Number of NLoS paths of eavesdroppers
SNR_dB = [-10:5:20]; % Signal-to-Noise Ratio
Monte_Carlo = 1e3; % Monte_Carlo Trials
A_Size = 3; % Aperture size of the transmit region normalized by the wavelength
BCD_Result = zeros(1,length(SNR_dB)); % The proposed BCD-based method
FPA_Result = zeros(1,length(SNR_dB)); % The fixed-posistion antenna (FPA)-based method
precision = 1e-4; % Convergence precision of the BCD method
precision_grid = 1000; % Resolution of the grid search
candidate = [(-A_Size/2):(A_Size/(precision_grid - 1)):(A_Size/2)]; % Candidate set of grid search
candidate = union(candidate, 1 / 2 * ([0:1:(M-1)]) - 1 / 2 * (M-1) / 2);
parfor Monte = [1:1:Monte_Carlo]  
    ele_Bob = pi*rand(Lk + 1,K);
    azi_Bob = pi*rand(Lk + 1,K);
    beta_Bob = (1/sqrt(2)*randn(K,Lk + 1) + sqrt(-1)*1/sqrt(2)*randn(K,Lk + 1)) .* (ones(K,1) * [1, 1/10 * ones(1,Lk)]);
    ele_Eve = pi*rand(Lj + 1,J);
    azi_Eve = pi*rand(Lj + 1,J);
    beta_Eve = (1/sqrt(2)*randn(J,Lj + 1) + sqrt(-1)*1/sqrt(2)*randn(J,Lj + 1)) .* (ones(J,1) * [1, 1/10 * ones(1,Lj)]);
    angle_user_Bob = cell(K,1);
    for bob_k = [1:1:K]
        angle_k = [sin(ele_Bob(:,bob_k)).*cos(azi_Bob(:,bob_k)),cos(ele_Bob(:,bob_k))];
        angle_user_Bob{bob_k} = angle_k;
    end
    angle_user_Eve = cell(J,1);
    for eve_j = [1:1:J]
        angle_j = [sin(ele_Eve(:,eve_j)).*cos(azi_Eve(:,eve_j)),cos(ele_Eve(:,eve_j))];
        angle_user_Eve{eve_j} = angle_j;
    end
    BCD_Result_Tmp = ones(1,length(SNR_dB));
    FPA_Result_Tmp = ones(1,length(SNR_dB));
    for SNR_Index = [1:1:length(SNR_dB)]
        [Monte, SNR_Index]
        %% The proposed BCD-based method
        snr = 10^(SNR_dB(SNR_Index)/10);
        alpha = ones(K,1); % Auxiliary Variable
        beta = ones(K,1); % Auxiliary Variable
        eta = ones(K,1); % Auxiliary Variable
        Tx = 1 / 2 * ([0:1:(M-1)]) - 1 / 2 * (M-1) / 2; % x-coordinates of the movable antennas
        Ty = zeros(1, M); % y-coordinates of the movable antennas
        T = [Tx;Ty]; % Locations of the movable antennas
        H = zeros(M,K); % Bobs' channel matrix
        G = zeros(M,J); % Eves' channel matrix
        g = snr * M * (Lj + 1) * sum(sum(abs(beta_Eve).^2)); 
        for bob_k = [1:1:K]
            angle_k = angle_user_Bob{bob_k};
            H(:,bob_k) = transpose(beta_Bob(bob_k,:) * exp(-1 * sqrt(-1) * 2 * pi * (angle_k * T)));
        end
        for eve_j = [1:1:J]
            angle_j = angle_user_Eve{eve_j};
            G(:,eve_j) = transpose(beta_Eve(eve_j,:) * exp(-1 * sqrt(-1) * 2 * pi * (angle_j * T)));
        end
        W = sqrt(snr) * H / sqrt(trace(H' * H)); % Initialization of the Digital Beamforming Matrix W using MRT
        b = ones(K, 1); % Initialization of the Rate Indicator b 
        Rs = 0;
        for bob_k = [1:1:K]
            SINR_k = abs(W(:,bob_k)'*H(:,bob_k))^2/(1 + sum(abs(W(:,setdiff([1:1:K],bob_k))'*H(:,bob_k)).^2));
            ESNR_k = sum(abs(W(:,bob_k)'*G).^2);
            Rk = log(1 + SINR_k) - log(1 + ESNR_k);
            Rs = Rs + (SINR_k > ESNR_k) * Rk;
            alpha(bob_k) = SINR_k; % Initialization of the Auxiliary Variable
            beta(bob_k) = (g - ESNR_k)/(1 + ESNR_k); % Initialization of the Auxiliary Variable
            eta(bob_k) = (H(:,bob_k)' * W(:,bob_k))/(1 + sum(abs(W'*H(:,bob_k)).^2)); % Initialization of the Auxiliary Variable
        end
        Sum_Rate = [0,Rs];
        while abs((Sum_Rate(end)-Sum_Rate(end-1))/Sum_Rate(end-1)) > precision
            %% Optimization of the Digital Beamforming Matrix W
            A = cell(K,1);
            a = cell(K,1);
            a1 = ones(M*K,1);
            v1 = ones(M*K,1);
            for bob_k = [1:1:K]
                a{bob_k} = b(bob_k) * (1 + alpha(bob_k)) * eta(bob_k) * H(:,bob_k);
                A{bob_k} = H * diag(b .* (1 + alpha) .* (abs(eta).^2)) * H' + b(bob_k) * (1 + beta(bob_k)) / (1 + g) * (G * G');
                [Uk, Vk, ~] = svd(A{bob_k});
                a1([1:1:M] + (bob_k-1) * M) = abs(Uk' * a{bob_k}).^2;
                v1([1:1:M] + (bob_k-1) * M) = diag(Vk);
            end
            if (abs(sum(a1./(v1.^2)))<snr)
                for bob_k = [1:1:K]
                    W(:,bob_k) = (A{bob_k})^(-1)*a{bob_k};
                end
            else
                a2 = 0; b2 = sqrt(sum(a1)/snr); % Bisection Search
                tol = 1e-13;
                m_bisection = 1 + round(round(log((b2-a2)/tol))/log(2));
                for k_bisection=1:m_bisection
                    p = (a2 + b2) / 2;    
                    if (real(sum(a1./((p+v1).^2)) - snr) * real(sum(a1./((b2+v1).^2)) - snr)) < 0
                            a2 = p;
                    else
                            b2 = p;
                    end
                    x = p;
                end
                lambda1 = x;
                for bob_k = [1:1:K]
                    W(:,bob_k) = (A{bob_k} + lambda1 * eye(M))^(-1) * a{bob_k};
                end
            end            
            %% Optimization of the Antenna Position Matrix T
            tmp_obj = 0;
            for bob_k = [1:1:K]
                tmp_obj = tmp_obj + b(bob_k) * ((1 + alpha(bob_k)) * (2 * real(eta(bob_k) * W(:,bob_k)' * H(:,bob_k)) - abs(eta(bob_k))^2 * sum(abs(W'*H(:,bob_k)).^2)) - (1 + beta(bob_k))/(1 + g) * sum(abs(W(:,bob_k)'*G).^2));
            end
            Element_Wise_Obj = [tmp_obj - 1, tmp_obj];
            while abs((Element_Wise_Obj(end)-Element_Wise_Obj(end-1))/Element_Wise_Obj(end-1)) > precision
                for m = [1:1:M]
                    %% Optimization of t_x
                    candidate_tmp = candidate;
                    for not_selected = setdiff([1:1:M], m)
                        distance = sqrt(sum(abs([candidate_tmp;T(2,m) * ones(1, length(candidate_tmp))] - T(:,not_selected) * ones(1, length(candidate_tmp))).^2));
                        candidate_tmp(find(distance < (1/2))) = [];
                    end
                    tmp_obj = 0;
                    for bob_k = [1:1:K]
                        tmp_obj1 = (1 + alpha(bob_k)) * (2 * real(eta(bob_k) * W(m,bob_k)' * (beta_Bob(bob_k,:) * exp(-1 * sqrt(-1) * 2 * pi * angle_user_Bob{bob_k} * [candidate_tmp;T(2,m) * ones(1, length(candidate_tmp))])))...
                        - abs(eta(bob_k))^2 * sum(abs((W(setdiff([1:1:M],m),:))' * H(setdiff([1:1:M],m),bob_k) * ones(1, length(candidate_tmp))...
                            + W(m,:)' * (beta_Bob(bob_k,:) * exp(-1 * sqrt(-1) * 2 * pi * angle_user_Bob{bob_k} * [candidate_tmp;T(2,m) * ones(1, length(candidate_tmp))]))).^2));
                        tmp_obj2 = 0; 
                        for eve_j = [1:1:J]
                            tmp_obj2 = tmp_obj2 + abs(W(setdiff([1:1:M],m),bob_k)' * G(setdiff([1:1:M],m),eve_j) + W(m,bob_k)' * (beta_Eve(eve_j,:) * exp(-1 * sqrt(-1) * 2 * pi * angle_user_Eve{eve_j} * [candidate_tmp;T(2,m) * ones(1, length(candidate_tmp))]))).^2;
                        end
                        tmp_obj2 = (1 + beta(bob_k))/(1 + g) * tmp_obj2;
                        tmp_obj = tmp_obj + b(bob_k) * (tmp_obj1 - tmp_obj2);
                    end
                    T(1,m) = candidate_tmp(find(tmp_obj == max(tmp_obj), 1));
                    %% Optimization of t_y
                    candidate_tmp = candidate;
                    for not_selected = setdiff([1:1:M], m)
                        distance = sqrt(sum(abs([T(1,m) * ones(1, length(candidate_tmp));candidate_tmp] - T(:,not_selected) * ones(1, length(candidate_tmp))).^2));
                        candidate_tmp(find(distance < (1/2))) = [];
                    end
                    tmp_obj = 0;
                    for bob_k = [1:1:K]
                        tmp_obj1 = (1 + alpha(bob_k)) * (2 * real(eta(bob_k) * W(m,bob_k)' * (beta_Bob(bob_k,:) * exp(-1 * sqrt(-1) * 2 * pi * angle_user_Bob{bob_k} * [T(1,m) * ones(1, length(candidate_tmp));candidate_tmp])))...
                        - abs(eta(bob_k))^2 * sum(abs((W(setdiff([1:1:M],m),:))' * H(setdiff([1:1:M],m),bob_k) * ones(1, length(candidate_tmp))...
                            + W(m,:)' * (beta_Bob(bob_k,:) * exp(-1 * sqrt(-1) * 2 * pi * angle_user_Bob{bob_k} * [T(1,m) * ones(1, length(candidate_tmp));candidate_tmp]))).^2));
                        tmp_obj2 = 0; 
                        for eve_j = [1:1:J]
                            tmp_obj2 = tmp_obj2 + abs(W(setdiff([1:1:M],m),bob_k)' * G(setdiff([1:1:M],m),eve_j) + W(m,bob_k)' * (beta_Eve(eve_j,:) * exp(-1 * sqrt(-1) * 2 * pi * angle_user_Eve{eve_j} * [T(1,m) * ones(1, length(candidate_tmp));candidate_tmp]))).^2;
                        end
                        tmp_obj2 = (1 + beta(bob_k))/(1 + g) * tmp_obj2;
                        tmp_obj = tmp_obj + b(bob_k) * (tmp_obj1 - tmp_obj2);
                    end
                    T(2,m) = candidate_tmp(find(tmp_obj == max(tmp_obj), 1));
                    for bob_k = [1:1:K]
                        angle_k = angle_user_Bob{bob_k};
                        H(:,bob_k) = transpose(beta_Bob(bob_k,:) * exp(-1 * sqrt(-1) * 2 * pi * (angle_k * T)));
                    end
                    for eve_j = [1:1:J]
                        angle_j = angle_user_Eve{eve_j};
                        G(:,eve_j) = transpose(beta_Eve(eve_j,:) * exp(-1 * sqrt(-1) * 2 * pi * (angle_j * T)));
                    end
                end
                tmp_obj = 0;
                for bob_k = [1:1:K]
                    tmp_obj = tmp_obj + b(bob_k) * ((1 + alpha(bob_k)) * (2 * real(eta(bob_k) * W(:,bob_k)' * H(:,bob_k)) - abs(eta(bob_k))^2 * sum(abs(W'*H(:,bob_k)).^2)) - (1 + beta(bob_k))/(1 + g) * sum(abs(W(:,bob_k)'*G).^2));
                end
                Element_Wise_Obj = [Element_Wise_Obj, tmp_obj];
            end            
            %% Optimization of the Rate Indicator b and the Auxiliary Variables
            Rs = 0;
            for bob_k = [1:1:K]
                SINR_k = abs(W(:,bob_k)'*H(:,bob_k))^2/(1 + sum(abs(W(:,setdiff([1:1:K],bob_k))'*H(:,bob_k)).^2));
                ESNR_k = sum(abs(W(:,bob_k)'*G).^2);
                k = bob_k;
                g_1_k = log(1 + alpha(bob_k)) - alpha(bob_k) + (1 + alpha(bob_k)) * (2 * real(eta(bob_k) * W(:,bob_k)' * H(:,bob_k)) - abs(eta(bob_k))^2 * (sum(abs(W'*H(:,bob_k)).^2) + 1));
                f_2_k = log(1 + beta(bob_k)) - beta(bob_k) + (1 + beta(bob_k)) / (1 + g) * (g - ESNR_k);
                b(bob_k) = ((g_1_k + f_2_k) > 0);
                alpha(bob_k) = SINR_k;
                beta(bob_k) = (g - ESNR_k)/(1 + ESNR_k);
                eta(bob_k) = (H(:,bob_k)' * W(:,bob_k))/(1 + sum(abs(W'*H(:,bob_k)).^2));
                Rk = log(1 + SINR_k) - log(1 + ESNR_k);
                Rs = Rs + (SINR_k > ESNR_k) * Rk;
            end
            Sum_Rate = [Sum_Rate,Rs];
        end
        BCD_Result_Tmp(SNR_Index) = Sum_Rate(end);
        %% The fixed-posistion antenna (FPA)-based method
        snr = 10^(SNR_dB(SNR_Index)/10);
        alpha = ones(K,1); % Auxiliary Variable
        beta = ones(K,1); % Auxiliary Variable
        eta = ones(K,1); % Auxiliary Variable
        Tx = 1 / 2 * ([0:1:(M-1)]) - 1 / 2 * (M-1) / 2; % x-coordinates of the movable antennas
        Ty = zeros(1, M); % y-coordinates of the movable antennas
        T = [Tx;Ty]; % Locations of the movable antennas
        H = zeros(M,K); % Bobs' channel matrix
        G = zeros(M,J); % Eves' channel matrix
        g = snr * M * (Lj + 1) * sum(sum(abs(beta_Eve).^2)); 
        for bob_k = [1:1:K]
            angle_k = angle_user_Bob{bob_k};
            H(:,bob_k) = transpose(beta_Bob(bob_k,:) * exp(-1 * sqrt(-1) * 2 * pi * (angle_k * T)));
        end
        for eve_j = [1:1:J]
            angle_j = angle_user_Eve{eve_j};
            G(:,eve_j) = transpose(beta_Eve(eve_j,:) * exp(-1 * sqrt(-1) * 2 * pi * (angle_j * T)));
        end
        W = sqrt(snr) * H / sqrt(trace(H' * H)); % Initialization of the Digital Beamforming Matrix W using MRT
        b = ones(K, 1); % Initialization of the Rate Indicator b 
        Rs = 0;
        for bob_k = [1:1:K]
            SINR_k = abs(W(:,bob_k)'*H(:,bob_k))^2/(1 + sum(abs(W(:,setdiff([1:1:K],bob_k))'*H(:,bob_k)).^2));
            ESNR_k = sum(abs(W(:,bob_k)'*G).^2);
            Rk = log(1 + SINR_k) - log(1 + ESNR_k);
            Rs = Rs + (SINR_k > ESNR_k) * Rk;
            alpha(bob_k) = SINR_k; % Initialization of the Auxiliary Variable
            beta(bob_k) = (g - ESNR_k)/(1 + ESNR_k); % Initialization of the Auxiliary Variable
            eta(bob_k) = (H(:,bob_k)' * W(:,bob_k))/(1 + sum(abs(W'*H(:,bob_k)).^2)); % Initialization of the Auxiliary Variable
        end
        Sum_Rate = [0,Rs];
        while abs((Sum_Rate(end)-Sum_Rate(end-1))/Sum_Rate(end-1)) > precision
            %% Optimization of the Digital Beamforming Matrix W
            A = cell(K,1);
            a = cell(K,1);
            a1 = ones(M*K,1);
            v1 = ones(M*K,1);
            for bob_k = [1:1:K]
                a{bob_k} = b(bob_k) * (1 + alpha(bob_k)) * eta(bob_k) * H(:,bob_k);
                A{bob_k} = H * diag(b .* (1 + alpha) .* (abs(eta).^2)) * H' + b(bob_k) * (1 + beta(bob_k)) / (1 + g) * (G * G');
                [Uk, Vk, ~] = svd(A{bob_k});
                a1([1:1:M] + (bob_k-1) * M) = abs(Uk' * a{bob_k}).^2;
                v1([1:1:M] + (bob_k-1) * M) = diag(Vk);
            end
            if (abs(sum(a1./(v1.^2)))<snr)
                for bob_k = [1:1:K]
                    W(:,bob_k) = (A{bob_k})^(-1)*a{bob_k};
                end
            else
                a2 = 0; b2 = sqrt(sum(a1)/snr); % Bisection Search
                tol = 1e-13;
                m_bisection = 1 + round(round(log((b2-a2)/tol))/log(2));
                for k_bisection=1:m_bisection
                    p = (a2 + b2) / 2;    
                    if (real(sum(a1./((p+v1).^2)) - snr) * real(sum(a1./((b2+v1).^2)) - snr)) < 0
                            a2 = p;
                    else
                            b2 = p;
                    end
                    x = p;
                end
                lambda1 = x;
                for bob_k = [1:1:K]
                    W(:,bob_k) = (A{bob_k} + lambda1 * eye(M))^(-1) * a{bob_k};
                end
            end            
            %% Optimization of the Rate Indicator b and the Auxiliary Variables
            Rs = 0;
            for bob_k = [1:1:K]
                SINR_k = abs(W(:,bob_k)'*H(:,bob_k))^2/(1 + sum(abs(W(:,setdiff([1:1:K],bob_k))'*H(:,bob_k)).^2));
                ESNR_k = sum(abs(W(:,bob_k)'*G).^2);
                k = bob_k;
                g_1_k = log(1 + alpha(bob_k)) - alpha(bob_k) + (1 + alpha(bob_k)) * (2 * real(eta(bob_k) * W(:,bob_k)' * H(:,bob_k)) - abs(eta(bob_k))^2 * (sum(abs(W'*H(:,bob_k)).^2) + 1));
                f_2_k = log(1 + beta(bob_k)) - beta(bob_k) + (1 + beta(bob_k)) / (1 + g) * (g - ESNR_k);
                b(bob_k) = ((g_1_k + f_2_k) > 0);
                alpha(bob_k) = SINR_k;
                beta(bob_k) = (g - ESNR_k)/(1 + ESNR_k);
                eta(bob_k) = (H(:,bob_k)' * W(:,bob_k))/(1 + sum(abs(W'*H(:,bob_k)).^2));
                Rk = log(1 + SINR_k) - log(1 + ESNR_k);
                Rs = Rs + (SINR_k > ESNR_k) * Rk;
            end
            Sum_Rate = [Sum_Rate,Rs];
        end
        FPA_Result_Tmp(SNR_Index) = Sum_Rate(end);
    end
    BCD_Result = BCD_Result + BCD_Result_Tmp;
    FPA_Result = FPA_Result + FPA_Result_Tmp;
end
plot(SNR_dB,BCD_Result/Monte_Carlo,'-o');hold on;
plot(SNR_dB,FPA_Result/Monte_Carlo,'-s');hold on;
