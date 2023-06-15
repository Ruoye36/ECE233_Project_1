clear;
%% Project 1
N=10; % Number of antennas of the base station | 64 | 10
M=10; % Number of antennas of each user | 16 | 10
K=1; % Number of users
d=2; % Number data streams required by each user | 6 | 2
Ns=K*d; % Number of total data streams
L=15; % Number of paths (paper suggests 15)
SNR=0:5:30;   % -10:2:6 | 0:5:30
tolerance=1e-6;
num_mc=100; %   Number of Monte Carlo trials
phase_list=[1 exp(1j*pi)];
P=10;    % TODO change values to see if any effect

% Complex channel
H=zeros(N,M);   % TODO N*M or M*N?

N_RF=Ns;    % Sec. IV assumes N_RF=Nt_RF=Nr_RF preserves generality
Nt_RF=N_RF;   % Number of transmit RF chains
Nr_RF=N_RF;   % Number of receive RF chains
% NOTE Nt_RF=Nr_RF=N_RF, the total number of RF chains is 2*N_RF which
% satisfies "if the number of RF chains is twice the total number of data streams, the
% hybrid beamforming structure can realize any fully digital beamformer exactly"

% Tx precoder and combiner
V_D=zeros(Nt_RF,Ns);    % (14)    
V_RF=zeros(N,Nt_RF);    % (11)
Vt=V_RF*V_D;
% Rx precoder and combiner
%W_RF=zeros(Ns,Nt_RF);
%W_D=zeros(Nt_RF,N);
%Wt=W_RF*W_D;


R_list=zeros(1,length(SNR));
R_quant_list=zeros(1,length(SNR));

for i_snr=1:length(SNR)
    
    SNR_lin=10^(SNR(i_snr)/10);
    sigma=sqrt(P/SNR_lin);
    R=0;
    R_quant=0;

    for i_mc=1:num_mc
        %i_mc
        H=zeros(M,N);
        for l=1:L
            alpha=sqrt(1/2)*(randn(1,1)+1j*randn(1,1));
            phi_r=2*pi*randn(1,1);
            phi_t=2*pi*randn(1,1);
            a_r=transpose(exp(1j*pi*(0:M-1)*sin(phi_r))/sqrt(M));
            a_t=transpose(exp(1j*pi*(0:N-1)*sin(phi_t))/sqrt(N));
            H=H+alpha*a_r*a_t';
        end
        H=sqrt(N*M/L)*H;
        F1=H'*H;
        
        
        %% Part 1: IV.B RF Precoder Design for N_RF=Ns

        % Necessary parameters
        V_RF=ones(N,N_RF);
        V_RF_quant=ones(N,N_RF);

        % Translating algorithm 1 to design V_RF
        not_converge_inf=1;
        not_converge_1bit=1;

        % Convergence loop for infinite phase shifter
        while not_converge_inf
            for j=1:N_RF
                V_RF_noj=V_RF;
                V_RF_noj(:,j)=[];
                C_j=eye(Nr_RF-1,Nr_RF-1)+(SNR_lin/N/N_RF)*V_RF_noj'*F1*V_RF_noj;
                G_j=(SNR_lin/N/N_RF)*F1-(SNR_lin/N/N_RF)^2*F1*V_RF_noj*pinv(C_j)*V_RF_noj'*F1;
                for i=1:N
                    eta_ij=0;
                    for l=1:N
                        if l~=i
                            eta_ij=eta_ij+G_j(i,l)*V_RF(l,j);
                        end
                    end
                    % Infinite phase shifter
                    if eta_ij==0
                        V_RF(i,j)=1;
                    else
                        V_RF(i,j)=eta_ij/abs(eta_ij);
                    end
                end
            end
            % Check convergence
            if abs(abs(V_RF).^2-ones(N,Nr_RF)) < tolerance*ones(N,Nr_RF) 
                not_converge_inf=0;
            end
        end
        
        % Convergence loop for 1-bit phase shifter
        while not_converge_1bit 
            for j=1:N_RF
                V_RF_noj=V_RF;
                V_RF_noj(:,j)=[];
                C_j=eye(Nr_RF-1,Nr_RF-1)+(SNR_lin/N/N_RF)*V_RF_noj'*F1*V_RF_noj;
                G_j=(SNR_lin/N/N_RF)*F1-(SNR_lin/N/N_RF)^2*F1*V_RF_noj*pinv(C_j)*V_RF_noj'*F1;
                for i=1:N
                    eta_ij=0;
                    for l=1:N
                        if l~=i
                            eta_ij=eta_ij+G_j(i,l)*V_RF(l,j);
                        end
                    end
                    % 1-bit resolution phase shifter
                    min_abs_sq=10000;
                    V_RF_quant_min=V_RF(i,j);
                    for i_phase=1:length(phase_list)
                        candidate=abs(phase_list(i_phase)-eta_ij/abs(eta_ij))^2;
                        if candidate<min_abs_sq
                            min_abs_sq=candidate;
                            V_RF_quant_min=phase_list(i_phase);
                        end
                    end
                    V_RF_quant(i,j)=V_RF_quant_min;
                end
            end
            % Check convergence
            if abs(abs(V_RF_quant).^2-ones(N,Nr_RF)) < tolerance*ones(N,Nr_RF) 
                not_converge_1bit=0;
            end
        end
        
    
    
        %% Part 2: IV.A. Digital Precoder Design for N_RF=Ns
        Q=V_RF'*V_RF;
        Q_quant=V_RF_quant'*V_RF_quant;
        Heff=H*V_RF;
        Heff_quant=H*V_RF_quant;
        % NOTE water-filling solution: refer to MIMO.pdf formula 7.11 and above
        [Left,Delta,U_e]=svd(Heff*pinv(sqrtm(Q)));
        [Left_quant,Delta_quant,U_e_quant]=svd(Heff_quant*pinv(sqrtm(Q_quant)));
        % REF https://scicoding.com/water-filling-algorithm-in-depth-explanation/#:~:text=Water%2Dfilling%20is%20a%20generic,in%20a%20technical%20sense%2C%20orthogonal
        % REF https://zhuanlan.zhihu.com/p/502453127
        Delta(N_RF+1:M,:)=[];
        Delta_quant(N_RF+1:M,:)=[];
        N0=sigma^2;    % NOTE MIMO.pdf one line below formula 7.12
        lambda_list=diag(Delta);
        lambda_quant_list=diag(Delta_quant);
        lambda_sq_list=lambda_list.^2;
        lambda_quant_sq_list=lambda_quant_list.^2;
        Nmin=min(Nt_RF,Nr_RF);
        mu=(P+sum(sigma^2/lambda_sq_list))/Nmin;
        mu_quant=(P+sum(sigma^2/lambda_quant_sq_list))/Nmin;
        P_star_list=max(mu-N0/lambda_sq_list,0);
        P_star_quant_list=max(mu_quant-N0/lambda_quant_sq_list,0);
        Gamma_e=diag(P_star_list);
        Gamma_e_quant=diag(P_star_quant_list);
        V_D=pinv(sqrtm(Q))*U_e*Gamma_e;
        V_D_quant=pinv(sqrtm(Q_quant))*U_e_quant*Gamma_e_quant;
    
    
    
        %% Part 3: Hybrid Combining Design for N_RF=NS
        
        % Necessary parameters
        W_RF=ones(M,N_RF);
        W_RF_quant=ones(M,N_RF);
        Vt=V_RF*V_D;
        Vt_quant=V_RF_quant*V_D_quant;
        F2=H*(Vt)*Vt'*H';
        F2_quant=H*(Vt_quant)*Vt_quant'*H';
    
        % Translating algorithm 1 to design W_RF for infinite phase shifter
        not_converge_inf=1;
        while not_converge_inf
            for j=1:N_RF
                W_RF_noj=W_RF;
                W_RF_noj(:,j)=[];
                big=W_RF_noj'*F2*W_RF_noj;
                C_j=eye(size(big))+(1/M/sigma^2)*big;
                G_j=(1/M/sigma^2)*F2-(1/M/sigma^2)^2*F2*W_RF_noj*pinv(C_j)*W_RF_noj'*F2;
                for i=1:M
                    eta_ij=0;
                    for l=1:M
                        if l~=i
                            eta_ij=eta_ij+G_j(i,l)*W_RF(l,j);
                        end
                    end
                    if eta_ij==0
                        W_RF(i,j)=1;
                    else
                        W_RF(i,j)=eta_ij/abs(eta_ij);
                    end
                end
            end
            % Check convergence
            if abs(abs(W_RF).^2-ones(M,N_RF)) < tolerance*ones(M,N_RF) 
                not_converge_inf=0;
            end
        end

        %  Translating algorithm 1 to design W_RF for 1-bit resolution phase shifter
        not_converge_1bit=1;
        while not_converge_1bit
            for j=1:N_RF
                W_RF_noj_quant=W_RF_quant;
                W_RF_noj_quant(:,j)=[];
                big_quant=W_RF_noj_quant'*F2_quant*W_RF_noj_quant;
                C_j_quant=eye(size(big_quant))+(1/M/sigma^2)*big_quant;
                G_j_quant=(1/M/sigma^2)*F2_quant-(1/M/sigma^2)^2*F2_quant*W_RF_noj_quant*pinv(C_j_quant)*W_RF_noj_quant'*F2_quant;
                for i=1:M
                    eta_ij_quant=0;
                    for l=1:M
                        if l~=i
                            eta_ij_quant=eta_ij_quant+G_j_quant(i,l)*W_RF_quant(l,j);
                        end
                    end
                    if eta_ij_quant==0
                        W_RF_quant(i,j)=1;
                    else
                        W_RF_quant(i,j)=eta_ij_quant/abs(eta_ij_quant);
                    end
                end
            end
            % Check convergence
            if abs(abs(W_RF_quant).^2-ones(M,N_RF)) < tolerance*ones(M,N_RF) 
                not_converge_1bit=0;
            end
        end
    
        % Design W_D
        J=W_RF'*H*(Vt)*Vt'*H'*W_RF+sigma^2*(W_RF')*W_RF;
        J_quant=W_RF_quant'*H*(Vt_quant)*Vt_quant'*H'*W_RF_quant+sigma^2*(W_RF_quant')*W_RF_quant;
        W_D=pinv(J)*W_RF'*H*Vt;
        W_D_quant=pinv(J_quant)*W_RF_quant'*H*Vt_quant;
        Wt=W_RF*W_D;
        Wt_quant=W_RF_quant*W_D_quant;

        % Spectral efficiency
        R=R+log2(det(eye(M,M)+(Wt*pinv(Wt'*Wt)*Wt'*H*(Vt)*Vt'*H')/sigma^2));
        R_quant=R_quant+log2(det(eye(M,M)+(Wt_quant*pinv(Wt_quant'*Wt_quant)*Wt_quant'*H*(Vt_quant)*Vt_quant'*H')/sigma^2));
        
    
        
        %% Part 4: Hybrid Beamforming Design for Ns<N_RF<2*Ns



    end

    R_list(i_snr)=R/num_mc;
    R_quant_list(i_snr)=R_quant/num_mc;

end

% !!! Change the parameters the same as Fig. 2 for this figure
figure;
plot(SNR,R_list,'-o');
grid on;
xlabel("SNR(dB)");
ylabel("Spectral efficiency (bits/s/Hz)");
title("SNR vs. Spectral Efficiency (Infinite phase shifter, 10*10 MIMO, single user, Ns=N_{RF}=2)");


% !!! Change the parameters the same as Fig. 3 for this figure
figure;
plot(SNR,R_quant_list,'-o');
grid on;
xlabel("SNR(dB)");
ylabel("Spectral efficiency (bits/s/Hz)");
title("SNR vs. Spectral Efficiency (1-bit phase shifter, 10*10 MIMO, single user, Ns=N_{RF}=2)");



