clear;
%% Project 1
N_RF=6; % Sec. IV assumes N_RF=Nt_RF=Nr_RF preserves generality
N=64; % Number of antennas of the base station
Nt_RF=N_RF;   % Number of transmit RF chains
M=16; % Number of antennas of each user
Nr_RF=N_RF;   % Number of receive RF chains
K=1; % Number of users
d=6; % Number data streams required by each user
Ns=K*d; % Number of total data streams
L=15; % Number of paths (paper suggests 15)
SNR=-10:2:6;
tolerance=1e-6;
num_mc=100; %   Number of Monte Carlo trials
phase_list=[1 exp(1j*pi)];


% Tx precoder and combiner
V_D=zeros(Nt_RF,Ns);    % (14)    
V_RF=zeros(N,Nt_RF);    % (11)
Vt=V_RF*V_D;
% Rx precoder and combiner
%W_RF=zeros(Ns,Nt_RF);
%W_D=zeros(Nt_RF,N);
%Wt=W_RF*W_D;

% Complex channel
H=zeros(N,M);   % TODO N*M or M*N?

N_RF=Ns;
% NOTE Nt_RF=Nr_RF=N_RF, the total number of RF chains is 2*N_RF which
% satisfies "if the number of RF chains is twice the total number of data streams, the
% hybrid beamforming structure can realize any fully digital beamformer exactly"

R_list=zeros(1,length(SNR));

for i_snr=1:length(SNR)
    
    SNR_lin=10^(SNR(i_snr)/10);
    sigma=sqrt(1/SNR_lin);
    R=0;

    for i_mc=1:num_mc
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
        not_converge_int=1;
        not_converge_1bit=1;
        while not_converge 
            for j=1:N_RF
                V_RF_noj=V_RF;
                V_RF_noj(:,j)=[];
                C_j=eye(Nr_RF-1,Nr_RF-1)+(SNR_lin/N/N_RF)*V_RF_noj'*F1*V_RF_noj;
                G_j=(SNR_lin/N/N_RF)*F1-(SNR_lin/N/N_RF)^2*F1*V_RF_noj*inv(C_j)*V_RF_noj'*F1;
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

                    % 1-bit resolution phase shifter
                    min_abs_sq=10000;
                    V_RF_quant_min=V_RF(i,j);
                    for i_phase=1:length(i_phase)
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
            if abs(abs(V_RF).^2-ones(N,Nr_RF)) < tolerance*ones(N,Nr_RF) 
                not_converge_inf=0;
            end
            if abs(abs(V_RF_quant).^2-ones(N,Nr_RF)) < tolerance*ones(N,Nr_RF) 
                not_converge_1bit=0;
            end
        end
        
    
    
        %% Part 2: IV.A. Digital Precoder Design for N_RF=Ns
        Q=V_RF'*V_RF;
        Heff=H*V_RF;
        [Left,Gamma_e,U_e]=svd(Heff*inv(sqrtm(Q)));
        %https://scicoding.com/water-filling-algorithm-in-depth-explanation/#:~:text=Water%2Dfilling%20is%20a%20generic,in%20a%20technical%20sense%2C%20orthogonal
        Gamma_e(7:16,:)=[];
        V_D=inv(sqrtm(Q))*U_e*Gamma_e;
    
    
    
        %% Part 3: Hybrid Combining Design for N_RF=NS
        
        % Necessary parameters
        W_RF=ones(M,N_RF);
        Vt=V_RF*V_D;
        F2=H*(Vt)*Vt'*H';
    
        % Translating algorithm 1 to design W_RF
        not_converge=1;
        while not_converge
            for j=1:N_RF
                W_RF_noj=W_RF;
                W_RF_noj(:,j)=[];
                big=W_RF_noj'*F2*W_RF_noj;
                C_j=eye(size(big))+(1/M/sigma^2)*big;
                G_j=(1/M/sigma^2)*F2-(1/M/sigma^2)^2*F2*W_RF_noj*inv(C_j)*W_RF_noj'*F2;
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
                not_converge=0;
            end
        end
        
    
    
        % Design W_D
        J=W_RF'*H*Vt*Vt'*H'*W_RF+sigma^2*W_RF'*W_RF;
        W_D=inv(J)*W_RF'*H*Vt;
        
        Wt=W_RF*W_D;
        % Spectral efficiency
        R=R+log2(det(eye(M,M)+(Wt*inv(Wt'*Wt)*Wt'*H*Vt*Vt'*H')/sigma^2));
        
    
        
        %% Part 4: Hybrid Beamforming Design for Ns<N_RF<2*Ns



    end

    R_list(i_snr)=R/num_mc;

end


figure;
plot(SNR,R_list);
xlabel("SNR(dB");
ylabel("Spectral efficiency (unit=?)");
title("SNR vs. Spectral Efficiency (64*16 MIMO, single user, Ns=N_{RF}=6)");



