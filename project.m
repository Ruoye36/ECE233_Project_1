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

% Tx precoder and combiner
V_D=zeros(Nt_RF,Ns);    % (14)    
V_RF=zeros(N,Nt_RF);    % (11)
Vt=V_RF*V_D;
% Rx precoder and combiner
W_RF=zeros(Ns,Nt_RF);
W_D=zeros(Nt_RF,N);
Wt=W_RF*W_D;

% Complex channel
H=zeros(N,M);   % TODO N*M or M*N?

N_RF=Ns;

for i_snr=1:length(SNR)
    
    SNR_lin=10^(SNR(i_snr)/10);
    sigma=sqrt(1/SNR_lin);

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

    % Translating algorithm 1
    not_converge=1;
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
                if eta_ij==0
                    V_RF(i,j)=1;
                else
                    V_RF(i,j)=eta_ij/abs(eta_ij);
                end
            end
        end
        % Check convergence
        if abs(abs(V_RF).^2-ones(N,Nr_RF)) < tolerance*ones(N,Nr_RF) 
            not_converge=0;
        end
    end
    % Spectral efficiency
    %R=log2(abs(eye(M,M)+(Wt*inv(Wt'*Wt)*Wt'*H*Vt*Vt'*H')/sigma^2)); % TODO
    
    %% Part 2: IV.A. Digital Precoder Design for N_RF=Ns
    Q=V_RF'*V_RF;
    Heff=H*V_RF;
    [Left,Gamma_e,U_e]=svd(Heff*inv(sqrtm(Q)));
    %https://scicoding.com/water-filling-algorithm-in-depth-explanation/#:~:text=Water%2Dfilling%20is%20a%20generic,in%20a%20technical%20sense%2C%20orthogonal
    Gamma_e(7:16,:)=[];
    V_D=inv(sqrtm(Q))*U_e*Gamma_e;
    




    %% Part 3: Hybrid Combining Design for N_RF=NS
    
    %% Part 4: Hybrid Beamforming Design for Ns<N_RF<2*Ns

end





