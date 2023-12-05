clc;
clear all;
%%%%%%%%%%%%************Model***************%%%%%%%%%%%  
% For         802.11n      802.11g     802.11b(DS)             
% Slots time  20e-6         20e-5            20e-6         
% SIFS        20e-6         10e-6           10e-6        
% DIFS        50 e-6        50e-6          50e-6          
% EIFS        92.6e-6      396e-6 ?         364e-6         
% PHY rate    6~64*1048576 1~2*1048576     1~2*1048576    
% Wmin        15           15              31             
% Wmin        15           15              31            

%Bit rates (data and phy)802.11g
%R_data = 54000000;          %(54 Mbps)
%R_phy  = 6000000;          %(6 Mbps)

%%%%% Input
%R_data = input('Data Rate: (in bps)\n');
%m = input('Maximum stage: \n');

R_data = 65000000;          %(12 Mbps)
R_phy  = 6000000;          %(6 Mbps from excel sheet)

% R_data_12 = 12000000;          %(12 Mbps)
% R_data_54 = 54000000;          %(54 Mbps)
%R_phy  = 6000000;          %(1 Mbps from excel sheet)

%Model Parameters (in MicroSec)
SIFS = 0.00001; %10e-06; 
slot_time = 0.00002; %20e-06; 
DIFS = 0.00005; %50e-06;
% delta= 0.000001; %1e-06;                %propagation delay(s)
delta= 0.0000020014; %1e-06;     %propagation delay(s) {d/c} from application scenario (largest distance to destination is 600.0m)
payload = 12000;            %bits (1500B)
mac_header = 272;           %bits (34B)
phy_header = 128;           %bits (16B)
ack_p = 112;                %bits (14B)
packet = payload + mac_header; %+ phy_header;
ACK = ack_p + phy_header;

%*********Data packet and ACK time
T_data = (packet/R_data)+ (phy_header/R_phy);
T_ack = (ACK/R_data)+((phy_header/R_phy));          %includes ack and header pckt

% T_data_12 = (packet/R_data_12)+ (phy_header/R_phy);
% T_ack_12 = (ACK/R_data_12)+((phy_header/R_phy)); 
% 
% T_data_54 = (packet/R_data_54)+ (phy_header/R_phy);
% T_ack_54 = (ACK/R_data_54)+((phy_header/R_phy)); 

%T_ack = ACK/R_phy;

%delay times (Basic access mechanism)
T_s = DIFS + T_data;               %transmitter active time
T_c = DIFS + T_data;       %collision during receiver active time

Trx=T_ack+2*delta+SIFS;

% T_s = DIFS + T_data + 2*delta + SIFS + T_ack;
% T_c = DIFS + T_data + delta;
T_ce = T_data + T_ack;
%*******************************************************
%Variable Number of Stations and Payload Size
N=0:8:200;
C=0.2*N;           %Assuming 10% of the nodes collide (i.e. 0.1x200)
%N=2:2:50;
N(1)=1; 
C(1)=1; 
%L_pck = 1500;                       %1500 in bytes

%parameters
W = 32;
m = 6;
Rtx=slot_time*(4*((2^0)*((W-1)/2)+(2^1)*((W-1)/2)+(2^2)*((W-1)/2)+(2^3)*((W-1)/2))); %Packets are sent in the 3rd tTx attempt
%Rtx=slot_time*((2^0)*(W/2)+(2^1)*((W-1)/2)+(2^2)*((W-1)/2)+(2^3)*((W-1)/2)+(2^4)*((W-1)/2));%
%Rtx=slot_time*((2^0)*(W/2)+(2^1)*((2*W)/2)+(2^2)*((2*W)/2)+(2^3)*((2*W)/2));  %+(2^4)*((2*W)/2));
E_s = delta + T_s+T_c + Rtx; %(crude approximation of the E[slot time])
lamda= 4;
mu = 1/E_s;
Q=lamda/mu;
%lamda= 0.25;%4;
%mu = 10;
%E_s = slot_time + T_s+T_c+T_ce; %+ P_tx*(1-P_suc)*T_c;
%Q=1./(lamda.*E_s);
%Q=0.25*E_s;
%Q=(lamda-E_s)./lamda;   %(Deterministic queue)
%Q=(1/lamda)./(1/E_s); 
%Q = (1-(1-(lamda./mu)));
%Q = (1-(1-exp(-(lamda.*E_s))));
%**************************************************************************
%Collision prob., Throughput and MAC delay calculated for different Num. of Users 

%*****Collision Probabilities*****#
P_fp=[];                    %Makama
% P_col_B=[];                 %Bianchi
% P_col_D=[];                 %Danesh
% P_col_F=[];                 %Felemban

%***Throughput********
S = [];                     %Done
% B_thr = [];                 %Bianchi throughput
% D_thr= [];                  %Danesh
% F_thr= [];                  %Felemban

 S_TrfMdl=[];               %Makama Simulation
 %r=[];
%******MAC Delay*********
Dmac = [];                  
% Dmac_B= [];                 %Bianchi
% Dmac_D= [];                 %Danesh
% Dmac_F= [];                 %Felemban

Dmac_sim = [];              %Makama Simulation
%******Packet drop time/Delay*********
%Dmac_drp = [];                  
% Dmac_Bdrp= [];                 %Bianchi
% Dmac_Ddrp= [];                 %Danesh
% Dmac_Fdrp= [];                 %Felemban

%*******Calculating prob to exceed a certain Delay*************
% P_suc_pvd=[];

%******************************************************
%iterate for every possible number of users
%for i=1:size(N, 2)
for k = 1:size(N, 2)        %Done
    if (N(k) == 1)
        p = 0.0;
        tau = 2/(W + 1);
      else
        inc = 0.000001;
        error_MAX = 0.000005;
        p= 0.0;
        tau1 = 0.0;
        tau2 = tau1 + 2*error_MAX;

        while abs(tau1 - tau2) > error_MAX
            p = p + inc;
%             tau1 = 2*Q*(1-2*p)*(1-(p).^(m+1))./(Q*((W*(1-p)*(1-(2*p).^(m+1))+ (1-2*p)*(1-(p).^(m+1))))+2*(1-Q)*(1-p)*(1-2*p));
            tau1 = 2*Q*P_cd*(1-2*p)*(1-(p).^(m+1))./(Q*((W*(1-p)*(1-(2*p).^(m+1))+ (1-2*p)*(2*P_cd-1)*(1-(p).^(m+1))))+2*(1-Q)*P_cd*(1-p)*(1-2*p));
            tau2 = 1 - (1-p).^(1./(N-1));
%             tau2 = 1 - (1-p)^(1/(N(1, k)-1));
            
        end
        tau = (tau1 + tau2)./2;          %avg value of both tau
        
    end
%*******************Probabilities***************

    %Probability of collision
    P_c = 1- (1-tau).^(N-1);
    %Probability of failure due to channel error
    P_ce = 0.05;          % (2%)
    %Probability of packet failure (P_f)
    P_f = P_c+P_ce-(P_c*P_ce);
    %Packet drop probability
    P_D = P_f.^(m+1);
    %Probability of Tx (probability that there is at least one transmission in a the considered time slot)
    P_tx = 1 - (1- tau).^N(1,k);
    %Probability of successful Tx (prob. that exactly one geophone Txmits)
    P_TxG = N(1,k).*tau.*(1- tau).^(N(1,k)-1);
    %Probability that a packet transmission ongoing on the channel is successful
    P_suc=P_TxG./P_tx;
    %P_suc=(N(1,k)*tau*(1- tau)^(N(1,k)-1))./(1 - (1- tau)^N(1,k));
    %Probability of empty slot time
    P_e = 1 - P_tx;
    %Probability of collision
    %P_c = P_tr - P_s;
    %Probability counter decrements
    P_cd = (1-tau).^(N(1, k)-1);             
    
%****************************************
%Average slot duration
%E_s = P_e*slot_time + P_suc*T_s + P_f*T_c;
%E_s = P_e.*slot_time + (P_suc.*P_tx).*T_s + (P_tx.*(1-P_suc)).*T_c;
%E_s = P_e*slot_time + P_suc*P_tx*(1-P_ce)*T_s + P_tx*(1-P_suc)*T_c+P_tx*P_suc*P_ce*T_ce; %include channel error prob.
%#####################################
% CW_avgS = 0.0;
%  for i=1:m
%      if j<m
%    %CW_avgS = mean(((1-P_f).*P_f./(1-P_D)).*((2^i)*(W-1)));
%    CW_avgS = CW_avgS+ ((1-P_f)./(1-P_D)).*((P_f.^i).*(((2^i)*(W-1))./2));
%      else
%          CW_avgS =CW_avgS+ ((1-P_f)./(1-P_D)).*((P_f.^m).*(((2^m)*(W-1))./2));
%        %CW_avgS = mean(((1-P_f).*P_f./(1-P_D)).*((2^m)*(W-1)));  
%      end
%  end
 for i=1:m
    for j= 1:i
     if j<m
   CW_avgS = (((1-P_f).*(P_f)./(1-P_D)).*((2^j)*(W-1)));
     else
       CW_avgS =(((1-P_f).*(P_f)./(1-P_D)).*((2^m)*(W-1)));  
     end
    end
 end
%###############################################
%CW_avgS = 167.5;        %calculated manually
%**********************Probabilities (pci,pcs and pcc)*************************

pci= (1-(1./CW_avgS)).^C;
pcs= C.*(1./CW_avgS).*(1-(1./CW_avgS)).^(C-1);
pcc= 1-pci-pcs;

%**********************Probabilities (pbi,pbs and pbc)*************************
pbi=(1-tau).^(N-1);
bs=nchoosek(199,1);             %199 = N-1          %**fix the problem wit N-1 here
pbs= bs.*tau.*(1-tau).^(N-2);
pbc=1-pbi-pbs;

%**********************Probabilities (pss and psi)*************************
pss=1/W;                   %1/CWmin
psi=1-pss;
%**********************Durations sigmaI, sigmaS and sigmaC*************************
%sigmaI= 1; %*slot_time;
sigmaI= slot_time;
sigmaS= (1/(1-pss))*T_s + (psi/(1-pss))*sigmaI;

%sumC=0.0;
for i=1:m
    for j=1:i
        if j<m
        %sumC = sumC + i.*pcc.^i;
        sumC = j.*(pcc.^j);
        else
%         sumC = sumC + m.*pcc.^m;
        sumC = m.*(pcc.^m);
        end
    end
end  
% for i=1:m
%     sumC = sumC + i.*pcc.^i;
% end  
sigmaC= sumC*T_c + (pcs./(1-pcc))*sigmaS + (pci./(1-pcc))*sigmaI;

%**Durations for entering backoff stages from previous bckoff and Tx stage deltaB and deltaT***
 %Probability that the counter decrements
Pcd=(1-tau).^(N-1);     %Pcd=pbi
del=sigmaI.*pbi+sigmaS.*pbs+sigmaC.*pbc;


deltaB=(1./Pcd).*(del);       

deltaT=Q*((1-(1./CW_avgS))).*del;
%********Durations a geophone spends in a backoff stage (delta)********
D= (1-tau).*deltaB + tau.*deltaT;

%******************************
% dmac_t2=0.0;
% for i=1:m
%     avg_CW =((2.^i).*(W-1)/2);
%     %avg_CW=167.5;
%  %dmac_t3=D.*(avg_CW*slot_time);
%      dmac_t3=D.*(avg_CW*slot_time);
%     
%     %dmac_t2=(1./(1-P_D)).*(dmac_t2 + (P_suc).*(T_s+ i.*T_c + dmac_t3));        %sum over i term of equation 40 in paper
%     dmac_t2=(1./(1-P_D)).*((P_suc).*(T_s+ i.*T_c + dmac_t3));
%  end

%Dmac = dmac_t2;              %Makama
%******************************
%#Working (Showing curve behavior)
dmac_t2=0.0;
for i=1:m
    %avg_CW=167.5;
    for j= 1:i
        if j<m
            avg_CW = ((2^j)*(W-1))/2;
             dmac_t3=D.*(avg_CW*slot_time);
             dmac_t2=(1./(1-P_D)).*(dmac_t2 + (P_suc.^j).*(T_s+ j.*T_c+ dmac_t3));   %(Working) %sum over i term of equation 40 in paper
             %dmac_t2=(1./(1-P_D)).*(dmac_t2 + (P_suc.^j).*(T_s+ j.*(T_c +T_ce)+ dmac_t3));        %sum over i term of equation 40 in paper
             %dmac_t2=(1./(1-P_D)).*(dmac_t2 + (P_suc).*(T_s+ i.*T_c + dmac_t3));        %sum over i term of equation 40 in paper
        else
            avg_CW = ((2^m)*(W-1))/2;
            dmac_t3=D.*(avg_CW*slot_time);
            dmac_t2=(1./(1-P_D)).*(dmac_t2 + (P_suc.^m).*(T_s+ m.*T_c + dmac_t3)); %(Working) 
            %dmac_t2=(1./(1-P_D)).*(dmac_t2 + (P_suc.^m).*(T_s+ m.*(T_c + T_ce)+ dmac_t3));
            %dmac_t2=(1./(1-P_D)).*(dmac_t2 + (P_suc).*(T_s+ m.*T_c + dmac_t3));        %sum over i term of equation 40 in paper
        end
        
    end
    
 end
%**********************************

% dmac_t2=0.0;
% for i=1:1:m
%     avg_CW=167.5;
%     dmac_t3=D.*(avg_CW*slot_time);
%     for j= 1:i
%         if j<m
%     dmac_t2=(1./(1-P_D)).*(dmac_t2 + (P_suc.^i).*(T_s+ i.*T_c + dmac_t3));        %sum over i term of equation 40 in paper
%         else
%             dmac_t2=(1./(1-P_D)).*(dmac_t2 + (P_suc.^m).*(T_s+ m.*T_c + dmac_t3));
%         end
%     end
%  end
%*************************
% dmac_t3=0;
% for i=1:1:m
%     avg_CW=0;
%     for j = 1:i
%         if j<m
%             avg_CW = ((2^j)*(W-1))/2;
%         else
%             avg_CW = ((2^m)*(W-1))/2;
%         end
%         dmac_t3=D.*(avg_CW*slot_time);
%         dmac_t2=(1./(1-P_D)).*(dmac_t2 + (P_suc).*(T_s+ i.*T_c + dmac_t3));
%     end
%             %sum over i term of equation 40 in paper
%  end
%*************************************
Dmac = dmac_t2;              %Makama
Dmac_sim_54 = ((Dmac).*20)./(N(1,k));
q_prob=lamda./(1./Dmac);%(lamda-Dmac)./lamda;
r=rand;
S_54 = (P_suc.*P_tx.*packet)./Dmac; 
 S_TrfMdl_54 = ((packet.*20).*N)./(5.0+((Dmac).*20));%(((4.75+Dmac).*20).*N); 
% S_TrfMdl = ((20*(packet).*N)./((4.75)+(Dmac.*20))); 
end
%################Plots#################
%%%%%%%%%%%%%54 Mbps
Thr_sim_54=[50526.31579 346559.7772	677242.7216	1013304.877	1348236.496	1681598.871	2015634.902	2350934.692	2686158.531	3021859.348	3355609.058	3691014.652	4025741.811	4360951.209	4693850.676	5029073.787	5364345.606	5698872.885	6033790.375	6362185.533	6695354.891	7030122.656	7364890.369	7698417.559	8033131.281	8367083.454
];
err_thr_54=[5.48645E-12 10506.02327	4271.679503	5227.245297	5451.643617	3925.490118	4107.878068	4852.857182	5483.956931	6137.418485	6497.778885	7101.621953	7612.324778	8306.373826	8595.238641	9183.696341	9795.942744	10482.71473	11155.62902	9227.628351	7003.115637	7353.434962	7703.683135	7887.303157	8230.175436	8637.158945
];
%OL=[50526.31579 404210.5263	808421.0526	1212631.579	1616842.105	2021052.632	2425263.158	2829473.684	3233684.211	3637894.737	4042105.263	4446315.789	4850526.316	5254736.842	5658947.368	6063157.895	6467368.421	6871578.947	7275789.474	7680000	8084210.526	8488421.053	8892631.579	9296842.105	9701052.632	10105263.16
%];
 %S_TrfMdl_54=[45937.7010630267 366946.770011481 732584.234122290 1096844.53133950 1459656.71685299 1820946.02307375 2180632.55899625 2538629.36122852 2894839.79028544 3249154.33212641 3601446.90875041 3951570.82084992 4299354.44218867 4644596.76419364 4987062.85622863 5326479.26838210 5662529.36501316 5994848.54328195 6323019.26458964 6646565.81008267 6964948.66488801 7277558.43963246 7583709.25193794 7882631.51504962 8173464.11614983 8455246.01460884];

 %####################################
 %******* Sim vs Analytical **********
% subplot(2,1,1);
% plot(N, OL, 'kp--');hold on;
% plot(N, S_TrfMdl, 'ks--');hold on;
% errorbar(N,Thr_sim_12,err_thr_12, 'k*--');hold on;
% title('Analytical Vs Simulation Throughput');xlabel('Number of geophones');ylabel('Throughput (bps)'); legend('Theoritical Offered Load','Analytical throughput (12 Mbps)','Simulated throughput (12 Mbps)');
% grid on;  
% subplot(2,1,2);
% plot(N, OL, 'kp--');hold on;
figure();
plot(N, S_TrfMdl_54,'s-','color', [0, 0.4470, 0.7410],'LineWidth', 1.5,'MarkerSize',6);hold on;
errorbar(N,Thr_sim_54,err_thr_54, '*-', 'color', [0.8500 0.3250 0.0980],'LineWidth', 1.5,'MarkerSize',6);hold on;
title('Analytical Vs Simulation Throughput');xlabel('Number of geophones');ylabel('Throughput (bps)');legend('Analytical','Simulation');%legend('Analytical throughput (54 Mbps)','Simulated throughput (54 Mbps)');
%legend('Theoritical Offered Load','Analytical throughput (54 Mbps)','Simulated throughput (54 Mbps)');
% grid on; 
%plot(N, S_TrfMdl, 'ks--');

figure();
plot(N, q_prob(:,:,end), 's-','color', [0, 0.4470, 0.7410],'LineWidth', 2,'MarkerSize',6);axis([0 200 0.00 1]);%hold on;
%######################################
% MAC_D_sim_54= [0.000315947	0.000311657	0.000315745	0.000322432	0.000325114	0.000327515	0.000331204	0.00033532	0.000341607	0.000346181	0.00034733	0.000349679	0.00035536	0.000359303	0.0003645	0.000373867	0.000377744	0.000384968	0.000389881	0.000395269	0.00040129	0.00040691	0.000411299	0.000420276	0.00042498	0.000431316
% ];
% 
% err_MAC_D_54= [6.56202E-16	9.77124E-06	6.5672E-06	1.04476E-05	1.07374E-05	8.15356E-06	9.89266E-06	8.83824E-06	1.42677E-05	1.21098E-05	1.41912E-05	1.3998E-05	1.42944E-05	1.20029E-05	1.00803E-05	1.10767E-05	1.33059E-05	1.39478E-05	1.28095E-05	1.71911E-05	1.70533E-05	1.25725E-05	1.34697E-05	1.50259E-05	1.32219E-05	1.30039E-05
% ];

%err_MAC_D= [0 7.07595E-05	4.45748E-05	6.15561E-05	5.91491E-05	5.27806E-05	6.35494E-05	4.88364E-05	8.72931E-05	8.35872E-05	8.97175E-05	0.000106468	0.000108018	0.000103798	9.77835E-05	0.000170373	0.000170857	0.000300435	0.000342006	0.000406274	0.000538557	0.00088955	0.00138065 0 0 0];

%******** Without sub-plot ***********
% figure();
% plot(N, Dmac_sim_54(:,:,end), 'ks-');axis([0 200 0.00 0.001]);hold on;
% errorbar(N,MAC_D_sim_54,err_MAC_D_54);
% title('Delay: Analytical Vs Simulation');xlabel('Number of geophones');ylabel('Delay (s)'); legend('Analytical','Simulation');

%######################################
% figure();
% % plot(N, OL, 'kp--');hold on;
%  plot(N, S_TrfMdl, 'ks--');%hold on;
%  title('Analytical Vs Simulation Throughput');xlabel('Number of geophones');ylabel('Throughput (bps)');
% % errorbar(N,Thr_sim_12,err_thr_12, 'k*--');hold on;
% % %plot(N, OL, 'kp--');hold on;
% % plot(N, S_TrfMdl_54, 'ko-');hold on;
% % errorbar(N,Thr_sim_54,err_thr_54, 'kd-');hold on;
% % 
% % title('Analytical Vs Simulation Throughput');xlabel('Number of geophones');ylabel('Throughput (bps)'); legend('Theoritical Offered Load','Analytical throughput (12 Mbps)','Simulated throughput (12 Mbps)','Analytical throughput (54 Mbps)','Simulated throughput (54 Mbps)');


% plot(N, OL, 'kp--');hold on;
% plot(N, S_TrfMdl_54, 'ks--');hold on;
% errorbar(N,Thr_sim_54,err_thr_54, 'k*--');
% title('Analytical Vs Simulation Throughput');xlabel('Number of geophones');ylabel('Throughput (bps)'); legend('Theoritical Offered Load','Analytical throughput','Simulated throughput');




%#####################################################################################
%Delay Plot
% N1=0:8:176;
N1=0:8:200;
C1=(0.2).*N1;  
N1(1)=1; 
C1(1)=1; 
for k1 = 1:size(N1, 2) 
    if (N1(k1) == 1)
        p1 = 0.0;
        tau_d = 2/(W + 1);
        
      else
        inc = 0.000001;
        error_MAX = 0.000005;
        p1= 0.0;
        tau1_d = 0.0;
        tau2_d = tau1_d + 2*error_MAX;

        while abs(tau1_d - tau2_d) > error_MAX
            p1 = p1 + inc;
%             tau1 = 2*Q*(1-2*p)*(1-(p).^(m+1))./(Q*((W*(1-p)*(1-(2*p).^(m+1))+ (1-2*p)*(1-(p).^(m+1))))+2*(1-Q)*(1-p)*(1-2*p));
            tau1_d = 2*Q*P_cd1*(1-2*p1)*(1-(p1).^(m+1))./(Q*((W*(1-p1)*(1-(2*p1).^(m+1))+ (1-2*p1)*(2*P_cd1-1)*(1-(p1).^(m+1))))+2*(1-Q)*P_cd1*(1-p1)*(1-2*p1));
            tau2_d = 1 - (1-p1).^(1./(N1-1));
            
        end
        tau_d = (tau1_d + tau2_d)./2;          %avg value of both tau
        
    end
%*******************Probabilities***************

    %Probability of collision
    P_c1 = 1- (1-tau_d).^(N1-1);
    %Probability of failure due to channel error
    P_ce1 = 0.05;          % (2%)
    %Probability of packet failure (P_f)
    P_f1 = P_c1+P_ce1-(P_c1*P_ce1);
    %Packet drop probability
    P_D1 = P_f1.^(m+1);
    %Probability of Tx (probability that there is at least one transmission in a the considered time slot)
    P_tx1 = 1 - (1- tau_d).^N1(1,k1);
    %Probability of successful Tx (prob. that exactly one geophone Txmits)
    P_TxG1 = N1(1,k1).*tau_d.*(1- tau_d).^(N1(1,k1)-1);
    %Probability that a packet transmission ongoing on the channel is successful
    % P_suc1=1-P_f1;
    P_suc1=P_TxG1./P_tx1;
    %P_suc=(N(1,k)*tau*(1- tau)^(N(1,k)-1))./(1 - (1- tau)^N(1,k));
    %Probability of empty slot time
    P_e1 = 1 - P_tx1;
    %Probability of collision
    %P_c = P_tr - P_s;
    %Probability counter decrements
    P_cd1 = (1-tau_d).^(N1(1, k1)-1);             
    
  %  E_s = P_e1.*slot_time + (P_suc1.*P_tx1).*T_s + (P_tx1.*(1-P_suc1)).*T_c;
%****************************************
%Average slot duration
%E_s = P_e*slot_time + P_suc*T_s + P_f*T_c;
%E_s = P_e.*slot_time + (P_suc.*P_tx).*T_s + (P_tx.*(1-P_suc)).*T_c;
%E_s = P_e*slot_time + P_suc*P_tx*(1-P_ce)*T_s + P_tx*(1-P_suc)*T_c+P_tx*P_suc*P_ce*T_ce; %include channel error prob.

%CW_avgS1 = 0.0;
 for i=1:m
    for j= 1:i
     if j<m
   CW_avgS1 = (((1-P_f1).*(P_f1)./(1-P_D1)).*((2^j)*(W-1)));
     else
       CW_avgS1 =(((1-P_f1).*(P_f1)./(1-P_D1)).*((2^m)*(W-1)));  
     end
    end
 end

%CW_avgS1 = 167.5;        %calculated manually
%**********************Probabilities (pci,pcs and pcc)*************************

pci1= (1-(1./CW_avgS1)).^C1;
pcs1= C1.*(1./CW_avgS1).*(1-(1./CW_avgS1)).^(C1-1);
pcc1= 1-pci1-pcs1;

%**********************Probabilities (pbi,pbs and pbc)*************************
pbi1=(1-tau_d).^(N1-1);
bs1=nchoosek(175,1);             %199 = N-1          %**fix the problem wit N-1 here
pbs1= bs1.*tau_d.*(1-tau_d).^(N1-2);
pbc1=1-pbi1-pbs1;

%**********************Probabilities (pss and psi)*************************
pss1=1/W;                   %1/CWmin
psi1=1-pss1;
%**********************Durations sigmaI, sigmaS and sigmaC*************************
%sigmaI= 1; %*slot_time;
sigmaI= slot_time;
sigmaS1= (1/(1-pss1))*T_s + (psi1/(1-pss1))*sigmaI;

%sumC1=0.0;
for i=1:m
    for j=1:i
        if j<m
        sumC1 = j.*(pcc1.^j);
        %sumC1 = i.*pcc1.^i;
        else
         sumC1 = m.*(pcc1.^m);
        %sumC1 = m.*pcc1.^m;
        end
    end
end  

% sumC1=0.0;
% for i=1:1:m
%     sumC1 = sumC1 + i.*pcc1.^i;
% end  
 sigmaC1= sumC1*T_c + (pcs1./(1-pcc1))*sigmaS1 + (pci1./(1-pcc1))*sigmaI;

%**Durations for entering backoff stages from previous bckoff and Tx stage deltaB and deltaT***
 %Probability that the counter decrements
 Pcd1 = (1-tau_d).^(N1-1); 
%Pcd1=(1-tau_d).^(N1-1);     %Pcd=pbi
del1=sigmaI.*pbi1+sigmaS1.*pbs1+sigmaC1.*pbc1;


deltaB1=(1./Pcd1).*(del1);       

deltaT1=Q*((1-(1./CW_avgS1))).*del1;
%********Durations a geophone spends in a backoff stage (delta)********
D1= (1-tau_d).*deltaB1 + tau_d.*deltaT1;
%*************************
% dmac_t21=0.0;
% for i=1:m
%     for j= 1:i
%         %i=i+1
%         if j<m
%              avg_CW1 = ((2^j)*(W-1))/2;
%              dmac_t31=D1.*(avg_CW1*slot_time);
%              dmac_t21=(1./(1-P_D1)).*((P_suc1.^i).*(T_s+ i.*T_c + dmac_t31));        %sum over i term of equation 40 in paper
%         else
%             avg_CW1 = ((2^m)*(W-1))/2;
%             dmac_t31=D1.*(avg_CW1*slot_time);
%             dmac_t21=(1./(1-P_D1)).*((P_suc1.^m).*(T_s+ m.*T_c + dmac_t31));
%         end
%     end
%  end
%************************************
dmac_t21=0.0;
% avg_CW1 = 167.5;
% dmac_t31=D1*(avg_CW1*slot_time);
for i=1:1:m
    for j= 1:i
        %i=i+1
        if j<m
             avg_CW1 = ((2^j)*(W-1))/2;
             dmac_t31=D1*(avg_CW1*slot_time);
             dmac_t21=(1./(1-P_D1)).*((dmac_t21+((P_suc1.^j).*(T_s+ j.*T_c + dmac_t31)))); 
              %sum over i term of equation 40 in paper
        else
            avg_CW1 = ((2^m)*(W-1))/2;
            dmac_t31=D1.*(avg_CW1*slot_time);
            dmac_t21=(1./(1-P_D1)).*((dmac_t21+((P_suc1.^m).*(T_s+ m.*T_c + dmac_t31)))); 
%             dmac_t31=D1.*(avg_CW1*slot_time);
%             dmac_t21=(1./(1-P_D1)).*((P_suc1.^m).*(T_s+ m.*T_c + dmac_t31));
        end    
        
    end
%     dmac_t31=D1.*(avg_CW1*slot_time);
%     dmac_t21=(1./(1-P_D1)).*((dmac_t21+((P_suc1).*(T_s+ i.*T_c + T_ce+ dmac_t31))));   
%     dmac_t21=(1./(1-P_D1)).*((dmac_t21+((P_suc1.^i).*(T_s+ i.*T_c + dmac_t31))));   
 end

%************************************
%###############
%working
% dmac_t21=0.0;
% for i=1:m
%     avg_CW1=900.5;
%     dmac_t31=D1.*(avg_CW1*slot_time);
%     for j= 1:i
%         if j<m
%         dmac_t21=(1./(1-P_D1)).*(dmac_t21 + (P_suc1.^i).*(T_s+ i.*T_c + dmac_t31));        %sum over i term of equation 40 in paper
%         else
%             dmac_t21=(1./(1-P_D1)).*(dmac_t21 + (P_suc1.^m).*(T_s+ m.*T_c + dmac_t31));
%         end
%     end
%  end
%#############
% dmac_t21=0.0;
% for i=1:m
%    % avg_CW=167.5;
%     
%  dmac_t31=D1.*(avg_CW*slot_time);
%     
%     %dmac_t21=(1./(1-P_D1)).*(dmac_t21 + (P_suc1).*(T_s+ i.*T_c + dmac_t31));        %sum over i term of equation 40 in paper
%     %dmac_t21=((1-P_D1).*(dmac_t21 + (P_suc1).*(T_s+ i.*T_c + dmac_t31)));
%  end

Dmac1 = dmac_t21;              %Makama
Dmac_sim = ((Dmac1).*20)./(N1(1,k1));
%Dmac_sim = (Dmac1.*20)./N1(1,k1);
end

%####### PLOTs #################
%****Original
% MAC_D_sim_54= [0.000315947	0.000311657	0.000315745	0.000322432	0.000325114	0.000327515	0.000331204	0.00033532	0.000341607	0.000346181	0.00034733	0.000349679	0.00035536	0.000359303	0.0003645	0.000373867	0.000377744	0.000384968	0.000389881	0.000395269	0.00040129	0.00040691	0.000411299	0.000420276	0.00042498	0.000431316
% ];
% 
% err_MAC_D_54= [6.56202E-16	9.77124E-06	6.5672E-06	1.04476E-05	1.07374E-05	8.15356E-06	9.89266E-06	8.83824E-06	1.42677E-05	1.21098E-05	1.41912E-05	1.3998E-05	1.42944E-05	1.20029E-05	1.00803E-05	1.10767E-05	1.33059E-05	1.39478E-05	1.28095E-05	1.71911E-05	1.70533E-05	1.25725E-05	1.34697E-05	1.50259E-05	1.32219E-05	1.30039E-05
% ];
%**
MAC_D_sim_54= [0.000	0.000311657	0.000315745	0.000322432	0.000325114	0.000327515	0.000331204	0.00033532	0.000341607	0.000346181	0.00034733	0.000349679	0.00035536	0.000359303	0.0003645	0.000373867	0.000377744	0.000384968	0.000389881	0.000395269	0.00040129	0.00040691	0.000411299	0.000420276	0.00042498	0.000431316
];

err_MAC_D_54= [0.00	9.77124E-06	6.5672E-06	1.04476E-05	1.07374E-05	8.15356E-06	9.89266E-06	8.83824E-06	1.42677E-05	1.21098E-05	1.41912E-05	1.3998E-05	1.42944E-05	1.20029E-05	1.00803E-05	1.10767E-05	1.33059E-05	1.39478E-05	1.28095E-05	1.71911E-05	1.70533E-05	1.25725E-05	1.34697E-05	1.50259E-05	1.32219E-05	1.30039E-05
];
%******** Without sub-plot ***********
% figure();
% plot(N1, Dmac_sim(:,:,end), 'ks-');axis([0 180 0.00 0.014]);hold on;
% errorbar(N1,MAC_D_sim,err_MAC_D);
% title('Delay: Analytical Vs Simulation');xlabel('Number of geophones');ylabel('Delay (s)'); legend('Analytical','Simulation');
%***************************************************
%%%%%%%%%%%%%%%%%
%Retransmission (%)
Rt=[0	0	0	0	0	0	0	0	0	0	0	0	0	0.105769231	0.178571429	0.291666667	0.34375	0.492647059	0.402777778	0.562091503	0.603125	0.601190476	0.602272727	0.885869565	0.755208333	0.975];
Rt_err=[0	0	0	0	0	0	0	0	0	0	0	0	0	0.082349892	0.14439841	0.356849533	0.33031844	0.366147147	0.293242595	0.558982362	0.438831048	0.227496275	0.299417155	0.373827682	0.288734141	0.30079508];

%******* Plot with collision prob. sub-plot *********************
figure();
%subplot(2,1,1);
plot(N1, Dmac_sim(:,:,end), 's-','color', [0, 0.4470, 0.7410],'LineWidth', 2,'MarkerSize',6);axis([0 200 0.00 0.01]);hold on;
errorbar(N1,MAC_D_sim_54,err_MAC_D_54, '*-','color', [0.8500, 0.3250, 0.0980],'LineWidth', 1.5,'MarkerSize',6 );
title('Delay: Analytical Vs Simulation');xlabel('Number of geophones');ylabel('Delay (s)'); legend('Analytical','Simulation');
% grid on; 
% subplot(2,1,2);
% plot(N1, P_f1(:,:,end), 's-','color', [0, 0.4470, 0.7410],'LineWidth', 1.5,'MarkerSize',6);axis([0 180 0.00 1.0]);
% title('Collision Probability');xlabel('Number of geophones');ylabel('Collision Probability'); %legend('Analytical','Simulation');
% grid on; 
figure();
%subplot(2,1,1);
%plot(N1, Dmac_sim(:,:,end), 's-','color', [0, 0.4470, 0.7410],'LineWidth', 2,'MarkerSize',6);axis([0 200 0.00 0.01]);hold on;
errorbar(N1,Rt,Rt_err, '*-','color', [0.8500, 0.3250, 0.0980],'LineWidth', 1.5,'MarkerSize',6 );
title('% Retransmission');xlabel('Number of geophones');ylabel('retransmission(%)'); %legend('Analytical','Simulation');


%Collision prob.
% figure();
% plot(N1, P_f1(:,:,end), 'ks-');%axis([0 180 ]);