%% Energy consumption model (ESP32 MCU)
clc;
clear all;

% %************Data Packet Transmission*************************
%************Data Packet Transmission*************************
N=120;                                      %Number of neighbours in this case number of nodes in the network

%******Time constants*******%
SIFS = 0.00001; %10e-06; 
slot_time = 0.00002; %20e-06; 
DIFS = 0.00005; %50e-06;
delta= 0.0000020014; %1e-06;     %propagation delay(s) {d/c} from application scenario (largest distance to destination is 600.0m)
Tsim=20;                                    %20s simulation time
TudpStopTime=15.75;                         % UDPApp stop time
sendinterval= 0.25;                         %packet send interval
Tremaining_sim = Tsim-TudpStopTime;         %time remaining after the UDP stopTime till end of simulation (Tsim-TstopTime_UDP)
P_ack=194;                                  %194 bits including PHY layer 
D_rate=65000000;                            %65 Mbps
T_ACK= SIFS+(P_ack/D_rate);                 %Ack packet reception time
Tr=SIFS+(P_ack/D_rate)+delta;
%*****power consumption value in diff. mode of ESP32 MCU operation*******%
Pt= 0.594;                                  %Tx in watts, 180mA @3.3v VDD
Pr=0.33;                                    %Rx power in watts, 100mA @3.3v VDD
Pidle=49.5e-3;                              %P idle in watts, 15mA @ 3.3v VDD
Pslp=2.64e-3;                               %P sleep mode in watts, 0.8mA @ 3.3v VDD 

%********************Neighbor Discovery**********************
T_NS=1.011365623272451e-04;                           % Neighbor Solicitation (NS) msg Tx time from DCF model in chpt 4 with 65Mbps rate (946 bits)
T_NA=1.011365623272451e-04;                           % Neighbor Advertisement (NA) msg (946 bits) Tx time from DCF model in chpt 4 with 65Mbps rate

%EC_ND= (Pt*T_NS)*(N-1) +(Pr*T_NA)*(N-1);    %EC from ND TXsion and Rxtion

EC_ND= Pt*T_NS + (N-1)*Pr*T_NS +  Pt*T_NA + (N-1)*Pr*T_NA; 
%********************RPL Control Message**********************
sendinterval_cntrl= 10; 

Ttx_DIO=9.161618841683301e-05; %9.92e-05;  %Time to transmit DIO cntrl msg=Time to receive DIO 
Trx_DIO=9.161618841683301e-05; %9.92e-05;

Tidlecntrldio = sendinterval_cntrl - (Ttx_DIO*119);

% EC_RPL= (Pr*Trx_DIO) + (Pt*Ttx_DIO)*(N-1) +(Pr*Trx_DIO)*(N-1) + Pidle*(Tidlecntrldio);    %DIO broadcast to all Ngbrs, receive rom all Ngbrs, and 1 reception from Root/gateway

EC_RPL= (Pt*Trx_DIO)+ (Pr*Trx_DIO)*(N-1) + (Pt*Ttx_DIO) + Pr*Ttx_DIO*(N-1) + Pidle*(Tidlecntrldio); 

%********************Data packets**********************
K=5/sendinterval;                            %No. of packets sent by each node over 5 second recording period
Ttx=0.000434383;                             %Tx time from DCF model in chpt 4 with 65Mbps rate
%Trx=3.902844102564103e-05;                   %Tr same as transmission time

Tidle= sendinterval - (Ttx+T_ACK);             %Time interval between end of a packet transmission and the beginning of the next


%EC_data=(Pt*Ttx)*20 + (Pr*T_ACK)*20 + (Pidle*Tidle)*20;
EC_data=((Pt*Ttx) + Pr*Ttx +Pt*T_ACK + (Pr*Tr) + (Pidle*Tidle))*20;

%*********TOTAL Energy Consumed Per Sweep***********%

EC_Total= EC_ND + EC_RPL + EC_data +(Pidle*Tremaining_sim);

%% Simulation Validation with 1 Shot (Energy Consumed per Shot per Geophone
% %%%
% EC_SimvsAnly =[0.966564279]; %0.9605];   %1st column 20 MHz non-overlapping, 2nd 40 MHz partially ovlpn
% err_chpt6sec2 = [0.008816741]; 
% %EC_Total=[];
% %Analytic= [0.972506400715762];
% 
% figure;
% ax = gca;
% 
% % Define the positions for each set of bars
% x = 1:2;
% 
% % Plot the bars
% bar(ax, x(1), EC_SimvsAnly);
% hold(ax, 'on');
% errorbar(ax, x(1), EC_SimvsAnly, err_chpt6sec2, 'LineStyle', 'none', 'Color', 'k');
% bar(ax, x(2), EC_Total);
% 
% % Set the axis limits and labels
% xlim(ax, [0 3]);
% xticks(ax, x);
% xticklabels(ax, {'Simulation', 'Analytical'});
% %xlabel(ax, 'Data Sets');
% ylabel(ax, 'EC (Joules)');

%%%***************************************
EC_SimvsAnly =[0.966564279]; %0.9605];   %1st column 20 MHz non-overlapping, 2nd 40 MHz partially ovlpn
err_chpt6sec2 = [0.008816741]; 
%EC_Total=[];
%Analytic= [0.972506400715762];

figure;
ax = gca;

% Define the positions for each set of bars
x = 1:2;

% Plot the bars
bar(ax, x(1), EC_Total);
hold(ax, 'on');
bar(ax, x(2), EC_SimvsAnly);
hold(ax, 'on');
errorbar(ax, x(2), EC_SimvsAnly, err_chpt6sec2, 'LineStyle', 'none', 'Color', 'k','LineWidth', 2);


% Set the axis limits and labels
xlim(ax, [0 3]);
xticks(ax, x);
xticklabels(ax, {'Analytical','Simulation'});
%xlabel(ax, 'Data Sets');
ylabel(ax, 'EC (Joules)');

%%%% Plotting EC in individual phases from Analytical and compare with
%%%% simulation
figure;
ax = gca;

x = 1:5;

bar(EC_ND);
hold(ax, 'on');
bar(ax, x(2), EC_RPL);
hold(ax, 'on');

bar(ax, x(3), EC_data);
hold(ax, 'on');

bar(ax, x(4), EC_Total);

% Plot the bars
bar(ax, x(5), EC_SimvsAnly);
hold(ax, 'on');
errorbar(ax, x(5), EC_SimvsAnly, err_chpt6sec2, 'LineStyle', 'none', 'Color', 'k','LineWidth', 2);
% Set the axis limits and labels
xlim(ax, [0 6]);
xticks(ax, x);
xticklabels(ax, {'ND', 'RPL','Data Tx','Total Analytical','Total Simulation'});
% xticklabels(ax, {'Network setup ND', 'Network setup RPL','Data Tx/node/Shot','Total Analytical','Total Simulation'});
%xlabel(ax, 'Data Sets');
ylabel(ax, 'EC (Joules)');


% figure();
% Analytic= [0.978363164207432];
% EC_SimvsAnly =[	0.966564279; 0.978363164207432]; %0.9605];   %1st column 20 MHz non-overlapping, 2nd 40 MHz partially ovlpn
% err_chpt6sec2 = [0.008816741; 0]; 
% b = bar(EC_SimvsAnly, 'grouped');
% xticklabels ({'Simulation', 'Analytical',});
% %xticklabels ({'Hop Count OF', 'Energy OF' , 'ETX OF', 'EC Analytical',});
% title=('Average Energy Consumed Per node Per Shot');
% %xlabel('PDR', 'FontSize',16,'FontWeight','bold');
% ylabel ('EC (Joules)','FontSize',16,'FontWeight','bold');
% %lg = legend ( 'Hop Count OF', 'Energy OF' , 'ETX OF', 'Analytical', 'AutoUpdate' , 'off' );
% %lg = legend ( 'Hop Count', 'ETX', 'Energy' , 'AutoUpdate' , 'off' );
% lg.Location = 'northwest' ;
% %For error bar plot
% hold on
% % Find the number of groups and the number of bars in each group
% [ngroups, nbars] = size(EC_SimvsAnly);
% % Calculate the width for each bar group
% groupwidth = min(0.8, nbars/(nbars + 1.5));
% % Set the position of each error bar in the centre of the main bar
% for i = 1:nbars
%     % Calculate center of each bar
%     x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
%    % errorbar(model_series,model_error);
%   errorbar(x,EC_SimvsAnly(:,i), err_chpt6sec2(:,i), 'k', 'linestyle', 'none');
% end
% %hold on
% %bar(Analytic);
% hold off

% %############################### OLD ###################################
% %Tx, Rx, and Idle times
% Tt=0.000434383;                             %Tx time from DCF model in chpt 4 with 65Mbps rate
% Tr=3.902844102564103e-05;                   %Tr same as transmission time
% 
% Tidle= sendinterval - (Tt+Tr);                   %Time interval between end of a packet transmission
%                                             %and the beginning of the next
% 
% %Energy consumed for a single packet transmission
% 
% Ep = (Pt*Tt) + (Pr*Tr) + (Pidle*Tidle) 
% 
% %Energy consumed per shot with 20 pckts sent excluding Control packets
% 
% Ep_pernode = (Pt*Tt)*20 + (Pr*Tr)*20 + (Pidle*Tidle)*20 +(Pidle*Tremaining_sim);     % (Pr*Tidle)--> (Pr=Plisten)--> Rx listening during idle period
%                                                     % (Pidle*Tidle)
% % %Energy consumed per shot with 20 pckts sent excluding Control packets
% %  Esht= (Ep*20)+(Pidle*Tremaining_sim);
% % 
% % 
%  %************Control(DIO/DAO) Packet Transmission*************************
% sendinterval_cntrl= 10;                           %default DIO interval in RPL is set to be 10 seconds
% Tdio=9.92e-05;	        %From calculations in file "Dio_Dao_TxtimeCalc.m"
% %Tdao=1.09207e-05;           %From calculations in file "Dio_Dao_TxtimeCalc.m"
% 
% % Tdio=1.10749e-05;	        %From calculations in file "Dio_Dao_TxtimeCalc.m"
% % Tdao=1.09207e-05;           %From calculations in file "Dio_Dao_TxtimeCalc.m"
% 
% Tidlecntrldio = sendinterval_cntrl - (Tdio*120);
% %Tidlecntrldao = sendinterval_cntrl - (Tdao*120);
% 
% % Tidlecntrldio = sendinterval_cntrl - Tdio;
% % Tidlecntrldao = sendinterval_cntrl - Tdao;
% 
% %Energy consumed for dio packet transmission per shot
% 
% % Ep_cntrl = (Pr*Tdio) + (Pt*Tdio)*120 + (Pidle*Tidlecntrldio);         %Energy consumed per shot due to dio Messages
% Ep_cntrl = 2*(Pr*Tdio) + (Pt*Tdio)*120*2 + (Pidle*Tidlecntrldio);       %Pr is for 1st dio received rom gateway
%                                                                         % Dio broadcasted to all nodes twice throughout simulation, 
%                                                                         % Tidlecontroldio is after the 1st dio txmsn
% 
% %Ep_cntrlsht = Pt*(Tdio+Tdao) + Pr*(Tdio+Tdao) + Pidle*(Tidlecntrldio);
%  %Pidle*(Tidlecntrldio+Tidlecntrldao);
% 
% 
% % ECpershot= Esht+Ep_cntrlsht;    %Total energy consumed per shot in joules 
% 
%  %************Control(DAO) Packet Transmission*************************
% Tdio=9.92e-05;	        %From calculations in file "Dio_Dao_TxtimeCalc.m"
% 
% Tidlecntrldio = sendinterval_cntrl - (Tdio*120);
% 
% %*******Total Energy Consumed in Joule per Shot*********%
% 
% ECpershot= Ep_pernode+Ep_cntrl;    %Total energy consumed per shot in joules 
% %############################### OLD END ###################################

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% Tsim=20;                                    %20s simulation time
% TudpStopTime=15.75;                         % UDPApp stop time
% sendinterval= 0.25;                         %packet send interval
% Tremaining_sim = Tsim-TudpStopTime;         %time remaining after the UDP stopTime till end of simulation (Tsim-TstopTime_UDP)
% %power consumption value in diff. mode of ESP32 MCU operation
% Pt= 0.594;                                  %Tx in watts, 180mA @3.3v VDD
% Pr=0.33;                                    %Rx power in watts, 100mA @3.3v VDD
% Pidle=49.5e-3;                              %P idle in watts, 15mA @ 3.3v VDD
% Pslp=2.64e-3;                               %P sleep mode in watts, 0.8mA @ 3.3v VDD 
% 
% %Tx, Rx, and Idle times
% Tt=0.000434383;                             %Tx time from DCF model in chpt 4 with 65Mbps rate
% Tr=Tt;                                      %Tr same as transmission time
% Tidle= sendinterval - Tt;                   %Time interval between end of a packet transmission
%                                             %and the beginning of the next
% 
% %Energy consumed for a single packet transmission
% 
% Ep = (Pt*Tt) + (Pr*Tr) + (Pidle*Tidle);
% 
% Ep_pernode = (Pt*Tt)*20 + (Pr*Tidle)*20 + (Pidle*Tidle)*20;     % (Pr*Tidle)--> (Pr=Plisten)--> Rx listening during idle period
%                                                     % (Pidle*Tidle)
% %Energy consumed per shot with 20 pckts sent excluding Control packets
%  Esht= (Ep*20)+(Pidle*Tremaining_sim);
% 
% 
%  %************Control(DIO/DAO) Packet Transmission*************************
% sendinterval_cntrl= 10;                           %default DIO interval in RPL is set to be 10 seconds
% Tdio=1.10749e-05;	        %From calculations in file "Dio_Dao_TxtimeCalc.m"
% Tdao=1.09207e-05;           %From calculations in file "Dio_Dao_TxtimeCalc.m"
% Tidlecntrldio = sendinterval_cntrl - Tdio;
% Tidlecntrldao = sendinterval_cntrl - Tdao;
% 
% %Energy consumed for a single packet transmission
% 
% Ep_cntrl = (Pt*Tdio) + (Pr*Tdio) + (Pidle*Tidlecntrldio);       %Energy consumed per shot due to dio Messages
% 
% Ep_cntrlsht = Pt*(Tdio+Tdao) + Pr*(Tdio+Tdao) + Pidle*(Tidlecntrldio);
%  %Pidle*(Tidlecntrldio+Tidlecntrldao);
% 
% 
% %Total Energy Consumed in Joule per Shot
% 
% ECpershot= Esht+Ep_cntrlsht;    %Total energy consumed per shot in joules 



% %% Simulation Validation with 1 Shot (Energy Consumed per Shot per Geophone
% 
% figure();
% EC_SimvsAnly =[1.01197138; 1.001442857;	1.066622912; 0.9714]; %0.9605];   %1st column 20 MHz non-overlapping, 2nd 40 MHz partially ovlpn
% err_chpt6sec2 = [0.035653179; 0.031773763;	0.036471177; 0]; 
% Analytic=[0.9605];
% b = bar(EC_SimvsAnly, 'grouped');
% xticklabels ({'Hop Count OF', 'Energy OF' , 'ETX OF', 'EC Analytical',});
% title=('Average Energy Consumed Per node Per Shot')
% %xlabel('PDR', 'FontSize',16,'FontWeight','bold');
% ylabel ('EC (Joules)','FontSize',16,'FontWeight','bold');
% %lg = legend ( 'Hop Count OF', 'Energy OF' , 'ETX OF', 'Analytical', 'AutoUpdate' , 'off' );
% %lg = legend ( 'Hop Count', 'ETX', 'Energy' , 'AutoUpdate' , 'off' );
% lg.Location = 'northwest' ;
% %For error bar plot
% hold on
% % Find the number of groups and the number of bars in each group
% [ngroups, nbars] = size(EC_SimvsAnly);
% % Calculate the width for each bar group
% groupwidth = min(0.8, nbars/(nbars + 1.5));
% % Set the position of each error bar in the centre of the main bar
% for i = 1:nbars
%     % Calculate center of each bar
%     x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
%    % errorbar(model_series,model_error);
%   errorbar(x,EC_SimvsAnly(:,i), err_chpt6sec2(:,i), 'k', 'linestyle', 'none');
% end
% %hold on
% %bar(Analytic);
% hold off

