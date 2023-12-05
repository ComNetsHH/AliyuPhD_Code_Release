 clc;
 clear all;

% Constants
numNodes = 120;              % Number of nodes
numRows = 6;                % Number of rows
numColumns = 20;            % Number of columns
rowSeparation = 200;        % Separation between rows in meters
columnSeparation = 25;      % Separation between nodes in a row in meters
transmittedPower = 5e-3;    % Transmitted power in watts (5mW)
transmittingAntennaGain = 2;    % Transmitting antenna gain in dBi
receivingAntennaGain = 2;       % Receiving antenna gain in dBi
transmittingAntennaHeight = 1;  % Transmitting antenna height in meters
receivingAntennaHeight = 1;     % Receiving antenna height in meters
frequency = 2.4e9;              % Frequency in Hz
speedOfLight = 3e8;             % Speed of light in m/s
thresholdPower = -85;           % Minimum received power threshold in dBm
thresholdPower_mW=0.000000003162;% In milliwatts
dist_Tx_Rx = 650;                %largest distance vertically to the furthest node
% Calculate wavelength
wavelength = speedOfLight / frequency;

% Calculate free-space path loss
freeSpacePathLoss = 20 * log10((4 * pi * rowSeparation *frequency) / speedOfLight); %in dB
freeSpacePathLoss_gw = 20 * log10((4 * pi * dist_Tx_Rx *frequency) / speedOfLight); %in dB
%freeSpacePathLoss = 20 * log10(1) + 20 * log10(frequency) + 20 * log10((4 * pi)/ speedOfLight);
%freeSpacePathLoss_gw = 20 * log10(dist_Tx_Rx) + 20 * log10(frequency) + 20 * log10(4 * pi / speedOfLight);
%freeSpacePathLoss = (4 * pi * rowSeparation) ^ 2 / wavelength ^ 2;
%freeSpacePathLoss_gw = (4 * pi * dist_Tx_Rx) ^ 2 / wavelength ^ 2;
% Calculate gateway position/center coordinates
gatewayX = ((numColumns-1) / 2) * columnSeparation;
gatewayY = ((numRows-1) / 2) * rowSeparation;

% gatewayX = ((numColumns) / 2) * columnSeparation;
% gatewayY = ((numRows) / 2) * rowSeparation;

% Calculate effective communication range of gateway
effectiveRange = sqrt((numColumns / 2)^2 + (numRows / 2)^2) * columnSeparation;

% Calculate received power of gateway
receivedPower_gw = (10 * log10(transmittedPower) + transmittingAntennaGain + receivingAntennaGain + 20 * log10(transmittingAntennaHeight) + 20 * log10(receivingAntennaHeight)) ...
- (10 * log10(effectiveRange^4) + freeSpacePathLoss_gw);

%********************************************************************    
    % Calculate the distances from the center for each node
    distances = zeros(numRows, numColumns);
    receivedPower = zeros(numRows, numColumns);
    freeSpacePathLoss_nodes=zeros(numRows, numColumns);

    for i = 1:numRows
        for j = 1:numColumns
        x = (j-1) * columnSeparation - gatewayX;
        y = (i-1) * rowSeparation - gatewayY;
        distances(i, j) = sqrt(x^2 + y^2);
        
        % x = (j-1) * columnSeparation;
        % y = (i-1) * rowSeparation;
        % distances(i, j) = sqrt((x - gatewayX)^2 + (y - gatewayY)^2);

        % freeSpacePathLoss_nodes(i, j) = (4 * pi * distances(i, j)) ^ 2 / wavelength ^ 2;
        freeSpacePathLoss_nodes(i, j) = 20 * log10((4 * pi * distances(i, j) *frequency) / speedOfLight); %20 * log10(distances(i, j)) + 20 * log10(frequency) + 20 * log10(4 * pi * frequency / speedOfLight);
        receivedPower(i,j) = (10 * log10(transmittedPower) + transmittingAntennaGain + receivingAntennaGain + 20 * log10(transmittingAntennaHeight) + 20 * log10(receivingAntennaHeight))... 
        - (10 * log10((distances(i, j)^4) * freeSpacePathLoss_nodes(i, j)));
        end
    end

% Calculate effective communication range of gateway based on minimum rcvd power
GW_eff_range = ((transmittedPower * 10^(transmittingAntennaGain/10) * 10^(receivingAntennaGain/10) * transmittingAntennaHeight^2 * receivingAntennaHeight^2) ...
            /((10^(receivedPower_gw/10) * 10^(freeSpacePathLoss/10))))^0.25;    %thrshldpower in mW


% Count the number of nodes falling within different gateway communication radius ranges
numNodesInRange = zeros(1, 4);
for i = 1:numRows
    for j = 1:numColumns
        dist = distances(i, j);
        if dist <= GW_eff_range
            numNodesInRange(1) = numNodesInRange(1) + 1;
        elseif dist > GW_eff_range && dist <= 2*GW_eff_range
            numNodesInRange(2) = numNodesInRange(2) + 1;
        elseif dist > 2*GW_eff_range && dist <= 3*GW_eff_range
            numNodesInRange(3) = numNodesInRange(3) + 1;
        elseif dist > 3*GW_eff_range && dist <= 4*GW_eff_range
            numNodesInRange(4) = numNodesInRange(4) + 1;
        end
    end
end

% Display the results
disp(['Number of nodes within R: ' num2str(numNodesInRange(1))]);
disp(['Number of nodes be tween R and 2R: ' num2str(numNodesInRange(2))]);
disp(['Number of nodes between 2R and 3R: ' num2str(numNodesInRange(3))]);
disp(['Number of nodes between 3R and 4R: ' num2str(numNodesInRange(4))]);

% % % HD_Tx5mW_2_4ghz_sim = Hop_5mW_simvsmodel{:,1};
% % % HD_Tx5mW_2_4ghz_model = Hop_5mW_simvsmodel{:,3};
% % 
% % % HD_Tx5mW_2_4ghz_sim = Hop_5mW_simvsmodel_530m{:,1};
% % % HD_Tx5mW_2_4ghz_model = Hop_5mW_simvsmodel_530m{:,5};

HD_Tx5mW_2_4ghz_sim = Hop_5mW_simvsmodel_650m{:,1};
HD_Tx5mW_2_4ghz_model = Hop_5mW_simvsmodel_650m{:,7};

%plot
distributionFitter(HD_Tx5mW_2_4ghz_model);hold on;
distributionFitter(HD_Tx5mW_2_4ghz_sim);
hold off