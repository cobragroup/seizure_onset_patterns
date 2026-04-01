% Author: Chris Fietkiewicz
% Reference: Christopher Fietkiewicz and Kenneth A. Loparo, �Analysis and Enhancements of a Prolific Macroscopic Model of Epilepsy,� Scientifica, vol. 2016, Article ID 3628247, 10 pages, 2016. doi:10.1155/2016/3628247
% Description: Original model used for Figure 2 in Fietkiewicz and Loparo 2016. It is based on the original model specified in Wendling F, % Bartolomei F, Bellanger JJ, Chauvel P. Epileptic fast activity can be
% explained by a model of impaired GABAergic dendritic inhibition. Eur J Neurosci, 2002;15(9):1499-1508.
function [eeg, N]=Wendl2005stimnoise(duration,aVals,bVals,gVals,S,N,fs,dt_scale,plotFlag)
% duration = [9 4 2 3 11 15];
% aVals = [3.5      4.6     7.7     8.7  ];
% bVals = [13.2     20.4    4.3     11.4 ];
% gVals = [10.76    11.48   15.1     2.1 ];

stepSize = 1/fs;
totalDuration = sum(duration)+length(duration)-1; %second part for 1 second long transitions
totalSteps = totalDuration / stepSize;

if or(S==0,isempty(S))
    S=zeros(1,totalSteps);
end

if isempty(N)
    N=randn(1,totalSteps);
elseif N==0
    N=zeros(1,totalSteps);
end

A = [];
B = [];
G = [];
for i = 1:length(duration)
    steps = duration(i) / stepSize;
    % segments
    aSegment = aVals(i) * ones(steps, 1);
    A = [A; aSegment];
    bSegment = bVals(i) * ones(steps, 1);
    B = [B; bSegment];
    gSegment = gVals(i) * ones(steps, 1);
    G = [G; gSegment];
    % transition
    if (i < length(duration) && length(duration)>1)
        if aVals(i) == aVals(i + 1)
            aTransition = aVals(i) * ones(1, 1/stepSize);
        else
            aStep = (aVals(i + 1) - aVals(i)) * stepSize; % Transitions are all 1 second
            aTransition = aVals(i) : aStep : aVals(i + 1);
        end
               
        if bVals(i) == bVals(i + 1)
            bTransition = bVals(i) * ones(1, 1/stepSize);
        else
            bStep = (bVals(i + 1) - bVals(i)) * stepSize; % Transitions are all 1 second
            bTransition = bVals(i) : bStep : bVals(i + 1);
        end
        
        if gVals(i) == gVals(i + 1)
            gTransition = gVals(i) * ones(1,  1/stepSize);
        else
            gStep = (gVals(i + 1) - gVals(i)) * stepSize; % Transitions are all 1 second
            gTransition = gVals(i) : gStep : gVals(i + 1);
        end
        A = [A; aTransition(1:1/stepSize)'];
        B = [B; bTransition(1:1/stepSize)'];
        G = [G; gTransition(1:1/stepSize)'];
    end
end

eeg = simulate(stepSize, totalDuration, A, B, G, S, N, dt_scale, 1, [0 0 0 0 0 0 0 0 0 0], plotFlag);


% Main simulation
function [eeg,epsp,ipsp_slow,ipsp_fast] = simulate(stepSize, endTime, A, B, G, S, N, dt_scale, transientTime, yInit, plotFlag)
steps = int32(endTime / stepSize);
transientSteps = int32(transientTime / stepSize);
eeg = zeros(steps, 1);
epsp = zeros(steps, 1);
ipsp_slow = zeros(steps, 1);
ipsp_fast = zeros(steps, 1);

y = yInit;

for i = 1:transientSteps
	y = takeStep(y, stepSize, A(1), B(1), G(1), S(i), N(i),dt_scale);
end

for i = 1:steps
	y = takeStep(y, stepSize, A(i), B(i), G(i), S(i), N(i),dt_scale);
	eeg(i) = y(2) - y(3) - y(4);
    epsp(i) = y(2);
    ipsp_slow(i) = y(3);
    ipsp_fast(i) = y(4);
end

if (plotFlag == 1)
    figure;
% 	plot(eeg);
% 	set(gcf, 'Position',[0 260 790 260]);
    t = stepSize * [0:length(eeg)-1];
    subplot(2,1,1);
    plot(t, eeg,t,S);
    xlabel('Time (sec)')
    ylabel('y-out (mV)')
    title('Adapted WendFietk Model')
    subplot(2,1,2);
    plot(t, epsp,t, ipsp_slow,t,ipsp_fast);
    ylabel('PSPs (mV)');
    legend('y1: EPSP','y2: IPSPslow','y3: IPSPfast')

end

% ***************************************************************************************************************
% Take a step
function yNew = takeStep(y, stepSize, A, B, G, S, N, dt_scale)

a = 100;%100;
b = 30; %30; %50;
g = 350;%350;
C = 135;
C1 = C;
C2 = 0.8 * C;
C3 = 0.25 * C;
C4 = 0.25 * C;
C5 = 0.3 * C; %0.1
C6 = 0.1 * C;
C7 = 0.8 * C;
MEAN = 90;
SIGMA = 30; 

yNew(1) = y(1) + y(6) * stepSize;
yNew(6) = y(6) + (A * a * Sig( S + y(2)-y(3)-y(4)) - 2 * a * y(6) - a * a * y(1)) * stepSize;
yNew(2) = y(2) + y(7) * stepSize;
if dt_scale
    yNew(7) = y(7) + (A * a * N *sqrt(stepSize) * sqrt(0.001) * SIGMA) + (A * a * (MEAN + C2 * Sig( S + C1* y(1))) - 2 * a * y(7) - a * a * y(2)) * stepSize;
else
    yNew(7) = y(7) + (A * a * (N * SIGMA + MEAN + C2 * Sig( S + C1* y(1))) - 2 * a * y(7) - a * a * y(2)) * stepSize;
end
yNew(3) = y(3) + y(8) * stepSize;
yNew(8) = y(8) + (B * b * (C4 * Sig( S + C3 * y(1)) ) - 2 * b * y(8) - b * b * y(3)) * stepSize;
yNew(4) = y(4) + y(9) * stepSize;
yNew(9) = y(9) + (G * g * C7 * Sig( S +(C5 * y(1) - C6 * y(5))) - 2 * g * y(9) - g * g * y(4)) * stepSize;
yNew(5) = y(5) + y(10) * stepSize;
yNew(10) = y(10) + (B * b * ( Sig( S + C3 * y(1)) ) - 2 * b * y(10) - b * b * y(5)) * stepSize;


% yNew(1) = y(1) + y(6) * stepSize;
% yNew(6) = y(6) + (A * a * Sig( y(2)-y(3)-y(4)) - 2 * a * y(6) - a * a * y(1)) * stepSize;
% yNew(2) = y(2) + y(7) * stepSize;
% if dt_scale
%     yNew(7) = y(7) + (A * a * N *sqrt(stepSize) * sqrt(0.001) * SIGMA) + (A * a * (MEAN + S + C2 * Sig( S + C1* y(1))) - 2 * a * y(7) - a * a * y(2)) * stepSize;
% else 
%    yNew(7) = y(7) + (A * a * randn * SIGMA + MEAN + S + C2 * Sig(C1* y(1))) - 2 * a * y(7) - a * a * y(2)) * stepSize;
% end
% yNew(3) = y(3) + y(8) * stepSize;
% yNew(8) = y(8) + (B * b * (C4 * Sig(C3 * y(1)) ) - 2 * b * y(8) - b * b * y(3)) * stepSize;
% yNew(4) = y(4) + y(9) * stepSize;
% yNew(9) = y(9) + (G * g * C7 * Sig((C5 * y(1) - C6 * y(5))) - 2 * g * y(9) - g * g * y(4)) * stepSize;
% yNew(5) = y(5) + y(10) * stepSize;
% yNew(10) = y(10) + (B * b * ( Sig(C3 * y(1)) ) - 2 * b * y(10) - b * b * y(5)) * stepSize;


% ***************************************************************************************************************
% Sigmoid function
function r = Sig(v)
r = 5 / (1.0 + exp(0.56 * (6.0 - v))); %2 x e0 = 5


