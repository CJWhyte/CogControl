%Cognitive Control Active Inference Agent; CBU application

%Christopher Whyte 23/11/2019

clear 
clc
rng('shuffle')

%% Simulation settings

%To simulate a reliabel cue z = 1
%To simulate an unreliable cue set z = 2

z = 1;

%% Level 1
%==========================================================================

% prior beliefs about initial states: D
%--------------------------------------------------------------------------

D{1} = [1 0 0 0]';% Stimuli {Cue, position 1, position 2, response}
D{2} = [1 1]';% Rule {1, 2}
D{3} = [1 0 0]';% Response mapping {null, left, right}

% probabilistic mapping from hidden states to outcomes: A
%--------------------------------------------------------------------------
%--Outcome modality 1: Stimulus presentation

for i = 1:2
    for j = 1:3
            A{1}(:,:,i,j) = [1 0 0 0; %cue
                             0 1 0 0; %position 1
                             0 0 1 0; %position 2
                             0 0 0 1];%response
    end                
end

if z == 1
    a = 1;
    b = 0;
elseif z == 2
    a = .6;
    b = .4;
end 

%--Outcome modality 2: Rule
for i = 1:2
    for j = 1:3
        A{2}(:,:,i,j) = [a .5 .5 .5;%Rule 1
                         b .5 .5 .5];%Rule 2
    end 
end


for i = 2
    for j = 1:3
        A{2}(:,:,i,j) = [b .5 .5 .5;%Rule 1
                         a .5 .5 .5];%Rule 2
    end
end



%--Outcome modality 3:Feedback
for i = 1:2
    for j = 1:3
        A{3}(:,:,i,j) = [0 0 0 0;%correct
                         0 0 0 0;%incorrect
                         1 1 1 1];%null
    end 
end

%Rule 1 - click left
A{3}(:,:,1,2) = [0 1 1 1;%correct
                 0 0 0 0;%incorrect
                 1 0 0 0];%null
%Rule 1 - click right
A{3}(:,:,1,3) = [0 0 0 0;%correct
                 0 1 1 1;%incorrect
                 1 0 0 0];%null
%Rule 2 - click left           
A{3}(:,:,2,2) = [0 0 0 0;%correct
                 0 1 1 1;%incorrect
                 1 0 0 0];%null
%Rule 2 - click right
A{3}(:,:,2,3) = [0 1 1 1;%correct
                 0 0 0 0;%incorrect
                 1 0 0 0];%null

a = A;

a{1} = a{1}*64;
a{2} = a{2}*64;
a{3} = a{3}*64;

% Transitions between states: B
%-------------------------------------------------------------------------

B{1} = [0 0 0 0;
        .5 0 0 0;
        .5 0 0 0; 
        0 1 1 1];

B{2} = eye(2,2);


B{3}(:,:,1) = [1 1 1;
               0 0 0;
               0 0 0];
B{3}(:,:,2) = [0 0 0;
               1 1 1;
               0 0 0];
B{3}(:,:,3) = [0 0 0;
               0 0 0;
               1 1 1];




% Policies
%--------------------------------------------------------------------------

T = 3; %number of timesteps
Nf = 3; %number of factors
Pi = 2; %number of policies
V = ones(T-1,Pi,Nf);

V(:,:,3) = [1 1;
            2 3];

% C matrices (outcome modality by timesteps)
%--------------------------------------------------------------------------

C{1} = zeros(4,5);

C{2} = zeros(2,5);

c = 2;   
i = -2;

C{3} = [c c c c c;%correct
        i i i i i;%incorrect
        0 0 0 0 0]; %null

% MDP Structure
%--------------------------------------------------------------------------

mdp.T = T;                      % number of updates
mdp.A = A;                      % observation model
mdp.a = a;
mdp.s = 1;
mdp.B = B;                      % transition probabilities
mdp.C = C;                      % preferred outcomes
mdp.D = D;                      % prior over initial states
mdp.V = V;                      % policies


MDP = spm_MDP_check(mdp);



%% Plot 
MDP = spm_MDP_VB_X(MDP);


spm_figure('GetWin','trial'); clf
spm_MDP_VB_trial(MDP(1));
spm_figure('GetWin','LFP'); clf
spm_MDP_VB_LFP(MDP(1),[],2);

