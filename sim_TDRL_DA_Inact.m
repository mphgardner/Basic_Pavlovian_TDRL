function results = sim_TDRL_DA_Inact(sim)

% This script simulates the experimental paradigms included in the paper by
%Maes et al. entitled:

% Causal evidence supporting the proposal that dopamine transients function 
% as a temporal difference prediction error. https://doi.org/10.1101/520965.


    %
    % USAGE: results = sim_TDRL(sim)
    %
    % INPUTS:
    %   sim - string specifying simulation:
    %           'Blocking' - blocking paradigm with inhibition of DA
    %           'Second_Order' - second order conditioning paradigm with inhibition
    %           of DA

    %
    % OUTPUTS:
    %   results - see linearTDRL.m
    %
    %Matt Gardner November 2018
    


%This provides a transfer function for determining Pavlovian conditioned responses
%from value
CR = @(V,max,Voff,ks) max./(1 + exp(-ks*(V - Voff)));
%parameters used:
Voff = .4;
ks = 3;
%Different maximal responding were used depending on whether the animals
%were connected to the patch fibers or not
MaxFree = 55;
MaxFiber = 40;


switch sim
    case 'Blocking'
        
        %Stimuli used in the Pavlovian Blocking Design. Each row of the arrays
        %represents a single time step. Each column represents whether an
        %element was present for the particular stimulus in a binary manner.
        %randTrials adds the ITI between trials
        A =  [1 0 0 0 0; 0 0 0 0 0];
        B =  [0 0 0 0 1; 0 0 0 0 0];
        AX = [1 1 0 0 0; 0 0 0 0 0];
        AY = [1 0 1 0 0; 0 0 0 0 0];
        BZ = [0 0 0 1 1; 0 0 0 0 0];                
        X =  [0 1 0 0 0; 0 0 0 0 0];
        Y =  [0 0 1 0 0; 0 0 0 0 0];
        Z =  [0 0 0 1 0; 0 0 0 0 0];

        %This sets the alphas for the SR and U learning:
        alpha = 0.05;
        
        R = struct();
            
        %stimuli are randomized within each execution of the TDRL model.
        %This is included for designs in which the order of stimulus
        %presentations affects the value.
        for i = 1:100
            
            %this randomizes trial order for each stage of conditioning
            S{1} = randTrials(48,A);%Trials: 48 8*6
            S{2} = randTrials(24,A,B);%Trials: 24 4*6
            S{3} = randTrials(24,AX,AY,BZ);%Trials: 24 4*6
            S{4} = randTrials(8,X,Y,Z);%Trials: 8 2*4
            
            %This gets a vector of the stage numbers. This is used incase
            %there are identical stimuli in each stage.
            stg = stagenum(S);
            
            M = cell2mat(S');
            
            %logicals of the first state of each stimulus
            L.a = fstate(A,M) & stg == 1;
            L.a2 = fstate(A,M) & stg == 2;
            L.b = fstate(B,M) & stg == 2;
            L.ax = fstate(AX,M) & stg == 3;
            L.ay = fstate(AY,M) & stg == 3;
            L.bz = fstate(BZ,M) & stg == 3;
            L.x = fstate(X,M) & stg == 4;
            L.y = fstate(Y,M) & stg == 4;
            L.z = fstate(Z,M) & stg == 4;
            
            FN = fieldnames(L);

            %rewarded states
            r = double(nstate(L.a,2) | nstate(L.a2,2) | nstate(L.ax,2) | nstate(L.ay,2) | nstate(L.bz,2));
            
            %This implements tonic inhibition during just A in AX
            %trials
            opto = -1*(L.ax) + -1*nstate(L.ay,-1);
            
            %Initial values of each of the stimuli. We set cross-modality
            %generalization at 0.2 and generalization of same modality
            %stimuli at 0.7
            ui = [0; 0.2; 0.2; 0.2; 0.7];
            
            results = linearTDRL_DA(M,r,ui,alpha,opto,'none');
            for k = 1:numel(FN)
                if strcmp(FN{k},'a')
                    R(1).(FN{k})(i,:) = CR([results(L.(FN{k})).V],MaxFree,Voff,ks);
                else    
                    R(1).(FN{k})(i,:) = CR([results(L.(FN{k})).V],MaxFiber,Voff,ks);
                end    
            end
            
            results = linearTDRL_DA(M,r,ui,alpha,opto,'tonic');
            for k = 1:numel(FN)
                R(2).(FN{k})(i,:) = CR([results(L.(FN{k})).V],MaxFiber,Voff,ks);
            end
            
            results = linearTDRL_DA(M,r,ui,alpha,opto,'V');
            for k = 1:numel(FN)
                R(3).(FN{k})(i,:) = CR([results(L.(FN{k})).V],MaxFiber,Voff,ks);
            end

        end
           
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
        %The following code creates a panel with line plots of each stage of conditioning
        %as a different column, as well as a final column with a bar plot
        %of the behavior during the extinction test
        
        
        %This determines the subplot locations for the panel
        plotwidths = [0 .2 .1 .2 .2];
        plotheights = [.15 .15 .15];
        spcx = [0.1 0.01 0.01 0.1];
        spcy = [0 0.15 0.15];
        positions = cell(3,4);
        
        for sx = 1:4
            for sy = 1:3
                
                positions{4-sy,sx} = [sum(spcx(1:sx)) + sum(plotwidths(1:sx)), sum(spcy(1:sy)) + sum(plotheights(1:sy)), plotwidths(sx + 1),plotheights(sy)];
            end
        end
        
        V = zeros(numel(R),3);
        
        for k = 1:numel(R)
            V(k,:) = mean([R(k).x(:,1) R(k).y(:,1) R(k).z(:,1)]);
        end
        
        figure
        set(gcf, 'Units', 'Normalized', 'OuterPosition', [.1, 0.1, .4, .8]);
        hold on
        models = {'Control', 'Error', 'V'};
        for m = 1:3
            
            subplot('Position',positions{m,4})
            b = bar(V(m,:));
            
            %legend({'X' 'Y' 'Z'},'FontSize',15,'Location','North');
            set(gca,'FontSize',11,'XLim', [.5 3.5],'YLim',[0 40],'XTickLabel', {'X' 'Y' 'Z'},'YTick', 0:10:50);
            xlabel(models{m},'FontSize',15)
            ylabel('CR','FontSize',15);
            %title('Blocking','FontSize',25,'FontWeight','Bold');
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            subplot('Position',positions{m,1})
            title(models{m})
            Resp = [];
            for i = 1:8
                range = (i-1)*6+1:i*6;
                Resp(i,:) = mean(mean(R(1).a(:,range),2));
            end
            plot(1:8,Resp','-ok')
            
            set(gca,'FontSize',11, 'XTick', 1:8, 'XLim',[0.5 8.5],'YLim',[0 50],'YTick', 0:10:50)
            ylabel('CR','FontSize',15);
            xlabel('Days','FontSize',15);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            subplot('Position',positions{m,2})
            title(models{m})
            Resp = [];
            for i = 1:4
                range = (i-1)*6+1:i*6;
                Resp(i,:) = mean([mean(R(m).a2(:,range),2) mean(R(m).b(:,range),2)]);
            end
            plot(9:12,Resp(:,1),'-ok')
            hold on
            plot(9:12,Resp(:,2),'-o','color',[0.5 0.5 0.5])
            set(gca,'FontSize',11, 'XTick', 9:12, 'XLim',[8.5 12.5],'YLim',[0 50],'YTick',[]);%0:10:50)
            xlabel('Days','FontSize',15);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            subplot('Position',positions{m,3})
            title(models{m})
            Resp = [];
            for i = 1:4
                range = (i-1)*6+1:i*6;
                Resp(i,:) = mean([mean(R(m).ax(:,range),2) mean(R(m).ay(:,range),2)+.5 mean(R(m).bz(:,range),2)]);
            end
            plot(1:4,Resp','-o')
            
            set(gca,'FontSize',11, 'XTick', 1:4, 'XLim',[0.5 4.5],'YLim',[0 50],'YTick',[]);
            %  ylabel('CR','FontSize',15);
            xlabel('Days','FontSize',15);
        end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
        
       
    case 'Second_Order'
        
        %randTrials adds the ITI between trials
        A =  [1 0 0; 1 0 0; 0 0 0];
        C_A = [0 1 0; 1 0 0; 0 0 0];
        D_A = [0 0 1; 1 0 0; 0 0 0];
        C = [0 1 0];
        D = [0 0 1];
        
        %This sets the alphas for the SR and U learning:
        alpha = [0.05];
        
        %stimuli are randomized within each execution of the TDSR model.
        for i = 1:100
            
            %this randomizes trial order for each stage of conditioning
            S{1} = randTrials(6,A);%48 8*6
            S{2} = randTrials(12,C_A,D_A);%12 2*6
            S{3} = randTrials(4,C,D);%8 2*4
            
            %This gets a vector of the stage numbers. This is used incase
            %there are identical stimuli in each stage.
            stg = stagenum(S);
            
            M = cell2mat(S');
            
            %logicals of the first state of each stimulus
            L.a = fstate(A,M) & stg == 1;
            L.c_a = fstate(C_A,M) & stg == 2;
            L.d_a = fstate(D_A,M) & stg == 2;
            L.ac = nstate(L.c_a,2) & stg == 2;
            L.ad = nstate(L.d_a,2) & stg == 2;
            L.c = fstate(C,M) & stg == 3;
            L.d = fstate(D,M) & stg == 3;
            
            FN = fieldnames(L);

            %rewarded states
            r = double(nstate(L.a,2));
            
            opto = -1*nstate(L.c_a,2) + -1*nstate(L.d_a,-1);
            
            %Rats had prior conditioning to A (see methods)
            ui = [0.95 0 0]';
            
            results = linearTDRL_DA(M,r,ui,alpha,opto,'none');
            for k = 1:numel(FN)
                R(1).(FN{k})(i,:) = CR([results(L.(FN{k})).V],MaxFiber,Voff,ks);
            end
            
            results = linearTDRL_DA(M,r,ui,alpha,opto,'tonic');
            for k = 1:numel(FN)
                R(2).(FN{k})(i,:) = CR([results(L.(FN{k})).V],MaxFiber,Voff,ks);
            end

            results = linearTDRL_DA(M,r,ui,alpha,opto,'V');
            for k = 1:numel(FN)
                R(3).(FN{k})(i,:) = CR([results(L.(FN{k})).V],MaxFiber,Voff,ks);
            end

        end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
        %The following code creates a panel with line plots of each stage of conditioning
        %as a different column, as well as a final column with a bar plot
        %of the behavior during the extinction test
        
        %this determines the locations of the subplots
        plotwidths = [0 .1 .2 .4];
        plotheights = [.15 .15 .15];
        spcx = [0.1 0.01 0.1];
        spcy = [0 0.15 0.15];
        
        positions = cell(3,3);
        for sx = 1:3
            for sy = 1:3
                
                positions{4-sy,sx} = [sum(spcx(1:sx)) + sum(plotwidths(1:sx)), sum(spcy(1:sy)) + sum(plotheights(1:sy)), plotwidths(sx + 1),plotheights(sy)];
            end
        end
        
        V = zeros(numel(R),2);
        for k = 1:numel(R)
            V(k,:) = mean([R(k).c(:,1) R(k).d(:,1)]);
        end
        
        figure
        set(gcf, 'Units', 'Normalized', 'OuterPosition', [.1, 0.1, .4, .8]);
        hold on
        models = {'Control', 'Error', 'V'};
        for m = 1:3
            
            subplot('Position',positions{m,3})
            b = bar(V(m,:));
            %legend({'X' 'Y' 'Z'},'FontSize',15,'Location','North');
            set(gca,'FontSize',11,'XLim', [.5 2.5],'YLim',[0 20],'XTickLabel', {'C' 'D'},'YTick', 0:10:20);
            xlabel(models{m},'FontSize',15)
            ylabel('CR','FontSize',15);
            %title('Blocking','FontSize',25,'FontWeight','Bold');
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            subplot('Position',positions{m,1})
            title(models{m})
            
            plot(mean(mean(R(m).a)),'-ok')
            set(gca,'FontSize',11, 'XTick', 1, 'XLim',[0.5 1.5],'YLim',[0 50],'YTick',0:10:50);
            ylabel('CR','FontSize',15);
            xlabel('Days','FontSize',15);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
            subplot('Position',positions{m,2})
            title(models{m})
            Resp = [];
            for i = 1:2
                range = (i-1)*6+1:i*6;
                Resp(i,:) = mean([mean(R(m).c_a(:,range),2) mean(R(m).d_a(:,range),2)+.5 mean(R(m).ac(:,range),2) mean(R(m).ad(:,range),2)]);
            end
            plot(1:2,Resp','-o')
            
            set(gca,'FontSize',11, 'XTick', 1:2, 'XLim',[0.5 2.5],'YLim',[0 50],'YTick',[]);%0:10:50)
            %ylabel('CR','FontSize',15);
            xlabel('Days','FontSize',15);
        end
        
end
end

function S1 = fstate(i,S)
%Provides a logical vector (N x 1) of the first state of 
%stimulus i by template matching i within the state matrix J. This script
%will check up to three rows of stimulus i within the state matrix J.
%Several rows of each stimulus are checked incase the first row or two of
%stimulus i are identical to another stimulus.
    S1 = false(size(S,1),1);
    
    if size(i,1) == 1
        S1 = all(S == i(1,:),2);
    elseif size(i,1) == 2
        S1 = all(S == i(1,:),2) & (all([S(2:end,:); zeros(1,size(i,2))] == i(2,:),2));
    elseif size(i,1) > 2  
        S1 = all(S == i(1,:),2) & (all([S(2:end,:); zeros(1,size(i,2))] == i(2,:),2)) & (all([S(3:end,:); zeros(2,size(i,2))] == i(3,:),2));
    end    
end

function Sn = nstate(X,n)
%Provides a logical vector of state n within a given stimulus by shifting
%the logicals of stimulus X by n values. For example, n = 2 provides an N x 1
%logical vector of the second row of the stimulus in X, where X is an N x 1
%logical vector. If n = -1, this returns the logicals for the prior state.

    if n >= 0
        Sn = [false(n-1,1); X(1:end-n+1)];
    else
        Sn = [X(-1*n+1:end); false(-1*n,1)];
    end
end

function Stage = stagenum(S)
%provides an N X 1 vector of the stage number of conditioning for
%the state matrix J. J is a cell array whose size is 1 X the number of
%stages.
    Stage = zeros(size(cell2mat(S'),1),1);
    
    ind_start = 1;
    for i = 1:numel(S)
        ind_end = ind_start + size(S{i},1) - 1;
        
        Stage(ind_start:ind_end) = i;
        ind_start = ind_end + 1;
        
    end    

end

function out = randTrials(n,varargin)
% This function randomizes presentations of stimuli for the sim_TDRL function

%INPUTS:
%n is the number of presentations for each stimulus sequence

%the varargin are the unique stimulus sequences to be included 
    %each stimulus is an [s X D] matrix. (s = # of states within a single
    %stimulus sequence, D = # of features 

%OUTPUT: [n*s*(nargin - 1) X D] of randomized stimulus presentations
    
%Note that zeros(1 X D) are added between each presentation

out = [];
ncue = nargin - 1;
curr = 0;
for i = 1:n

    
    A = randperm(ncue,ncue);
    
    for j = 1:ncue
        
        str = size(varargin{A(j)},1);
        
        for k = 1:str
            curr = curr + 1;
            out(curr, :) = varargin{A(j)}(k,:);
        end
        curr = curr + 1;
        out(curr, :) = zeros(1, size(varargin{A(j)},2));
    end
end    
end    
