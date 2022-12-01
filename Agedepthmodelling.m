%% --- INTRO ---
% Lithology based age depth modelling
% s.j.schouten - Universitat Bern - 2020 orignal

% Code last modified 14/06/2022
% Compatibility with any data signal

%% --- PREPARING DATA INTO SHEETS ---
% Running an age depth model based on signal values
% Make sure that your excel file is in the right format with the following
% order: 
% Data_TAB: Depth - proxy 1, proxy 2 etc. 
% Tiepoints_TAB - Point information, Depth(compatible with excel sheet), Age(BP)
% Checkpoints_TAB - Point information, Depth(compatible with excel sheet), Age(BP)

% Model forces the extrapolation trough the tiepoints not the checkpoints,
% The checkpoints just get plotted for reference

%% --- HANDLING CORING GAPS ---
% When encountering a coring gap one should manually put the thickness of
% the gap in the depth collumn and calculate the average signal values
% (typically an mix of average values before and after the gap)
% that can be expected during that coring gap.

%% --- HANDLING HIATUSSUS ---
% When one wants to model a hiatus, one can manually insert a collumn
% between which samples a hiatus is expected. Fill in a fictional 
% thickness (guess) and calculate the average weight values during that
% interval using excel. 

%% --- RUNNING AS MATLAB SCRIPT FUNCTION---
% function [out, fig, f] = Agedepthmodelling(inputfilename,inputfiletabname,exponent,timefix,depthfix,fixnames,tiepnts_polx,tiepnts_poly,polnames,choice,method,endpoint);

%% --- TYPICAL INPUT OF FUNCTION ---
% inputfilename = 'Excel Retournemer 2.xlsx'
% inputfiletabname = 'Retournemer_2020_LOI_composite'
% exponent = [0.6:0.05:0.9];
% 
% depthfix = [2905.5,2859.5,2854.5,2596.5]; % depthfix = [2905.5,2895.5,2859.5,2854.5,2801.5,2596.5];
% timefix = [15842,12950,12900,9250]; % timefix = [15842,14542,12950,12900,11650,9250];
% fixnames = {'Pleniglacial end','GI-1 start','LST','LST','Start Holocene','Tilia'};
% 
% tiepnts_polx = [2762.5,2762.5,2595,2595];
% tiepnts_poly = [10800,11000,9400,9200];
% polnames = {'Corylus','Tilia'};
% endpoint = 2392.5;
% choice = 0.8;
% method = 1;

inputfilename = "Input_agedepthmodel.xlsx"              % Specify path to file of Excel
inputfiletabname = ["Data","Tiepoints","Checkpoints"]   % Specify sheet name
exponent = [0:0.1:2];                                   % Exponent range

input.data = flip(xlsread(inputfilename,inputfiletabname(1)));

[~,~, input.tpnts] = xlsread(inputfilename,inputfiletabname(2));
tp.names = (input.tpnts(:,1));
tp.depth = cell2mat(input.tpnts(:,2));
tp.time = cell2mat(input.tpnts(:,3));

[~,~, input.cpnts] = xlsread(inputfilename,inputfiletabname(3));
cp.names = (input.cpnts(:,1));
cp.depth = cell2mat(input.cpnts(:,2));
cp.time = cell2mat(input.cpnts(:,3));

signames = {'K^{s} (more K is les MAR)','Ti^{s} (more Ti is less MAR)','K^{-s} (more K is more MAR)','Ti^{-s} (more Ti is more MAR)'}

endpoint = 400.2; 
choice = -0.8; 
method = 4; 

%% Calculating and generating the input values

inp.middepth = input.data(:,1);
inp.samplewidth = -diff(inp.middepth(:,1));
inp.data = input.data(:,(2:size(input.data,2)));

ma = (size(inp.data,2))*2

for i = 1:(size(inp.data,2))
    inp.sig{i} = inp.data(:,i);
end

d_1 = tp.depth(1:(length(tp.depth)-1));
t_1 = tp.time(1:(length(tp.time)-1));

d_2 = tp.depth(2:length(tp.depth));
t_2 = tp.time(2:length(tp.time));

%% Modelling the Age-depth relation for the full model

fig = figure,

for m = 1:ma

    subplot(2,ma,m)
    
    for s = exponent

        if m <= 2; select = inp.sig{m}.^(s); title(signames{m});
        elseif m >= 2; select = inp.sig{m-2}.^(-s); title(signames{m});
        end

        % Defining the top and bottom points of the model
        
        x1 = find(inp.middepth == d_1(1));
        x2 = find(inp.middepth == d_2(length(d_2)));
        x3 = find(inp.middepth == endpoint);
        
        % Calculating the mu value 
        
        deviding = (t_1(1)-t_2(length(d_2)))/nansum(select(x1:x2).*inp.samplewidth(x1:x2));
        
        % Calculating the Age
        
        for n = x1:x3
            if n == x1
                Age(x1) = t_1(1);
            else
                Age(n) = Age(n-1) - (select(n)*inp.samplewidth(n))*deviding;
            end
        end

        % Plotting the model with special issues like choice and exponent 1
        
        if s == choice; plot(inp.middepth(x1:x3),Age(x1:x3),'k','LineWidth',2);grid;hold on,
        elseif s == 1; plot(inp.middepth(x1:x3),Age(x1:x3),'--k','LineWidth',1);grid;hold on,
        else; plot(inp.middepth(x1:x3),Age(x1:x3));grid; hold on,
        end
    
        if s == min(exponent)
           Diffcomp.mini{m} = Age;
        elseif s == max(exponent)
           Diffcomp.maxi{m} = Age;
        end

    end

    % Plotting the tie and fix points
    
    hold on, scatter(tp.depth,tp.time,50,'xk')
    for b = 1:length(tp.time); hold on, text(tp.depth(b)-4,tp.time(b),tp.names(b));end

    y = 0; 
    for v = 1:length(cp.depth)
        hold on, plot(cp.depth(v),cp.time(v),'o'); 
        hold on, text(cp.depth(v),cp.time(v),cp.names(v));
        y = y+1; 
    end
    
    % Calculating a nice neat looking 10 of covered depth limitation factor
    % for figure plotting
    
    dlimfact = ((max(tp.depth)-endpoint)/100)*10;
    tlimfact = ((max(tp.time)-7500)/100)*10; 
    xlim([endpoint-dlimfact,max(tp.depth)+dlimfact]);ylim([7500-tlimfact,max(tp.time)+tlimfact]);
   
    xlabel('Depth(cm)');ylabel('Age (cal yr BP)');

    Diffcomp.full{m} = Diffcomp.maxi{m} - Diffcomp.mini{m};
    
end

lgd = legend(string(exponent));title(lgd,'^{x} (exponent)');

%% Segmentbased age modelling

for m = 1:ma

    subplot(2,ma,m+4)

    for s = exponent
        
        for h = 1:length(d_1)

            % Defining all the methodologies

        if m <= 2; select = inp.sig{m}.^(s); title(signames{m});
        elseif m >= 2; select = inp.sig{m-2}.^(-s); title(signames{m});
        end
%             if m == 1; select = inp.LOI.^(-s); title('LOI^{-x} More organic is higher sed. rate');
%             elseif m == 2; select = inp.LOI.^(s); title('LOI^x More organic is lower sed. rate');
%             elseif m == 3; select = inp.DBD.^(s); title('\rho_{dry}^x High dry bulk density is lower sed. rate');
%             elseif m == 4; select = inp.ratio.^(s); title('( \rho_{inorganic}/\rho_{organic} )^x More organic is higher sed. rate');
%             end

            % Defining the top and bottom points of the model
        
            x1 = find(inp.middepth == d_1(h));
            if h == 1; targ = x1; end
            x2 = find(inp.middepth == d_2(h));
            if h == length(d_1); x3 = find(inp.middepth == endpoint); end
            
            % Calculating the mu value
            
            deviding(h) = (t_1(h)-t_2(h))/sum(select(x1:x2).*inp.samplewidth(x1:x2,1));

            % Calculating the Age
            if h == length(d_1)
                for n = x1:x3
                    if n == x1
                        Age(x1) = t_1(h);
                    else
                        Age(n) = Age(n-1) - (select(n)*inp.samplewidth((n),1))*deviding(h);
                    end
                end
            else
                for n = x1:x2
                    if n == x1
                        Age(x1) = t_1(h);
                    else
                        Age(n) = Age(n-1) - (select(n)*inp.samplewidth((n),1))*deviding(h);
                    end
                end
            end
        end

        % Plotting the model with special issues like choice and exponent 1
        
        if s == choice; plot(inp.middepth(targ:x3),Age(targ:x3),'k','LineWidth',2);grid;hold on,
        elseif s == 1; plot(inp.middepth(targ:x3),Age(targ:x3),'--k','LineWidth',1);grid;hold on,
        else plot(inp.middepth(targ:x3),Age(targ:x3));grid; hold on,
        end
        
        if s == max(exponent)
           Diffcomp.partmax{m} = Age;
        elseif s == min(exponent)
           Diffcomp.partmin{m} = Age;
        end
    end

    % Plotting the tie and fix points
    
    hold on, scatter(tp.depth,tp.time,50,'xk')
    for b = 1:length(tp.time); hold on, text(tp.depth(b)-4,tp.time(b),tp.names(b));end

    y = 0; 
    for v = 1:length(cp.depth)
        hold on, plot(cp.depth(v),cp.time(v),'.'); 
        hold on, text(cp.depth(v),cp.time(v),cp.names(v));
        y = y+1; 
    end
    
    % Calculating a nice neat looking 10 of covered depth limitation factor
    % for figure plotting
    
    dlimfact = ((max(tp.depth)-endpoint)/100)*10;
    tlimfact = ((max(tp.time)-7500)/100)*10; 
    xlim([endpoint-dlimfact,max(tp.depth)+dlimfact]);ylim([7500-tlimfact,max(tp.time)+tlimfact]);
   
    xlabel('Depth(cm)');ylabel('Age (cal yr BP)');
    Diffcomp.part{m} = Diffcomp.partmax{m} - Diffcomp.partmin{m};
    
end

lgd = legend(string(exponent));title(lgd,'^{x} (exponent)');

figure,
for m = 1:ma
subplot(2,1,1),plot(inp.middepth(1:x3),abs(Diffcomp.full{m})),hold on; grid; legend('LOI^x','LOI^{-x}','\rho_{dry bulk}','( \rho_{inorganic}/\rho_{organic} )^x'); xlabel('Depth(cm)'); ylabel('\delta Age'); title('Full residuals')
subplot(2,1,2), plot(inp.middepth(targ:x3),abs(Diffcomp.part{m}(targ:x3))),hold on; grid; legend('LOI^x','LOI^{-x}','\rho_{dry bulk}','( \rho_{inorganic}/\rho_{organic} )^x'); xlabel('Depth(cm)'); ylabel('\delta Age'); title('Partial residuals')
end

%% Single age depth model plotting and sedimentation rate calculation

s = choice;

f = figure, subplot(2,1,1);

    for s = choice
        
        for h = 1:length(d_1)

            % Defining all the methodologies

           select = inp.sig{method-2}.^(-s); title(signames{method});
%             if method == 1; select = inp.LOI.^(-s); title('LOI^{-x} More organic is higher sed. rate');
%             elseif method == 2; select = inp.LOI.^(s); title('LOI^x More organic is lower sed. rate');
%             elseif method == 3; select = inp.DBD.^(s); title('\rho_{dry}^x High dry bulk density is lower sed. rate');
%             elseif method == 4; select = inp.ratio.^(s); title('( \rho_{inorganic}/\rho_{organic} )^x More organic is higher sed. rate');
%             end

            % Defining the top and bottom points of the model
        
            x1 = find(inp.middepth == d_1(h));
            if h == 1; targ = x1; end
            x2 = find(inp.middepth == d_2(h));
            if h == length(d_1); x3 = find(inp.middepth == endpoint); end
            
            % Calculating the mu value
            
            deviding(h) = (t_1(h)-t_2(h))/sum(select(x1:x2).*inp.samplewidth(x1:x2,1));

            % Calculating the Age
            if h == length(d_1)
                for n = x1:x3
                    if n == x1
                        Age(x1) = t_1(h);
                    else
                        sedr(n) = ((select(n)*deviding(h)).^-1)*10;
                        Age(n) = Age(n-1) - (select(n)*inp.samplewidth((n),1))*deviding(h);
                    end
                end
            else
                for n = x1:x2
                    if n == x1
                        Age(x1) = t_1(h);
                    else
                        sedr(n) = ((select(n)*deviding(h)).^-1)*10;
                        Age(n) = Age(n-1) - (select(n)*inp.samplewidth((n),1))*deviding(h);
                    end
                end
            end
        end

        % Plotting the model with special issues like choice and exponent 1

        x = inp.middepth(targ:x3)
        xx = [x,x]
        yy = [Diffcomp.partmax{3}(targ:x3);Diffcomp.partmin{3}(targ:x3)]
        fill(xx',yy,'k','EdgeColor',[0.7 0.6 0.8],'LineWidth',3); hold on,
            
        yyaxis right,plot(inp.middepth(targ:x3),sedr(targ:x3),'--k','LineWidth',1)
        if s == choice; yyaxis left, plot(inp.middepth(targ:x3),Age(targ:x3),'k','LineWidth',2);grid;hold on,
        elseif s == 1; yyaxis left, plot(inp.middepth(targ:x3),Age(targ:x3),'--k','LineWidth',1);grid;hold on,
        else yyaxis left, plot(inp.middepth(targ:x3),Age(targ:x3));grid; hold on,
        end
        
    end

    % Plotting the tie and fix points
    
    hold on, scatter(tp.depth,tp.time,50,'xk')
    for b = 1:length(tp.time); hold on, text(tp.depth(b)-4,tp.time(b),tp.names(b));end

    y = 0; 
    for v = 1:length(cp.depth)
        hold on, plot(cp.depth(v),cp.time(v),'o'); 
        hold on, text(cp.depth(v),cp.time(v),cp.names(v));
        y = y+1; 
    end

% Calculating a nice neat looking 10 of covered depth limitation factor
% for figure plotting

dlimfact = ((max(tp.depth)-endpoint)/100)*10;
tlimfact = ((max(tp.time)-7500)/100)*10;
xlim([endpoint-dlimfact,max(tp.depth)+dlimfact]);ylim([7500-tlimfact,max(tp.time)+tlimfact]);

xlabel('Depth(cm)');ylabel('Age (cal yr BP)');

% lgd = legend('Age depth','','','tie points','Corylus estimate','Tilia estimate','Sedimentation rates');title(lgd,'Legend');

subplot(2,1,2), area(inp.middepth(targ:x3),abs(Diffcomp.part{method}(targ:x3))/2,'FaceColor','red','EdgeColor','k','FaceAlpha',.7),grid; legend('LOI^x','LOI^{-x}','\rho_{dry bulk}','( \rho_{inorganic}/\rho_{organic} )^x'); xlabel('Depth(cm)'); ylabel('\delta Age'); title('Section split residuals')
xlim([endpoint-dlimfact,max(tp.depth)+dlimfact]);

% Defining the output

out.Depth = inp.middepth(targ:x3);
out.Age = Age(targ:x3)';
out.sedr = sedr(targ:x3);
for j = 1:length(inp.sig)
    out.sig{j} = inp.sig{j}(targ:x3);
end
out.difffull = abs(Diffcomp.full{method})/2
out.diffpart = abs(Diffcomp.part{method}(targ:x3))/2

% end
