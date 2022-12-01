% Running an age depth model based on Loss on Ignition values
% Make sure that your excel file is in the right format with the following
% order: Sample thickness - Crusable nr. - Lab ID - Bot depth - Top depth -
% Mid depth - Weight crusable - Wet weight - Dry weight - Organic weight

% When encountering a coring gap one should manually put the thickness of
% the gap in the leftmost collumn and calculate the average weight values
% that can be expected during that coring gap.

% When one wants to model a hiatus, one can manually insert a collumn
% between which samples a hiatus is expected. Fill in a fictional 
% thickness (guess) and calculate the average weight values during that
% interval using excel. 

function [out, fig, f] = Agedepthmodelling(inputfilename,inputfiletabname,exponent,timefix,depthfix,fixnames,tiepnts_polx,tiepnts_poly,polnames,choice,method,endpoint);

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

%% Calculating and generating the input values

dat = xlsread(inputfilename,inputfiletabname);   %'Excel Timmelteich.xlsx', 'TT-11 composite'
inp.middepth = dat(:,6);
inp.weights = dat(:,(11:14));
inp.DBD = (inp.weights(:,2))./dat(:,1);
inp.CD = (inp.weights(:,4))./dat(:,1);
inp.KlassD = (inp.weights(:,3))./dat(:,1);
inp.ratio = inp.KlassD./inp.CD;
inp.LOI = ((inp.weights(:,4))./(inp.weights(:,2)))*100;

d_1 = depthfix(1:(length(depthfix)-1));
t_1 = timefix(1:(length(timefix)-1));

d_2 = depthfix(2:length(depthfix));
t_2 = timefix(2:length(timefix));

%% Modelling the Age-depth relation for the full model

fig = figure,

for m = 1:4

    subplot(2,4,m)

    for s = exponent
        
        % Defining all the methodologies
        
        if m == 1; select = inp.LOI.^(-s); title('LOI^{-x} More organic is higher sed. rate');
        elseif m == 2; select = inp.LOI.^(s); title('LOI^x More organic is lower sed. rate');
        elseif m == 3; select = inp.DBD.^(s); title('\rho_{dry}^x High dry bulk density is lower sed. rate');
        elseif m == 4; select = inp.ratio.^(s); title('( \rho_{inorganic}/\rho_{organic} )^x More organic is higher sed. rate');
        end
        
        % Defining the top and bottom points of the model
        
        x1 = find(inp.middepth == d_1(1));
        x2 = find(inp.middepth == d_2(length(d_2)));
        x3 = find(inp.middepth == endpoint);
        
        % Calculating the mu value 
        
        deviding = (t_1(1)-t_2(length(d_2)))/nansum(select(x1:x2).*dat(x1:x2,1));
        
        % Calculating the Age
        
        for n = x1:x3
            if n == x1
                Age(x1) = t_1(1);
            else
                Age(n) = Age(n-1) - (select(n)*dat((n),1))*deviding;
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
    
    hold on, scatter(depthfix,timefix,50,'xk')
    for b = 1:length(timefix); hold on, text(depthfix(b)-4,timefix(b),fixnames(b));end

    y = 0; 
    for v = 1:length(tiepnts_polx)/2
        hold on, plot(tiepnts_polx((v+y):(v+y+1)),tiepnts_poly((v+y):(v+y+1)),'-^','LineWidth',2); 
        hold on, text(tiepnts_polx((v+y):(v+y+1)),tiepnts_poly((v+y):(v+y+1)),polnames(v));
        y = y+1; 
    end
    
    % Calculating a nice neat looking 10 of covered depth limitation factor
    % for figure plotting
    
    dlimfact = ((max(depthfix)-endpoint)/100)*10;
    tlimfact = ((max(timefix)-7500)/100)*10; 
    xlim([endpoint-dlimfact,max(depthfix)+dlimfact]);ylim([7500-tlimfact,max(timefix)+tlimfact]);
   
    xlabel('Depth(cm)');ylabel('Age (cal yr BP)');

    Diffcomp.full{m} = Diffcomp.maxi{m} - Diffcomp.mini{m};
    
end

lgd = legend(string(exponent));title(lgd,'^{x} (exponent)');

%% Segmentbased age modelling

for m = 1:4

    subplot(2,4,m+4)

    for s = exponent
        
        for h = 1:length(d_1)

            % Defining all the methodologies

            if m == 1; select = inp.LOI.^(-s); title('LOI^{-x} More organic is higher sed. rate');
            elseif m == 2; select = inp.LOI.^(s); title('LOI^x More organic is lower sed. rate');
            elseif m == 3; select = inp.DBD.^(s); title('\rho_{dry}^x High dry bulk density is lower sed. rate');
            elseif m == 4; select = inp.ratio.^(s); title('( \rho_{inorganic}/\rho_{organic} )^x More organic is higher sed. rate');
            end

            % Defining the top and bottom points of the model
        
            x1 = find(inp.middepth == d_1(h));
            if h == 1; targ = x1; end
            x2 = find(inp.middepth == d_2(h));
            if h == length(d_1); x3 = find(inp.middepth == endpoint); end
            
            % Calculating the mu value
            
            deviding(h) = (t_1(h)-t_2(h))/sum(select(x1:x2).*dat(x1:x2,1));

            % Calculating the Age
            if h == length(d_1)
                for n = x1:x3
                    if n == x1
                        Age(x1) = t_1(h);
                    else
                        Age(n) = Age(n-1) - (select(n)*dat((n),1))*deviding(h);
                    end
                end
            else
                for n = x1:x2
                    if n == x1
                        Age(x1) = t_1(h);
                    else
                        Age(n) = Age(n-1) - (select(n)*dat((n),1))*deviding(h);
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
    
    hold on, scatter(depthfix,timefix,50,'xk')
    for b = 1:length(timefix); hold on, text(depthfix(b)-4,timefix(b),fixnames(b));end

    y = 0; 
    for v = 1:length(tiepnts_polx)/2
        hold on, plot(tiepnts_polx((v+y):(v+y+1)),tiepnts_poly((v+y):(v+y+1)),'-^','LineWidth',2); 
        hold on, text(tiepnts_polx((v+y):(v+y+1)),tiepnts_poly((v+y):(v+y+1)),polnames(v));
        y = y+1; 
    end
    
    % Calculating a nice neat looking 10 of covered depth limitation factor
    % for figure plotting
    
    dlimfact = ((max(depthfix)-endpoint)/100)*10;
    tlimfact = ((max(timefix)-7500)/100)*10; 
    xlim([endpoint-dlimfact,max(depthfix)+dlimfact]);ylim([7500-tlimfact,max(timefix)+tlimfact]);
   
    xlabel('Depth(cm)');ylabel('Age (cal yr BP)');
    Diffcomp.part{m} = Diffcomp.partmax{m} - Diffcomp.partmin{m};
    
end

lgd = legend(string(exponent));title(lgd,'^{x} (exponent)');

figure,
for m = 1:4
subplot(2,1,1),plot(inp.middepth(1:508),abs(Diffcomp.full{m})), hold on,; grid; legend('LOI^x','LOI^{-x}','\rho_{dry bulk}','( \rho_{inorganic}/\rho_{organic} )^x'); xlabel('Depth(cm)'); ylabel('\delta Age'); title('Full residuals')
subplot(2,1,2), plot(inp.middepth(targ:x3),abs(Diffcomp.part{m}(targ:x3))),hold on,; grid; legend('LOI^x','LOI^{-x}','\rho_{dry bulk}','( \rho_{inorganic}/\rho_{organic} )^x'); xlabel('Depth(cm)'); ylabel('\delta Age'); title('Partial residuals')
end

%% Single age depth model plotting and sedimentation rate calculation

s = choice;

f = figure, subplot(2,1,1),;

    for s = choice
        
        for h = 1:length(d_1)

            % Defining all the methodologies

            if method == 1; select = inp.LOI.^(-s); title('LOI^{-x} More organic is higher sed. rate');
            elseif method == 2; select = inp.LOI.^(s); title('LOI^x More organic is lower sed. rate');
            elseif method == 3; select = inp.DBD.^(s); title('\rho_{dry}^x High dry bulk density is lower sed. rate');
            elseif method == 4; select = inp.ratio.^(s); title('( \rho_{inorganic}/\rho_{organic} )^x More organic is higher sed. rate');
            end

            % Defining the top and bottom points of the model
        
            x1 = find(inp.middepth == d_1(h));
            if h == 1; targ = x1; end
            x2 = find(inp.middepth == d_2(h));
            if h == length(d_1); x3 = find(inp.middepth == endpoint); end
            
            % Calculating the mu value
            
            deviding(h) = (t_1(h)-t_2(h))/sum(select(x1:x2).*dat(x1:x2,1));

            % Calculating the Age
            if h == length(d_1)
                for n = x1:x3
                    if n == x1
                        Age(x1) = t_1(h);
                    else
                        sedr(n) = ((select(n)*deviding(h)).^-1)*10;
                        Age(n) = Age(n-1) - (select(n)*dat((n),1))*deviding(h);
                    end
                end
            else
                for n = x1:x2
                    if n == x1
                        Age(x1) = t_1(h);
                    else
                        sedr(n) = ((select(n)*deviding(h)).^-1)*10;
                        Age(n) = Age(n-1) - (select(n)*dat((n),1))*deviding(h);
                    end
                end
            end
        end

        % Plotting the model with special issues like choice and exponent 1

        x = inp.middepth(targ:x3)
        xx = [x,x]
        yy = [Diffcomp.partmax{1}(targ:x3);Diffcomp.partmin{1}(targ:x3)]
        fill(xx',yy,'k','EdgeColor',[0.7 0.6 0.8],'LineWidth',3); hold on,
            
        yyaxis right,plot(inp.middepth(targ:x3),sedr(targ:x3),'--k','LineWidth',1)
        if s == choice; yyaxis left, plot(inp.middepth(targ:x3),Age(targ:x3),'k','LineWidth',2);grid;hold on,
        elseif s == 1; yyaxis left, plot(inp.middepth(targ:x3),Age(targ:x3),'--k','LineWidth',1);grid;hold on,
        else yyaxis left, plot(inp.middepth(targ:x3),Age(targ:x3));grid; hold on,
        end
        
    end

    % Plotting the tie and fix points
    
    hold on, scatter(depthfix,timefix,50,'xk')
    for b = 1:length(timefix); hold on, text(depthfix(b)-4,timefix(b),fixnames(b));end

    y = 0; 
    for v = 1:length(tiepnts_polx)/2
        hold on, plot(tiepnts_polx((v+y):(v+y+1)),tiepnts_poly((v+y):(v+y+1)),'-^','LineWidth',2); 
        hold on, text(tiepnts_polx((v+y):(v+y+1)),tiepnts_poly((v+y):(v+y+1)),polnames(v));
        y = y+1; 
    end

% Calculating a nice neat looking 10 of covered depth limitation factor
% for figure plotting

dlimfact = ((max(depthfix)-endpoint)/100)*10;
tlimfact = ((max(timefix)-7500)/100)*10;
xlim([endpoint-dlimfact,max(depthfix)+dlimfact]);ylim([7500-tlimfact,max(timefix)+tlimfact]);

xlabel('Depth(cm)');ylabel('Age (cal yr BP)');

% lgd = legend('Age depth','','','tie points','Corylus estimate','Tilia estimate','Sedimentation rates');title(lgd,'Legend');

subplot(2,1,2), area(inp.middepth(targ:x3),abs(Diffcomp.part{method}(targ:x3))/2,'FaceColor','red','EdgeColor','k','FaceAlpha',.7),grid; legend('LOI^x','LOI^{-x}','\rho_{dry bulk}','( \rho_{inorganic}/\rho_{organic} )^x'); xlabel('Depth(cm)'); ylabel('\delta Age'); title('Section split residuals')
xlim([endpoint-dlimfact,max(depthfix)+dlimfact]);

% Defining the output

out.Depth = inp.middepth(targ:x3);
out.Age = Age(targ:x3)';
out.sedr = sedr(targ:x3);
out.LOI = inp.LOI(targ:x3);
out.difffull = abs(Diffcomp.full{method})/2
out.diffpart = abs(Diffcomp.part{method}(targ:x3))/2

end
