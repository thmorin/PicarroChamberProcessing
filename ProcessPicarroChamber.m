function [HM_CO2, NRMSE_CO2, HM_CH4,NRMSE_CH4,BubbleFlux]=ProcessPicarroChamber(meas_date,start_time,end_time,L_data_sheet,Start_offset,End_offset,ChamberVolume,ChamberArea, plotYN,Site,Start_temp,End_temp,~,~)
% ProcessChamber - Written by Tim Morin, 7/5/2022
% Added ebullition and HM fit - Erin Hassett and Tim Morin 7/26/2022
% 
% meas_date - provide as a datetime, because that is how Picarro will read the time it must compare it to
% nb_start_time - datetime for start of chamber measurement, recorded in notebook with stopwatch
% nb_end_time - datetime for end of chamber measurement, recorded in notebook with stopwatch
% Start_temp - degrees C - Temp at start of measurement period
% End_temp - degrees C - Temp at end of measurent period
% start_time- epoch time recorded by picarro (Seconds since 1/1/1970)
% stop_time- epoch time recorded by picarro
% Start_offset - How many points to ignore at the start of the measurement period
% End_offset - How many points to ignore at the end of the measurement period
% Offset - Time measurement offset to GMT= unused currently

%nb_start_time=start_time+hours(Offset); %Replaced with the time from the SPF_lite file
%nb_end_time=end_time+hours(Offset); %Replaced with the time from the SPF_lite file
%nb_start_time=start_time-seconds(Shift_back); %irrelevant if you use the SPF_lite file. Same for line below
%nb_end_time=end_time-seconds(Shift_back);

%% Import raw Picarro data
laserdata=ImportPicarroFile(L_data_sheet);

%% filter to relevant day
dayFilt=laserdata.DATE==meas_date;
laserdata.TIME(~dayFilt)=NaT;

%% Remove bad diagnostic data. NOTE: Picarro alarm doesn't necessarily seem to indicate bad data
% bad=laserdata.ALARM_STATUS~=0;
% laserdata.CO2_dry(bad | ~dayFilt)=nan;
% laserdata.CH4_dry(bad | ~dayFilt)=nan;

%% Filter to relevant time of measurment
% tfilt = laserdata.TIME>=start_time & laserdata.TIME<=end_time;
% TIME=datenum(laserdata.TIME(tfilt));
% TIME=(TIME-TIME(1))*24*3600;    % convert from days to seconds
tfilt = laserdata.EPOCH_TIME>=start_time & laserdata.EPOCH_TIME<=end_time;%in between start time and end time
TIME=laserdata.EPOCH_TIME(tfilt);
TIME=TIME-TIME(1); %if error: Index exceeds the number of array elements (0). Make sure data_sheet is correct in metadata file. 

CO2=laserdata.CO2_dry(tfilt);
CH4=laserdata.CH4_dry(tfilt);

%% Correct for chamber volume and area
ChamberVolume=ChamberVolume/100/100/100;%originally in cm3
ChamberArea=ChamberArea/100/100; %conversion from cm2

%% Ideal gas law correction
TMP=(linspace(Start_temp,End_temp,length(TIME))+273.15)';   %Conversion to K, puts into vector
%R=8.205*10^(-5); %m3 atm K-1 mol-1
%CO2=CO2/R./TMP;
%CH4=CH4/R./TMP;

%% Trim out any data at beginning or end of measurement period that were affected by chamber settling
[~,IndSt]=min(abs(TIME-(Start_offset+1)));
[~,IndEn]=min(abs(TIME-(End_offset+1))); 

T_TIME=TIME(IndSt:(length(TIME)-IndEn));
% T_H2O =H2O (Start_offset+1:(length(TIME)-End_offset));
T_CO2 =CO2 (IndSt:(length(TIME)-IndEn));
T_CH4 =CH4 (IndSt:(length(TIME)-IndEn)); 

if length(T_CH4) < 50 %if you only had 50 points for a sample, don't count it

        HM_CH4 = nan; %or 0?
        %% Too few points to run regression on
        HM_CO2=nan;
        BubbleFlux=nan;
        NRMSE_CH4= nan;
        NRMSE_CO2=nan;
        %Bubble_RMSE=nan;
else

%% %Bubbles
%P = 99100;% laserdata.Pressure %Pa conversion from 991 hPa; add into function
P = 99300; %August
%P= 99600; %september
%P= 99400;  %October
R = 8.314; % m3 Pa K-1 mol-1

%EbulFlux = [];
%HM=[]; %hutchison mosier
%BubThr = 0.59; %July 2022 (ppm/s)
BubThr = 0.078; %August 2022
%BubThr=0.078; %September 2022
%BubThr = 0.05; %October 2022

CHDif = diff(T_CH4);  %creates matrix of changes in methane concentration as the sample progresses; 
                      %goes through SFP and calculates the diff of methane conc within start and stop times
fprintf('Chamber %s had max diff of %0.03f\n',Site,max(CHDif));
BubbleAdj = (diff(T_CH4)./ diff( T_TIME) > BubThr) | (diff(T_CH4)./diff( T_TIME) < -BubThr); %makes 1 or 0 if above or below threshold (| = or);creates a matrix of when bubbles cross the threshold in either (+) or (-) direction
%Bubble = diff(T_CH4) > BubThr; %only selecting non-negative bubbles - NOT USED. Allowing for negative bubbles to avoid reading instrument errors as bubbles; creates a matrix or bubbles that only cross threshold in (+) direction

CH4Add = zeros(size(CHDif));  %matrix with zeros + same size as CHDif
CH4Add(BubbleAdj) = CHDif(BubbleAdj);%pairs bubbleadjust with CHDif; when they have same paired event, 
                                       %it puts it into CH4add; (i.e. when
                                       %Bubble Adj =1; pulls value from
                                       %CHDiff)= all methane above thresh
TotBubble = cumsum(CH4Add);%adds up the total amount of methane attributed to bubbling events for the sample
T_CH4_NoBubbles = [T_CH4(1); T_CH4(2:end) - TotBubble]; %new matrix; T_CH4 start and stop; starts at 1 and goes from 2:end.
                                                        %creates a matrix of methane exluding CH4 from bubbling events= total diffusive flux
%%%%%[~,Ind] = max(T_CH4); %make new index; m=max methane for the whole run without bubbles 
Ind=length(T_CH4); %= it returns the row number with the largest methane in T_CH4
   % Ind is used to find the peak of methane for a sample,in order to deal with fitting issues...
    %typically if theres a high bubble, methane will then sink back to a lower concentration, making the fit really screwy. 
    %Instead they made Ind so they could stop the sample at it's peak, avoiding artifical sinks  
   
 if Ind == 1  %if our highest [methane] is at the first data point in the sample, then there was no flux due to bubbles if not, 
               %then we run calculations to get concentrations for the bubbles
    BubbleFlux = 0;    
 else
    BubbleFlux = (TotBubble(Ind-1)/(T_TIME(Ind) - T_TIME(1))) * ...%added T_TIME (Ind -1) instead of (Ind)= may need to change back
                  ((P.*ChamberVolume)./ (R*TMP(Ind).*ChamberArea)); % umol m-2 s-1  ... %take Ind (max flux value)and subtract 1 because cumsum makes one less value (leaves out first one)
%for all the ones that aren't 1, goes down index and uses start and stop times with the bubbles, 
%takes bubble events and when it pairs with the index, it takes the methane for those values and then runs through the equation
%edited so that there is only one value output: TMP (Ind)    
% TotBubble=TotBubble(1:Ind); %TIM CHANGED THIS LINE - MAKE SURE HE DIDN'T BREAK ANYTHING AFTER THIS
 end
 
%Bubbling = sum(diff(T_CH4(1:Ind)) > BubThr); % Number of bubbling events; %Currently unused
    %EbulFlux = [EbulFlux; BubbleFlux(j) Bubbling ChCodeEach(j) Month Day Patch Plot Rep];    
    
% figuring out if the chamber is sinking methane
    
ini_T_CH4 = mean(T_CH4(1:5));       %sets the mean of the first 5 points
end_T_CH4 = mean(T_CH4(end-5:end)); %sets mean of last 5 points

% if our first 5 is greater than last 5, we have a real sink, it runs for the entirety of the sample    
 if ini_T_CH4 > end_T_CH4 
%         
    Time_HM = T_TIME(1:end); %[s]
    T_CH4_HM = T_CH4_NoBubbles(1:end);% [ppm]
%     
%if methane is not sinking but rises to a max later in the sample, it only uses data up to this max,
% this is to prevent artificial sinking that occurs after the max
else                       
    Time_HM = T_TIME(1:Ind); %[s]
    T_CH4_HM = T_CH4_NoBubbles(1:Ind);% [ppm]
 end
%%  Nonlinear fit: hutch-mosier fit
%METHANE: Calculating diffusive fluxes with Hutchinson and Mosier one-dimension fit

%CH4 without outliers
T_CH4_rm=T_CH4_HM;
[~,TF] = rmoutliers(detrend(T_CH4_rm), "mean");
T_CH4_rm(TF==1)=nan;

Start_co_CH4 = nanmean(T_CH4_rm(1:5)); % Start value
Start_cs_CH4 = max(T_CH4_rm);       % Saturation value
Start_k_CH4 = 0.0035; 

x0=[Start_cs_CH4,Start_co_CH4,Start_k_CH4];
x = lsqcurvefit(@myfun2,x0,Time_HM(~isnan(T_CH4_rm)),T_CH4_rm(~isnan(T_CH4_rm)));
cs_CH4=x(1);
co_CH4=x(2);
k_CH4=x(3);%shape of curve

HM_CH4 = mean((-k_CH4.*(co_CH4 - cs_CH4).*P.*ChamberVolume)./(R.*TMP.*ChamberArea)); %changed to -K not K! 9/26/22 incorrect first derivative math HM equation in paper

%CO2 without outliers
T_CO2_rm=T_CO2;
[~,TF] = rmoutliers(detrend(T_CO2), "mean");
T_CO2_rm(TF==1)=nan;

Start_co_CO2 = mean(T_CO2_rm(1:5),'omitnan');    
Start_cs_CO2 = 0; 
Start_k_CO2 = 0; 

z0=[Start_cs_CO2,Start_co_CO2,Start_k_CO2];
z = lsqcurvefit(@myfun3,z0,T_TIME(~isnan(T_CO2_rm)),T_CO2_rm(~isnan(T_CO2_rm)));
cs_CO2=z(1);
co_CO2=z(2);
k_CO2=z(3);

HM_CO2 =mean((-k_CO2.*(co_CO2 -cs_CO2).*P.*ChamberVolume)./(R.*TMP.*ChamberArea)); 

% Slope = K*(Co-Cs); (Slope*(P*V)/(R*Temp*A)); [umol m-2 s-1]
% see Pendersen et al.2000 (https://onlinelibrary.wiley.com/doi/full/10.1111/j.1365-2389.2010.01291.x) and Kutzbach et al. 2007(https://www.biogeosciences.net/4/1005/2007/)

%% Optionally plot the chamber raw data and the fit
if plotYN==1
f1=figure();
f1.Units='Inches';
f1.Position=[2 2 4 4];

%% Useful lines to plot what you're doing with
subplot(2,1,1);
y=myfun2(x,Time_HM);
y2=myfun3(z,T_TIME);

%%RMSE without outliers
%CH4 RMSE
err= y - T_CH4_rm;
squareError= err.^2;
sumSqaureError= sum(squareError, 'omitnan'); %added ~isnan= doesn't work
n=sum(~isnan(T_CH4_rm));
RMSE_CH4=sqrt(sumSqaureError/n);
NRMSE_CH4=RMSE_CH4/mean(T_CH4_rm, 'omitnan')*100; 

%CO2 RMSE
err_CO2= y2 - T_CO2_rm;
squareError_CO2= err_CO2.^2;
sumSqaureError_CO2= sum(squareError_CO2,'omitnan'); % keeps returning NANs
n_CO2=sum(~isnan(T_CO2_rm));
RMSE_CO2=sqrt(sumSqaureError_CO2/n_CO2);
NRMSE_CO2=RMSE_CO2/mean(T_CO2_rm, 'omitnan')*100; %added omitnan

%%%%%%%%%%%%%%%%%%% PLOTS
%CH4 plots without outliers
hold on;
plot(Time_HM,y,'-r','MarkerSize',15); %this gives the curve fit with HM line (should pair with no bubbles plot)
plot(Time_HM,T_CH4_rm,'.c','MarkerSize',10); %This is the plot with no bubbles 
errorbar(Time_HM,T_CH4_rm,RMSE_CH4*ones(length(T_CH4_rm),1));
title(Site);
ylabel('CH_4 (ppm)');
legend ({'H-M Fit','Diffusive Flux'},... %
     'FontSize',18,'TextColor','black' );
 
%Orginal CH4 plot
% hold on;
% subplot(2,1,2);
% %plot(T_TIME(2:end),TotBubble,'g*','MarkerSize',10); %Total Bubbles
% plot(TIME,CH4,'.b','MarkerSize',10); %This is the diffusive plot with bubbles 
% title(Site);
% ylabel('CH_4 (ppm)');
% errorbar(Time_HM,T_CH4_HM,RMSE_CH4*ones(length(T_CH4_HM),1));
% legend ({'Diffusive Flux'},... %
%     'FontSize',12,'TextColor','black' );

%CO2 without outliers
subplot(2,1,2);
hold on;
plot(T_TIME,y2,'-r','MarkerSize',15); %CO2 HM line
plot(T_TIME,T_CO2_rm,'.c','MarkerSize',10);
ylabel('CO_2 (\mumol m^{-3})');
xlabel('Time (s)');
errorbar(T_TIME,T_CO2_rm,RMSE_CO2*ones(length(T_TIME),1));
legend ({'H-M Fit','No Outliers',},... 
    'FontSize',18,'TextColor','black' );

end
end
end
function F =myfun2(x,xdata)
% F=cs+(co-cs)*exp(-k*x)
F=x(1)+(x(2)-x(1)).*exp(-x(3).*xdata);
end
function F2 =myfun3(z,xdata)
 % F=cs+(co-cs)*exp(-k*x)
 F2=z(1)+(z(2)-z(1)).*exp(-z(3).*xdata);
end


%%%%Old code
%Time_HM = T_TIME(1:Ind); %[s]= not needed anymore because of above
%condition %previous line 143
%T_CH4_HM = T_CH4_NoBubbles(1:Ind);% [ppm] = not needed anymore because of above condition

% CH4
% Start_co_CH4 = mean(T_CH4_HM(1:5)); % Start value
% Start_cs_CH4 = max(T_CH4_HM);       % Saturation value
% Start_k_CH4 = 0.0035;               % arbitrary value, chosen from looking at previopus samples, changing it can help/hurt the fit
%                                     % this is what has worked best for Jorge in the past, you can change, he says it usually works well with this valu
% x0=[Start_cs_CH4,Start_co_CH4,Start_k_CH4];
% x = lsqcurvefit(@myfun2,x0,Time_HM,T_CH4_HM);
% cs_CH4=x(1);
% co_CH4=x(2);
% k_CH4=x(3);%shape of curve

%CARBON DIOXIDE
% Start_co_CO2 = mean(T_CO2(1:5));    
% Start_cs_CO2 = 0; 
% Start_k_CO2 = 0; 
% 
% z0=[Start_cs_CO2,Start_co_CO2,Start_k_CO2];
% z = lsqcurvefit(@myfun3,z0,T_TIME(~isnan(T_CO2)),T_CO2(~isnan(T_CO2)));
% cs_CO2=z(1);
% co_CO2=z(2);
% k_CO2=z(3);

%y_orig=myfun2([Start_cs_CH4,Start_co_CH4,Start_k_CH4],Time_HM);
%x = lsqcurvefit(@myfun2,x0,Time_HM,T_CH4_HM);
% 
%CH4 RMSE
% err= y - T_CH4_HM;
% squareError= err.^2;
% sumSqaureError= sum(squareError);
% n=sum(~isnan(T_CH4_HM));
% RMSE_CH4=sqrt(sumSqaureError/n);
% NRMSE_CH4=RMSE_CH4/mean(T_CH4_HM)*100; 

% CO2 RMSE
% err_CO2= y2 - T_CO2;
% squareError_CO2= err_CO2.^2;
% sumSqaureError_CO2= sum(squareError_CO2);
% n_CO2=sum(~isnan(T_CO2));
% RMSE_CO2=sqrt(sumSqaureError_CO2/n_CO2);
% NRMSE_CO2=RMSE_CO2/mean(T_CO2)*100;
% 
% %HM Fit CH4====this doesn't take out outliers
% hold on;
% plot(Time_HM,y,'-r','MarkerSize',15); %this gives the curve fit with HM line (should pair with no bubbles plot)
% plot(Time_HM,T_CH4_HM,'.c','MarkerSize',10); %This is the plot with no bubbles 
% errorbar(Time_HM,T_CH4_HM,RMSE_CH4*ones(length(T_CH4_HM),1));
% title(Site);
% ylabel('CH_4 (ppm)');
% legend ({'H-M Fit','Diffusive Flux'},... %
%     'FontSize',12,'TextColor','blue' );

% % CH4 original plot
% subplot(2,1,2);
% hold on;
% plot(T_TIME(2:end),TotBubble,'g*','MarkerSize',10); %Total Bubbles
% plot(TIME,CH4,'.b','MarkerSize',10); %This is the diffusive plot with bubbles 
% title(Site);
% ylabel('CH_4 (ppm)');
% %errorbar(Time_HM,T_CH4_HM,RMSE_CH4*ones(length(T_CH4_HM),1));
% legend ({'Diffusive Flux'},... %
%     'FontSize',12,'TextColor','blue' );

%CO2 plots
% %original CO2
% subplot(2,1,2);
% hold on;
% plot(T_TIME,y2,'-r','MarkerSize',15); %CO2 HM line
% plot(T_TIME,T_CO2,'.c','MarkerSize',10);
% ylabel('CO_2 (\mumol m^{-3})');
% xlabel('Time (s)');
% errorbar(T_TIME,T_CO2,RMSE_CO2*ones(length(T_CO2),1));
% legend ({'H-M Fit','CO2',},... 
%     'FontSize',12,'TextColor','blue' );