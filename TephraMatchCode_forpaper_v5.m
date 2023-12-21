load tephramatch_basic
% Begin assuming the following variables are loaded into matlab:
% AgeFamCore1(2): Age Model iterations, each column should be a different age
% output sorted by depth
% depthsCore1(2): the depths that correspond to each row in the variable
% AgeFamCore 1(2)

t = 25; %Set t to be the maximum age difference between the correlated tephra (or any correlated depths)
Core1tephradepth = [108, 131, 192.5, 202]; %Create a variable that defines the correlated depths in Core1
Core2tephradepth = [13, 36.5, 178.5, 200]; %Create a variable that defines the correlated depths in Core2
%% Matching age model iterations, runtimes can be long
goodij = [0 0 0];
h = waitbar(0, 'Please Wait...');
for i = 1:size(AgeFamCore1,2)
    for j = 1:size(AgeFamCore2,2)
        
        Core1tephraage = interp1(depthsCore1, AgeFamCore1(:,i), Core1tephradepth);
        Core2tephraage = interp1(depthsCore2, AgeFamCore2(:,j), Core2tephradepth);
        tephraagediff = abs(Core1tephraage-Core2tephraage);
        if any(tephraagediff>t)
            
        else
            goodij=[goodij ; i j max(tephraagediff)];
        end
    end
    waitbar(i/size(AgeFamCore1,2),h)
end
delete(h)

MatchGood=goodij(2:end,:);
Core1unique=unique(MatchGood(:,1));
Core2unique=unique(MatchGood(:,2));
save tephramatch.mat
%% Exporting Ages
depthCore1_out=0:0.5:300; %Set the depth resolution and length of the core you want to export
depthCore2_out=0:0.5:450; %Set the depth resolution and length of the core you want to export

for i=1:size(AgeFamCore1,2)
    AgeFamCore1_5mm(:,i)= interp1(depthsCore1, AgeFamCore1(:,i), depthCore1_out);
end
for i=1:size(AgeFamCore2,2)
    AgeFamCore2_5mm(:,i)= interp1(depthsCore2, AgeFamCore2(:,i), depthCore2_out);
end

AgeFamCore1_Out = [depthCore1_out' quantile(AgeFamCore1_5mm,[.5 .025 .975],2)];
AgeFamCore2_Out = [depthCore2_out' quantile(AgeFamCore2_5mm,[.5 .025 .975],2)];

for i=1:size(Core1unique)
    Core1_Out(:,i)= interp1(depthsCore1, AgeFamCore1(:,Core1unique(i)), depthCore1_out);
end

Core1_OutSorted=sort(Core1_Out, 2);

Core1_medAge=Core1_OutSorted(:,round(0.5*size(Core1unique,1)));
Core1_cl025=Core1_OutSorted(:,round(0.025*size(Core1unique,1)));
Core1_cl975=Core1_OutSorted(:,round(0.975*size(Core1unique,1)));

depthCore1_out2=transpose(depthCore1_out);
Core1Age_Out=[depthCore1_out2, Core1_medAge, Core1_cl025, Core1_cl975];
csvwrite('Core1Ages.csv', Core1Age_Out)

for i=1:size(Core2unique)
    Core2_Out(:,i)= interp1(depthsCore2, AgeFamCore2(:,Core2unique(i)), depthCore2_out);
end

Core2_OutSorted=sort(Core2_Out, 2);

Core2_medAge=Core2_OutSorted(:,round(0.5*size(Core2unique,1)));
Core2_cl025=Core2_OutSorted(:,round(0.025*size(Core2unique,1)));
Core2_cl975=Core2_OutSorted(:,round(0.975*size(Core2unique,1)));

depthCore2_out2=transpose(depthCore2_out);
Core2Age_Out=[depthCore2_out2, Core2_medAge, Core2_cl025, Core2_cl975];
csvwrite('Core2Ages.csv', Core2Age_Out)
save AgesOut.mat Core1Age_Out Core2Age_Out

%plot age models
figure(1)
clf
subplot(221);
h11 = plot(AgeFamCore1_Out(:,1),AgeFamCore1_Out(:,2:4),'k')
hold on
h12 = plot(Core1Age_Out(:,1),Core1Age_Out(:,2:4),'r')
xlabel('Depth (cm)')
ylabel('Age (BP)')
title('Core1 Models')
legend([h11(1) h12(1)],'Original Model','Tephra limited model','location','NorthWest')

subplot(222);
h21 = plot(AgeFamCore2_Out(:,1),AgeFamCore2_Out(:,2:4),'k')
hold on
h22 = plot(Core2Age_Out(:,1),Core2Age_Out(:,2:4),'r')
xlabel('Depth (cm)')
ylabel('Age (BP)')
title('Core2 Models')
legend([h21(1) h22(1)],'Original Model','Tephra limited model','location','NorthWest')


subplot(223);
plot(AgeFamCore1_Out(:,1),AgeFamCore1_Out(:,2:4)- Core1Age_Out(:,2:4))
xlabel('Depth (cm)')
ylabel('Age (BP) Difference')
title('Core1 Model change')
legend('Median','0.975 quantile','0.025 quantile','location','SouthWest')

subplot(224);
plot(AgeFamCore2_Out(:,1),AgeFamCore2_Out(:,2:4)- Core2Age_Out(:,2:4))
xlabel('Depth (cm)')
ylabel('Age (BP) Difference')
title('Core2 Model change')
legend('Median','0.975 quantile','0.025 quantile','location','NorthEast')

%% Calculating tephra age estimates using all iterations from both cores that meet the matching criteria
% "Core1tephradepth(1)*2+1" represents the row number
% corresponding to the depth of the first tephra, assuming you are using a
% depth resolution of 0.5. If you are using a depth resolution other than 0.5, you will need to edit
% the code below to ensure that the ages are extracted from the correct row
% in Core1(2)_Out
Core1T1=Core1_Out(Core1tephradepth(1)*2+1,:);
Core2T1=Core2_Out(Core2tephradepth(1)*2+1,:);
T1=horzcat(Core1T1, Core2T1);

Core1T2=Core1_Out(Core1tephradepth(2)*2+1,:);
Core2T2=Core2_Out(Core2tephradepth(2)*2+1,:);
T2=horzcat(Core1T2, Core2T2);

Core1T3=Core1_Out(Core1tephradepth(3)*2+1,:);
Core2T3=Core2_Out(Core2tephradepth(3)*2+1,:);
T3=horzcat(Core1T3, Core2T3);

Core1T4=Core1_Out(Core1tephradepth(4)*2+1,:);
Core2T4=Core2_Out(Core2tephradepth(4)*2+1,:);
T4=horzcat(Core1T4, Core2T4);

% This example uses four correlated tephra, you may need to add or remove
% the above code segments depending on how many tiepoints you are using

Tephra_all_ages=vertcat(T1, T2, T3, T4);
Tephra_ages_sorted=sort(Tephra_all_ages, 2);
Tephra_medAge=Tephra_ages_sorted(:,round(0.5*length(Tephra_ages_sorted)));
Tephra_cl025=Tephra_ages_sorted(:,round(0.025*length(Tephra_ages_sorted)));
Tephra_cl975=Tephra_ages_sorted(:,round(0.975*length(Tephra_ages_sorted)));
Tephra_out=[Tephra_medAge,Tephra_cl025,Tephra_cl975];
csvwrite('Tephra_out.csv', Tephra_out)

% This will produce a spreasheet in which each row includes a median age estimate
% for each of the correlated depths, as well as the 95% confidenceintervals