sessionDate = '04_27_2018';
session = '1';
directory = ['\\130.49.229.252\gharbawie\Lab\Gilligan\All Data\Gilligan_',sessionDate, '\'];
load([directory,'Phasespace\Results\Gilligan_', sessionDate, '_', session, '.mat']);
data = sortedPSData.PhasespaceData.SmallSphere.allTrialData;
f1 = gca(figure(1));
hold on;
f2 = gca(figure(2));
hold on;
for t = 1:size(data,1)
    data = sortedPSData.PhasespaceData.SmallSphere.allTrialData{t,4};
    lift = sortedPSData.PhasespaceData.SmallSphere.allTrialData{t,5};
    frameRange = sortedPSData.PhasespaceData.SmallSphere.trialDataForEachLED.FrameRange;
    frameRange = [frameRange(:,3) frameRange(:,4)-1];
    LEDs = sortedPSData.LEDs;
    dp = 0;
    LED = [];
    for d = 1:size(LEDs,1)
        ind = data.LED == LEDs(d);
        indG = lift.LED == LEDs(d);
        if(any(ind))
            dp = dp+1;
            frames = data.FrameNumber(ind);
            droppedFrames = ~ismember(frameRange(t,1):frameRange(t,2),frames)';
            LED(droppedFrames,1,dp) = NaN;
            LED(droppedFrames,2,dp) = NaN;
            if(any(~droppedFrames))
                LED(~droppedFrames,1,dp) = data.x(ind);
                LED(~droppedFrames,2,dp) = data.z(ind);
            end
            LEDG(dp,1) = lift.x(find(indG,1));
            LEDG(dp,2) = lift.z(find(indG,1));
        end
    end
    combs = nchoosek(1:size(LED,3), 2);
    gangle = [];
    angles = [];
    for c = 1:size(combs,1)
        angles(:,c) = atan2(LED(:,1,combs(c,1))-LED(:,1,combs(c,2)),...
            LED(:,2,combs(c,1))-LED(:,2,combs(c,2)));
        gangle(:,c) = atan2(LEDG(combs(c,1),1)-LEDG(combs(c,2),1),...
            LEDG(combs(c,1),2)-LEDG(combs(c,2),2));
        angles(angles(:,c)<0,c) = angles(angles(:,c)<0,c) + 2*pi;
        gangle(gangle(:,c)<0,c) = gangle(gangle(:,c)<0,c) + 2*pi;
        angles(:,c) = 180/pi * angles(:,c);
        gangle(:,c) = 180/pi * gangle(:,c);
    end
    gangle = repmat(gangle,size(angles,1),1);
    difference = angles - gangle;
    difference = abs(sum(difference,2));
    [val, ind] = sort(difference);
    histVals(t,:) = ind(1:10)/size(difference,1);
    %plot(abs(sum(difference,2)));
    plot(f1,(1:size(difference,1))/size(difference,1),difference/max(difference));
end
histogram(f2, reshape(histVals,t*10,1),30);