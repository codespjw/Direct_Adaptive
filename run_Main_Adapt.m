clear
rng(0)


period = 348;
uffwdPreloaded.ma = zeros(1,period);
uffwdPreloaded.vcm = zeros(1,period);
%%

for actuators = 1:2
    if ~exist('allHarmonics','var')
        allHarmonics = 1:58;
        actuatorType = 'vcm';
        
    else
        allHarmonics = 59:173;
        actuatorType = 'ma';
    end

    clear region
    region{1} = allHarmonics(1:floor(end/2));
    region{2} = allHarmonics(floor(end/2)+1:end);
    
    t = 1;
    i=1;
    while i <= numel(region)
        harm = region{i};
        fprintf('Current region: %i - %i\n',harm(1), harm(end));
        out = MAIN_Adapt(harm,uffwdPreloaded,actuatorType);
        if (out.terminated == 1)
            fprintf('Problematic harmonic: %i\n',harm(out.maxidx))
            region{i} = harm(1:ceil(end/2));
            if i < numel(region)
                region{i+1} = [harm(floor(end/2)+1:end) region{i+1}];
            else
                region{i+1} = harm(floor(end/2)+1:end);
            end
        else
            for j = 1:3, beep; pause(0.5); end
            uffwd{i} = out.uffwd;
            uffwdPreloaded.(actuatorType) = uffwd{i};
            i = i+1;
        end
        for j = 1:numel(region)
            region{j}([1,end])
        end
        col = true;
        img = ones(1,period/2-1);
        for j = 1:numel(region)
            idx = region{j};
            img(1,idx) = col;
            col = ~col;
        end
        IMG.(actuatorType)(t,:) = img;
        t = t+1;
    end
end

plotPartitions
