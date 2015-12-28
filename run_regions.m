clc;clear;close all;
projectPath = pwd;
addpath(genpath(projectPath));
period = 348;
uffwdPreloaded.ma = zeros(1,period);
uffwdPreloaded.vcm = zeros(1,period);
t = 1;
for actuators = 1:2
    if actuators == 1
        region = {1:29,30:43,44:47,48:58};
        actuatorType = 'vcm';
    else
        region = {59:73,74:97,98:116,117:120,121:122,123:145,146:173};
        actuatorType = 'ma';
    end

    i=1;
    while i <= numel(region)
        harm = region{i};
        fprintf('Current region: %i - %i\n',harm(1), harm(end));
        out = MAIN_Adapt(harm,uffwdPreloaded,actuatorType);
        uffwd{i} = out.uffwd;
        uffwdPreloaded.(actuatorType) = uffwd{i};
        i = i+1;
        U{t} = out.uffwd;
        t = t+1;
    end
end


figurename('uffwd: 1'); 
x1 = U{1};
plot([1:period], x1);
grid on
xlabel('Step','Interpreter','latex');
ylabel('Feedforward Control','Interpreter','latex');
tmp1 = min(x1); tmp1 = tmp1 * (1-sign(tmp1)*0.1);
tmp2 = max(x1); tmp2 = tmp2 * (1+sign(tmp2)*0.1);
tmp3 = max(-tmp1,tmp2);
ylim(tmp3*[-1 1])
xlim([1,period]);
saveImgPdf(6,3,['./Figures/uffwd-1-29']);

figurename('uffwd: last'); 
x1 = U{end};
plot([1:period], x1);
grid on
xlabel('Step','Interpreter','latex');
ylabel('Feedforward Control','Interpreter','latex');
tmp1 = min(x1); tmp1 = tmp1 * (1-sign(tmp1)*0.1);
tmp2 = max(x1); tmp2 = tmp2 * (1+sign(tmp2)*0.1);
tmp3 = max(-tmp1,tmp2);
ylim(tmp3*[-1 1])
xlim([1,period]);
saveImgPdf(6,3,['./Figures/uffwd-146-173']);
