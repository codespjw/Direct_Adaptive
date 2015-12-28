% Novel directive adaptive control for attenuating periodic disturbances
% clear all;
function [out] = MAIN_Adapt(harmCancel,uffwdPreloaded,actuatorType)
close all;
rng(0)
%%
period = 348;
compOn = 1;
m = 40e-4;
sampleTime = 1;
Nt = 10e4;
nrroOn = 1;
saturateLimit = 1.25;
% harmCancel = 1:1:173;
rroHarmonics = 1:173;
delFHz = 0.001;
igain = 4e-5;
ugain = 25;
na = 17;
nah = 3;
nbh = 3;
% actuatorType = 'ma';
plotFlag = 1;
maxValidationFail = 6;

% Lowpass filter coeff on estimated phase and amplitude
pargain = 5e-4;
normTrhLpCoef = 2e-3;

if ~exist('uffwdPreloaded','var')
    uffwdPreloaded = zeros(period,1);
end

%% Set Frequencies
nharm = length(harmCancel);
tR = zeros(2*nharm,1);

wrad = 2*pi * harmCancel'/period;
delWrad = 2*pi*delFHz*sampleTime;
%%
if strcmpi(actuatorType, 'ma')
    g = load('ma_1_173');
    sys = tf(g.sys);
    
    g = load('vcm_1_173');
    outtmp = lsim(g.vcmModel,repmat(uffwdPreloaded.vcm(:),50,1));
    otherActuatorOutput = outtmp(end-period+1:end,1);
    uffwdPreloaded = uffwdPreloaded.ma;
    
elseif strcmpi(actuatorType, 'vcm')
    g = load('vcm_1_173');
    sys = tf(g.vcmModel);
    
    g = load('ma_1_173');
    outtmp = lsim(g.sys,repmat(uffwdPreloaded.ma(:),50,1));
    otherActuatorOutput = outtmp(end-period+1:end,1);
    uffwdPreloaded = uffwdPreloaded.vcm;
end


sys.Ts = sampleTime;
sys = balred(sys,na, 'Elimination', 'Truncate');
% figurename('Sys Bode');
% bode(sys);
% iopzmap(sys)

%% Project poles to the inside of unit disk
eps = 1e-2;
projRadius = 1 - eps;
[z, p, k] = zpkdata(sys,'v');
idx = abs(p)>projRadius;
p(idx) = p(idx)./abs(p(idx))*projRadius;
sys = zpk(z,p,k,sys.Ts);

%%

[num,den]=tfdata(sys,'v');
As = den(2:end);
relDegree = find(num~=0,1) - 1;
B = num(relDegree+1:end);


%%
tA = -As';
tB = B';
tC = [1,1]';
tE = tC(2:end,1);

na = length(tA);
pA = zeros(na,1);
pAh = zeros(nah,1);

nb = length(tB);
pB = zeros(nb,1);
pBh = zeros(nbh,1);

pD = zeros(nb,1);

nr = length(tR);
pR = zeros(nr,1);

ne = length(tE);
pE = zeros(ne,1);

tAh = zeros(nah,1);
tBh = zeros(nbh,1);
tRh = tR*0;
tEh = tE*0;

%%

ma = m*ones(Nt,1)*10e2*1; 
mb = m*ones(Nt,1)*10e2*1; 
mr = m*20e-0;
me = m;

%% Initialize Parameters
TA = zeros(size(tAh,1),Nt);
TB = zeros(size(tBh,1),Nt);
TR = zeros(size(tRh,1),Nt);
U = zeros(Nt,1);
Utot = zeros(Nt,1);
yf = 0;
Y = zeros(Nt,1);
E = zeros(Nt,1);
b1 = zeros(Nt,1);
tRd = tRh*0;
magAbs = 1;
phsRad = 0;
normTrh = 0.2;
normTR = zeros(Nt,1);
e = 0;
cnt1 = 0;

%%
zinv = exp(-1i*wrad);
zvec = zeros(nharm,nbh);
for i = 1:nbh
    zvec(:,i) = zinv.^i;
end
%% Noise model
g = load('noiseModel');
tmp = g.tmp;
% tmp = NoiseModel(0);
rro = tmp.rro;
noiseModel = tmp.noiseModel;
w = lsim(noiseModel, randn(Nt,1))*logical(nrroOn);
% save noiseModel tmp


%% RRO: Separate RRO contents
x = repmat(rro,20,1);
wtmp = (1:period/2-1)/period*pi*2;
[y,mag,phs] = fftfreq(x,wtmp);
rroFreq = zeros(period/2-1,period);
for i = 1:length(mag)
    rroFreq(i,:) = mag(i)*sin((0:period-1)*wtmp(i)+phs(i));
end
% figure; plot(rro-sum(rroFreq)')
rro = sum(rroFreq(rroHarmonics,:))';
% plotfft(repmat(rro,10,1),1/period/10,1);
% plotfft(repmat(rroCancel,10,1),1/period/10,1,'r');


%% Decreasing gain
a = 2e3/3;
b = a;
%     sc = 1./(1+exp(a-b*[1:Nt]/Nt));
sc = a./(b+(1:Nt));
id = sc < 1e-2;
sc(id) = sc(id)*0;
%     figure; plot(sc)
    

%% Filter for input signal
tmpgain = 1;
filtOrder = 10;
bpFilt = designfilt('bandpassiir','FilterOrder',filtOrder, ...
         'HalfPowerFrequency1',max([1,harmCancel(1)-2]),'HalfPowerFrequency2',min([harmCancel(end)+2, period/2-1]), ...
         'SampleRate',period);
utmp = randn(Nt,1);
urand = filter(bpFilt,utmp);
urand = urand/std(urand)*tmpgain;
[a,b,c,d]=ss(bpFilt);
bpss = ss(a,b,c,d,-1);
[bf,df] = tfdata(bpss,'v');
af = df(2:end);
% fvtool(bpFilt)

% Response of filter
[filtMagAbs,filtPhsDeg] = bode(bpss,wrad); 
filtPhsRad = mod(deg2rad(filtPhsDeg(:)),2*pi);
filtMagAbs = filtMagAbs(:);

% rro = filter(bpFilt, repmat(rro,10,1));
% rro = rro(end-period+1:end);
puf = zeros(length(bf),1);
pyf = zeros(length(af),1);

puf2 = zeros(length(bf),1);
pyf2 = zeros(length(af),1);

a = struct();
ufil = bpss;
yfil = bpss;
a= addstate_LTISIM(a,ufil,zeros(filtOrder,1),'ss');
a= addstate_LTISIM(a,yfil,zeros(filtOrder,1),'ss');
opt.storeState = 0;
opt.storeOutput = 0;
    
%% Validation for early termination
maxTrhNorm = 0;
validationCounter = 0;
out.terminated = 0;
%%
%  tAh = [-1.1287    0.2609    0.6425]';
%  tBh = [-0.0419   -0.0253   -0.0114]';
%     scalar = 1/(1+exp(a-b*k));
%%
tic
for k = 1:Nt
    krot = mod(k,period)+1;
    
    % ------ Generate Y ------
    y = tA'*pA + tB'*pB + tB'*pD ;
    nrro = w(k);
    r = rro(krot);
    meas = y+nrro+r+ otherActuatorOutput(krot);
    if meas > 2, break; end
    
    
    % ------ Update parameters ------
    yh = tAh'*pAh + tBh'*pBh + tRh' * pR ;
    if isnan(yh), break; end
    [measf,a]=updatestate_LTISIM2(a,'yfil',meas,opt);
    e = measf - yh;
    ma(k) = ma(k)*sc(k);
    tmp = pAh'*pAh + pBh'*pBh;
    tAh = tAh + ma(k)/(1+tmp) * pAh * e;
    if ~mod(k,3)
        tAh = fstab([1;-tAh],1);
        tAh = -tAh(2:end,:);
    end
    mb(k) = mb(k)*sc(k);
    tBh = tBh + mb(k)/(1+tmp) * pBh * e;
    at = 1e4;
    mrt = mr*((at)/(at+k));
    tRh = tRh + mrt/(1+pR'*pR) * pR * e;
    
    % ------ Invert B_hat(q^-1) ------
    % tmp = tf(tBh',[1 zeros(1,na)],sampleTime); [magAbs,phsDeg] = bode(tmp,wrad); phsRad = deg2rad(phsDeg);
    resp = zvec*tBh;
    magAbs = (1-pargain)*magAbs + pargain*abs(resp);
    phsRad = (1-pargain)*phsRad + pargain*angle(resp);
    phsTmp = wrad*(k-1)-phsRad-filtPhsRad;
    magTmp = magAbs.*filtMagAbs;
    sD = (-1./[magTmp; magTmp]).*[sin(phsTmp); cos(phsTmp)];
    
    
    % ------ Compensation signal ------
    dold = d;
    dcGain = inf;
    tRd = (1-igain/dcGain)*tRd + igain * tRh;
    d = tRd'*sD;
    if isnan(d)
        d = dold;
    end
    if compOn == false
        d = 0;
    end
    d = d + uffwdPreloaded(krot);


    % ------ Symmetric sinusoidal signal for sys-ID ------
    pR = [sin(wrad*k); cos(wrad*k)];
    at = 1e4;
    gain = ugain*(at/(at+k));
    normTrh = (1-normTrhLpCoef)*normTrh+normTrhLpCoef*norm(tRh)*gain;
    gain = normTrh;%(normTrh > 0.002) * normTrh;
    udel = cos(delWrad*k)*sum(pR(1:end/2))*2* gain;


    % ------ Total excitation signal for sys-ID ------
    u = udel;
%     u = udel;% + randn(1)/urandGain;
%     u = urand(k);

    % ------ Saturate total injection signal ------
    u = min(u, saturateLimit-d);
    u = max(u, -saturateLimit-d);


    [uf,a]=updatestate_LTISIM2(a,'ufil',u,opt);

    % ------ Update feature vectors ------
    pA = [y;pA(1:end-1,:)];
    pAh = [measf;pAh(1:end-1,:)];
    pB = [u; pB(1:end-1,:)];
    pBh = [uf; pBh(1:end-1,:)];
    pD = [d; pD(1:end-1,:)];
    
    %--------- Early Termination ---------
    if ~mod(k,1000)
        maxTrhNormOld = maxTrhNorm;
        [maxTrhNorm, maxidx] = max(abs(tRh(1:end/2)+1i*tRh(end/2+1:end)));
        if maxTrhNorm > maxTrhNormOld
            validationCounter = validationCounter + 1;
            if validationCounter > maxValidationFail 
                out.terminated = 1;
                out.maxidx = maxidx;
                break;
            end
        else
            validationCounter = 0;
        end
    end
    
    %--------- Log variables --------
    b1(k) = r;
    normTR(k) = normTrh;
    Y(k)=meas;
    E(k) = e;
    U(k) = u;
    Utot(k) = u+d;
    %D(k)=d;
    uffwd(krot) = d;
    if ~mod(k,100)
        cnt1 = cnt1+1;
        TR(:,cnt1) = tRh;
        TA(:,cnt1) = tAh;
        TB(:,cnt1) = tBh;
    end
    
    if ~mod(k,ceil(Nt/10))
        fprintf('step: %i   failed validation: %i\n',k, validationCounter);
    end
            
end
toc

%% Remove zero tail of data
TR = TR(:,1:cnt1);
TA = TA(:,1:cnt1);
TB = TB(:,1:cnt1);

%% Visualize results
tpi = 400e3; 
tpNano = 25.4/tpi*1e6;
mkdir('Figures');
tscale = 1;

% --------- Position Error (nm) -------------
figurename('PES');
plot([1:length(Y)]*tscale,Y*tpNano); 
ylabel('Position Error (nm)','Interpreter','latex');
xlabel('Step','Interpreter','latex');
grid on;
saveImgPdf(6,3,['./Figures/pes-compensated-' num2str(harmCancel(1)),'-',num2str(harmCancel(end)) ]);

% --------- Estimated A and B coefficients -------------
figurename('A,B');
plot(linspace(1,Nt*tscale,size(TA,2)),[-TA' TB'],'linewidth',2);
hold on; grid on;
for i = 1:size(TA,1)
    legend_str{i} = ['a_',num2str(i)];
end
for i = 1:size(TB,1)
    legend_str{i+size(TA,1)} = ['b_',num2str(i)];
end
legend(legend_str)
xlim([1,Nt*tscale])
ylim([min([-TA(:); TB(:)])-0.1,max(([-TA(:); TB(:)]))+0.1])
xlabel('Step','Interpreter','latex'); 
ylabel('Estimated A and B coefficients','Interpreter','latex');
saveImgPdf(6,3,['./Figures/AB-' num2str(harmCancel(1)),'-',num2str(harmCancel(end)) ]);

% --------- R (\hat{\theta}_M) coefficients -------------
figurename('R');
plot(linspace(1,Nt*tscale,size(TR,2)),TR','linewidth',1);
hold on
grid on
ylabel('$\hat{\theta}_M$ coefficients','Interpreter','latex');
xlabel('Step','Interpreter','latex');
xlim([1,Nt*tscale]);
tmp1 = min(TR(:)); tmp1 = tmp1 * (1-sign(tmp1)*0.1);
tmp2 = max(TR(:)); tmp2 = tmp2 * (1+sign(tmp2)*0.1);
ylim([tmp1 tmp2])
saveImgPdf(6,3,['./Figures/thetaM-' num2str(harmCancel(1)),'-',num2str(harmCancel(end)) ]);

% --------- Estimation Error -------------
time = linspace(1,Nt*tscale,size(E,1));
x = E*tpNano;
figurename('E');
t1 = time(1:end);
x1 = x(1:end);
subplot 121
plot(t1,x1,'linewidth',1);
grid on
ylabel('$\epsilon$: Estimation error (nm)','Interpreter','latex');
xlabel('Step','Interpreter','latex');
xlim([t1(1),t1(end)]);
tmp1 = min(x1); tmp1 = tmp1 * (1-sign(tmp1)*0.1);
tmp2 = max(x1); tmp2 = tmp2 * (1+sign(tmp2)*0.1);
tmp3 = max(-tmp1,tmp2);
ylim(tmp3*[-1 1])

tnmp = 0.3;
t1 = time(end*tnmp:end);
x1 = x(end*tnmp:end);
subplot 122
plot(t1,x1,'linewidth',1);
grid on
% ylabel('$\epsilon$: Estimation error (nm)','Interpreter','latex');
xlabel('Step','Interpreter','latex');
xlim([t1(1),t1(end)]);
tmp1 = min(x1); tmp1 = tmp1 * (1-sign(tmp1)*0.1);
tmp2 = max(x1); tmp2 = tmp2 * (1+sign(tmp2)*0.1);
tmp3 = max(-tmp1,tmp2);
% ylim(tmp3*[-1 1])
saveImgPdf(6,3,['./Figures/epsilon-' num2str(harmCancel(1)),'-',num2str(harmCancel(end)) ]);


% ---------- PES spectrum all frequencies ------------
g = load('pes_hgst_raw');
pes = g.pes(1:end/5);
h = figurename('PES and NRRO Spectrum');
plotfft(pes*tpNano,1/period,period*10,h,'g',1.5);
plotfft(Y(end*0.6:end)*tpNano,1/period,period*10,h,'r',1.5);
legend('Adaptive Controller Off', 'Adaptive Controller On', 'Orientation','Horizontal')
box on
ylabel('Amplitude Spectrum (nm)','Interpreter','latex')
xlabel('Harmonic','Interpreter','latex');
saveImgPdf(6,3,['./Figures/spectrum-allfreq-' num2str(harmCancel(1)),'-',num2str(harmCancel(end)) ]);

% ---------- Spectrum only at harmonics ------------
nrev = floor(Nt/period);
n = w(1:nrev*period);
r = repmat(rro, nrev,1);
noi = n + r;
[y, mag, phs] = fftfreq(pes*tpNano, [1:period/2-1]/period*2*pi);
[y2, mag2, phs2] = fftfreq(Y(end*0.6:end)*tpNano, [1:period/2-1]/period*2*pi);
h = figurename('Spectrum - comp and uncomp - bars');
bar([1:period/2-1],mag,'g','EdgeColor','g');
hold on
bar([1:period/2-1],mag2,'r','EdgeColor','r');
ylabel('Amplitude Spectrum (nm)','Interpreter','latex');
xlabel('Harmonic','Interpreter','latex');
grid on; xlim([1,period/2-1]);
legend('Adaptive Controller Off', 'Adaptive Controller On', 'Orientation','Horizontal')
saveImgPdf(6,3,['./Figures/spectrum-bars-',num2str(harmCancel(1)),'-',num2str(harmCancel(end))]);

% ---------- Excitation signal gain ------------
figurename('Norm of tRh');
semilogy(linspace(1,Nt*tscale,size(normTR,1)),normTR,'linewidth',2); 
xlabel('Step','Interpreter','latex')
ylabel('$\alpha_u$: Excitation gain','Interpreter','latex')
grid on;
saveImgPdf(6,3,['./Figures/alphaU-' num2str(harmCancel(1)),'-',num2str(harmCancel(end)) ]);


% ---------- All injenction signals ------------
figurename('Injection Signals');
plot([Utot U]); grid on
legend('Total: excitation + control','Excitation')
xlabel('Step','Interpreter','latex')
ylabel('Injection signals','Interpreter','latex')
saveImgPdf(6,3,['./Figures/injection-' num2str(harmCancel(1)),'-',num2str(harmCancel(end)) ]);


Ainv = tf(1,[1;-tA]',sampleTime);
yuncomp = lsim(sys,randn(Nt,1));
yuncomp = yuncomp + lsim(Ainv,b1);
yuncomp = yuncomp(end/2:end);

ystd = std(Y(end*0.7:end));
fprintf('std(Y) = %f\n',ystd);
ystdMin = norm(sys);
fprintf('min std(Y) = %f\n',ystdMin);
yuncompStd = std(yuncomp);
fprintf('std(Y_uncomp) = %f\n',yuncompStd);


%% Bode options
Ts = 1/period;
sys.Ts = Ts; 
bodeOp = bodeoptions();             bodeOp.FreqScale = 'linear';
bodeOp.PhaseWrapping = 'on';        bodeOp.FreqUnits = 'Hz';
bodeOp.XLim = [0.01 1/Ts/2]*1.1;      bodeOp.PhaseMatching = 'on';
bodeOp.PhaseMatchingFreq = (1/Ts)*0.3;   bodeOp.grid = 'on';

sysh = tf([0 tBh' zeros(1,nah-nbh)],[1; -tAh]',Ts);
% -- FIG: Bodes -- 
h = figurename('Bode of estimated system'); 
customizedBode(sys,bodeOp,0,h);
customizedBode(sysh,bodeOp,0,h);
for i = 1:2
    subplot(2,1,i)
    yl = ylim; 
    ylim(yl);
    xpatch = [harmCancel(1),harmCancel(end),harmCancel(end),harmCancel(1)];
    ypatch = [yl(1) yl(1) yl(2) yl(2)];
    p=patch(xpatch,ypatch,'r','LineStyle','none');
    set(p,'FaceAlpha',0.1);
end
saveImgPdf(6,3,['./Figures/bode-',num2str(harmCancel(1)),'-',num2str(harmCancel(end))]);

% -- FIG: Diff of bodes -- 
h = figurename('Diff in Bode of estimated system'); 
customizedBode(sys/sysh,bodeOp,0,h);
for i = 1:2
    subplot(2,1,i)
    yl = ylim; 
    ylim(yl);
    xpatch = [harmCancel(1),harmCancel(end),harmCancel(end),harmCancel(1)];
    ypatch = [yl(1) yl(1) yl(2) yl(2)];
    p=patch(xpatch,ypatch,'r','LineStyle','none');
    set(p,'FaceAlpha',0.1);
end
saveImgPdf(6,3,['./Figures/bode-diff-',num2str(harmCancel(1)),'-',num2str(harmCancel(end))]);


fprintf('Stability: %i\n',isstable(sysh));


%% Output argument
out.uffwd = uffwd;
% keyboard