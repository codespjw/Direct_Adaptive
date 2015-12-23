function out = NoiseModel(plotFlag)
%% Creates a transfer function such that its response to white noise
%  matches NRRO in frequency domain.

%% Parameters
na = 13;
nc = na;
period = 348;
tpi = 400e3; 
tpNano = 25.4/tpi*1e6;
%% Load PES data
g = load('pes_hgst_raw');
pes = g.pes(1:end/5);
tmp = rmvRro(pes,period);
ynrro = tmp.nrro;
rro = tmp.rro;
N = length(ynrro);

if plotFlag
    h = figurename('PES and NRRO Spectrum');
    plotfft(pes*tpNano,1/period,period*10,h,'r');
    plotfft(ynrro*tpNano,1/period,period*10,h);
    legend('Raw PES','Non-Repeatable Runout (NRRO)')
    box on
    ylabel('Amplitude Spectrum (nm)')
    xlabel('Harmonic');
    saveImgPdf(6,3,[mfilename,'_spectrum_rawPes_nrro']);
end
%% Model noise by an ARMA model: y = C/A e where e is white noise

% Filter the PES
yf = idfilt(ynrro,[1/(period/2) 1]*pi,1);

% Fit ARMA model
Opt = armaxOptions; 
if plotFlag
    Opt.Display = 'full';
end
model = armax(yf,[na nc], Opt);
[a,~,c]=polydata(model);
noiseModel = tf(c,a,-1);

% Adjust gain
noiseModel = noiseModel/norm(noiseModel)*std(yf);


%% Compare spectrum of NRRO with noise model response
ypredict = lsim(noiseModel,randn(N,1));

if plotFlag
    h = figurename('Actual and modeled NRRO spectrums');
    plotfft(ynrro*tpNano,1/period,period*10,h);
    plotfft(ypredict*tpNano,1/period,period*10,h,'r');
    legend('Actual spectrum','Modeled spectrum')
    box on
    ylabel('Amplitude Spectrum (nm)')
    xlabel('Harmonic');
    saveImgPdf(6,3,[mfilename,'_spectrum_nrro_actual_model']);
end

%% Output argument
out.noiseModel = noiseModel;
out.rro = rro;
