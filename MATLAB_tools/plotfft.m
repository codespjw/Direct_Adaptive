function [y fn w] = plotfft(x,Ts,nfft,fn,color,linewidth)
%% plotfft(x, 1/period, figureNumber)
%
% figureNumber: the figure number on which the results should be plotted
%               0 if plots are not required
%               [] for generating a new figure

if ~exist('Ts','var')
    Ts = 1;
end

if exist('color','var')==0, color = 'b'; end

Fs = 1/Ts;

if ~exist('nfft','var')
    nfft = Fs*1;
end
if nargin < 6
    linewidth = 1;
end
n = length(x);
xFft = fft(x,nfft)/min(n,nfft)*2;
n = length(xFft);
y = xFft(1:floor(n/2));
y = abs(y);
w = linspace(0,Fs/2*(1-1/length(y)),length(y));
if ~exist('color','var')
    fn = figure;
elseif fn == 0
    return;
else
    figure(fn)
end

% subplot(211)
hold all;
plot(w(:),y(:),color,'linewidth',linewidth)
grid on

id = find( y>0);
xlim([0 Fs/2])
% subplot(212)
% hold all
% plot(id/length(y)*Fs/2,mag2db(y(id)) )
% grid on

% figure
% semilogx(linspace(0,Fs/2,length(y)),mag2db(abs(y)) )

