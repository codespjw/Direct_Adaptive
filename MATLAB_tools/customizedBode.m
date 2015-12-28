function FigHandle = customizedBode(g,bodeOp,onlyMag,FigHandle)
%%
if nargin < 3
    onlyMag = 1;
end
if nargin < 4
    yw = 500;
    xw = 300;
    set(0,'units','pixels');  
    %Obtains this pixel information
    Pix_SS = get(0,'screensize');
    FigHandle = figure('Position', [Pix_SS(3)/4, Pix_SS(4)/3, yw, xw]);
end
%%
[magout,phase,w] = bode(g, bodeOp);
w = w/2/pi;
if ~onlyMag
    subplot(211);
end
plot(w,20*log10(squeeze(magout)),'linewidth',2);
hold on
xlim(bodeOp.XLim{1})
if ~strcmp(bodeOp.YLimMode,'auto')
    ylim(bodeOp.YLim{1})
end
line(1/g.Ts/2*[1 1],ylim,'color','k','linewidth',2)
ylabel('Magnitude (dB)')
grid on; 
if onlyMag
    xlabel('Harmonic')
    return;
end

subplot(212)
phase = mod(phase,360);
phase(phase>180) = phase(phase>180)-360;
plot(w,squeeze(phase),'linewidth',2);
hold on
set(gca,'YTick',[-180:90:180])
grid on; 
xlim(bodeOp.XLim{1})
ylim([-180 180])
line(1/g.Ts/2*[1 1],ylim,'color','k','linewidth',2)
xlabel('Harmonic')
ylabel('Phase (deg)')