function [y,p] = updatestate_LTISIM2(p,sysName,u,opt)
%% [y,p] = updatestate_LTISIM2(p,sysName,u,opt)

% if ~exist('opt','var')
%     opt.storeState = 1;
%     opt.storeOutput = 1;
% end

sysStruct = p.(sysName);
x = sysStruct.x;

if strcmpi(sysStruct.type, 'ss')
    xNext = sysStruct.a * x(:,end) + sysStruct.b * u;
    y = sysStruct.c * x(:,end) + sysStruct.d * u;

elseif strcmpi(sysStruct.type, 'iir')
    denDeg = sysStruct.denDeg;
    Unew = [u; x(denDeg+1:end-1, end)]; % [u(k) u(k-1) ... u(k-denDeg)]
    Y = x(1:denDeg,end); % [-y(k-1) -y(k-2) ... -y(k-denDeg)]
    y = [Y; Unew]' * sysStruct.w;
    Ynew = [-y; Y(1:denDeg-1)];
    xNext = [Ynew ; Unew];

elseif strcmpi(sysStruct.type, 'fir')
    y = x(:,end)' * sysStruct.w;
    xNext = [u; x(1:end-1, end) ];
end
    


if opt.storeState ~= 0
    sysStruct.x = [x xNext];
else
    sysStruct.x = xNext;
end

if opt.storeOutput ~= 0
    sysStruct.y = [sysStruct.y y];
else
    sysStruct.y = y;
end

p.(sysName) = sysStruct;

