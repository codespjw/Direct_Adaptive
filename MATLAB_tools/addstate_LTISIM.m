function p = addstate_LTISIM(p,sys,x0, type)
%% p = addstate_LTISIM(p,sys,x0, type)
if nargin < 4
    type = 'ss';
end

if ( (nargin <=2) && strcmp(type,'ss') )
    if isempty(sys.A)%Static gain
        x0 = 0;
    else
        x0 = zeros(size(sys.A,1),1);
    end
end
    
if strcmpi(type, 'ss')
    if ( (sys.InputDelay + sys.OutputDelay  ~= 0) || ( ~isempty(sys.InternalDelay) && (sys.InternalDelay ~=0) ) )
        error('The input system should have zero Input/Output/Internal delay. If the system has delays, they should be modeled by adding delayed states in the state vector.');
    end
    if isempty(sys.a), sys.a = 0; end
    if isempty(sys.b), sys.b = 0; end
    if isempty(sys.c), sys.c = 0; end
    if isempty(sys.d), sys.d = 0; end
    p.(inputname(2)) = struct('x',x0,'y',[],'a',sys.a,'b',sys.b,'c',sys.c,'d',sys.d,'Ts',sys.Ts, 'type', type);
    
elseif strcmpi(type, 'iir')
    
    [b a] = tfdata(sys, 'v');
    denDeg = length(pole(sys));
    numDeg = length(zero(sys));
    relDegree = denDeg - numDeg;
%     b = b(2:end); %Strictly Causal systems
    a = a(2:end); %remove 1 from a
    w = [a b];
    if isempty(x0)
        x0 = zeros(length(w),1);
    end
    Ts = sys.Ts;
    
    p.(inputname(2)) = struct('x', x0, 'y',[],'w',[a b]', 'numDeg', numDeg, 'denDeg', denDeg,'Ts',Ts, 'type',type,'W',[]);
elseif strcmpi(type, 'FIR')
    if ( ~isvector(sys) )
        error( 'Second argument, in2, should be a vector such that in2(i) gives the coefficient for u(k-i)');
    end
    sys = reshape(sys,[],1);
    
    p.(inputname(2)) = struct('x', x0, 'y',[],'w',sys, 'type',type,'W',[]);
end