function [history,phase,driver] = netdev(A,parameter,varargin)

r = parameter(1);
g = parameter(2);
e = parameter(3);

N    = length(A);
step = 10000;
Ct     = 0;
driver = [];
initialstate = rand(N,1);
initialphase = zeros(N,1);

out = sum(A,1)';   in = sum(A,2);   tot = out-in ;

j = 1;
while j <= nargin-2
    switch ischar(varargin{j})
        case 0
            error('No define classification of property for %dth value\n',j)
        case 1
            switch varargin{j}
                case 'step'
                    step = varargin{j+1};
                    j = j+2;
                case 'initialstate'
                    initialstate = varargin{j+1};
                    j = j+2;
                case 'initialphase'
                    initialphase = varargin{j+1};
                    j = j+2;
                case 'controlT'
                    Ct = sort(varargin{j+1}) ;
                    j = j+2;
                case 'driver'
                    driver_strategy = varargin{j+1};
                    switch ischar(driver_strategy)
                        case 0
                            driver = driver_strategy;
                            j = j+2;
                        case 1
                            controlrate = varargin{j+2};
                            switch driver_strategy
                                case 'rand'
                                    [s,Index]   = sort(rand(N,1),'descend');
                                    driver      = sort(Index(1:ceil(controlrate*N)));
                                case 'OID'
                                    [O_I,Index] = sort(out-in,'descend');
                                    driver      = sort(Index(1:ceil(controlrate*N)));
                                case 'OD'
                                    [OD,Index]  = sort(tot,'descend');
                                    driver      = sort(Index(1:ceil(controlrate*N)));
                                case 'water'
                                    water = out ;
                                    for i = 1 : 10000
                                        water = A*(water./out) ;
                                    end
                                    [W,Index] = sort(water,'descend') ;
                                    driver    = sort(Index(1:ceil(controlrate*N))) ;
                            end
                            j = j+3;
                    end
            end
    end
end
%===============================================================
w = A.*g;
M = (w-diag(sum(w,2)))./(in*ones(1,N));
M(isnan(sum(M,2)),:) = A(isnan(sum(M,2)),:);
M = M+eye(N);
Cphase  = initialphase;
space   = initialstate;
phase   = zeros(N,step);
history = zeros(N,step);
%========== Evolution ==========================
t = 1 ;
for iii = 1:step
    fx = r*space.*(1-space);
    %========== control ========================
    if t <= length(Ct) && iii == Ct(t)
        c = Cphase;
        t = t+1;
    else
        c = zeros(N,1) ;
    end
    
    fxc = (r+c*e).*space.*(1-space);
    fx(driver) = fxc(driver);
    
    history(:,iii) = M*fx;
    phase(:,iii)   = sign(history(:,iii)-space) ;
    space  = history(:,iii) ;
    Cphase = phase(:,iii) ;
end

end