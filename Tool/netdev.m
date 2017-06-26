function [history,delta,c] = netdev(adj,parameter,varargin)

r = parameter(:,1);
g = parameter(:,2);
e = parameter(:,3);

A = adj';
N = length(A);
step = 10000;
Ct     = 0;
driver = [];
initialstate = rand(N,1);
initialdelta = zeros(N,1);

L = 1 - 1/(r-e) ;
R = 1 - 1/(r+e) ;
out = sum(A,1)' ;   in = sum(A,2) ;
in_zero    = (in == 0) ;
in_nonzero = (in ~= 0) ;

j = 1;
while j <= nargin-2
    switch ischar(varargin{j})
        case 0
            error('No define classification of property for %dth value\n',j)
        case 1
            switch varargin{j}
                case 'L'
                    L = varargin{j+1} ;
                    j = j+2;
                case 'R'
                    R = varargin{j+1} ;
                    j = j+2;
                case 'threshold'
                    threshold = varargin{j+1} ;
                    j = j+2;
                case 'step'
                    step = varargin{j+1};
                    j = j+2;
                case 'initialstate'
                    initialstate = varargin{j+1};
                    j = j+2;
                case 'initialdelta'
                    initialdelta = varargin{j+1};
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
                            switch driver_strategy
                                case 'rand'
                                    numCN = varargin{j+2};
                                    [~,Index] = sort(rand(N,1),'descend');
                                    driver    = Index(1:numCN);
                                    j = j + 3 ;
                                case 'dynamic'
                                    evalf  = varargin{j+2} ;
                                    numCN  = varargin{j+3} ;
                                    j = j + 4 ;
                            end
                    end
                otherwise
                    error('" %s " is not defined.\n',varargin{j})
            end
    end
end
%===============================================================
x  = initialstate ;
dx = initialdelta ;
ph = sign(dx) ;
c  = false(N,step) ;
delta   = zeros(N,step) ;
history = zeros(N,step) ;
%========== Evolution ==========================
t = 1 ;
distance = dist(A,ceil(0.5*(N-N^0.5))) ; %%%
for iii = 1:step
    %========== control ========================
    if t <= length(Ct) && iii == Ct(t)
        if strcmp(driver_strategy,'dynamic')
            switch evalf
                case 'laplacian'
                    D = abs((A*x)./in - x) ;
                case 'unknown'
                    D = abs((A*(x.*ph))./in - x) ;
                case 'tent'
                    D = (x<L).*x + (x>R).*(L/(1-R)*(1-x)) ;
%                     D = D.*(distance <= ceil(1/10*(iii-Ct(1)))) ; %%%
                case '+0-'
%                     threshold = 0.026 ;
                    D = abs(dx) ;
                    D = (1-D).*(D>threshold) ;
            end
            [D_sort,Index] = sort(D,'descend') ;
            ND = sum(D_sort>0) ;
            switch ND > numCN
                case 0
                    driver = Index(1:ND) ;
                case 1
                    driver = Index(1:numCN) ;
            end
        end
        c(driver,iii) = true ;
        t = t + 1 ;
    end
    
    fx = (r+c(:,iii).*ph.*e).*x.*(1-x);
    history(:,iii) = (1-g.*in_nonzero).*fx + g./(in+in_zero).*(A*fx) ;
    delta(:,iii)   = history(:,iii) - x ;
    
    x  = history(:,iii) ;
    dx = delta(:,iii) ;
    ph = sign(dx) ;
end

end