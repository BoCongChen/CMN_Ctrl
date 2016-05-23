function [history,driver] = netdev(A,parameter,varargin)

r = parameter(1);
g = parameter(2);
e = parameter(3);

N    = length(A);
preT = 1000;
step = 10000;
driver_strategy = [];
initialstate = rand(N,1);

% if nargin > 2
%     input = cell(nargin-2,1);
%     for i = 1:nargin-2
%         eval(['input(',num2str(i),') = {v',num2str(i),'};'])
%     end
% end

j = 1;
while j <= nargin-2
    switch ischar(varargin{j})
        case 0
            error('No define classification of property for %dth value\n',j)
        case 1
            switch varargin{j}
                case 'preT'
                    preT = varargin{j+1};
                    j = j+2;
                case 'step'
                    step = varargin{j+1};
                    j = j+2;
                case 'initialstate'
                    initialstate = varargin{j+1};
                    j = j+2;
                case 'driver'
                    driver_strategy = varargin{j+1};
                    switch ischar(driver_strategy)
                        case 0
                            j = j+2;
                        case 1
                            controlrate = varargin{j+2};
                            j = j+3;
                    end
            end
    end
end

%========== Decide control node ================================
out = sum(A,1)';   in = sum(A,2);   tot = out-in;
switch ischar(driver_strategy)
    case 0
        driver = driver_strategy;
    case 1
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
end
%===============================================================
w = A.*g;
M = (w-diag(sum(w,2)))./(in*ones(1,N));
M(isnan(sum(M,2)),:) = A(isnan(sum(M,2)),:);
M = M+eye(N);
space   = initialstate;
history = zeros(N,step);
%========== Evolution ==========================
for iii = 1:preT
    fx    = r*space.*(1-space);
    space = M*fx;
    history(:,iii) = space;
end
%==================================================================================================
for iii=iii+1:step
    fx = r*space.*(1-space);
    %========== control ========================
    c  = history(:,iii-1)-history(:,iii-2);
    cc = abs(c) > 0;
    c  = sign(cc.*c);
    
    fxc = (r+c*e).*space.*(1-space);
    fx(driver) = fxc(driver);
    
    space = M*fx;
    history(:,iii) = space;
end

end