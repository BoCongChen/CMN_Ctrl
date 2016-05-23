function adj = netgenerate(N,property,varargin)

switch property
    case 'ring'
        ee  = eye(N) ;
        adj = cat(1,ee(2:N,:),ee(1,:))+cat(1,ee(N,:),ee(1:N-1,:)) ;
    case 'line'
        ee  = eye(N) ;
        adj = cat(1,ee(2:N,:),zeros(1,N))+cat(1,zeros(1,N),ee(1:N-1,:)) ;
    case 'rand'
        p   = varargin{1} ;
        adj = floor(rand(N)+p) ;
        adj = adj - diag(diag(adj)) ;
    case 'rand_2'
        p   =  varargin{1};
        adj = floor(triu(rand(N)+p,1)) ;
        adj = adj + adj';
    case '1-dir'
        adj = circshift(eye(N),[1 0]) ;
        adj(1,end) = 0 ;
    case '1-dir_ring'
        adj = circshift(eye(N),[1 0]) ;
    case 'square'
        x = N(1) ;
        y = N(2) ;
        adj = zeros(x*y) ;
        B = reshape(1:x*y,[x,y]) ;
        C = [B B B ; B B B ; B B B] ;
        for n = 1 : x*y
            i = rem(n,x)  + x ;
            j = ceil(n/x) + y ;
            v = [C(i+1,j) C(i,j+1) C(i-1,j) C(i,j-1)] ;
            adj(v,n) = 1 ;
        end
end

end