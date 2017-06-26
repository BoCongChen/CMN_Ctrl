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
        p   = varargin{1} ;
        adj = floor(triu(rand(N)+p,1)) ;
        adj = adj + adj';
    case 'BA_model'
        adj = false(N) ;
        adj([3 N+3 2*N+1 2*N+2]) = true ;
        v = 3 ;
        while v < N
            d = sum(adj,2) ;
            p = d/v ;
            for link = 1 : 2
                i = 0 ;
                number = rand*sum(p) ;
                while number > 0
                    i = i + 1 ;
                    number = number - p(i) ;
                end
                adj(i,v+1) = true ; adj(v+1,i) = true ;
                p(i) = 0 ;
            end
            v = v + 1 ;
        end
    case 'radiation'
        adj = false(N) ;
        adj(1,2:end) = true ;
        adj(2:end,1) = true ;
    case '1-dir'
        adj = circshift(eye(N),[1 0]) ;
        adj(1,end) = 0 ;
    case '1-dir_ring'
        adj = circshift(eye(N),[1 0]) ;
    case 'square'
        x = N(1) ;
        y = N(2) ;
        adj = false(x*y) ;
        B = reshape(1:x*y,[x,y])' ;
        C = [B B B ; B B B ; B B B] ;
        for n = 1 : x*y
            i = ceil(n/x) + y ;
            j = rem(n,x)  + x ;
            v = [C(i+1,j) C(i,j+1) C(i-1,j) C(i,j-1)] ;
            adj(v,n) = true ;
        end
end

switch length(adj) > 70
    case 0
        adj = logical(adj) ;
    case 1
        adj = sparse(adj) ;
end

end