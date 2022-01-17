function x_dn = FFT_derivative_3d(x,delta,dt,dim)
% x     : M-N-P matrix
% delta : the cell size
% dt    : horizontal translation
% dim   : dimension derivative is to be taken along
% x_dn  : output derivative matrix


%check input matrix has 3 dimensions
total_dim = size(size(x),2);
if total_dim ~= 3
    error('Input Matrix x does not have 3 dimensions');
end

%permute matrix x so dimension to have derivative calculated is the first
%dimension
% permute_arg = 1:3;
% inrange = permute_arg == dim;
% permute_arg(inrange)=1;
% permute_arg(1)=dim;
% x = permute(x,permute_arg);

%pad with zeros
len = size(x,dim);
if mod(size(x,dim),2)==1
    if dim==1
      xt = zeros(len+1,size(x,2),size(x,3));
      xt(1:len,:,:)=x;
      x = xt;
    elseif dim==2
      xt = zeros(size(x,1),len+1,size(x,3));
      xt(:,1:len,:)=x;
      x = xt;   
    elseif dim==3
      xt = zeros(size(x,1),size(x,2),len+1);
      xt(:,:,1:len)=x;
      x = xt;    
    end    
end

%calculate terms for derivative calculation
n = size(x,dim)
X = fft(x,[],dim);
P = delta*n;

k_array = [(0:(n/2)),(n/2+1:n-1)-n];
G_array = (exp(1i*k_array*2*pi()*dt/P)).*(1i*k_array*2*pi()/P);

save G_array G_array X;
%dt/delta 
%G_array(n/2+1)*delta
%calculate derivative and translate by dt storing in X_d
X_d1 = zeros(size(X));
if dim==1
if mod(n,2) == 0
% $$$     for k=0:(n/2)
% $$$         X_d1(k+1,:,:) = X(k+1,:,:)*exp(1i*k*2*pi()*dt/P)*1i*k*2*pi()/P;
% $$$     end
% $$$     for k=(n/2+1):(n-1)
% $$$         X_d1(k+1,:,:) = -X(k+1,:,:)*exp(-1i*(n-k)*2*pi()*dt/P)*1i*(n-k)*2*pi()/P;
% $$$     end
         X_d = repmat(G_array',[1 size(X,2) size(X,3)]).*X;
    else
    error('n is not even')
    end
%         if ~isequal(X_d1,X_d)
%         error('dim1')
%         end
elseif dim==2
    if mod(n,2) == 0
% $$$     for k=0:(n/2)
% $$$         X_d1(:,k+1,:) = X(:,k+1,:)*exp(1i*k*2*pi()*dt/P)*1i*k*2*pi()/P;
% $$$     end
% $$$     for k=(n/2+1):(n-1)
% $$$         X_d1(:,k+1,:) = -X(:,k+1,:)*exp(-1i*(n-k)*2*pi()*dt/P)*1i*(n-k)*2*pi()/P;
% $$$     end
     X_d = repmat(G_array,[size(X,1) 1 size(X,3)]).*X;
    else
    error('n is not even')
    end
%     if ~isequal(X_d1,X_d)
%         error('dim2')
%     end
elseif dim==3
    if mod(n,2) == 0
% $$$     for k=0:(n/2)
% $$$         X_d(:,:,k+1) = X(:,:,k+1)*exp(1i*k*2*pi()*dt/P)*1i*k*2*pi()/P;
% $$$     end
% $$$     for k=(n/2+1):(n-1)
% $$$         X_d(:,:,k+1) = -X(:,:,k+1)*exp(-1i*(n-k)*2*pi()*dt/P)*1i*(n-k)*2*pi()/P;
% $$$     end
    X_d = repmat(permute(G_array,[1 3 2]),[size(X,1) size(X,2) 1]).*X;
    else
    error('n is not even')
    end
%     if ~isequal(X_d1,X_d)
%         error('dim3')
%     end
end

%convert out of frequency domain and remove excess row if present
x_dn2 = real(ifft(X_d,[],dim));
if dim==1
x_dn = x_dn2(1:len,:,:);
elseif dim==2
x_dn = x_dn2(:,1:len,:);
elseif dim==3
x_dn = x_dn2(:,:,1:len);    
end 
