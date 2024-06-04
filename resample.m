% input = state system
%   S_hat:      3xN
% output = state system resampled
%   S:          3xN
function S = resample(S_hat)

N = size(S_hat, 2);
resampleMethod = 1; % 0 vanilla, 1 systematic

%vanilla
if resampleMethod == 0
    
    CDF=cumsum(S_hat(3,:));
    S = zeros(3,N);
    
    for i=1:N
        index = find(CDF >= rand(), 1, 'first');
        S(:,i) = S_hat(:,index);
        S(3,i) = 1/N;
    end
    
%systematic
elseif resampleMethod == 1
    CDF=cumsum(S_hat(3,:));
    S = zeros(3,N);
    
    r = rand()*(1/N);

    for j=1:N
        index = find(CDF>=r+(j-1)/N, 1, 'first');
        S(:,j) = S_hat(:,index);
        S(3,j) = 1/N;
    end
end

end







