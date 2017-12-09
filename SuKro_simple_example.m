addpath('./misc/','./DEMO_Image_Denoising/toolbox/ompbox/')

N = 64;     % Data dimension
K = 2000;   % Number of training samples
M = N*4;    % Number of atoms

% Data to approximate
Y = randn(N,K); % random data
Y = Y./repmat(sqrt(sum(Y.^2,1)),size(Y,1),1); % normalize data columns (optional)

% Initial dictionary
D = randn(N,M); D = D./repmat(sqrt(sum(D.^2,1)),size(D,1),1); % random
%D = odctndict([sqrt(N) sqrt(N)],M,2); % 2D-ODCT

%% Parameters for SuKro
params.iternum = 20;        % Number of iterations
params.memusage = 'high';   % Memory usage
params.initdict = D;        % Initial dictionary
params.data = Y;            % Data to approximate
params.u = 500;             % ADMM coefficient
params.alpha = 2.15;        % Regularization parameter
% Sparse coding parameters
params.codemode = 'sparsity';   % 'error' or 'sparsity'
params.Edata = 1e-1;            % when using codemode = 'error'
params.Tdata = ceil(size(D,1)/10);     % when using codemode = 'sparsity'
% Dimensions of the subdictionaries
params.kro_dims.N1 = sqrt(N); params.kro_dims.N2 = sqrt(N);
params.kro_dims.M1 = sqrt(M); params.kro_dims.M2 = sqrt(M);

% Checking subdictionary dimensions
assert( params.kro_dims.N1*params.kro_dims.N2==N && ...
        params.kro_dims.M1*params.kro_dims.M2==M, ...
        'N (resp. M) should be equal to N1*N2 (resp. M1*M2)')
assert( round(params.kro_dims.N1)==params.kro_dims.N1 && ...
        round(params.kro_dims.N2)==params.kro_dims.N2 && ...
        round(params.kro_dims.M1)==params.kro_dims.M1 && ...
        round(params.kro_dims.M2)==params.kro_dims.M2,   ...
        'N1,N2,M1 and M2 should all be integers')

%% Running SuKro algorithm
flops_dense = 2*size(Y,1)*size(Y,2);

[D, D_not_normalized, X] = sum_separable_dict_learn(params);

%% Results
% Verifying the number of separable terms on the learned dictionary
D_reord = reord(D_not_normalized,params.kro_dims.N1,params.kro_dims.N2,params.kro_dims.M1,params.kro_dims.M2);
fprintf('Learned dictionary has %d separable terms.\n',rank(D_reord,norm(D_reord,'fro')*2e-7));