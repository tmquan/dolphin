clc; clear all; close all;
addpath(genpath(pwd)); 
%% Size of x:   [y, x, k, n]
load x;

%% Declare in-place functions
soft_thresh     = @(x, th) sign(x)*max(abs(x)-th, 0);
vec             = @(x) x(:);
norm2           = @(x) abs(vec(x)'*vec(x));
norm1           = @(x) sum(abs(vec(x)));


%% Check the input vector
b = squeeze(x);
%figure; imagesc(b(:,:,1)); axis square off; colormap gray;

%% Declare the data size per dimension
H = 100;     % Number of dimension y
W = 100;     % Number of dimension x
K = 25;      % Number of filters
N = 10;      % Number of training images
M = 11;      % Filter dimension

%% Declare the size of each array
size_xx = [H, W, 1, N];     % x     :   Original image
size_xh = [H, W, 1, N];     % xhat  :   Fourier of x
size_dd = [H, W, K, 1];     % dk    :   a filter
size_dh = [H, W, K, 1];
size_ss = [H, W, K, 1];
size_sh = [H, W, K, 1];
size_zz = [H, W, K, 1];
size_zh = [H, W, K, 1];
size_tt = [H, W, K, 1];
size_th = [H, W, K, 1];


%% Construct the masking matrix C, a weird trick of shifting
C = ones(M, M);
C = padarray(C, [ceil(0.5*(H-M)),  ceil(0.5*(W-M))],  0, 'pre'); 
C = padarray(C, [floor(0.5*(H-M)), floor(0.5*(W-M))], 0, 'pos'); 
C = circshift(C, [H/2, W/2]);

% nnz(C)
% figure; surf(C);
% figure; imagesc(C); axis square off; colormap gray;
C = C>0;
C = repmat(C, [1, 1, K]);
% C = reshape(C, [H, W, K, 1]);

%% Control the random number generator, make it repeatable
rng(2016);

%% Convert to single precision for faster computation
x = single(x);

dd = rand(size_dd); 
dd = single(dd);

ss = rand(size_ss);
ss = single(ss);

zz = rand(size_zz);
zz = single(zz);

tt = rand(size_tt);
tt = single(tt);

lambda_ss = rand(size_ss);
lambda_ss = single(lambda_ss);

lambda_tt = rand(size_tt);
lambda_tt = single(lambda_tt);

%% Perform the Fourier Transform of the corresponding variables
xf = fft2(x);

df = fft2(dd);
sf = fft2(ss);

zf = fft2(zz);
tf = fft2(tt);

lambda_sf = fft2(lambda_ss);
lambda_tf = fft2(lambda_tt);

%% Constant variables go here
mu_max      = 1e+5;
mu_s        = 1e-2;
mu_t        = 1e-2;
tau         = 1.05;
tol         = 1e-3;
beta        = 1;
max_iter    = 500;


%% Solver goes here
for iter = 1:max_iter
    %% Subproblem zz
    %% Subproblem tt
    %% Subproblem dd
    %% Subproblem ss
    %% Multiplier update
    %% Penalty update
    %% Collect the cost
    %% Visualize the filter
end