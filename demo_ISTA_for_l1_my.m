%% Reference
% https://people.rennes.inria.fr/Cedric.Herzet/Cedric.Herzet/Sparse_Seminar/Entrees/2012/11/12_A_Fast_Iterative_Shrinkage-Thresholding_Algorithmfor_Linear_Inverse_Problems_(A._Beck,_M._Teboulle)_files/Breck_2009.pdf

%% COST FUNCTION
% x^* = argmin_x { 1/2 * || A(X) - Y ||_2^2 + lambda * || X ||_1 }
%
% x^k+1 = threshold(x^k - 1/L*AT(A(x^k)) - Y), lambda/L)

%%
clear ;
close all;
home;

%% GPU Processing
% If there is GPU device on your board, 
% then isgpu is true. Otherwise, it is false.
bgpu    = false;
bfig    = true;

%%  SYSTEM SETTING
is_real = 1;
A       = @(x, mask) ifdct_wrapping(x, is_real).*mask;
A2       = @(x) ifdct_wrapping(x, is_real);
AT      = @(y, mask) fdct_wrapping(y.*mask, is_real);
AINV    = @(y) fdct_wrapping(y, is_real);

%% DATA GENERATION
data0 = load('real_801_150_expand.mat');
data0 = squeeze(data0.data);
mask = load('real_801_150_regular_missing_0.5_mask.mat');
mask = squeeze(mask.data);
[h, w] = size(mask);
data0_miss = load('real_801_150_regular_missing_0.5.mat');
data0_miss = squeeze(data0_miss.data);

pad_size = 0;
st=1;ed=w-0;
w0 = ed-st+1;
data0 = padarray(data0(:, st:ed), [0,pad_size], 'symmetric');
data0_miss = padarray(data0_miss(:, st:ed), [0,pad_size], 'symmetric');
mask = padarray(mask(:, st:ed), [0,pad_size], 'symmetric');

x       = fdct_wrapping(data0, is_real);
p       = data0_miss;
x_full  = AINV(p);

%% LOW-DOSE SINOGRAM GENERATION

x_low   = AINV(p);

%% NEWTON METHOD INITIALIZATION
LAMBDA  = [0.001, 0.0001, 0.00001];
T       = [1e0, 1e0, 1e0];

y       = p;
x0      = fdct_wrapping(zeros(size(data0)), is_real);
niter   = 202;

L1              = @(x) norm_my_l1(x);
L2              = @(x) power(norm(x, 'fro'), 2);
COST.equation   = '1/2 * || A(X) - Y ||_2^2 + lambda * || X ||_1';
COST.function	= @(x, mask, lambda) 1/2 * L2(A(x, mask) - y) + lambda * L1(x);
COST.function1	= @(x, mask, lambda) 1/2 * L2(A(x, mask) - y);
COST.function2	= @(x, mask, lambda) lambda * L1(x);

%% RUN NEWTON METHOD
if bgpu
    y  = gpuArray(y);
    x0 = gpuArray(x0);
end

 [x_ista, obj]	= ISTA_my(A, A2, AT, x0, y, mask, LAMBDA, 1./T, niter, COST, bfig, data0);

 %%
 y_est = A2(x_ista);
 
 figure, imagesc(y_est(:,pad_size+1:pad_size+w0));
 y_est = y_est(:,pad_size+1:pad_size+w0);
 save curvelet2.mat y_est
 
img = fdct_wrapping_dispcoef(fdct_wrapping(data0_miss, is_real));
hh=figure; imagesc(abs(img)); title('input'); saveas(hh,'input-curvelet.png')
img = fdct_wrapping_dispcoef(fdct_wrapping(data0, is_real));
hh=figure; imagesc(abs(img)); title('target'); saveas(hh,'target-curvelet.png')
img = fdct_wrapping_dispcoef(fdct_wrapping(y_est, is_real));
hh=figure; imagesc(abs(img)); title('output'); saveas(hh,'input-curvelet.png')
