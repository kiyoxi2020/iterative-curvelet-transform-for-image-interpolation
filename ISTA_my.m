%% Reference
% https://people.rennes.inria.fr/Cedric.Herzet/Cedric.Herzet/Sparse_Seminar/Entrees/2012/11/12_A_Fast_Iterative_Shrinkage-Thresholding_Algorithmfor_Linear_Inverse_Problems_(A._Beck,_M._Teboulle)_files/Breck_2009.pdf

%% COST FUNCTION
% x^* = argmin_x { 1/2 * || A(X) - Y ||_2^2 + lambda * || X ||_1 }
%
% x^k+1 = threshold(x^k - 1/L*AT(A(x^k)) - Y), lambda/L)

%%
function [x, obj]  = ISTA_my(A, A2, AT, x, b, mask, LAMBDA_list, L_list, n, COST, bfig, data0)

if (nargin < 9)
    bfig = false;
end

if (nargin < 8 || isempty(COST))
    COST.function	= @(x) (0);
    COST.equation	= [];
end

if (nargin < 7)
    n   = 1e2;
end

obj     = zeros(n, 1);
obj1     = zeros(n, 1);
obj2     = zeros(n, 1);
snr_list      = zeros(n, 1);
count = 1;
LAMBDA = LAMBDA_list(count);
L = L_list(count);
for i = 1:n
    if i==200
        LAMBDA = LAMBDA_list(2);
    elseif i == 201
        LAMBDA = LAMBDA_list(3);
    end
    x_new = update_my(x, 1/L, AT(A(x, mask) - b, mask));
    x       = threshold_my(x_new, LAMBDA);
%     x       = x - 1/L*AT(A(x) - b);
    obj(i)  = COST.function(x, mask, LAMBDA);
    obj1(i)  = COST.function1(x, mask, LAMBDA);
    obj2(i)  = COST.function2(x, mask, LAMBDA);
    
    if (bfig)
        figure(1); 
        y_est = A2(x);
        snr0 = snr(data0, y_est-data0);
        snr_list(i) = snr0;
        subplot(2,3,1);colormap gray; imagesc(y_est);
        title([num2str(i) ' / ' num2str(n), ', snr=', num2str(snr0)]);
        
         img = fdct_wrapping_dispcoef(x);
        subplot(2,3,2); imagesc(img);
        
        subplot(2,3,3); semilogy(obj, '*-');  title(COST.equation);  xlabel('# of iteration');   ylabel('Objective'); 
                                            xlim([1, n]);   grid on; grid minor;
                                            
        subplot(2,3,4); semilogy(obj1, '*-r');  title(COST.equation);  xlabel('# of iteration');   ylabel('Objective'); 
                                            xlim([1, n]);   grid on; grid minor;
        hold on; semilogy(obj2, '*-b');  title(COST.equation);  xlabel('# of iteration');   ylabel('Objective'); 
                                            xlim([1, n]);   grid on; grid minor;
%         legend('obj1', 'obj2');
        subplot(2,3,5); semilogy(snr_list, '*-');  title('snr');  xlabel('# of iteration');   ylabel('snr'); 
                                            xlim([1, n]);   grid on; grid minor;
        drawnow();
    end
end

x   = gather(x);
obj = gather(obj);

end
