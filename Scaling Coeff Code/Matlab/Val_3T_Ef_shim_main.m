function [] = Val_3T_Ef_shim_main()
y_grid = [2];
r_pwr = [0.5];
corr_all = zeros(length(y_grid),length(r_pwr));
std_err = corr_all;
max_err = corr_all;
for ii=1:1:length(y_grid)
    for jj = 1:1:length(r_pwr)
        [corr, err_r] = Val_3T_Ef_shim('3T-III-empty-quad',y_grid(ii),r_pwr(jj),'l','on');
        corr_all(ii,jj) = corr;
        std_err(ii,jj) = std(err_r);
        max_err(ii,jj) = max(abs(err_r));
    end
end
disp(corr_all);
disp(std_err);
disp(max_err);