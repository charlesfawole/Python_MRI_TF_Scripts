function [cN,TF,TF_result] = TFscale(p_TF,cN,mN,mode,TF)
switch p_TF
    case 'T'
        switch mode
            case 'abs'
                r = sqrt(sum(cN.*mN)/sum(cN.^2)); %%% scaling coefficient
            case 'rel'
                r = sqrt(sum(cN./mN)/sum((cN./mN).^2));
            case 'max'
                [~,idx] = max(mN);
                r = sqrt(mN(idx)/cN(idx));
            case 'mm' % minimize maximum error
                cN_1 = min(cN./mN);
                cN_2 = max(cN./mN);
                r = sqrt(2/(cN_1+cN_2));
            otherwise
                [~,idx] = max(mN);
                r = sqrt(mN(idx)/cN(idx));
        end
        TF.y = TF.y.*r;
        cN = r^2*cN;
    case 'V'
        switch mode
            case 'abs'
                r = sum(cN.*mN)/sum(cN.^2); %%% scaling coefficient
            case 'rel'
                r = sum(cN./mN)/sum((cN./mN).^2);
            case 'max'
                [~,idx] = max(mN);
                r = mN(idx)/cN(idx);
            case 'mm' % minimize maximum error
                cN_1 = min(cN./mN);
                cN_2 = max(cN./mN);
                r = 2/(cN_1+cN_2);
            otherwise
                [~,idx] = max(mN);
                r = mN(idx)/cN(idx);
        end
        TF.y = TF.y.*r;
        cN = r*cN;
end

TF.y = interp1(TF.x,TF.y,(0.01:0.01:0.43),'linear','extrap');
TF.x = 0.01:0.01:0.43;
TF.mag = abs(TF.y);
TF.phase = angle(TF.y);
TF.phase = TF.phase./pi.*180;
result_file = 'E:\human model\TF\test.xlsx';
TF_result = [TF.mag.',TF.phase.'];
xlswrite(result_file,TF_result,strcat('A1',':','B',num2str(length(TF.phase))));

figure
subplot(2,1,1)
plot(TF.x*100,TF.mag,'b')
title('scaled transfer function magnitude')
subplot(2,1,2)
plot(TF.x*100,TF.phase,'b');
title('scaled transfer function phase')
end