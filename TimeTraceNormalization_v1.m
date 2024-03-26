function [nn, base] = TimeTraceNormalization_v1(data, type, NormStartTime)
% NormStartTime은 Fl type에만 적용됨

switch type
    case 'abs'
        for j = 1 : data.size.freq
            for i = 1 : data.size.iter
                for k = 1 : data.size.rep
                    % base(k, i, j) = mean(data.abs.raw(data.baselinerange,k,i,j));
                    base(k, i, j) = mean(data.abs.raw(NormStartTime*1e3/data.dt:19*1e3/data.dt,k,i,j));
                    nn(:,k,i,j) = 1-(data.abs.raw(:,k,i,j)/base(k,i,j));
                end
            end
        end

    case 'pfm'
        for j = 1 : data.size.freq
            for i = 1 : data.size.iter
                for k = 1 : data.size.rep
                    % base(k, i, j) = mean(data.abs.pfm(data.baselinerange,k,i,j));
                    base(k, i, j) = mean(data.abs.pfm(NormStartTime*1e3/data.dt:19*1e3/data.dt,k,i,j));
                    nn(:,k,i,j) = 1-(data.abs.pfm(:,k,i,j)/base(k,i,j));
                end
            end
        end

    case 'fl'
        for j = 1 : data.size.freq
            for i = 1 : data.size.iter
                for k = 1 : data.size.rep
                    base(k, i, j) = mean(data.fl.raw(NormStartTime*1e3/data.dt:end,k,i,j));
                    nn(:,k, i,j) = (data.fl.raw(:,k, i,j)-base(k,i,j)); % time trace 2D image를 그려보고 8 ms으로 정햇음...
                end
            end
        end
    otherwise
        error('try abs,pfm, or fl')
end

end