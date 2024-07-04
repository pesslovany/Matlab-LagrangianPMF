function evalRes(nx, k, annes_PMFout, annes_PFout, rmsePMF, rmsePF, ...
    astdPMF11, astdPMF22, astdPF11, astdPF22, tocPMFavg, tocPFavg, ...
    x, measMean2, measVar2, xEst, varEst)

    % Normalize outputs
    annes_PMFout = annes_PMFout / (nx * k);
    annes_PFout = annes_PFout / (nx * k);
    rmsePMFout = mean(rmsePMF, 2);
    rmsePFout = mean(rmsePF, 2);
    tocPMFavgOut = mean(tocPMFavg, 2);
    tocPFavgOut = mean(tocPFavg, 2);

    % Create table T2
    T2 = table([mean(rmsePMFout(1:2)), mean(rmsePFout(1:2))]', ...
              [mean([astdPMF11 astdPMF22]), mean([astdPF11 astdPF22])]', ...
              [mean(annes_PMFout), mean(annes_PFout)]', ...
              [mean(tocPMFavgOut), mean(tocPFavgOut)]', ...
              'VariableNames', {'RMSE POS','ASTD POS','ANNES','TIME'}, ...
              'RowName', {'PMF','PF bootstrap'})

    % Plot figures
    endTime_act = k - 1;
    figure

    if nx == 3
        plotTitles = {'pos-x', 'pos-y', 'bias'};
        plotUnits = {'pos-x [m]', 'pos-y [m]', 'bias [m]'};
        sp = [3,1]; % subplot indices
    elseif nx == 4
        plotTitles = {'pos-x', 'pos-y', 'vel-x', 'vel-y'};
        plotUnits = {'pos-x [m]', 'pos-y [m]', 'vel-x [m/s]', 'vel-y [m/s]'};
        sp = [2,2]; % subplot indices
    end

    for idx = 1:nx
        subplot(sp(1), sp(2), idx)
        plot(x(idx, 1:endTime_act) - measMean2(idx, 1:endTime_act), 'b', 'LineWidth', 2)
        hold on
        plot(sqrt(squeeze(measVar2(idx, idx, 1:endTime_act))), 'b--', 'LineWidth', 2)
        plot(-sqrt(squeeze(measVar2(idx, idx, 1:endTime_act))), 'b--', 'LineWidth', 2)
        plot(x(idx, 1:endTime_act) - xEst(idx, 1:endTime_act), 'r', 'LineWidth', 2)
        plot(sqrt(varEst(idx, 1:endTime_act)), 'r--', 'LineWidth', 2)
        plot(-sqrt(varEst(idx, 1:endTime_act)), 'r--', 'LineWidth', 2)
        xlabel('time [s]')
        ylabel(plotUnits{idx})
        title(plotTitles{idx})
        if idx == 1
            legend('PMF err', 'PMF +std', 'PMF -std', 'PF err', 'PF +std', 'PF -std')
        end
    end
end
