%% ===================== INITIALIZATION =====================
clc; clear; close all;

% 常数
hbar = 1.055e-34;   % J·s
me   = 9.109e-31;   % kg
eV2J = 1.602e-19;   % J
um2m = 1e6;         % μm^-1 -> m^-1

% 统一系数 (单位: eV·μm^2)
m_coef = (hbar^2 * um2m^2) / (2 * me * eV2J);

%% 文件路径设置
filebox = 'E:\课题组资料\5-汇报类\课题\BA2PbI4\ARPL\z5aa13new';

% 获取所有 _ekcorrection.dat 文件
ek_files = dir(fullfile(filebox,'*_ekcorrection.dat'));
% 获取所有 _lpbfit.dat 文件
lpb_files = dir(fullfile(filebox,'*_lpbfit.dat'));

%% ===================== PROCESS =====================
for ii = 1:length(ek_files)
    temp_ek = fullfile(filebox, ek_files(ii).name);
    temp_lpb = fullfile(filebox, lpb_files(ii).name);
    
    % ================= 导入校正后的能量-k数据 =================
    data = readmatrix(temp_ek);
    k_corr = data(1,2:end);     % 校正后的k
    E = data(2:end,1);          % 能量
    I = data(2:end,2:end);      % 强度矩阵
    
    % 校正后的k-E图
    mask = abs(k_corr) <= 10;
    figure;
    pcolor(k_corr(mask), E, I(:,mask));
    shading interp;   % 平滑显示，没有格子线
    colormap(magma);
    colorbar;
    xlabel('k (\mum^{-1})');
    ylabel('Energy (eV)');
    title('Corrected k-E Map');
    caxis([0,1]);  % 可根据需要调整
    
    hold on;
    
    % ================= 导入LPB拟合参数 =================
    param_cell = readcell(temp_lpb,'Delimiter',',');
    % 将字段名和值对应
    params = struct();
    for jj = 1:size(param_cell,1)
        params.(param_cell{jj,1}) = cell2mat(param_cell(jj,2));
    end
    
    E_C0    = params.E_C0+0.03;
    alpha   = params.alpha+0.00045;
    E_X     = 2.42;%params.E_X;
    gfactor = 0.228;
    
    % ================= 计算LPB和UPB曲线 =================
    E_C = E_C0 + alpha*(k_corr).^2;
    
    E_LPB = (E_C + E_X)/2 - 0.5*sqrt((E_C - E_X).^2 + gfactor^2);
    E_UPB = (E_C + E_X)/2 + 0.5*sqrt((E_C - E_X).^2 + gfactor^2);
    
    % ================= 绘制LPB/UPB =================
    plot(k_corr, E_LPB, 'r-', 'LineWidth', 2);
    plot(k_corr, E_UPB, 'b-', 'LineWidth', 2);
    plot(k_corr, E_C, 'y--', 'LineWidth', 2);
    yline(E_X,'g--', 'LineWidth', 2);
    legend('Reflectivity','LPB','UPB');
    
    grid on;
    hold off;
    
    % ================= 保存图像 =================
    [~, name, ~] = fileparts(temp_ek);
    saveas(gcf, fullfile(filebox, [name,'_LPB_UPB.png']));
    
    lam_c = 1240./E_C;
    lam_x = 1240./E_X;
    lam_l = 1240./E_LPB;
    lam_u = 1240./E_UPB;

    %% ================= 拟合上下能支 =================
    Elp_real = nan(size(k_corr));
    Eup_real = nan(size(k_corr));
    
    for kk = 1:length(k_corr)
        % 只处理 |k| <= 10 的点
        if abs(k_corr(kk)) > 10
            continue;
        end
        
        % 找 LPB 附近最小值
        idx_range = find(E >= E_LPB(kk)-0.05 & E <= E_LPB(kk)+0.05);
        if ~isempty(idx_range)
            [~, min_idx] = min(I(idx_range, kk));
            Elp_real(kk) = E(idx_range(min_idx));
        end
        
        % 找 UPB 附近最小值
        idx_range = find(E >= E_UPB(kk)-0.05 & E <= E_UPB(kk)+0.05);
        if ~isempty(idx_range)
            [~, min_idx] = min(I(idx_range, kk));
            Eup_real(kk) = E(idx_range(min_idx));
        end
    end
    
    % 只绘制 k<=10 的点
    mask_k = abs(k_corr) <= 10;
    
    figure; 
    pcolor(k_corr(mask_k), E, I(:,mask_k));
    shading interp; colormap(magma); colorbar;
    xlabel('k (\mum^{-1})'); ylabel('Energy (eV)'); title('Corrected k-E Map');
    hold on;
    plot(k_corr(mask_k), Elp_real(mask_k),'r.-','LineWidth',1.5);
    plot(k_corr(mask_k), Eup_real(mask_k),'b.-','LineWidth',1.5);
    plot(k_corr(mask_k), E_C(mask_k),'y--','LineWidth',2);
    yline(E_X,'g--','LineWidth',2);
    legend('Reflectivity','LPB_{real}','UPB_{real}','Cavity','Exciton');
    hold off;

    %% ================= 联合拟合（不含线宽） =================
    % 注意：此处去掉激子/腔线宽的影响，使用实值强耦合公式
    gamma_ex = 0.0517;
    gamma_ph = 0.0318;

    % 构造 Elp_real 和 Eup_real 已在上面计算好（含 NaN）
    % 仅取有效数据 (k<=10 & 非 NaN)
    mask_k = abs(k_corr) <= 10;
    valid_idx = find(mask_k & ~isnan(Elp_real) & ~isnan(Eup_real));
    if isempty(valid_idx)
        warning('没有满足 |k|<=10 且 Elp/Eup 都有效的数据，跳过拟合。');
    else
        k_fit = k_corr(valid_idx);
        Elp_fit = Elp_real(valid_idx);
        Eup_fit = Eup_real(valid_idx);
    
        % 联合残差函数（不含线宽，实值）
        joint_res_lw = @(p) [ ...
            real( 0.5*(p(4) + (p(1)+p(2)*k_fit.^2) - 1i*(gamma_ph+gamma_ex)) ...
                 - sqrt( p(3)^2 + 0.25*((p(4)-(p(1)+p(2)*k_fit.^2)) + 1i*(gamma_ph-gamma_ex)).^2 ) ) - Elp_fit; ...
            real( 0.5*(p(4) + (p(1)+p(2)*k_fit.^2) - 1i*(gamma_ph+gamma_ex)) ...
                 + sqrt( p(3)^2 + 0.25*((p(4)-(p(1)+p(2)*k_fit.^2)) + 1i*(gamma_ph-gamma_ex)).^2 ) ) - Eup_fit ...
        ];
    
        % 初始猜测 [E_C0, alpha, Omega_R]
        p0 = [E_C0, alpha, gfactor, E_X];
    
        if exist('lsqnonlin','file')
            opts = optimoptions('lsqnonlin','Display','iter','MaxFunctionEvaluations',1e4);
            [p_fit_lw, ~] = lsqnonlin(joint_res_lw, p0, [], [], opts);
        else
            sqsum = @(p) sum(joint_res_lw(p).^2);
            p_fit_lw = fminsearch(sqsum, p0);
        end    

        % 输出拟合结果
        fprintf('\n=== Joint fit (no linewidth) result for file: %s ===\n', ek_files(ii).name);
        fprintf('E_C0  = %.6f eV\n', p_fit_lw(1));
        fprintf('alpha = %.8f eV·μm^2\n', p_fit_lw(2));
        fprintf('gfactor = %.6f eV\n', p_fit_lw(3));
    
        % 计算拟合曲线

        % 带线宽拟合曲线
        ELPB_fit_lw = real( 0.5*(p_fit_lw(4) + (p_fit_lw(1)+p_fit_lw(2)*k_fit.^2) - 1i*(gamma_ph+gamma_ex)) ...
                     - sqrt( p_fit_lw(3)^2 + 0.25*((p_fit_lw(4)-(p_fit_lw(1)+p_fit_lw(2)*k_fit.^2)) + 1i*(gamma_ph-gamma_ex)).^2 ) );
        EUPB_fit_lw = real( 0.5*(p_fit_lw(4) + (p_fit_lw(1)+p_fit_lw(2)*k_fit.^2) - 1i*(gamma_ph+gamma_ex)) ...
                     + sqrt( p_fit_lw(3)^2 + 0.25*((p_fit_lw(4)-(p_fit_lw(1)+p_fit_lw(2)*k_fit.^2)) + 1i*(gamma_ph-gamma_ex)).^2 ) );
        EPH_fit = p_fit_lw(1) + p_fit_lw(2)*k_fit.^2;

        % 绘图：实测点 + 联合拟合曲线
        figure;
        pcolor(k_corr(mask_k), E, I(:,mask_k));
        shading interp; colormap(magma); colorbar;
        xlabel('k (\mum^{-1})'); ylabel('Energy (eV)'); title('LPB/UPB Real + Joint Fit (no linewidth)');
        hold on;
        plot(k_fit, ELPB_fit_lw,'r-','LineWidth',2);
        plot(k_fit, EUPB_fit_lw,'b-','LineWidth',2);
        plot(k_fit, Elp_fit,'ro','MarkerSize',6,'MarkerFaceColor','r');
        plot(k_fit, Eup_fit,'bo','MarkerSize',6,'MarkerFaceColor','b');
        plot(k_fit, p_fit_lw(1) + p_fit_lw(2)*k_fit.^2,'y--','LineWidth',1.5); % cavity dispersion
        yline(E_X,'g--','LineWidth',2);
        legend('Reflectivity','LPB Fit','UPB Fit','LPB Real','UPB Real','Cavity','Exciton','Location','northwest');
        hold off;
    
        % 保存拟合图（可选）
        [~, name, ~] = fileparts(temp_ek);
        saveas(gcf, fullfile(filebox, [name,'_LPB_UPB_jointfit_nolw.png']));
    end


    %% ================= 导出结果 =================
    % 拼接表头
    header = {'k (um^-1)', 'ELP_real (eV)', 'EUP_real (eV)', ...
              'ELP_fit (eV)', 'EUP_fit (eV)', 'EPH_fit (eV)'};
    
    % 拼接数据（只导出有效 k）
    export_data = [k_fit(:), Elp_fit(:), Eup_fit(:), ELPB_fit_lw(:), EUPB_fit_lw(:), EPH_fit(:)];
    
    % 保存为 dat 文件
    out_file = fullfile(filebox, [name,'_LPB_UPB_results.dat']);
    fid = fopen(out_file, 'w');
    fprintf(fid, '%s\t%s\t%s\t%s\t%s\t%s\n', header{:});
    fclose(fid);
    dlmwrite(out_file, export_data, '-append', 'delimiter', '\t', 'precision', 6);
    
    fprintf('结果已导出到: %s\n', out_file);
end
