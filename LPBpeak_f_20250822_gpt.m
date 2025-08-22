%% ===================== INITIALIZATION =====================
clc; clear; close all;

%% 文件路径
filebox = 'G:\2024年数据\20250114 ARPL BA2PbI4(E2)\NEWfit';
file = search_folder(filebox, 'dat');   % 搜索dat文件
n_max = length(file);

% 能量区间设置（eV）
E_min = 1240/540;  % 例子，可修改
E_max = 1240/520;  % 例子，可修改

%% ==================== 数据结构初始化 ======================
% 创建结构数组存储文件信息：
% name - 文件名
% data - 原始角分辨光谱矩阵
% ref  - 反射光谱矩阵（白光背景）
for ii = 1:n_max
    cfsj(ii).name = cell2mat(file(ii));  % 转换单元格为字符串
end

%% ===================== PROCESS EACH FILE =====================
for ii = 1:n_max
    temp_file = cfsj(ii).name;
    
    if contains(temp_file, 'ekout') && contains(temp_file, '535') % 读取能量域数据
        data = readmatrix(temp_file);
        
        % 第一行是能量，第一列是 k
        k = data(1,2:end);
        E = data(2:end,1);
        I = data(2:end,2:end);
        
        % 预分配峰值信息
        peak_center = zeros(length(k),1);
        peak_width  = zeros(length(k),1);
        peak_height = zeros(length(k),1);
        
        % 循环每个 k 点
        for jj = 1:length(k)
            I_k = I(:,jj);  % 当前 k 点的反射谱
            
            % 限定能量区间
            idx_range = E >= E_min & E <= E_max;
            E_sub = E(idx_range);
            I_sub = I_k(idx_range);
            
            % 找峰值（基础版：直接找最大值）
            [Imin, imin] = min(I_sub);
            E_peak = E_sub(imin);
            
            % 保存峰值信息
            peak_center(jj) = E_peak;
            peak_width(jj)  = NaN;  % 后续可用洛伦兹拟合计算
            peak_height(jj) = Imin;
        end
        
        %% 可视化
        figure;
        plot(k, peak_center, 'o-');
        xlabel('k (\mum^{-1})');
        ylabel('Peak Energy (eV)');
        title(['Peak energy vs k']);

        % 原始 k 和峰值能量
        k_all = k;                  % k 数据向量
        E_peak_all = peak_center;    % 对应峰值能量
        
        % 筛选 k 范围 [-9, 9]
        mask = (k_all >= -6.3) & (k_all <= 6.3);
        k_fit = k_all(mask);
        E_peak_fit = E_peak_all(mask);
        
        % 调用拟合函数
        params = fit_LPB(k_fit(:), E_peak_fit(:));
        
        % 导出拟合参数
        disp(params);

        % ===== 保存参数到dat =====
        % 提取字段名和数值
        fields = fieldnames(params);
        values = struct2cell(params);
        
        % 拼成 cell，第一列字段名，第二列数值
        param_cell = [fields num2cell(cell2mat(values))];
        gang_position = regexp(temp_file,'\_');
        file_lpbout = [temp_file(1:gang_position(end)-1), '_lpbfit.dat'];
        
        % 保存为 dat
        writecell(param_cell, file_lpbout, 'Delimiter','tab');
        disp(['拟合参数已保存: ', file_lpbout]);

    end
end

function params = fit_LPB(k_data, E_peak)
    % 拟合下能支极化激元（带k偏移）
    % k_data: 波矢 (μm^-1)
    % E_peak: 峰值能量 (eV)
    %
    % 输出 params 结构体:
    %   E_C0    - 空腔模能量 (eV)
    %   alpha   - 光子色散系数
    %   E_X     - 激子能量 (eV)
    %   Omega_R - 拉比分裂 (eV)
    %   k_offset- k 偏移 (μm^-1)

    % 下能支模型函数 (带 k offset)
    LPB_fun = @(p, k) (p(1) + p(2)*(k - p(5)).^2 + p(3))/2 + ...
                       (-1/2) * sqrt((p(3) - (p(1) + p(2)*(k - p(5)).^2)).^2 + p(4)^2);
    % p = [E_C0, alpha, E_X, Omega_R, k_offset]

    % 初始猜测参数
    E_C0_0    = min(E_peak);
    alpha_0   = 0.0008;  % 光子色散初值，可根据实验调
    Omega0    = 0.228;    % 拉比劈裂初值 (eV)
    E_X0      = max(E_peak) + Omega0/2;
    k_offset0 = 0;      % k 偏移初值

    p0 = [E_C0_0, alpha_0, E_X0, Omega0, k_offset0];

    % 设置拟合上下界
    lb = [2, 0, 2, 0, -5];
    ub = [3, 0.005, 3, Inf, 5];

    % 使用 lsqcurvefit 拟合
    opts = optimset('Display','off');
    pFit = lsqcurvefit(LPB_fun, p0, k_data, E_peak, lb, ub, opts);

    % 输出拟合参数
    params.E_C0    = pFit(1);
    params.alpha   = pFit(2);
    params.E_X     = pFit(3);
    params.Omega_R = pFit(4);
    params.k_offset= pFit(5);

    % 绘图对比
    figure;
    scatter(k_data, E_peak, 'b','filled'); hold on;
    k_fit = linspace(min(k_data), max(k_data), 200);
    plot(k_fit, LPB_fun(pFit, k_fit), 'r-', 'LineWidth', 2);
    xlabel('k (\mum^{-1})'); ylabel('Energy (eV)');
    title('Lower Polariton Branch Fit with k_{offset}');
    legend('Experimental Peaks', 'LPB Fit');
    grid on;
end
