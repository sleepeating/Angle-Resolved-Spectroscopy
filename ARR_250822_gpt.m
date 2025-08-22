%% ===================== INITIALIZATION =====================
%% Clear workspace
clc;        % Clear command window
clear;      % Clear workspace variables
close all;  % Close all figures

%% 实验数据文件路径设置
filebox = ['G:\2024年数据\20250114 ARPL BA2PbI4(E2)\NEWfit'];
file = search_folder(filebox,'txt');  % 搜索指定路径下的txt文件
n_max = length(file);  % 获取文件总数

%% ====================== 初始参数配置 ======================
xsize = 1340;   % X轴维度尺寸（根据原始数据维度设置）
ysize = 400;    % Y轴维度尺寸
wavmin = 495;   % 波长最小值(nm)
wavmax = 575;   % 波长最大值(nm)
enemin = 1240./wavmax;   % 能量最小值(eV)
enemax = 1240./wavmin;   % 能量最大值(eV)
NA = 0.9;       % 数值孔径
nair = 1;       % 空气折射率
theta = asind(NA/nair);  % 最大接收角（度）

%% ==================== 数据结构初始化 ======================
% 创建结构数组存储文件信息：
% name - 文件名
% data - 原始角分辨光谱矩阵
% ref  - 反射光谱矩阵（白光背景）
for ii = 1:n_max
    cfsj(ii).name = cell2mat(file(ii));  % 转换单元格为字符串
    cfsj(ii).data = [];
    cfsj(ii).ref  = [];
end

%% ==================== 白光背景数据处理（改进版） ====================
for ii = 1:n_max
    % 筛选符合特征的白光基准文件（名称包含'WHITE'、'150'和'580'）
    if contains(cfsj(ii).name,'WHITE') 
        if contains(cfsj(ii).name,'600') && contains(cfsj(ii).name,'535')
            
            % 数据导入与维度重塑
            arplfile = importdata(cfsj(ii).name);
            wavelength = reshape(arplfile(:,1),[xsize,ysize]);  % 波长矩阵
            pixellength = reshape(arplfile(:,2),[xsize,ysize]); % 像素坐标矩阵
            intensity = reshape(arplfile(:,3),[xsize,ysize]);   % 强度矩阵

            % 角度校准参数提取（上下边界位置）
            pixeldata = linspace(1,xsize,xsize);
            key_pixel = [100,400,700,1000,1300]; % 光谱矩阵中的上下边界索引
            figure;hold on
            h = plot(pixellength(1,:), intensity(key_pixel,:));  % 绘制多条曲线
            % 创建对应的 legend 标签
            legend_labels = arrayfun(@(x) ['Pixel = ', num2str(x)], key_pixel, 'UniformOutput', false);
            % 设置 legend
            legend(h, legend_labels);
            xlabel('Pixel coordinate');
            ylabel('Intensity (a.u.)');
            title('Key pixel intensity profiles');
            for tt = 1:length(key_pixel)
                param = dual_logistic(pixellength(tt,:), intensity(key_pixel(tt),:));
                pixel_l_manual(tt) = param.x1;
                pixel_r_manual(tt) = param.x2;
            end
            tz = 0; lt = 0; rt = 0;
            pixel_l_manual = pixel_l_manual-tz+lt;
            pixel_r_manual = pixel_r_manual+tz-rt;
            
            %% ---------------- 多项式拟合边界随波长变化 ----------------
            pixel_l_poly = polyfit(key_pixel, pixel_l_manual, 2);
            pixel_r_poly = polyfit(key_pixel, pixel_r_manual, 2);
            pixel_l_fit = polyval(pixel_l_poly, pixeldata);
            pixel_r_fit = polyval(pixel_r_poly, pixeldata);
            figure;
            plot (pixeldata,pixel_l_fit,pixeldata,pixel_r_fit);
            
            %% ---------------- 波长依赖的像素到角度转换 ----------------
            thetalength = zeros(size(pixellength));
            for jj = 1:xsize
                center_pixel = (pixel_l_fit(jj) + pixel_r_fit(jj))/2;
                pixel2theta_jj = theta./((pixel_r_fit(jj) - pixel_l_fit(jj))/2);               
                thetalength(jj,:) = (pixellength(jj,:) - center_pixel) * pixel2theta_jj;
            end
            
            %% ---------------- 插值到规则网格 ----------------
            toffset = 5;
            thetaout = linspace(-theta-toffset, theta+toffset, 501);  % 角度网格
            waveout  = wavelength(:,1);               % 波长基准
            [thetaout, waveout] = meshgrid(thetaout, waveout);
            
            F = scatteredInterpolant(thetalength(:), wavelength(:), intensity(:));
            white_bg = F(thetaout, waveout);  % 生成背景白光基准数据
            
            %% ---------------- 可视化 ----------------
            figure;
            pcolor(thetaout, waveout, white_bg);
            shading interp;
            colormap magma; colorbar;
            xlabel('Angle (°)'); ylabel('Wavelength (nm)');
            title(['Calibrated white background']);
            hold on;  % 保持当前图像
            % theta = 50° 的竖直线
            theta_line = theta;  % 角度
            plot([theta_line theta_line], [min(waveout(:)) max(waveout(:))], 'w--', 'LineWidth', 1.5);
            plot([-theta_line -theta_line], [min(waveout(:)) max(waveout(:))], 'w--', 'LineWidth', 1.5);
        end
    end
end

%% ==================== 样品数据处理 ========================
for ii = 1:n_max
    % 筛选样品数据文件（排除白光文件，包含'E2'和'600'）
    if isempty(findstr(cfsj(ii).name,'WHITE'))   
      if ~isempty(findstr(cfsj(ii).name,'600')) && ~isempty(findstr(cfsj(ii).name,'S'))
        % 数据导入与插值处理
        arplfile = importdata(cfsj(ii).name);
        F = scatteredInterpolant(thetalength(:),wavelength(:),arplfile(:,3));
        cfsj(ii).data = F(thetaout,waveout);  % 生成规则网格数据
        cfsj(ii).ref = cfsj(ii).data./white_bg;
        cfsjtref = cfsj(ii).ref;

        % 反射光谱计算与可视化
        figure
        subplot(1,2,1);
        pcolor(thetaout,waveout,cfsj(ii).data);
        shading interp
        ylim([wavmin,wavmax]);  % 设置波长显示范围
        set(gcf,'Colormap',magma);
        xlabel('Angle (°)','FontName','Arial','FontSize',10);
        ylabel('Wavelength (nm)','FontName','Arial','FontSize',10);
        set(gca, 'Fontname', 'Arial', 'Fontsize', 10);
        set(gca,'Position',[.17 .17 .3 .73]);  % 调整子图位置
        subplot(1,2,2);
        pcolor(thetaout,waveout,cfsj(ii).ref);
        shading interp
        ylim([wavmin,wavmax]);
        caxis([-1,1]);  % 设置色标范围
        set(gcf,'Colormap',magma);
        xlabel('Angle (°)','FontName','Arial','FontSize',10);
        ylabel('Wavelength (nm)','FontName','Arial','FontSize',10);
        set(gca, 'Fontname', 'Arial', 'Fontsize', 10);
        set(gca,'Position',[.57 .17 .3 .73]);

        temp_outfile = cfsj(ii).name;
        dot_position = regexp(temp_outfile,'\.');
        print([temp_outfile(1:dot_position(end)-1),'_ref.tif'], '-dtiffn','-r600');


        % 波矢转换：k = 2π/λ * sinθ
        kpout = 1000.*2.*pi./waveout.*sind(thetaout);
        
        % ---------------- 生成规则波矢-波长网格 ----------------
        k_min = -15;      % μm^-1
        k_max = 15;       % μm^-1
        k_res = 300;
        wave_res = 400;
        k_new    = linspace(k_min, k_max, k_res);        
        wave_new = linspace(wavmin, wavmax, wave_res);
        [kq_mesh, waveq_mesh] = meshgrid(k_new, wave_new);
        
        % ---------------- 插值到新规则网格 ----------------
        Fk = scatteredInterpolant(kpout(:), waveout(:), cfsjtref(:));
        cfsjrefk = Fk(kq_mesh, waveq_mesh);
        
        % ---------------- 后处理（填充NaN） ----------------
        nan_mask = isnan(cfsjrefk);
        cfsjrefk(nan_mask) = 0;
        
        % ---------------- 绘图 ----------------
        figure;
        pcolor(kq_mesh, waveq_mesh, cfsjrefk);
        shading interp;
        colormap("magma"); colorbar;
        xlabel('k (μm^{-1})');
        ylabel('Wavelength (nm)');
        line_position = regexp(temp_outfile,'\\');
        title(['Reflection k2wave Map - ', temp_outfile(line_position(end)+1:end)]);
        caxis([0.6,1]);  % 设置色标范围

        % ---------------- 导出二维矩阵 ----------------
        % --- (1) angle-wave 数据 (arrdata)
        arrdata = [0, thetaout(1,:); waveout(:,1), cfsjtref];
        file_arr = [temp_outfile(1:dot_position(end)-1),'_arrdata.dat'];
        writematrix(arrdata, file_arr);

        % --- (2) k-wave 数据 (kmapdata)
        kmapdata = [0, kq_mesh(1,:); waveq_mesh(:,1), cfsjrefk];
        file_kmap = [temp_outfile(1:dot_position(end)-1),'_kmapdata.dat'];
        writematrix(kmapdata, file_kmap);

        % 可选保存图片
        print([temp_outfile(1:dot_position(end)-1),'_k_lambda_ref.tif'], '-dtiffn','-r600');
      end
    end
end



function params = dual_logistic(x, y)
    % 双逻辑函数拟合
    % 调用方法：
    % params = dual_logistic(x, y);
    % 返回参数结构体: params.A, params.c, params.x1, params.x2, params.k1, params.k2
    
    % 双逻辑函数形式
    dualLogisticFun = @(p, x) p(2) + ...
        p(1) ./ (1 + exp(-p(5) * (x - p(3)))) - ...
        p(1) ./ (1 + exp(-p(6) * (x - p(4))));
    
    % 初始参数估计
    A0 = (max(y) - min(y))/2;   % 初始幅值
    c0 = min(y);                % 初始基线
    [~, idxMax] = max(y);       % 最大值附近的点
    [~, idxMin] = min(y);       % 最小值附近的点
    x1_0 = x(idxMin);           % 上升点猜测
    x2_0 = x(idxMax);           % 下降点猜测
    k1_0 = 1;                   % 初始斜率
    k2_0 = 1;                   % 初始斜率
    
    p0 = [A0, c0, x1_0, x2_0, k1_0, k2_0];
    
    % 使用 lsqcurvefit 拟合
    opts = optimset('Display', 'off');
    lb = [0, -Inf, min(x), min(x), 0, 0];  % 参数下限
    ub = [Inf, Inf, max(x), max(x), Inf, Inf]; % 参数上限
    
    pFit = lsqcurvefit(dualLogisticFun, p0, x, y, lb, ub, opts);
    
    % 输出参数
    params.A = pFit(1);
    params.c = pFit(2);
    params.x1 = pFit(3);
    params.x2 = pFit(4);
    params.k1 = pFit(5);
    params.k2 = pFit(6);
    
    % 绘图对比
    figure;
    scatter(x, y, 'b','filled'); hold on;
    xx = linspace(min(x), max(x), 500);
    plot(xx, dualLogisticFun(pFit, xx), 'r-', 'LineWidth', 2);
    legend('Data','Dual Logistic Fit');
    xlabel('x'); ylabel('y');
    title('Dual Logistic Fit');
end