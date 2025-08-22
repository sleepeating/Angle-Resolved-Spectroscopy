%% ===================== INITIALIZATION =====================
%% Clear workspace
clc;        % Clear command window
clear;      % Clear workspace variables
close all;  % Close all figures

%% 文件路径
filebox = ['G:\2024年数据\20250114 ARPL BA2PbI4(E2)\NEWfit'];
file = search_folder(filebox, 'dat');
n_max = length(file);

for ii = 1:n_max
    cfsj(ii).name = cell2mat(file(ii));  % 转换单元格为字符串
end

for ii = 1:n_max
    % 读取数据
    temp_file = cfsj(ii).name;
    if contains(temp_file,'kmap')
        data = readmatrix(temp_file);  % 假设第一列波长（nm），第一行是k
        
        k = data(1,2:end);           % k 值
        lambda = data(2:end,1);      % 波长 (nm)
        I_k = data(2:end,2:end);  % 强度矩阵
        
        %% 绘制波长域k-map
        figure;
        imagesc(k, lambda, I_k);
        set(gca,'YDir','normal');
        xlabel('k (\mum^{-1})');
        ylabel('Wavelength (nm)');
        title(['k-map (Wavelength)']);
        caxis([0,1]);
        colormap("magma"); colorbar;
        
        %% 波长 -> 能量转换 (eV)
        E = 1240 ./ lambda;  % eV
        
        %% 绘制能量域k-map
        figure;
        imagesc(k, E, I_k);
        set(gca,'YDir','normal');
        xlabel('k (\mum^{-1})');
        ylabel('Energy (eV)');
        title(['k-map (Energy)']);
        caxis([0,1]);
        colormap("magma"); colorbar;

        %% 保存能量域数据
        % 第一行是能量，第一列是 k，后续是强度矩阵
        ekmap = [0, k; E, I_k];
        gang_position = regexp(temp_file,'\_');
        file_ekout = [temp_file(1:gang_position(end)-1), '_ekout.dat'];
        writematrix(ekmap, file_ekout);
        
        disp(['已保存能量域数据: ', file_ekout]);
    end
end
