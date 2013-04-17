%data analysis for 32nt64 kl3

all_traj=1100:40:10000;

t_size=128;
m_range=8:1:120;

sep={30};
kl3_range={8:1:22};

pl=0.157962;
ps=0.454173;

max_size=length(all_traj);

kaon00 = cell(max_size,1);
kaon01 = cell(max_size,1);
pion00 = cell(max_size,1);
pion01 = cell(max_size,1);
kl3_01g0 = cell(max_size,1);
kl3_01g3 = cell(max_size,1);
kl3_10g0 = cell(max_size,1);
kl3_10g3 = cell(max_size,1);

% import data
n=0;
for traj = all_traj
    trajs=int2str(traj);
    if ~exist(['results.tavg/kaon-00WW.', trajs],'file')
        break
    end

    n=n+1;
    kaon00{n}=importdata(['results.tavg/kaon-00WW.', trajs]);
    kaon10{n}=importdata(['results.tavg/kaon-10WW.', trajs]);
    pion00{n}=importdata(['results.tavg/pion-00WW.', trajs]);
    pion01{n}=importdata(['results.tavg/pion-01WW.', trajs]);
    kl3_01g0{n}=importdata(['results.tavg/kl3-01-g0.',trajs]);
    kl3_01g3{n}=importdata(['results.tavg/kl3-01-g3.',trajs]);
    kl3_10g0{n}=importdata(['results.tavg/kl3-10-g0.',trajs]);
    kl3_10g3{n}=importdata(['results.tavg/kl3-10-g3.',trajs]);
end

m_len=length(m_range);
kl3_len=0;
for i = 1:1:length(kl3_range)
    kl3_len = kl3_len + length(kl3_range{i});
end

% build necessary vectors/matrices
V = cell(n, 1); %% covariance matrix
C = zeros(n, 4 * m_len + 4 * kl3_len);

% construct C matrix
for i = 1:1:n
    % C(i,:) = [pion00{i}(m_range, 2)', ...
    %           pion01{i}(m_range, 2)', ...
    %           kaon00{i}(m_range, 2)', ...
    %           kaon10{i}(m_range, 2)', ...
    %           kl3_01g3{i}(128 * sep{1} + kl3_range{1}, 3)', ...
    %           kl3_01g0{i}(128 * sep{1} + kl3_range{1}, 4)', ...
    %           kl3_10g3{i}(128 * sep{1} + kl3_range{1}, 3)', ...
    %           kl3_10g0{i}(128 * sep{1} + kl3_range{1}, 4)', ...
    %           kl3_01g3{i}(128 * sep{2} + kl3_range{2}, 3)', ...
    %           kl3_01g0{i}(128 * sep{2} + kl3_range{2}, 4)', ...
    %           kl3_10g3{i}(128 * sep{2} + kl3_range{2}, 3)', ...
    %           kl3_10g0{i}(128 * sep{2} + kl3_range{2}, 4)'];

    C(i,1:1:4 * m_len) = [
        pion00{i}(m_range, 2)', ...
        pion01{i}(m_range, 2)', ...
        kaon00{i}(m_range, 2)', ...
        kaon10{i}(m_range, 2)' ];

    shift = 4 * m_len;
    for j = 1:1:length(sep)
        C(i, shift + (1:1:4 * length(kl3_range{j}))) = [
            kl3_01g3{i}(128 * sep{j} + kl3_range{j}, 3)', ...
            kl3_01g0{i}(128 * sep{j} + kl3_range{j}, 4)', ...
            kl3_10g3{i}(128 * sep{j} + kl3_range{j}, 3)', ...
            kl3_10g0{i}(128 * sep{j} + kl3_range{j}, 4)'];
        shift = shift + 4 * length(kl3_range{j});
    end
end

clear pion00 pion01 kaon00 kaon01 
clear kl3_01g0 kl3_01g3 kl3_10g0 kl3_10g3

%% return value is [chi_sq, mk, mp, Zk00, Zk10, Zp00, Zp01, fp, fm]
ret = jackknife(@(C) chi_sq_fit(C, pl, ps, sep, t_size, m_range, kl3_range), C);
vmean = mean(ret);
vsd = (n - 1) / sqrt(n) * std(ret);
