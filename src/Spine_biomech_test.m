close all; 
clear; 
clc;


%%
path = ('sample_data');
cd (path);
addpath('code');

%%
Spine_static_raw = ezc3dRead('SpineVST_Buiding_ST.c3d'); % read all info from the c3d file
fs = Spine_static_raw.parameters.POINT.RATE.DATA; % sampling frequency

Marker_data_raw = Spine_static_raw.data.points;
Mrker_label = Spine_static_raw.parameters.POINT.LABELS.DATA;

MrkerCell = cell(63,1);
for m = 1:63
    MrkerCell{m} = squeeze(Marker_data_raw(:,m,:))';  
end

% filter the marker data
fc = 6;          % <-- cutoff (Hz), try 6-10 for marker traj
order = 4;
[b,a] = butter(order, fc/(fs/2), 'low');

MrkerCell_f = MrkerCell;  % same size
for m = 1:numel(MrkerCell)
    X = MrkerCell{m};              % N x 3
    MrkerCell_f{m} = filtfilt(b,a, X);   % filters each column
end




%% ---------- 1) Build SpineCell (L / S / R) and reorder anatomically ----------
labels = string(Mrker_label(:));     % 63x1
isL = endsWith(labels,"L");
isR = endsWith(labels,"R");
isS = ~isL & ~isR;

base = labels;
base(isL | isR) = extractBefore(base(isL | isR), strlength(base(isL | isR))); % drop last char

% Canonical anatomy order you want:
anatOrder = ["C6","C7", ...
             "T1","T2","T3","T4","T5","T6","T7","T8","T9","T10","T11","T12", ...
             "L1","L2","L3","L4","L5", ...
             "S1","S2","S3","S4","S5"];

% Keep only those that actually exist in your dataset (spinous labels drive the "levels" list)
levels_in_data = base(isS);
levels_in_data = unique(levels_in_data,'stable');

% Reorder to anatOrder (and remove any that don't exist)
levels = anatOrder(ismember(anatOrder, levels_in_data));
nLevels = numel(levels);

SpineCell = cell(nLevels,3);   % col1=L, col2=S, col3=R
for k = 1:nLevels
    lev = levels(k);

    idxL = find(labels == lev + "L", 1, "first");
    idxS = find(labels == lev,        1, "first");
    idxR = find(labels == lev + "R",  1, "first");

    if ~isempty(idxL), SpineCell{k,1} = MrkerCell_f{idxL}; end
    if ~isempty(idxS), SpineCell{k,2} = MrkerCell_f{idxS}; end
    if ~isempty(idxR), SpineCell{k,3} = MrkerCell_f{idxR}; end
end

SpineLevels = cellstr(levels);

% quick missing report
fprintf("Missing L: %d, S: %d, R: %d\n", ...
    sum(cellfun(@isempty,SpineCell(:,1))), ...
    sum(cellfun(@isempty,SpineCell(:,2))), ...
    sum(cellfun(@isempty,SpineCell(:,3))) );


%% ---------- 2) Plot ONE frame for inspection (all markers + labels) ----------
frameIdx = 1;   % change this to inspect another frame

figure('Color','w'); hold on; grid on; axis equal;
xlabel('X'); ylabel('Y'); zlabel('Z');
title(sprintf('Spine markers @ frame %d', frameIdx));

for k = 1:nLevels
    lev = string(SpineLevels{k});

    % Left
    if ~isempty(SpineCell{k,1})
        p = SpineCell{k,1}(frameIdx,:);     % 1x3
        plotAndLabel(p, lev+"L");
    end

    % Spinous
    if ~isempty(SpineCell{k,2})
        p = SpineCell{k,2}(frameIdx,:);
        plotAndLabel(p, lev);
    end

    % Right
    if ~isempty(SpineCell{k,3})
        p = SpineCell{k,3}(frameIdx,:);
        plotAndLabel(p, lev+"R");
    end
end

view(3);





%% Inputs you already have:
% SpineCell : nLevels x 3 cell  [L, S(mid), R], each cell is N x 3
% SpineLevels : nLevels x 1 cellstr or string (e.g., "C7","T1"...)

% frameIdx = 1;  % pick a frame to inspect (static)

% --- Assemble raw column points at this frame (nLevels x 3) ---
P_L = nan(numel(SpineLevels),3);
P_M = nan(numel(SpineLevels),3);
P_R = nan(numel(SpineLevels),3);

for k = 1:numel(SpineLevels)
    if ~isempty(SpineCell{k,1}), P_L(k,:) = SpineCell{k,1}(frameIdx,:); end
    if ~isempty(SpineCell{k,2}), P_M(k,:) = SpineCell{k,2}(frameIdx,:); end
    if ~isempty(SpineCell{k,3}), P_R(k,:) = SpineCell{k,3}(frameIdx,:); end
end

% keep only rows where all 3 columns exist (recommended for spline surface)
valid = all(isfinite(P_L),2) & all(isfinite(P_M),2) & all(isfinite(P_R),2);
P_L = P_L(valid,:);  P_M = P_M(valid,:);  P_R = P_R(valid,:);
levels_valid = string(SpineLevels(valid));
nLev = size(P_M,1);

% s: parameter = spine level index (1..nLev)
s = (1:nLev);

nPieces = 6;   % six segments (as paper)
kOrder  = 4;   % cubic spline => order 4

% Left column splines
spLx = spap2(nPieces, kOrder, s, P_L(:,1)'); 
spLy = spap2(nPieces, kOrder, s, P_L(:,2)'); 
spLz = spap2(nPieces, kOrder, s, P_L(:,3)');

% Mid column splines
spMx = spap2(nPieces, kOrder, s, P_M(:,1)'); 
spMy = spap2(nPieces, kOrder, s, P_M(:,2)'); 
spMz = spap2(nPieces, kOrder, s, P_M(:,3)');

% Right column splines
spRx = spap2(nPieces, kOrder, s, P_R(:,1)'); 
spRy = spap2(nPieces, kOrder, s, P_R(:,2)'); 
spRz = spap2(nPieces, kOrder, s, P_R(:,3)');

% Dense evaluation
u = linspace(1, nLev, 100);
S_L = [fnval(spLx,u); fnval(spLy,u); fnval(spLz,u)]';
S_M = [fnval(spMx,u); fnval(spMy,u); fnval(spMz,u)]';
S_R = [fnval(spRx,u); fnval(spRy,u); fnval(spRz,u)]';


% --- Index each vertebral level to nearest midline spline point (paper step) ---
u_idx = zeros(nLev,1);
for k = 1:nLev
    d = vecnorm(S_M - P_M(k,:), 2, 2);
    [~,u_idx(k)] = min(d);
end

% --- Build LCS at each level ---
% Rotation matrix columns will be [ML AP SI] (x=ML, y=AP, z=SI)
R = nan(3,3,nLev);
O = nan(nLev,3);

% derivatives of midline spline for tangent (SI)
dspMx = fnder(spMx); dspMy = fnder(spMy); dspMz = fnder(spMz);

for k = 1:nLev
    uk = u(u_idx(k));

    Om = [fnval(spMx,uk), fnval(spMy,uk), fnval(spMz,uk)];
    Ol = [fnval(spLx,uk), fnval(spLy,uk), fnval(spLz,uk)];
    Or = [fnval(spRx,uk), fnval(spRy,uk), fnval(spRz,uk)];

    % SI = tangent to midline
    t = [fnval(dspMx,uk), fnval(dspMy,uk), fnval(dspMz,uk)];
    SI = t / norm(t);

    % LR vector (right - left): defines frontal plane with SI
    LR = (Or - Ol);
    LR = LR / norm(LR);

    % AP = perpendicular to SI and LR
    AP = cross(SI, LR);
    AP = AP / norm(AP);

    % ML = perpendicular to SI and AP (right-handed)
    ML = cross(AP, SI);
    ML = ML / norm(ML);

    % store
    O(k,:) = Om;
    R(:,:,k) = [ML(:), AP(:), SI(:)];

    % enforce right-handed (det=+1)
    if det(R(:,:,k)) < 0
        R(:,1,k) = -R(:,1,k); % flip ML if needed
    end
end


%%
figure('Color','w'); hold on; grid on; axis equal;
plot3(S_L(:,1),S_L(:,2),S_L(:,3),'k-','LineWidth',1);
plot3(S_M(:,1),S_M(:,2),S_M(:,3),'g-','LineWidth',2);
plot3(S_R(:,1),S_R(:,2),S_R(:,3),'k-','LineWidth',1);

% raw points
plot3(P_L(:,1),P_L(:,2),P_L(:,3),'ko','MarkerSize',4);
plot3(P_M(:,1),P_M(:,2),P_M(:,3),'ko','MarkerSize',4);
plot3(P_R(:,1),P_R(:,2),P_R(:,3),'ko','MarkerSize',4);

% axes at each level (scale controls arrow length)
axLen = 30;  % change depending on your units (mm)
for k = 1:nLev
    T = eye(4);
    T(1:3,1:3) = R(:,:,k);
    T(1:3,4)   = O(k,:)';
    showcs(T, axLen); % uses your showcs.m
    text(O(k,1),O(k,2),O(k,3), " "+levels_valid(k), 'FontSize',8, 'Interpreter','none');
end
xlabel('X'); ylabel('Y'); zlabel('Z');
view(3);
title(sprintf('Spine LCS inspection @ frame %d', frameIdx));


%%
% helper: Cardan XYZ (intrinsic) from rotation matrix
cardanXYZ = @(Rrel) deal( ...
    atan2( Rrel(3,2), Rrel(3,3) ), ...                         % alpha (X)
    -asin( Rrel(3,1) ), ...                                    % beta  (Y)
    atan2( Rrel(2,1), Rrel(1,1) ) );                            % gamma (Z)

% Example: thoracic = C7 relative to T12, lumbar = T12 relative to S1
% You must know which rows correspond to C7, T12, S1 in levels_valid.
iC7  = find(levels_valid=="C7",1);
iT12 = find(levels_valid=="T12",1);
iS1  = find(levels_valid=="S1",1);

R_C7  = R(:,:,iC7);
R_T12 = R(:,:,iT12);
R_S1  = R(:,:,iS1);

R_thor = R_T12' * R_C7;   % C7 relative to T12 (prox=T12, dist=C7)
R_lumb = R_S1'  * R_T12;  % T12 relative to S1

[aFE,bLB,cAT] = cardanXYZ(R_thor);
thoracic_deg = rad2deg([aFE,bLB,cAT]);

[aFE,bLB,cAT] = cardanXYZ(R_lumb);
lumbar_deg   = rad2deg([aFE,bLB,cAT]);

fprintf('Thoracic [FE LB AT] deg = %.2f, %.2f, %.2f\n', thoracic_deg);
fprintf('Lumbar   [FE LB AT] deg = %.2f, %.2f, %.2f\n', lumbar_deg);



%% functions
function plotAndLabel(xyz, txt)
    plot3(xyz(1), xyz(2), xyz(3), '.', 'MarkerSize', 18); 
    text(xyz(1), xyz(2), xyz(3), [' ' txt], ...
        'FontSize', 9, 'Interpreter', 'none');
end





