


%% ===== Load dynamic trial =====
dyn = ezc3dRead('Spine_Extention01.c3d');

fs = dyn.parameters.POINT.RATE.DATA;
labels_raw = dyn.parameters.POINT.LABELS.DATA;
labels = string(labels_raw(:));

pts = dyn.data.points;                 % ezc3d: 4 x nMarkers x nFrames (usually)
nMarkers = numel(labels);
nFrames  = size(pts,3);

% Build marker cell: each is nFrames x 3
MrkerCell = cell(nMarkers,1);
for m = 1:nMarkers
    MrkerCell{m} = squeeze(pts(1:3,m,:))';    % nFrames x 3
end

% Filter
fc = 6; order = 4;
[b,a] = butter(order, fc/(fs/2), 'low');

MrkerCell_f = MrkerCell;
for m = 1:nMarkers
    MrkerCell_f{m} = filtfilt(b,a, MrkerCell{m});
end

%% ===== Build SpineCell in anatomical order =====
isL = endsWith(labels,"L");
isR = endsWith(labels,"R");
isS = ~isL & ~isR;

base = labels;
base(isL | isR) = extractBefore(base(isL | isR), strlength(base(isL | isR)));

anatOrder = ["C6","C7", ...
             "T1","T2","T3","T4","T5","T6","T7","T8","T9","T10","T11","T12", ...
             "L1","L2","L3","L4","L5", ...
             "S1","S2","S3","S4","S5"];

levels_in_data = unique(base(isS),'stable');
levels = anatOrder(ismember(anatOrder, levels_in_data));
SpineLevels = cellstr(levels);
nLevels = numel(levels);

SpineCell = cell(nLevels,3); % [L S R], each: nFrames x 3
for k = 1:nLevels
    lev = levels(k);
    idxL = find(labels == lev+"L", 1);
    idxS = find(labels == lev,     1);
    idxR = find(labels == lev+"R", 1);

    if ~isempty(idxL), SpineCell{k,1} = MrkerCell_f{idxL}; end
    if ~isempty(idxS), SpineCell{k,2} = MrkerCell_f{idxS}; end
    if ~isempty(idxR), SpineCell{k,3} = MrkerCell_f{idxR}; end
end


%% ===== Parameters for spline/LCS =====
nPieces = 6;      % 6 segments
kOrder  = 4;      % cubic
nDense  = 100;    % spline points
uDense  = [];     % will be set after we know nLev

% Choose which levels you actually want for LCS/angles (paper often C7..S1)
useMin = "C7"; useMax = "S1";
useMask = string(SpineLevels) >= useMin & string(SpineLevels) <= useMax; % works if labels are consistent; if not, set manually
% If that string compare is sketchy in your MATLAB, do this instead:
% useMask = ismember(string(SpineLevels), ["C7","T1","T2",...,"S1"]);

useIdx_all = find(useMask);
SpineLevels_use = string(SpineLevels(useIdx_all));

% Pre-allocate outputs
% R_all: 3x3xnLev_usexnFrames, O_all: nLev_use x 3 x nFrames
% If memory is a concern, store only angles instead (recommended).
storeR = false;

% Angles over time
thoracic_deg = nan(nFrames,3);   % [FE LB AT]
lumbar_deg   = nan(nFrames,3);   % [FE LB AT]

% Intersegmental (between adjacent levels in the "use" list)
nLev_use = numel(useIdx_all);
FE_int = nan(nFrames, nLev_use-1);
LB_int = nan(nFrames, nLev_use-1);
AT_int = nan(nFrames, nLev_use-1);

% Helper: Cardan XYZ (X=ML, Y=AP, Z=SI)
cardanXYZ = @(Rrel) deal( ...
    atan2( Rrel(3,2), Rrel(3,3) ), ...   % X
    -asin( Rrel(3,1) ), ...              % Y
    atan2( Rrel(2,1), Rrel(1,1) ) );     % Z

%% ===== Main loop over frames =====
for f = 1:nFrames

    % --- Assemble raw column points at this frame ---
    P_L = nan(nLevels,3);
    P_M = nan(nLevels,3);
    P_R = nan(nLevels,3);

    for k = 1:nLevels
        if ~isempty(SpineCell{k,1}), P_L(k,:) = SpineCell{k,1}(f,:); end
        if ~isempty(SpineCell{k,2}), P_M(k,:) = SpineCell{k,2}(f,:); end
        if ~isempty(SpineCell{k,3}), P_R(k,:) = SpineCell{k,3}(f,:); end
    end

    valid = all(isfinite(P_L),2) & all(isfinite(P_M),2) & all(isfinite(P_R),2);
    P_Lv = P_L(valid,:);  P_Mv = P_M(valid,:);  P_Rv = P_R(valid,:);
    levels_valid = string(SpineLevels(valid));
    nLev = size(P_Mv,1);

    if nLev < 6
        % not enough levels to do 6 segments reliably
        continue;
    end

    s = (1:nLev);
    nPieces_eff = min(nPieces, nLev-1);  % safety
    uDense = linspace(1, nLev, nDense);

    % --- Fit splines (spap2) ---
    spLx = spap2(nPieces_eff, kOrder, s, P_Lv(:,1)');  spLy = spap2(nPieces_eff, kOrder, s, P_Lv(:,2)');  spLz = spap2(nPieces_eff, kOrder, s, P_Lv(:,3)');
    spMx = spap2(nPieces_eff, kOrder, s, P_Mv(:,1)');  spMy = spap2(nPieces_eff, kOrder, s, P_Mv(:,2)');  spMz = spap2(nPieces_eff, kOrder, s, P_Mv(:,3)');
    spRx = spap2(nPieces_eff, kOrder, s, P_Rv(:,1)');  spRy = spap2(nPieces_eff, kOrder, s, P_Rv(:,2)');  spRz = spap2(nPieces_eff, kOrder, s, P_Rv(:,3)');

    % Dense midline for nearest-point indexing
    S_M = [fnval(spMx,uDense); fnval(spMy,uDense); fnval(spMz,uDense)]';

    % --- Index each vertebral level to nearest midline spline point ---
    u_idx = zeros(nLev,1);
    for k = 1:nLev
        d = vecnorm(S_M - P_Mv(k,:), 2, 2);
        [~,u_idx(k)] = min(d);
    end

    % --- Build LCS at each valid level ---
    Rv = nan(3,3,nLev);
    Ov = nan(nLev,3);

    dspMx = fnder(spMx); dspMy = fnder(spMy); dspMz = fnder(spMz);

    for k = 1:nLev
        uk = uDense(u_idx(k));

        Om = [fnval(spMx,uk), fnval(spMy,uk), fnval(spMz,uk)];
        Ol = [fnval(spLx,uk), fnval(spLy,uk), fnval(spLz,uk)];
        Or = [fnval(spRx,uk), fnval(spRy,uk), fnval(spRz,uk)];

        % SI tangent
        t = [fnval(dspMx,uk), fnval(dspMy,uk), fnval(dspMz,uk)];
        SI = t / norm(t);

        % LR vector (R-L)
        LR = (Or - Ol);
        LR = LR / norm(LR);

        % AP and ML
        AP = cross(SI, LR);  AP = AP / norm(AP);
        ML = cross(AP, SI);  ML = ML / norm(ML);

        Ov(k,:) = Om;
        Rv(:,:,k) = [ML(:), AP(:), SI(:)];

        if det(Rv(:,:,k)) < 0
            Rv(:,1,k) = -Rv(:,1,k);
        end
    end

    % --- Extract only the levels you care about (C7..S1 etc.) ---
    % Map “use levels” into the current frame’s valid levels list.
    mapUse = nan(nLev_use,1);
    for k = 1:nLev_use
        mapUse(k) = find(levels_valid == SpineLevels_use(k), 1);
    end
    if any(isnan(mapUse)), continue; end

    R_use = Rv(:,:,mapUse);

    % --- Thoracic and lumbar angles ---
    iC7  = find(SpineLevels_use=="C7",1);
    iT12 = find(SpineLevels_use=="T12",1);
    iS1  = find(SpineLevels_use=="S1",1);

    if ~isempty(iC7) && ~isempty(iT12)
        R_thor = R_use(:,:,iT12)' * R_use(:,:,iC7);  % C7 relative to T12
        [a,b,c] = cardanXYZ(R_thor);
        thoracic_deg(f,:) = rad2deg([a,b,c]);
    end

    if ~isempty(iT12) && ~isempty(iS1)
        R_lumb = R_use(:,:,iS1)' * R_use(:,:,iT12);  % T12 relative to S1
        [a,b,c] = cardanXYZ(R_lumb);
        lumbar_deg(f,:) = rad2deg([a,b,c]);
    end

    % --- Intersegmental angles (adjacent) ---
    for k = 1:nLev_use-1
        Rrel = R_use(:,:,k+1)' * R_use(:,:,k);
        [a,b,c] = cardanXYZ(Rrel);
        FE_int(f,k) = rad2deg(a);
        LB_int(f,k) = rad2deg(b);
        AT_int(f,k) = rad2deg(c);
    end
end

%%
f = 100;
fprintf('Frame %d thoracic [FE LB AT] = %.2f, %.2f, %.2f\n', f, thoracic_deg(f,:));
fprintf('Frame %d lumbar   [FE LB AT] = %.2f, %.2f, %.2f\n', f, lumbar_deg(f,:));


%% ===== Animation settings =====
frameStep = 20;         % speed
axLen = 30;            % axis length
showLevels = ["C7","T12","L1","L3","L5","S1"];  % change

figure('Color','w'); hold on; grid on; axis equal; view(3);
xlabel('X'); ylabel('Y'); zlabel('Z');
title('Spine_Extention01: splines + LCS animation');

% Pre-create plot objects
hPts = plot3(nan,nan,nan,'k.','MarkerSize',14);
hSL  = plot3(nan,nan,nan,'k-','LineWidth',1);
hSM  = plot3(nan,nan,nan,'g-','LineWidth',2);
hSR  = plot3(nan,nan,nan,'k-','LineWidth',1);

% quiver handles for X(ML), Y(AP), Z(SI)
hqX = gobjects(numel(showLevels),1);
hqY = gobjects(numel(showLevels),1);
hqZ = gobjects(numel(showLevels),1);
htL = gobjects(numel(showLevels),1);

for i = 1:numel(showLevels)
    hqX(i) = quiver3(nan,nan,nan,nan,nan,nan,0,'LineWidth',2);
    hqY(i) = quiver3(nan,nan,nan,nan,nan,nan,0,'LineWidth',2);
    hqZ(i) = quiver3(nan,nan,nan,nan,nan,nan,0,'LineWidth',2);
    htL(i) = text(nan,nan,nan,'','FontSize',9,'Interpreter','none');
end

for f = 1:frameStep:nFrames

    % ---- Build per-frame points (fast enough) ----
    P_L = nan(nLevels,3); P_M = nan(nLevels,3); P_R = nan(nLevels,3);
    for k = 1:nLevels
        if ~isempty(SpineCell{k,1}), P_L(k,:) = SpineCell{k,1}(f,:); end
        if ~isempty(SpineCell{k,2}), P_M(k,:) = SpineCell{k,2}(f,:); end
        if ~isempty(SpineCell{k,3}), P_R(k,:) = SpineCell{k,3}(f,:); end
    end
    valid = all(isfinite(P_L),2) & all(isfinite(P_M),2) & all(isfinite(P_R),2);
    P_Lv = P_L(valid,:); P_Mv = P_M(valid,:); P_Rv = P_R(valid,:);
    levels_valid = string(SpineLevels(valid));
    nLev = size(P_Mv,1);
    if nLev < 6, continue; end

    s = (1:nLev);
    nPieces_eff = min(6, nLev-1);
    uDense = linspace(1, nLev, 100);

    % splines
    spLx = spap2(nPieces_eff, 4, s, P_Lv(:,1)'); spLy = spap2(nPieces_eff, 4, s, P_Lv(:,2)'); spLz = spap2(nPieces_eff, 4, s, P_Lv(:,3)');
    spMx = spap2(nPieces_eff, 4, s, P_Mv(:,1)'); spMy = spap2(nPieces_eff, 4, s, P_Mv(:,2)'); spMz = spap2(nPieces_eff, 4, s, P_Mv(:,3)');
    spRx = spap2(nPieces_eff, 4, s, P_Rv(:,1)'); spRy = spap2(nPieces_eff, 4, s, P_Rv(:,2)'); spRz = spap2(nPieces_eff, 4, s, P_Rv(:,3)');

    S_L = [fnval(spLx,uDense); fnval(spLy,uDense); fnval(spLz,uDense)]';
    S_M = [fnval(spMx,uDense); fnval(spMy,uDense); fnval(spMz,uDense)]';
    S_R = [fnval(spRx,uDense); fnval(spRy,uDense); fnval(spRz,uDense)]';

    % update plots
    set(hSL,'XData',S_L(:,1),'YData',S_L(:,2),'ZData',S_L(:,3));
    set(hSM,'XData',S_M(:,1),'YData',S_M(:,2),'ZData',S_M(:,3));
    set(hSR,'XData',S_R(:,1),'YData',S_R(:,2),'ZData',S_R(:,3));

    % markers (all markers that exist in labels list)
    % (If you only want spine markers, plot P_L/P_M/P_R instead.)
    PM_all = [P_L(isfinite(P_L(:,1)),:); P_M(isfinite(P_M(:,1)),:); P_R(isfinite(P_R(:,1)),:)];
    set(hPts,'XData',PM_all(:,1),'YData',PM_all(:,2),'ZData',PM_all(:,3));

    % Build LCS only at showLevels
    dspMx = fnder(spMx); dspMy = fnder(spMy); dspMz = fnder(spMz);

    for i = 1:numel(showLevels)
        lev = showLevels(i);
        kLev = find(levels_valid == lev, 1);
        if isempty(kLev)
            set(hqX(i),'XData',nan,'YData',nan,'ZData',nan,'UData',nan,'VData',nan,'WData',nan);
            set(hqY(i),'XData',nan,'YData',nan,'ZData',nan,'UData',nan,'VData',nan,'WData',nan);
            set(hqZ(i),'XData',nan,'YData',nan,'ZData',nan,'UData',nan,'VData',nan,'WData',nan);
            set(htL(i),'Position',[nan nan nan],'String','');
            continue;
        end

        % nearest dense point for this level
        d = vecnorm(S_M - P_Mv(kLev,:), 2, 2);
        [~,ii] = min(d);
        uk = uDense(ii);

        Om = [fnval(spMx,uk), fnval(spMy,uk), fnval(spMz,uk)];
        Ol = [fnval(spLx,uk), fnval(spLy,uk), fnval(spLz,uk)];
        Or = [fnval(spRx,uk), fnval(spRy,uk), fnval(spRz,uk)];

        t  = [fnval(dspMx,uk), fnval(dspMy,uk), fnval(dspMz,uk)];
        SI = t / norm(t);
        LR = (Or - Ol); LR = LR / norm(LR);
        AP = cross(SI, LR); AP = AP / norm(AP);
        ML = cross(AP, SI); ML = ML / norm(ML);

        % quivers (X=ML, Y=AP, Z=SI)
        set(hqX(i),'XData',Om(1),'YData',Om(2),'ZData',Om(3),'UData',axLen*ML(1),'VData',axLen*ML(2),'WData',axLen*ML(3));
        set(hqY(i),'XData',Om(1),'YData',Om(2),'ZData',Om(3),'UData',axLen*AP(1),'VData',axLen*AP(2),'WData',axLen*AP(3));
        set(hqZ(i),'XData',Om(1),'YData',Om(2),'ZData',Om(3),'UData',axLen*SI(1),'VData',axLen*SI(2),'WData',axLen*SI(3));

        set(htL(i),'Position',Om,'String'," "+lev);
    end

    title(sprintf('Spine_Extention01: frame %d / %d', f, nFrames));
    drawnow;
end
