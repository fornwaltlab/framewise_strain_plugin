classdef FramewiseDENSEanalysis < plugins.DENSEanalysisPlugin
    % FramewiseDENSEanalysis - A DENSEanalysis plugin
    %
    %   Plugin to compute strains for each frame independently
    %
    %   h/t to Jonathan Suever and Gregory Wehner for writing the initial
    %   implementation of this function
    %
    % Copyright (c) 2016, Cardiac Imaging Technology Lab

    methods
        function validate(~, data)
            % validate - Check if the plugin can run.
            %
            %   Performs validation to ensure that the state of the program
            %   is correct to be able to run the plugin.
            %
            % USAGE:
            %   FramewiseDENSEanalysis.validate(data)
            %
            % INPUTS:
            %   data:   Object, DENSEdata object containing all underlying
            %           data from the DENSEanalysis program.

            % Assert that image data base been loaded
            assert(~isempty(data.seq), ...
                'You must load imaging data into DENSEanalysis first.')
            assert(~isempty(data.spl),'You must Run Analysis first.')
        end

        function run(~, data)
            % run - Method executed when user selects the plugin
            %
            % USAGE:
            %   FramewiseDENSEanalysis.run(data)
            %
            % INPUTS:
            %   data:   Object, DENSEdata object containing all underlying
            %           data from the DENSEanalysis program.

            api.ResampleMethod = data.spl.ResampleMethod;
            api.Mag = data.spl.Mag;
            api.res = 1;
            api.Type = data.spl.ROIType;
            api.XunwrapPx = data.spl.Xunwrap * data.spl.Multipliers(1);
            api.YunwrapPx = data.spl.Yunwrap * data.spl.Multipliers(2);
            api.ZunwrapPx = data.spl.Zunwrap * data.spl.Multipliers(3);
            api.xrng = data.spl.jrng;
            api.yrng = data.spl.irng;
            api.Resolution = 1;
            api.FramesForAnalysis = data.spl.frrng;
            api.ResampleMethod = data.spl.ResampleMethod;
            api.PositionA = [];
            api.PositionB = [];

            %'PositionA',
            %'PositionB',
            %'Nmodel',
            %'Nseg',
            %'Clockwise',

            strains = [];
            FV = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%
            if data.spl.TemporalOrder>0 %if spline data have temporal smoothing
                didx = data.UIDtoIndexDENSE(data.spl.DENSEUID);
                ridx = data.UIDtoIndexROI(data.spl.ROIUID);
                imdata = data.imagedata(didx);
                cndata = data.contourdata(ridx);
                options = struct('FramesForAnalysis',cndata.ValidFrames);
                tagsA = {'ResampleMethod','SpatialSmoothing','TemporalOrder'};
                tagsB = {'Xseed','Yseed','Zseed'};
                for ti = 1:numel(tagsA)
                    options.(tagsA{ti}) = data.spl.(tagsA{ti});
                end
                options.TemporalOrder = -1; %override value

                dlast = UIDtoIndexDENSE(data,data.spl.DENSEUID);
                rlast = UIDtoIndexROI(data,data.spl.ROIUID);
                if isequal(dlast,didx) && isequal(rlast,ridx)
                    for ti = 1:numel(tagsB)
                        options.(tagsB{ti}) = data.spl.(tagsB{ti});
                    end
                    options.FramesForAnalysis = data.spl.frrng;
                end


                % 'UnwrapConnectivity' option based on ROI type
                if any(strcmpi(data.spl.ROIType,{'open','closed'}));
                    options.UnwrapConnectivity = 8;
                end
                spdata = DENSEspline(imdata,cndata,options,...
                'SeedFrame',data.spl.Xseed(1,3),'OptionsPanel',true);
                %see if I can pass seed points as varargin to make this
                %easier

                api.spldx = spdata.spldx;
                api.spldy = spdata.spldy;
                api.spldz = spdata.spldz;
            else
                api.spldx = data.spl.spldx;
                api.spldy = data.spl.spldy;
                api.spldz = data.spl.spldz;
             end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % These are the contours at each frame drawn by the user
            C       = data.spl.Contour;
            C0      = data.spl.Contour;
            frrng   = data.spl.frrng;
            method  = data.spl.ResampleMethod;
            xrng    = data.spl.jrng;
            yrng    = data.spl.irng;

            % Prepare this structure to use as input to patch calculator
            inputs = api;

            % Image grid generation
            xi = xrng(1):xrng(2);
            yi = yrng(1):yrng(2);
            %    [Xi, Yi] = meshgrid(xi, yi);

            % C0 contains the configuration of the contours at time=0 using each
            % contour as a reference geometry
            for fr = frrng(1):frrng(2)

                % current point locations & current-to-origin deformation from
                % Eulerian phase data
                Xphase = api.XunwrapPx(:,:,fr);
                Yphase = api.YunwrapPx(:,:,fr);
                tmp = find(~isnan(Xphase));
                [tmpYpts,tmpXpts] = ind2sub(size(Xphase),tmp);
                pts  = [tmpXpts,tmpYpts];
                dpts = [-Xphase(tmp),-Yphase(tmp)];

                % THIN PLATE SPLINE
                if isequal(method,'tpaps')

                    % spline current tissue configuration to
                    % current-to-origin deformation
                    st = tpaps_v3(pts',dpts',1);

                    % project contours back to origin
                    for k = 1:size(C,2)
                        dxy = fnval(st,C0{fr,k}');
                        C0{fr,k} = C0{fr,k} + dxy';
                    end

                % GRIDFIT
                else

                    % resample current-to-origin deformation on
                    % uniform grid of current positions
                    sm = 0.1;
                    dxfr = gridfit(pts(:,1),pts(:,2),dpts(:,1),xi,yi,...
                    'extend','always','smoothness',sm);
                    dyfr = gridfit(pts(:,1),pts(:,2),dpts(:,2),xi,yi,...
                    'extend','always','smoothness',sm);

                    % spline current tissue configuration to
                    % current-to-origin deformation
                    ppx = spapi({2,2},{yi,xi},dxfr);
                    ppy = spapi({2,2},{yi,xi},dyfr);

                    % project contours back to origin
                    for k = 1:size(C,2)
                        dxk = fnval(ppx,C0{fr,k}(:,[2 1])');
                        dyk = fnval(ppy,C0{fr,k}(:,[2 1])');
                        C0{fr,k} = C0{fr,k} + [dxk(:),dyk(:)];
                    end

                end

            % If the user doesn't provide any position information, bring up dialog
                if isempty(api.PositionA) && fr == 1 && ~strcmp(api.Type,'curve')
                    inputs.CardiacModelPanel = true;
                else
                    inputs.CardiacModelPanel = false;
                end

                % Update the resting contour to the current frame resting contour
                inputs.RestingContour   = C0(fr,:);
                inputs.Type             = api.Type;
                inputs.Mag              = api.Mag;

                %generate mask for contours
                [X,Y] = meshgrid(xi,yi);
                mask = data.spl.MaskFcn(X,Y,inputs.RestingContour);

            % Evaluate and store as inputs in case the user had to use cardiacmodel
                inputs = spl2patch(inputs);
                inputs.mask=mask;
                % Store the resting meshes in a structure
                FV = cat(1, FV, inputs);

                % Now call patch strain
                straininputs = struct('vertices',   inputs.vertices,...
                                      'faces',      inputs.faces,...
                                      'times',      fr,...
                                      'spldx',      api.spldx,...
                                      'spldy',      api.spldy,...
                                      'spldz',      api.spldz,...
                                      'Orientation',inputs.orientation);

                if any(strcmpi(api.Type,{'open','closed'}))
                    str = contourstrain(straininputs);
                else
                    str = patchstrain(straininputs);
                end
                strains = cat(1,strains,str);

            end

            strains(1).Nmodel = inputs.Nmodel;
            strains(1).PositionA = inputs.PositionA;
            strains(1).PositionB = inputs.PositionB;
            strains(1).Clockwise = inputs.Clockwise;
            exportFcn(data,strains,FV);

        end
    end
end

function exportFcn(data, strain,fv)
%copied from AnalysisViewer to maintain consistency of output

    startpath = pwd;
    obj.hdata = data; %for handle consistency
    obj.straindata = strain;


    % DENSE & ROI indices
    duid = obj.hdata.spl.DENSEUID;
    didx = obj.hdata.UIDtoIndexDENSE(duid);

    ruid = obj.hdata.spl.ROIUID;
    ridx = obj.hdata.UIDtoIndexROI(ruid);

    % file name
    header = sprintf('%s_%s',...
        obj.hdata.dns(didx).Name,...
        obj.hdata.roi(ridx).Name);

    expr = '[\\\/\?\%\:\*\"\<\>\|]';
    header = regexprep(header,expr,'_');

    file = fullfile(startpath,[header '.mat']);
    cnt = 0;
    while isfile(file)
        cnt = cnt+1;
        file = fullfile(startpath,sprintf('%s (%d).mat',header,cnt));
    end

    % allow user to change file name
    [uifile,uipath] = uiputfile('*.mat',[],file);
    if isequal(uifile,0)
        file = []; %#ok<NASGU>
        return;
    end

    % check extension
    file = fullfile(uipath,uifile);
    [~,~,e] = fileparts(file);
    if ~isequal(e,'.mat')
        file = [file, '.mat'];
    end

    % start waitbartimer
    hwait = waitbartimer;
    cleanupObj = onCleanup(@()delete(hwait(isvalid(hwait))));
    hwait.String = 'Saving MAT file...';
    hwait.WindowStyle = 'modal';
    hwait.AllowClose  = false;
    start(hwait);
    drawnow

    % magnitude/phase indices
    midx = obj.hdata.dns(didx).MagIndex;
    pidx = obj.hdata.dns(didx).PhaIndex;
    sidx = [midx; pidx];

    % image information
    tags = {'DENSEType','Multipliers','Mag'};
    tfsingle = logical([0 0 0]);
    if ~isnan(pidx(1))
        tags = cat(2,tags,{'Xwrap','Xunwrap'});
        tfsingle = [tfsingle, true, true];
    end
    if ~isnan(pidx(2))
        tags = cat(2,tags,{'Ywrap','Yunwrap'});
        tfsingle = [tfsingle, true, true];
    end
    if ~isnan(pidx(3))
        tags = cat(2,tags,{'Zwrap','Zunwrap'});
        tfsingle = [tfsingle, true, true];
    end

    ImageInfo = struct;
    for ti = 1:numel(tags)
        tag = tags{ti};
        if tfsingle(ti)
            ImageInfo.(tag) = single(obj.hdata.spl.(tag));
        else
            ImageInfo.(tag) = obj.hdata.spl.(tag);
        end
    end

    % ROI information
    tags = {'ROIType','RestingContour','Contour'};
    ROIInfo = struct;
    for ti = 1:numel(tags)
        ROIInfo.(tags{ti}) = obj.hdata.spl.(tags{ti});
    end

    % analysis information
    tags = {'ResampleMethod','ResampleDistance','SpatialSmoothing',...
        'TemporalOrder','Xseed','Yseed','Zseed'};
    AnalysisInfo = struct;
    for ti = 1:numel(tags)
        AnalysisInfo.(tags{ti}) = obj.hdata.spl.(tags{ti});
    end

    % frame range analyized
    frrng  = obj.hdata.spl.frrng;
    frames = frrng(1):frrng(2);
    AnalysisInfo.FramesForAnalysis = frrng;

    % segment model & orientation
    if any(strcmpi(obj.hdata.spl.ROIType,{'sa','la'}))
        AnalysisInfo.Nmodel    = obj.straindata.Nmodel;
        AnalysisInfo.PositionA = obj.straindata.PositionA;
        AnalysisInfo.PositionB = obj.straindata.PositionB;
        AnalysisInfo.Clockwise = obj.straindata.Clockwise;
    end

    % DENSE group information
    DENSEInfo = obj.hdata.dns(didx);

    % original sequence header information
    SequenceInfo = repmat(struct,[2 3]);
    for k = 1:6
        if ~isnan(sidx(k))
            tags = fieldnames(obj.hdata.seq(sidx(k)));
            for ti = 1:numel(tags)
                tag = tags{ti};
                SequenceInfo(k).(tag) = ...
                    obj.hdata.seq(sidx(k)).(tag);
            end
        end
    end

    % Displacement Info
    Isz = size(obj.hdata.spl.Xwrap(:,:,1));
    Nfr = size(obj.hdata.spl.Xwrap,3);

    x = 1:Isz(2);
    y = 1:Isz(1);
    [X,Y] = meshgrid(x,y,0);

    mask0 = obj.hdata.spl.MaskFcn(...
        X,Y,obj.hdata.spl.RestingContour);
    Npts = sum(mask0(:));

    DisplacementInfo = struct(...
        'X',    X(mask0),...
        'Y',    Y(mask0),...
        'dX',   NaN([Npts,Nfr]),...
        'dY',   NaN([Npts,Nfr]),...
        'dZ',   NaN([Npts,Nfr]));

    pts = [Y(:),X(:),zeros(size(X(:)))];
    pts = pts(mask0,:)';

    for fr = frames
        pts(3,:) = fr;
        DisplacementInfo.dX(:,fr) = fnvalmod(obj.hdata.spl.spldx,pts);
        DisplacementInfo.dY(:,fr) = fnvalmod(obj.hdata.spl.spldy,pts);
        DisplacementInfo.dZ(:,fr) = fnvalmod(obj.hdata.spl.spldz,pts);
    end

    % short-axis angle
    % for 6-segment model, [0,60)=anterior, [60,120)=anteroseptal, etc.
    % for 4-segment model, [0,90)=anterior, [90 180)=septal, etc.
    if strcmpi(obj.hdata.spl.ROIType,'sa')
        origin  = obj.straindata.PositionA;
        posB    = obj.straindata.PositionB;
        flag_clockwise = obj.straindata.Clockwise;

        theta0 = atan2(posB(2)-origin(2),posB(1)-origin(1));
        theta  = atan2(Y(mask0)-origin(2),X(mask0)-origin(1)) - theta0;
        if ~flag_clockwise, theta = -theta; end

        theta(theta<0) = theta(theta<0) + 2*pi;
        theta(theta>=2*pi) = 0;

        DisplacementInfo.Angle = theta(:);
    end

    % export file
    AnalysisInstanceUID = dicomuid;

    out.ImageInfo = ImageInfo;
    out.ROIInfo = ROIInfo;
    out.DisplacementInfo = DisplacementInfo;
    out.AnalysisInfo = AnalysisInfo;
    out.DENSEInfo = DENSEInfo;
    out.SequenceInfo = SequenceInfo;
    out.StrainInfo = obj.straindata;
    out.MaskInfo = fv;
    out.AnalysisInstanceUID = AnalysisInstanceUID; %#ok<STRNU>

    save(file, '-struct', 'out')
end

