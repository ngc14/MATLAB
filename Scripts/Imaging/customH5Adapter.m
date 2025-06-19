%% STEP 1: INHERIT FROM ADAPTER CLASS
classdef customH5Adapter < images.blocked.Adapter
    properties(Access=public)
        File (1,1) string
        Ds (1,1) string
        Size (1,:)
        IOBlockSize (1,:)
    end
    properties(Access=private)
        FID
        FILEA
    end
    properties(Hidden)
        Info struct
    end
    properties(Constant, Hidden)
        classNames = ["H5T_INTEGER", "H5T_FLOAT","H5T_STRING","H5T_BITFIELD", ...
            "H5T_OPAQUE","H5T_COMPOUND","H5T_ENUM", "H5T_VLEN","H5T_ARRAY"];
    end
    methods
        function obj = customH5Adapter(ds)
            obj.Ds = "/"+ds;
        end
        % Define the openToRead method
        function openToRead(obj,source)
            filea = H5P.create('H5P_FILE_ACCESS');
            H5P.set_fapl_sec2(filea);
            H5P.set_fclose_degree(filea,'H5F_CLOSE_STRONG');
            H5P.set_libver_bounds(filea,'H5F_LIBVER_LATEST','H5F_LIBVER_LATEST');
            obj.FID = H5F.open(source,"H5F_ACC_RDONLY|H5F_ACC_SWMR_READ",filea);
            obj.File = source;
            obj.FILEA = filea;
        end
        % Define the getInfo method
        function info = getInfo(obj)
            if(H5L.exists(obj.FID,obj.Ds,"H5P_DEFAULT"))
                openedDS = H5D.open(obj.FID,obj.Ds);
                openedSP = H5D.get_space(openedDS);
                [~,obj.Size] = H5S.get_simple_extent_dims(openedSP);
                obj.Size = fliplr(obj.Size);
                % read all trials of a given frame
                obj.IOBlockSize = [obj.Size(1:2) 1 obj.Size(end)];
                typID = H5T.get_class(H5D.get_type(openedDS));
                info.Datatype = "double";%obj.classNames(typID+1);
                info.IOBlockSize = obj.IOBlockSize;
                info.InitialValue = ones(1,numel(obj.Size));
                info.Size = obj.Size;
                H5S.close(openedSP);
                H5D.close(openedDS);
            else
                info = [];
            end
        end

        % Define the getIOBlock method
        function block = getIOBlock(obj,f,level)
            oF = H5D.open(obj.FID,obj.Ds);
            H5D.refresh(oF);
            dOF = H5D.get_space(oF);
            H5S.select_hyperslab(dOF,"H5S_SELECT_SET",fliplr(f-1),[],[],...
                fliplr(obj.IOBlockSize));
            block = H5D.read(oF,'H5ML_DEFAULT',H5S.create_simple(numel(obj.Size),...
                fliplr(obj.IOBlockSize),fliplr(obj.Size)),dOF,'H5P_DEFAULT');
            H5S.close(dOF);
            H5D.close(oF);
        end
        function close(obj)
            H5P.close(obj.FILEA);
            H5F.close(obj.FID);
            waitfor(H5F.get_obj_count(obj.FID,1));
        end
    end
end
