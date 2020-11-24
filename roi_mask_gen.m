function [mask,mask_id] = roi_mask_gen(name,pix_length,hmap,map,pix2m,isroof,output_path_save,allow_custom,force_mask)
%ROI_MASK_GEN runs a MATLAB App to generate a "mask" for useless regions in
% the height map, based on user input.
    outputname = output_path_save+name+"_"+pix_length+"pix_mask.mat";
    if ~isfile(outputname) || force_mask
        if allow_custom
            % Create output mask!
            obj = roi_app_test(name,hmap,map,pix2m,isroof,outputname,allow_custom);
            waitfor(obj);
        else
            % Get mask automatically
            roi_app_noUI(name,hmap,map,pix2m,isroof,outputname,allow_custom);
        end
    end
    load(outputname,'mask','mask_id');
end