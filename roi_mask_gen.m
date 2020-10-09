function [mask,mask_id] = roi_mask_gen(name,pix_length,hmap,map,output_path_save)
%ROI_MASK_GEN runs a MATLAB App to generate a "mask" for useless regions in
% the height map, based on user input.
    outputname = output_path_save+name+"_"+pix_length+"pix_mask.mat";
    if ~isfile(outputname)
        % Create output mask!
        obj = roi_app_test(name,hmap,map,outputname);
        waitfor(obj);
    end
    load(outputname,'mask','mask_id');
end