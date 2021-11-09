import HDF5

function get_reflection_data(file::String)
    fid = HDF5.h5open(file, "r")
    angle = HDF5.read(fid["angle"])
    energy = HDF5.read(fid["energy"])
    reflectivity = HDF5.read(fid["reflectivity/ir"])
    close(fid)

    keydata = cat(angle, reflectivity, dims=2)
    keydata = cat([0.0 energy'], keydata, dims=1)

    return keydata
end