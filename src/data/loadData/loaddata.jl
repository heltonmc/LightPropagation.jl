
#filename1 = "/home/heltonmc/Desktop/DTOF.asc"
#filename = "/home/heltonmc/Desktop/IRF.asc"


t, IRF, DTOF = load_asc_data(filename, filename1)

function load_asc_data(IRFfilename, DTOFfilename)

    #read .asc data and skip first 10 lines
    data_array = readdlm(IRFfilename,skipstart=10)
    data_array1 = readdlm(DTOFfilename,skipstart=10)

    #separate counts and time columns (ignore last line)
    RtIRF = data_array[1:end-1,2]
    RtIRF = RtIRF./maximum(RtIRF)
    tIRF = data_array[1:end-1,1]

    RtDTOF = data_array1[1:end-1,2]
    RtDTOF = RtDTOF./maximum(RtDTOF)
    tDTOF = data_array1[1:end-1,1]

    @assert tDTOF == tIRF

    return tDTOF, RtIRF, RtDTOF
end




