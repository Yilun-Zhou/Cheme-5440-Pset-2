function calculate_transcription_kinetics(t::Float64, x::Array{Float64,1}, problem::Dict{String,Any})::Float64

    # initailzie -
    transcription_rate = 0.0

    # alias -
    mRNA = x[1]
    G = x[2]
    σ70 = x[3]

    # get parameters from problem dictionary -
    toml_dictionary = TOML.parsefile(path_to_parameters_file)["biophysical_constants"]
    KX = problem["transcription_saturation_constant"]
    tau_X = problem["transcription_time_constant"]
    Vmax_X = problem["maximum_transcription_velocity"]
    RNAPII_concentration =toml_dictionary["RNAPII_concentration"]
    transcription_time_constant =toml_dictionary["transcription_time_constant"]
    transcription_saturation_constant =toml_dictionary["transcription_saturation_constant"]
    # TODO: compute the transcription rate -
    transcription_rate = Vmax_X*(RNAPII_concentration/(transcription_time_constant*transcription_saturation_constant+(1+transcription_time_constant)*RNAPII_concentration))

    # return -
    return transcription_rate
end

function calculate_mRNA_degradation_kinetics(t::Float64, x::Array{Float64,1}, problem::Dict{String,Any})::Float64

    # alias -
    mRNA = x[1]
    G = x[2]
    σ70 = x[3]

    # get the mRNA degradation constant -
    mRNA_degradation_constant = problem["mRNA_degradation_constant"]

    # return -
    return (mRNA_degradation_constant*mRNA)
end
