@with_kw struct fitDTOF
   
    #required inputs
    t::Union{AbstractRange{Float64}, AbstractVector{Float64}}
    DTOF::Array{Float64,1}
    IRF::Array{Float64,1}

    nmed::Float64 = 1.5 
	ndet::Float64 = 1.51
	lambda::Float64 = 750.0
     
    #semi-infinite and slab
    ρ::Union{Float64, Missing} = missing

    #slab
    s::Union{Float64, Missing} = missing

    #parallelepiped
    rd::Union{Array{Float64,1}, Missing} = missing
    rs::Union{Array{Float64,1}, Missing} = missing
    L::Union{Array{Float64,1}, Missing} = missing
end



@with_kw struct fitg2
   
    #required inputs
    τ::Union{AbstractRange{Float64}, AbstractVector{Float64}}
    g2::Array{Float64,1}
     
    #semi-infinite
    ρ::Union{Float64, Missing}
end




@with_kw struct DTOF_fitparams
    model::Function = TPSF_DA_semiinf_refl          # model to be fitted against
	initparams::Array{Float64,1} = [0.2, 10.0]      # initial starting point of LsqFit
	lb::Array{Float64,1} = [0.001, 1.0]             # lowerbound of fit for [mua, musp]
	ub::Array{Float64,1} = [1.0, 20.0]              # upperbound of fit for [mua, musp]
	risefactor::Float64  = 0.9			            # percent of peak to fit on rising edge
	tailfactor::Float64	 = 0.1			            # percent of peak to fit on falling tail
end

@with_kw struct g2_fitparams
    model::Function = g2_DA_SI_R                    # model to be fitted against
	initparams::Array{Float64,1} = [2.e-8, 0.5]     # initial starting point of LsqFit
	lb::Array{Float64,1} = [1.e-12, 0.0]            # lowerbound of fit for [BFI, β]
	ub::Array{Float64,1} = [5.e-6, 1.0]             # upperbound of fit for [BFI, β]
end

struct fitresult
    #Fit inputs
    xdata::Array{Float64,1}     #independent variable
    ydata::Array{Float64,1}     #dependent variable
    yerrors::Array{Float64,1}   #uncertainties in dependent variable
	model::Function             #model to be fitted against

    #Fit results
    residuals::Array{Float64,1}  	 #weighted residuals
    param::Array{Float64,1}    		#best fit parameters for model
	marginerror::Array{Float64,1}        #product of standard error and the critical value of each parameter at a certain significance level (default is 5%) from t-distribution.
    perrors::Array{Float64,1}            #uncertainties in best fit parameters
    dof::Int64                  		#number of degrees of freedom in fit
    chisq::Float64             		 #chi^2 value of fit
    probability::Float64      		  #probability of chi^2 value
end
