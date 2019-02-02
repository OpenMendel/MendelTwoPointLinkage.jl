__precompile__()

"""
This module orchestrates a two-point linkage analysis.
"""
module MendelTwoPointLinkage
#
# Required OpenMendel packages and modules.
#
using MendelBase
# namely: DataStructures, ModelConstruction,
# ElstonStewartPreparation, ElstonStewartEvaluation
using MendelSearch
#
# Required external modules.
#
using CSV
using DataFrames
using LinearAlgebra

export TwoPointLinkage

"""
This is the wrapper function for the Two-Point Linkage analysis option.
"""
function TwoPointLinkage(control_file = ""; args...)

  TWO_POINT_LINKAGE_VERSION :: VersionNumber = v"0.5.0"
  #
  # Print the logo. Store the initial directory.
  #
  print(" \n \n")
  println("     Welcome to OpenMendel's")
  println("Two-Point Linkage analysis option")
  println("         version ", TWO_POINT_LINKAGE_VERSION)
  print(" \n \n")
  println("Reading the data.\n")
  initial_directory = pwd()
  #
  # The user specifies the analysis to perform via a set of keywords.
  # Start the keywords at their default values.
  #
  keyword = set_keyword_defaults!(Dict{AbstractString, Any}())
  #
  # Keywords unique to this analysis should be first defined here
  # by setting their default values using the format:
  # keyword["some_keyword_name"] = default_value
  #
  keyword["gender_neutral"] = true
  keyword["lod_score_table"] = "Lod_Score_Frame.txt"
  #
  # Process the run-time user-specified keywords that will control the analysis.
  # This will also initialize the random number generator.
  #
  process_keywords!(keyword, control_file, args)
  #
  # Check that the correct analysis option was specified.
  #
  lc_analysis_option = lowercase(keyword["analysis_option"])
  if (lc_analysis_option != "" &&
      lc_analysis_option != "twopointlinkage")
     throw(ArgumentError(
       "An incorrect analysis option was specified.\n \n"))
  end
  keyword["analysis_option"] = "TwoPointLinkage"
  #
  # Read the genetic data from the external files named in the keywords.
  #
  (pedigree, person, nuclear_family, locus, snpdata,
    locus_frame, phenotype_frame, pedigree_frame, snp_definition_frame) =
    read_external_data_files(keyword)
  #
  # Execute the specified analysis.
  #
  println(" \nAnalyzing the data.\n")
  execution_error = two_point_linkage_option(pedigree, person, nuclear_family,
    locus, locus_frame, phenotype_frame, pedigree_frame, keyword)
  if execution_error
    println(" \n \nERROR: Mendel terminated prematurely!\n")
  else
    println(" \n \nMendel's analysis is finished.\n")
  end
  #
  # Finish up by closing, and thus flushing, any output files.
  # Return to the initial directory.
  #
  close(keyword["output_unit"])
  cd(initial_directory)
  return nothing
end # function TwoPointLinkage

"""
This function maps a trait locus by linkage.
"""
function two_point_linkage_option(pedigree::Pedigree, person::Person,
  nuclear_family::NuclearFamily, locus::Locus, locus_frame::DataFrame,
  phenotype_frame::DataFrame, pedigree_frame::DataFrame,
  keyword::Dict{AbstractString, Any})

  io = keyword["output_unit"]
  #
  # Problem formulation checks.
  #
  if locus.loci <= 1
    println("Error: This option requires at least two loci.")
    println(io, "Error: This option requires at least two loci.")
    return execution_error = true
  end
  trait_locus = 0
  for l = 1:locus.loci
    if locus.name[l] == keyword["trait"]
      trait_locus = l
      break
    end
  end
  if trait_locus == 0
    println("Error: No match to the trait among the possible loci.")
    println(io, "Error:  No match to the trait among the possible loci.")
    return execution_error = true
  end
  #
  # Prepare for likelihood analysis.
  #
  keyword["eliminate_genotypes"] = true
  keyword["lump_alleles"] = true
  locus.model_loci = 2
  locus.model_locus = [trait_locus, 0]
  if keyword["travel"] == "grid"
    keyword["parameters"] = 1
    keyword["points"] = 9
  else
    if keyword["gender_neutral"]
      keyword["parameters"] = 1
    else
      keyword["parameters"] = 2
    end
  end
  keyword["goal"] = "maximize"
  #
  # Define the parameter data structure.
  #
  parameter = set_parameter_defaults(keyword)
  parameter =
    initialize_optimization_two_point_linkage!(locus, parameter, keyword)
  #
  # Define a lod score data frame.
  #
  lodscore_frame = DataFrame(Trait = AbstractString[],
    Marker = AbstractString[], XXtheta = Float64[],
    XYtheta = Float64[], LodScore = Float64[])
  #
  # Execute two-point linkage analysis over all trait-marker pairs.
  #
  for loc = 1:locus.loci
    if loc == trait_locus; continue; end
    locus.model_locus[2] = loc
    parameter.title =
      "Two point linkage analysis between " * locus.name[trait_locus] *
      " and " * locus.name[loc]
    locus.theta = model_recombination_fractions(locus, keyword)
    #
    # Fetch the instructions for conducting the Elston-Stewart algorithm.
    #
    (instruction, elston_stewart_count) = orchestrate_likelihood(pedigree,
      person, nuclear_family, locus, keyword)
    if any(elston_stewart_count .>  keyword["complexity_threshold"])
      println(io, "Marker $locus.name[loc] was skipped because one or more ",
        "pedigrees exceeds the complexity threshold.")
      continue
    end
    #
    # Pass the variables to search for maximum likelihood estimation.
    #
    function fun(par)
      copyto!(parameter.par, par)
      f = elston_stewart_loglikelihood(penetrance_two_point_linkage,
        prior_two_point_linkage, transmission_two_point_linkage,
        pedigree, person, locus, parameter, instruction, keyword)
      return (f, nothing, nothing)
    end # function fun
    (best_par, best_value) = mendel_search(fun, parameter)
    #
    # Insert the results in a data frame.
    #
    if parameter.travel == "search"
      lod = parameter.function_value[end] - parameter.function_value[1]
      lod = log10(exp(1.0)) * lod
      if parameter.parameters == 1
        push!(lodscore_frame, [locus.name[trait_locus], locus.name[loc],
        parameter.par[1], parameter.par[1], lod])
      else
        push!(lodscore_frame, [locus.name[trait_locus], locus.name[loc],
        parameter.par[1], parameter.par[2], lod])
      end
    else
      for i = 1:parameter.points
        lod = parameter.function_value[i] - parameter.function_value[1]
        lod = log10(exp(1.0)) * lod
        push!(lodscore_frame, [locus.name[trait_locus], locus.name[loc],
        parameter.grid[i, 1], parameter.grid[i, 1], lod])
      end
    end
  end
  lod_table_file = string(keyword["lod_score_table"])
  CSV.write(lod_table_file, lodscore_frame;
    writeheader = true, delim = keyword["output_field_separator"],
    missingstring = keyword["output_missing_value"])
  show(lodscore_frame)
  return execution_error = false
end # function two_point_linkage_option

"""
Supply a penetrance for individual i.
"""
function penetrance_two_point_linkage(person::Person, locus::Locus,
  multi_genotype::Matrix{Int}, par::Vector{Float64},
  keyword::Dict{AbstractString, Any}, start::Int, finish::Int, i::Int)

  pen = 1.0
  for l = start:finish
    allele1 = multi_genotype[1, l]
    allele2 = multi_genotype[2, l]
    loc = locus.model_locus[l]
    p = 1.0 # for reduced penetrance let p depend on loc
    pen = p * pen
  end
  return pen
end # function penetrance_two_point_linkage

"""
Supply a prior probability for founder i.
"""
function prior_two_point_linkage(person::Person, locus::Locus,
  multi_genotype::Matrix{Int}, par::Vector{Float64},
  keyword::Dict{AbstractString, Any}, start::Int, finish::Int, i::Int)

  prior_prob = 1.0
  for l = start:finish
    loc = locus.model_locus[l]
    allele = multi_genotype[1, l]
    frequency = dot(vec(person.admixture[i, :]),
                    vec(locus.frequency[loc][:, allele]))
    prior_prob = prior_prob * frequency
    if !locus.xlinked[loc] || !person.male[i]
      allele = multi_genotype[2, l]
      frequency = dot(vec(person.admixture[i, :]),
                      vec(locus.frequency[loc][:, allele]))
      prior_prob = prior_prob * frequency
    end
  end
  return prior_prob
end # function prior_two_point_linkage

"""
Supply the transmission probability that a parent i with a particular
genotype transmits a particular gamete to his or her child j.
"""
function transmission_two_point_linkage(person::Person, locus::Locus,
  gamete::Vector{Int}, multi_genotype::Matrix{Int}, par::Vector{Float64},
  keyword::Dict{AbstractString, Any}, start::Int, finish::Int, i::Int, j::Int)
  #
  # For male to male inheritance at an x-linked locus,
  # set the transmission probability equal to 1.
  #
  loc = locus.model_locus[start]
  xlinked = locus.xlinked[loc]
  if xlinked && person.male[i] && person.male[j]
    return 1.0
  end
  #
  # Equate recombination fractions to parameters for two-point linkage
  # analysis.
  #
  if length(par) == 2
    locus.theta[:, 1] = par
  elseif length(par) == 1
    locus.theta[:, 1] .= par[1]
  end
  #
  # Store an indicator of the sex of the parent.
  #
  if person.male[i]
    i_sex = 2
  else
    i_sex = 1
  end
  #
  # Reduce the computations by considering only the heterozygous loci.
  # Use Trow's formula to express the recombination fraction
  # between two heterozygous loci in terms of the recombination
  # fractions between the adjacent loci that separate them.
  # Set the logical variable found to true when the first heterozygous
  # parental locus is found. Phase records the phase of the most
  # recent heterozygous parental locus.
  #
  trans = 1.0
  found = false
  phase = true
  r = 0.5
  for l = start:finish
    match1 = multi_genotype[1, l] == gamete[l]
    match2 = multi_genotype[2, l] == gamete[l]
    #
    # Check whether either the first or second parental gene at
    # the current locus matches the gamete gene at this locus.
    # If not, then return with 0 for the transmission probability.
    #
    if !match1 && !match2
      return 0.0
    end
    #
    # Check whether the current locus is heterozygous.
    #
    if match1 != match2
      if found
        if phase == match1
          trans = trans * (0.5 + r)
        else
          trans = trans * (0.5 - r)
        end
      else
        found = true
        if start == 1 || start == finish
          trans = 0.5
        else
          trans = 1.0
        end
      end
      phase = match1
      r = 0.5
    end
    if found && l < finish
      r = r * (1.0 - 2.0 * locus.theta[i_sex, l])
    end
  end
  if !found; trans = 1.0; end
  return trans
end # function transmission_two_point_linkage

"""
Initialize the optimization problem.
"""
function initialize_optimization_two_point_linkage!(locus::Locus,
  parameter::Parameter, keyword::Dict{AbstractString, Any})
  #
  # Initialize, bound, and name the parameters.
  #
  fill!(parameter.name, "theta")
  if parameter.travel == "grid"
    parameter.grid[:, 1] = [0.5, 0.4, 0.3, 0.2, 0.15, 0.1, 0.05, 0.01, 0.001]
  else
    fill!(parameter.par, 0.5 - 1e-5)
    fill!(parameter.min, 1e-5)
    fill!(parameter.max, 0.5)
    if parameter.parameters == 2
      parameter.name[1] = "xxtheta"
      parameter.name[2] = "xytheta"
    end
  end
  return parameter
end # function initialize_optimization_two_point_linkage!
#
# Method to obtain path to this package's data files
# so they can be used in the documentation and testing routines.
# For example, datadir("Control file.txt") will return
# "/path/to/package/data/Control file.txt"
#
datadir(parts...) = joinpath(@__DIR__, "..", "data", parts...)

end # module MendelTwoPointLinkage
