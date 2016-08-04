"""
This module orchestrates a two-point linkage analysis.
"""
module MendelTwoPointLinkage
#
# Required OpenMendel packages and modules.
#
using MendelBase
# using DataStructures                  # Now in MendelBase.
# using ModelConstruction               # Now in MendelBase.
# using ElstonStewartPreparation        # Now in MendelBase.
# using ElstonStewartEvaluation         # Now in MendelBase.
using Search
using SearchSetup
#
# Required external modules.
#
using DataFrames                        # From package DataFrames.

export TwoPointLinkage

"""
This is the wrapper function for the Two-Point Linkage analysis option.
"""
function TwoPointLinkage(control_file = ""; args...)

  const TWO_POINT_LINKAGE_VERSION :: VersionNumber = v"0.1.0"
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
  keyword = set_keyword_defaults!(Dict{ASCIIString, Any}())
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
  keyword::Dict{ASCIIString, Any})

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
  parameter = initialize_optimization(locus, parameter, keyword)
  #
  # Define a lod score data frame.
  #
  lodscore_frame = DataFrame(Trait = ASCIIString[], Marker = ASCIIString[],
    XXtheta = Float64[], XYtheta = Float64[], LodScore = Float64[])
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
    # Pass the variables to optimize for maximum likelihood estimation.
    #
    function fun(par)
      copy!(parameter.par, par)
      f = elston_stewart_loglikelihood(pedigree, person, locus, parameter, 
        instruction, keyword)
      return (f, nothing, nothing)
    end # function fun
    (best_par, best_value) = optimize(fun, parameter)
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
  writetable(keyword["lod_score_table"], lodscore_frame)
  show(lodscore_frame)
  return execution_error = false
end # function two_point_linkage_option

end # module MendelTwoPointLinkage

