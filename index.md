### Overview
Mendel Two Point Linkage is a component of the umbrella [OpenMendel](https://openmendel.github.io) project. This analysis option maps a trait locus using linkage analysis.

<!--- ### Appropriate Problems and Data Sets
Genetic distance is proportional to the expected number of recombination events per meiosis separating two loci.
... --->

### Installation
*Note: The three OpenMendel packages (1) [SnpArrays](https://openmendel.github.io/SnpArrays.jl/latest/), (2) [Search](https://openmendel.github.io/Search.jl), and (3) [MendelBase](https://openmendel.github.io/MendelBase.jl) must be installed before any other OpenMendel package will run. It is easiest if these three packages are installed in the above order and before any other OpenMendel package.*

Within Julia, use the package manager to install MendelTwoPointLinkage:

    Pkg.clone("https://github.com/OpenMendel/MendelTwoPointLinkage.jl.git")

This package supports Julia v0.4 and v0.5.

### Input Files
The Mendel Two Point Linkage analysis package uses the following input files. Example input files can be found in the [docs]( https://github.com/OpenMendel/MendelTwoPointLinkage.jl/tree/master/docs) subfolder of the Mendel Two Point Linkage project. (An analysis won't always need every file type below.)

* [Control File](#control-file): Specifies the names of your data input and output files and any optional parameters (*keywords*) for the analysis. (For a list of common keywords, see [Keywords Table](#keywords-table)).
* [Locus File](https://openmendel.github.io/MendelBase.jl/#locus-file): Names and describes the genetic loci in your data.
* [Pedigree File](https://openmendel.github.io/MendelBase.jl/#pedigree-file): Gives information about your individuals, such as name, sex, family structure, and ancestry.
* [Phenotype File](https://openmendel.github.io/MendelBase.jl/#phenotype-file): Lists the available phenotypes.
* [SNP Definition File](#snp-definition-file): Defines your SNPs with information such as SNP name, chromosome, position, allele names, allele frequencies.
* [SNP Data File](https://openmendel.github.io/MendelBase.jl/#snp-data-file): Holds the genotypes for your data set. Must be a standard binary PLINK BED file in SNP major format. If you have a SNP data file you must have a SNP definition file.

### Control file<a id="control-file"></a>
The Control file is a text file consisting of keywords and their assigned values. The format of the Control file is:

	Keyword = Keyword_Value(s)

Below is an example of a simple Control file to run Two Point Linkage:

	#
	# Input and Output files.
	#
	locus_file = two-point linkage LocusFrame.txt
	pedigree_file = two-point linkage PedigreeFrame.txt
	phenotype_file = two-point linkage PhenotypeFrame.txt
	output_file = two-point linkage Output.txt
	lod-score-table = two-point linkage LOD Table Output.txt
	#
	# Analysis parameters for Two-Point Linkage option.
	#
	disease_status = RADIN
	gender-neutral = true
	standard_errors = true
	travel = grid

In the example above, there are nine keywords. The first three keywords specify input files: *two-point linkage LocusFrame.txt*, *two-point linkage PedigreeFrame.txt*, *two-point linkage PhenotypeFrame.txt*. The next two keywords specify output files with results of the analysis: *two-point linkage Output.txt* and *two-point linkage LOD Table Output.txt*. The last four keywords specify analysis parameters: *disease_status*, *gender-neutral*, *standard_errors* and *travel*. The text after the '=' are the keyword values.

### Keywords<a id="keywords-table"></a>
This is a list of OpenMendel keywords specific to Two Point Linkage. A list of OpenMendel keywords common to most analysis package can be found [here](https://openmendel.github.io/MendelBase.jl/#keywords-table). The names of keywords are *not* case sensitive. (The keyword values *may* be case sensitive.)

 Keyword          |   Default Value    | Allowed Values |  Short Description       
----------------  |  ----------------  |  ------------- |  ----------------
   gender_neutral | true               |   true, false  | Forces equal recombination fractions
   goal           |  maximize          |                |  
   lod_score_table|Lod_Score_Frame.txt | User defined output file name  |  Creates a lod score table output file
   output_unit    |                    |                |  
   parameters     |  1                 |                |  
   points         |   9                |                |  
   travel         |  grid              |                |  Mode of sampling parameter space

### Data Files
Two Point Linkage requires a [Control file](https://openmendel.github.io/MendelBase.jl/#control-file), and a [Pedigree file](https://openmendel.github.io/MendelBase.jl/#pedigree-file). Genotype data can be included in the Pedigree file, in which case a [Locus file](https://openmendel.github.io/MendelBase.jl/#locus-file) is required. Alternatively, genotype data can be provided in a [SNP data file](https://openmendel.github.io/MendelBase.jl/#snp-data-file), in which case a [SNP Definition File](https://openmendel.github.io/MendelBase.jl/#snp-definition-file) is required. OpenMendel will also accept [PLINK format](http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml) FAM and BIM files. Details on the format and contents of the Control and data files can be found on the [MendelBase](https://openmendel.github.io/MendelBase.jl) documentation page. There are example data files in the Two Point Linkage [docs](https://github.com/OpenMendel/MendelTwoPointLinkage.jl/tree/master/docs) folder.

### Running the Analysis

To run this analysis package, first launch Julia. Then load the package with the command:

     julia> using MendelTwoPointLinkage

Next, if necessary, change to the directory containing your files, for example,

     julia> cd("~/path/to/data/files/")

Finally, to run the analysis using the parameters in the control file Control_file.txt use the command:

     julia> TwoPointLinkage("Control_file.txt")

*Note: The package is called* MendelTwoPointLinkage *but the analysis function is called simply* TwoPointLinkage.

<!--- ### Interpreting the results
 ... --->

### Citation

If you use this analysis package in your research, please cite the following reference in the resulting publications:

*Lange K, Papp JC, Sinsheimer JS, Sripracha R, Zhou H, Sobel EM (2013) Mendel: The Swiss army knife of genetic analysis programs. Bioinformatics 29:1568-1570.*

<!--- ### Contributing
We welcome contributions to this Open Source project. To contribute, follow this procedure ... --->

### Acknowledgments

This project is supported by the National Institutes of Health under NIGMS awards R01GM053275 and R25GM103774 and NHGRI award R01HG006139.
