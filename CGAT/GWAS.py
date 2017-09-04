#########################################################################
#########################################################################
# Classes for handling genome-wide association input and output files, ##
# analysis and qc programs, and post-hoc analyses                      ##
#########################################################################
#########################################################################

import CGAT.Experiment as E
import CGAT.IOTools as IOTools
import numpy as np
import pandas as pd
import pandas.io.sql as pdsql
import re
import random
import os
import subprocess
import rpy2.robjects as ro
from rpy2.robjects import r as R
from rpy2.robjects import pandas2ri as py2ri
# set matplotlib non-interactive backend to Agg to
# allow running on cluster
import collections
import sqlite3 as sql
from math import *
import scipy.stats as stats


class FileGroup(object):
    '''
    An object for holding, formatting and processing files for genome-wide
    association analysis including compressed and binary files

    File types supported:
    * plink - .ped and .map files
    * plink binary - .bim, .fam. and .bed files
    * variant call format - .vcf and .bcf (including gzipped vcf)
    * Oxford format - .gen or .bgen with matched sample text file (must
                      be .sample)
    * GRM_binary - genetic relationship matrix calculated in an appropriate
      program in binary format.  File suffixes are *.grm.bin, *.grm.N.bin
      and *.grmid
    * GRM_gz - previously calcualted gzip compressed GRM, file suffixes
      are *.grm.gz and *.grm.id

    Phenotypes are assumed to be contained in the relevant files, if not
    then an additional phenotypes files can be included using the
    `phenotypes` argument.  Covariate files (if different from the phenotypes
    file) can also be included in the instantiation of a :FileGroup:
    object using the `covarite_files` argument.

    Only the `files` and `file_format` arguments are required.

    Genotype data are assumed to be raw genotype calls.  This can be modified
    using the `genotype_format` argument upon instantiation.  Values allowed
    are:
    * calls - standard bi-allelic genotype calls, i.e. AA, AB, BB
    * imputed_call - discrete genotype calls from imputed data,
                     essentially treated the same as ``calls``
    * genotype_prob - posterior probabilities for each genotype class,
                      i.e. 0.88 0.07 0.05 corresponding to homozygote
                      reference, heterozygote then homozygote rare allele.
    '''

    # Defaults for file formats
    ped_file = None
    map_file = None
    bim_file = None
    fam_file = None
    bed_file = None
    sample_file = None
    gen_file = None
    bgen_file = None
    vcf_file = None
    bcf_file = None

    def __init__(self, files, file_format, phenotypes=None,
                 genotype_format="calls", covariate_files=None):

        self.files = files
        self.file_format = file_format
        self.pheno_file = phenotypes
        self.genotype_format = genotype_format
        self.covariate_files = covariate_files
        self.set_file_prefix(files)

    def set_file_prefix(self, infiles):
        '''Get file prefixes from input files.  These are used across all
        file formats, e.g. myfile.bed, myfile.bim, myfile.fam name=myfile.
        Only use periods, '.' to denote file suffixes. use hyphens and
        underscores for separating file names.

        Set these to the appropriate attributes.
        '''

        file_prefixes = set()

        for f in infiles:
            # get all input file prefixes
            if len(f.split("/")) > 1:
                g = f.split("/")[-1]
                fdir = f.split("/")[:-1]
                fdir = "/".join(fdir)
                ffile = fdir + "/" + g.split(".")[0]
                file_prefixes.add(ffile)
            else:
                file_prefixes.add(f.split(".")[0])

        # if only prefix then use this for all data files
        if len(file_prefixes) == 1:
            self.name = [xf for xf in file_prefixes][0]
        else:
            # if there are multiple prefixes then use separate
            # flags for file inputs
            self.name = None

        # define file types by their suffix instead
        if self.file_format == "plink":
            self.ped_file = [pf for pf in infiles if re.search(".ped",
                                                               pf)][0]
            self.map_file = [mf for mf in infiles if re.search(".map",
                                                               mf)][0]

            # check files exist (i.e. are not the default None values)
            try:
                assert self.ped_file
            except AssertionError:
                raise ValueError(".ped file is missing, please "
                                 "specify")
            try:
                assert self.map_file
            except AssertionError:
                raise ValueError(".map file is missing, please "
                                 "specify")

        elif self.file_format == "plink_binary":
            self.fam_file = [ff for ff in infiles if re.search(".fam",
                                                               ff)][0]
            self.bim_file = [fb for fb in infiles if re.search(".bim",
                                                               fb)][0]
            self.bed_file = [bf for bf in infiles if re.search(".bed",
                                                               bf)][0]
            # check files exist (i.e. are not the default None values)
            try:
                assert self.fam_file
            except AssertionError:
                raise ValueError(".fam file is missing, please "
                                 "specify")
            try:
                assert self.bim_file
            except AssertionError:
                raise ValueError(".bim file is missing, please "
                                 "specify")
            try:
                assert self.bed_file
            except AssertionError:
                raise ValueError(".bed file is missing, please "
                                 "specify")

        elif self.file_format == "oxford":
            self.gen_file = [gf for gf in infiles if re.search(".gen",
                                                               gf)][0]
            self.sample_file = [sf for sf in infiles if re.search(".sample",
                                                                  sf)][0]
            # check files exist (i.e. are not the default None values)
            try:
                assert self.gen_file
            except AssertionError:
                raise ValueError(".gen file missing, please "
                                 "specify")
            try:
                assert self.sample_file
            except AssertionError:
                raise ValueError(".sample file missing, please "
                                 "specify")

        elif self.file_format == "oxford_binary":
            self.bgen_file = [bg for bg in infiles if re.search(".bgen",
                                                                bg)][0]
            self.sample_file = [sf for sf in infiles if re.search(".sample",
                                                                  sf)][0]
            # check files exist (i.e. are not the default None values)
            try:
                assert self.bgen_file
            except AssertionError:
                raise ValueError(".bgen file is missing, please "
                                 "specify")
            try:
                assert self.sample_file
            except AssertionError:
                raise ValueError(".sample file is missing, please "
                                 "specify")

        elif self.file_format == "vcf":
            self.vcf_file = [vf for vf in infiles if re.search(".vcf",
                                                               vf)][0]
            # check files exist (i.e. are not the default None values)
            try:
                assert self.vcf_file
            except AssertionError:
                raise ValueError(".vcf file is missing, please "
                                 "specify")

        elif self.file_format == "bcf":
            self.bcf_file = [bv for bv in infiles if re.search(".bcf",
                                                               bv)][0]
            # check files exist (i.e. are not the default None values)
            try:
                assert self.bcf_file
            except AssertionError:
                raise ValueError(".bcf file is missing, please "
                                 "specify")

        elif self.file_format == "GRM_binary":
            self.id_file = [ig for ig in infiles if re.search(".grm.id",
                                                              ig)][0]
            self.n_file = [gn for gn in infiles if re.search(".grm.N.bin",
                                                             gn)][0]
            self.bin_file = [gb for gb in infiles if re.search(".grm.bin",
                                                               gb)][0]
            # check files exits
            try:
                assert self.id_file
            except AssertionError:
                raise ValueError("GRM ids file is missing, please "
                                 "specify")
            try:
                assert self.n_file
            except AssertionError:
                raise ValueError("grm.N file is missing, please "
                                 "specify")
            try:
                assert self.bin_file
            except AssertionError:
                VaueError("GRM genotype is missing, please "
                          "specify")

        elif self.file_format == "GRM_plink":
            self.id_file = [ig for ig in infiles if re.search(".rel.id",
                                                              ig)][0]
            self.rel_file = [gn for gn in infiles if re.search(".rel.N.bin",
                                                               gn)][0]
            # check files exits
            try:
                assert self.id_file
            except AssertionError:
                raise ValueError("GRM ids file is missing, please "
                                 "specify")
            try:
                assert self.rel_file
            except AssertionError:
                raise ValueError("rel.N file is missing, please "
                                 "specify")

    def set_phenotype(self, pheno_file=None, pheno=1):
        '''
        Set the phenotype for a set of individuals
        using an external phenotypes file.

        Default is to use the (n+2)th column, designated
        as pheno 1.
        '''

        if type(pheno) == int:
            pheno = str(pheno)
        elif type(pheno) == str:
            pass
        else:
            raise AttributeError("Type of pheno unknown. "
                                 "Must be str or int.")
        self.pheno_file = pheno_file
        self.pheno = pheno


class GWASProgram(object):
    '''
    A base level object for programs designed to perform genome-wide
    association analysis and operate on genome-wide genotyping data.

    [INSERT PROPER DOCSTRING - see style guide]
    '''

    def __init__(self, executable=None, required_format=None):
        self.executable = executable
        self.require_format = required_format

    def program_call(self, infiles, outfile):
        '''build a statement to perform genome-wide
        analysis using infiles
        '''

        return ""

    def postprocess(self, infiles, outfile):
        '''collect and process output files from
        program - format for Result class objects'''

        return ""

    def build(self, infiles, outfile):
        '''run analysis program'''

        cmd_program = self.program_call(infile, outfile)
        cmd_postprocess = self.postprocess(infiles, outfile)

        if cmd_postprocess:
            cmd_postprocess = cmd_postprocess.strip().endswith(";")
            assert cmd_postprocess
        else:
            pass

        statement = " checkpoint; ".join((cmd_program,
                                          cmd_postprocess))

        return statement


class GCTA(GWASProgram):
    '''
    GCTA is designed for computing genetic relationship matrices, linear
    mixed model analyses and phenotype estimation/prediction.
    It can also perform SNP-wise GWAS.

    Files MUST be in Plink binary format
    '''

    def __init__(self, files, options=None, settings=None,
                 design=None):
        self.infiles = files
        self.options = options
        self.settings = settings
        self.design = design
        self.executable = "gcta64"
        self.statement = {}
        self.filters = []

    def program_call(self, infiles, outfile):
        '''build GCTA call statement on infiles'''

        statement = []
        statement.append(self.executable)

        if infiles.name:
            inputs = self._build_single_file_input(infiles,
                                                   infiles.file_format)
            statement.append(inputs)
        else:
            raise AttributeError("Files must be in binary plink format "
                                 "or as a GRM to use GCTA.  Please "
                                 "convert and try again.")

        if infiles.pheno_file:
            statement.append(" --pheno %s --mpheno %s " % (infiles.pheno_file,
                                                           infiles.pheno))
        else:
            pass

        self.statement["program"] = " ".join(statement)

    def _build_single_file_input(self, infiles, file_format):
        '''internal function only. Use it to construct the
        file input flags with --file, --bfile or --data
        '''

        statement = None

        if file_format == "plink":
            statement = " --file %s " % infiles.name
        elif file_format == "plink_binary":
            statement = " --bfile %s " % infiles.name
        elif file_format == "oxford" or file_format == "oxford_binary":
            statement = " --data %s" % infiles.name
        elif file_format == "GRM_binary" or file_format == "GRM_plink":
            statement = " --grm %s " % infiles.name
        else:
            raise AttributeError("file format is not defined or recognised."
                                 "Please define the input corectly when "
                                 "instantiating a FileGroup object")

        return statement

    def PCA(self, n_pcs="20"):
        '''
        Perform PCA analysis on previosly generated GRM, output the number n
        principal componets, default = 20
        '''

        self._run_tasks(pca=n_pcs)

    def apply_filters(self, filter_type, filter_value):
        '''
        * chromosome - exclude all variants not on the specified chromosome(s).
          [str/list]
        * autosome_number - for non-human species, the number of chromosomes to
          be considered autosomes
        * exclude_snps - text file list of variant IDs to exclude from analysis
          [file]
        * extract - text file list of variant IDs to include in analysis,
          ignores all others. [file]
        * min_allele_frequency - only include SNPs with cohort/case allele
          frequency above this threshold. [float]
        * max_allele_frequency - include all SNPs with a MAF equal to or below
          this value. [float]
        '''

        if filter_type == "chromosome":
            self._construct_filters(chromosome=filter_value)
        elif filter_type == "autosome_number":
            self._construct_filters(autosome_number=filter_value)
        elif filter_type == "exclude_snps":
            self._construct_filters(exclude_snps=filter_value)
        elif filter_type == "extract":
            self._construct_filters(extract=filter_value)
        elif filter_type == "min_allele_frequency":
            self._construct_filters(min_allele_frequency=filter_value)
        elif filter_type == "max_allele_frequency":
            self._construct_filters(max_allele_frequency=filter_value)
        elif filter_type == "keep":
            self._construct_filters(keep=filter_value)
        elif filter_type == "remove":
            self._construct_filters(remove=filter_value)

    def _construct_filters(self, **kwargs):
        '''
        Add filter to each GCTA run.

        The filters accepted are defined below.  These are input as keyword
        arguments supported by this function.

        * min_allele_frequency - only include SNPs with cohort/case allele
          frequency above this threshold. [float]
        * max_allele_frequency - include all SNPs with a MAF equal to or below
          this value. [float]
        * keep - keep individuals with matching individual and family IDs.
          [file]
        * remove - remove all individuals with matching individual and family
          IDs. [file]
        * extract - text file list of variant IDs to include in analysis,
          ignores all others. [file]
        * exclude - text file list of variant IDs to exclude from analysis.
          [file]
        * chromosome - exclude all variants not on the specified chromosome(s).
          [str/list]
        * autosome - exclude all non-place and non-autosomal variants.
          [boolean]
        * covariates_file - specify the covariates file with family and
          individual IDs in the first two columns.  Covariates are in the
          (n+2)th column. Only used in conjunction with `covariate_filter`.
          [file]
        * covariate_filter - covariate columns value to filter on.  Can be
          used with non-numeric values to filter out individuals with
          covariate =/= `covariate_filter` value. [str/int/float]
        * covariate_column - column number to apply filtering to if more
          than one covariate in the file. [int]
        * update_gender - provide gender information in a separate text
          file. [file]
        * grm_threshold - remove one of a pair of individuals with
          estimated relatedness greater than this value.
        * ld_significance - p-value threshold for regression test
          of LD significance
        * genotype_call - GenCall score cut-off for calling raw
          genotypes into Plink PED format
        * meta_pval - p-value threshold cut-off for conditional
          and joint genome-wide analysis
        * cojo_window - distance in kb beyond wich SNPs this
          distance apart are assumed to be in linkage equilibrium
        * cojo_collinear - multiple regression R^2 on selected SNPs
          value above which the testing SNP will not be selected.
        * cojo_inflation - adjust COJO analysis test statistics
          for genomic control. [boolean]
        * reml_iterations - maximum number of iterations to use
          during reml analysis.  Default is 100. [int]
        '''

        statement = []

        # map of keyword arguments recognised to Plink2 filtering flags
        filter_map = {"min_allele_frequency": " --maf %s ",
                      "max_allele_frequency": " --max-maf %s ",
                      "keep": " --keep %s ",
                      "remove": " --remove %s ",
                      "extract": " --extract %s ",
                      "exclude": " --exclude %s ",
                      "chromosome": " --chr %s ",
                      "autosome": " --autosome ",
                      "autosome_number": " --autosome-num %s ",
                      "grm_threshold": " --grm-cutoff %s ",
                      "ld_significance": " --ls-sig %s ",
                      "genotype_call": " --gencall %s ",
                      "meta_pval": " --cojo-p %s ",
                      "cojo_window": " --cojo-wind %s ",
                      "cojo_collinear": " --cojo-collinear %s ",
                      "cojo_inflation": " --cojo-gc ",
                      "reml_iterations": " --reml-maxit %s "}

        # compile all filters together, checking for dependencies.
        # use a mapping dictionary to extract the relevant flags and
        # combinations to use.
        filters = []
        filter_dict = {}
        for key, value in kwargs.items():
            filter_dict[key] = value

        for each in filter_dict.keys():
            try:
                assert filter_map[each]
                # check for data type <- behaviour is type dependent
                if type(filter_dict[each]) == 'bool':
                    filters.append(filter_map[each])
                else:
                    filter_val = filter_dict[each]
                    filters.append(filter_map[each] % filter_val)

            except KeyError:
                E.warn("%s filter not recognised, please see "
                       "documentation for allowed filters" % each)

        self.filters.append(" ".join(filters))
        self.statement["filters"] = " ".join(self.filters)

    def mixed_model(self, lmm_method, grm=None, qcovar=None,
                    dcovar=None):
        '''
        Run a linear mixed model with the GRM used to model
        random effects of an estimated genetic relationshi
        between individuals
        '''

        # add the mlm flag to the statement
        self._run_tasks(lmm=lmm_method)

        # construct the rest of mlm statement
        statement = []
        if qcovar:
            statement.append(" --qcovar %s " % qcovar)
        else:
            pass

        if dcovar:
            statement.append(" --covar %s " % dcovar)
        else:
            pass

        try:
            statement.append(" --grm %s " % grm)
        except ValueError:
            E.warn("No GRM has been provided, the GRM ")

        self.statement["mlm"] = " ".join(statement)

    def reml_analysis(self, method, parameters, prevalence=None,
                      qcovariates=None, discrete_covar=None):
        '''
        Use REML to estimate the proportion of phenotypic variance
        explained by the estimated genetic relationship between
        individuals.

        Arguments
        ---------
        method: string
          GCTA method to use for REML estimation of h2.  Includes:
          * snpBLUP - calculate the SNP BLUPs from the genotype
            data and the estimated total genetic value/ breeding value
          * fixed_cor -
          * priors - provide initial priors for the variance components
            estimation
          * unconstrained - allow variance estimates to fall outside
            of the normal parameter space, bounded [0, ).
          * GxE - estimate the contribution of GxE with covariates
            to the phenotype variance
          * BLUP_EBV - output individual total genetic effect/breeding
            values

        '''

        statement = []

        try:
            params = parameters.split(",")
            if len(params) == 1:
                params = params[0]
            else:
                pass
        except AttributeError:
            params = parameters

        self._run_tasks(parameter=params,
                        greml=method)

        if prevalence:
            statement.append(" --prevalence %0.3f " % prevalence)
        else:
            pass

        if qcovariates:
            statement.append(" --qcovar %s " % qcovariates)
        else:
            pass

        if discrete_covar:
            statement.append(" --covar %s " % discrete_covar)
        else:
            pass

        self.statement["reml"] = " ".join(statement)

    def _run_tasks(self, parameter=None, **kwargs):
        '''
        The principal functions of GCTA revolve around GRM estimation
        and variance components analysis, such as REML estimation of
        heritability and variance components, BLUP and phenotype prediciton.

        It can also be used to do PCA and conditional and joint GWAS.

        Tasks
        -----
        * pca - perform principal components analysis on a GRM
        * greml - perform restricted maximum likelihood analysis
          for estimation of variance components
        * estimate_ld - estimate the linkage disequilibrium structure
          over the genomic regions specified
        * simulate_gwas - simulate genome-wide association data based
          on observed genotype data
        * cojo - conditional and joint genome-wide association
          analysis across SNPs and covariates
        * bivariate_reml - perform GREML on two traits, either both
          binary, both quantitative or one of each
        * lmm - perform a linear mixed model based association analysis
        '''

        statement = []

        # set up a dictionary of recognised tasks with key word argument
        # values as further dictionaries. Use the parameter argument
        # to pass arguments by value to string formatting

        # put all of the other tasks as options in the calling function

        task_map = {"pca": " --pca %s ",
                    "greml": {"standard": " --reml ",
                              "priors": " --reml --reml-priors %s ",
                              "reml_algorithm": " --reml --reml-alg %s ",
                              "unconstrained": " --reml --reml-no-constrain ",
                              "GxE": " --reml --gxe %s ",
                              "LRT": " --reml --reml-lrt %s ",
                              "BLUP_EBV": " --reml --reml-pred-rand ",
                              "snpBLUP": " --blup-snp %s "},
                    "estimate_ld": " --ld %s ",
                    "simulate_gwas": {"quantitative": " --simu-qt ",
                                      "case_control": " --simu-cc %s %s "},
                    "cojo": {"stepwise": " --cojo-file %s --cojo-slct ",
                             "no_selection": " --cojo-file %s --cojo-joint ",
                             "snp_conditional": " --cojo-file %s --cojo-cond %s "},
                    "bivariate_reml": {"standard": " --reml-bivar %s ",
                                       "no_residual": " --reml-bivar %s --reml-bivar-nocove ",
                                       "fixed_cor": " --reml-bivar %s --reml-bivar-lrt-rg %s "},
                    "lmm": {"standard": " --mlma ",
                            "loco": " --mlma-loco ",
                            "no_covar": " --mlma-no-adj-covar "},
                    "remove_relations": {"cutoff": " --grm-cutoff %s "}}

        for task, value in kwargs.items():
            # check for PCA first as it is not nested in task_map
            if task == "pca":
                try:
                    state = task_map[task] % value
                    statement.append(state)
                except TypeError:
                    statement.append(task_map[task])
                statement.append
            # LD estimation is likewise not nested
            elif task == "estimate_ld":
                try:
                    state = task_map[task] % value
                    statement.append(state)
                except TypeError:
                    raise IOError("no SNP file list detected")
            elif task != "parameter":
                try:
                    # sub_task is a nested dictionary
                    sub_task = task_map[task]
                    try:
                        assert sub_task[value]
                        try:
                            # some tasks do not contain task values for the
                            # parameter argument - catch these with the TypeError
                            # exception
                            statement.append(sub_task[value] % parameter)

                            # the default for parameter is None, check this is appropriate
                            if not parameter:
                                E.warn("Parameter value is set to NoneType. "
                                       "Please check this is an appropriate value "
                                       "to pass for this task")
                            else:
                                pass
                        except TypeError:
                            statement.append(sub_task[value])
                    except KeyError:
                        raise KeyError("% Task not recognised, see docs for details of "
                                       "recognised tasks" % task)
                except KeyError:
                    raise KeyError("Task not recognised, see docs for details of "
                                   "recognised tasks")
            else:
                pass

        self.statement["tasks"] = " ".join(statement)

    def genetic_relationship_matrix(self, compression="binary", metric=None,
                                    shape="square", options=None):
        '''
        Calculate the estimated genetic relationship matrix from
        genotyping data

        * estimate_grm - estimate the realized genetic relationship
          matrix between individuals from genotyping data
        '''

        mapf = {"binary": " --make-grm-bin ",
                "gzip": " --make-grm-gz ",
                "no_compress": " --make-grm ",
                "X_chr": " --make-grm-chr ",
                "X_chr_gz": " --make-grm-gz ",
                "inbreeding": " --ibc "}

        if options == "X_chr":
            if compression == "gz":
                state = mapf["X_chr_gz"]
            else:
                state = mapf["X_chr"]
        elif options == "inbreding":
            state = mapf["inbreeding"]
        else:
            pass

        # check compression is compatible
        if compression == "gz":
            state = mapf["gzip"]
        elif compression == "bin":
            state = mapf["binary"]
        elif compression is None and not options:
            state = mapf["no_compress"]

        self.statement["matrix"] = state

    def build_statement(self, infiles, outfile, threads=None,
                        memory=None, parallel=None):
        '''
        Build statement and execute from components
        '''

        statement = []
        exec_state = self.executable
        # calls to function add to the self.statement dictionary
        try:
            statement.append(self.statement["program"])
        except KeyError:
            raise AttributeError("Input files and format not detected")

        try:
            statement.append(self.statement["filters"])
        except KeyError:
            pass

        try:
            statement.append(self.statement["tasks"])
        except KeyError:
            pass

        try:
            statement.append(self.statement["matrix"])
        except KeyError:
            pass

        try:
            statement.append(self.statement["mlm"])
        except KeyError:
            pass

        try:
            statement.append(self.statement["reml"])
        except KeyError:
            pass

        if threads:
            statement.append(" --thread-num %i " % threads)
        else:
            pass

        # add output flag
        statement.append(" --out %s " % outfile)

        os.system(" ".join(statement))


class Plink2(GWASProgram):
    '''
    Run various Plink functions and analysis, including file processing, GRM
    calculation, PCA and other GWA tasks

    Require Plink v1.9 to be in the users PATH variable as ``plink2`` to
    distinguish it from Plink v1.07.
    '''

    def __init__(self, files, options=None,
                 settings=None, design=None):
        self.infiles = files
        self.options = options
        self.settings = settings
        self.design = design
        self.executable = "plink2"
        self.statement = {}
        self.filters = []

    def program_call(self, infiles, outfile):
        ''' build Plink call statement on infiles'''

        statement = []
        statement.append(self.executable)

        if infiles.name:
            inputs = self. _build_single_file_input(infiles,
                                                    infiles.file_format)
            statement.append(inputs)

        else:
            inputs = self._build_multiple_file_input(infiles,
                                                     infiles.file_format)
            statement.append(inputs)

        # check for the presence of an additional phenotypes file
        try:
            if infiles.pheno_file:
                statement.append(" --pheno %s --mpheno %s " % (infiles.pheno_file,
                                                               infiles.pheno))
            else:
                pass
        except AttributeError:
            pass

        self.statement["program"] = " ".join(statement)

    def hamming_matrix(self, shape, compression, options):
        '''
        Calculate genomic pair-wise distance matrix between
        individuals using Hamming distance across all variants
        '''

        # check shape is compatible
        if not shape:
            shape = "triangle"
        elif shape in ["square", "square0", "triangle"]:
            pass
        else:
            raise ValueError("matrix shape %s not recognised."
                             "Valid options are square, square0, "
                             "and triangle." % shape)

        # check compression is compatible
        if compression in ["gz", "bin", "bin4"]:
            pass
        else:
            raise ValueError("compression %s not recognised. Accepted "
                             "formats are gz, bin and bin4." % compression)

        if options:
            state = self._matrices(matrix_type="hamming", shape=shape,
                                   compression=compression, options=options)
        else:
            state = self._matrices(matrix_type="hamming", shape=shape,
                                   compression=compression)

        self.statement["matrix"] = state

    def ibs_matrix(self, shape, compression, options):
        '''
        Calculate genomic pair-wise similarity matrix between
        individuals using proportion of IBS alleles
       '''

        # check shape is compatible
        if shape in ["square", "square0", "triangle"]:
            pass
        else:
            raise ValueError("matrix shape %s not recognised."
                             "Valid options are square, square0, "
                             "and triangle." % shape)

        # check compression is compatible
        if compression in ["gz", "bin", "bin4"]:
            pass
        else:
            raise ValueError("compression %s not recognised. Accepted "
                             "formats are gz, bin and bin4." % compression)

        if options:
            state = self._matrices(matrix_type="ibs", shape=shape,
                                   compression=compression, options=options)
        else:
            state = self._matrices(matrix_type="ibs", shape=shape,
                                   compression=compression)

        self.statement["matrix"] = state

    def genome_matrix(self, shape, compression, options):
        '''
        Calculate genomic pair-wise distance matrix between
        individuals using 1 - proportion of IBS alleles
       '''

        # check shape is compatible
        if shape in ["square", "square0", "triangle"]:
            pass
        else:
            raise ValueError("matrix shape %s not recognised."
                             "Valid options are square, square0, "
                             "and triangle." % shape)

        # check compression is compatible
        if compression in ["gz", "bin", "bin4"]:
            pass
        else:
            raise ValueError("compression %s not recognised. Accepted "
                             "formats are gz, bin and bin4." % compression)

        if options:
            state = self._matrices(matrix_type="genomic", shape=shape,
                                   compression=compression, options=options)
        else:
            state = self._matrices(matrix_type="genomic", shape=shape,
                                   compression=compression)

        self.statement["matrix"] = state

    def genetic_relationship_matrix(self, shape, compression, metric,
                                    options=None):
        '''
        Calculate genomic pair-wise distance matrix between
        individuals using proportion of IBS alleles
        Requires the use of the Plink2 parallelisation to run with large
        cohorts of patients
        '''

        # check shape is compatible
        if shape in ["square", "square0", "triangle"]:
            pass
        else:
            raise ValueError("matrix shape %s not recognised."
                             "Valid options are square, square0, "
                             "and triangle." % shape)

        # check compression is compatible
        if compression in ["gz", "bin", "bin4"]:
            pass
        else:
            raise ValueError("compression %s not recognised. Accepted "
                             "formats are gz, bin and bin4." % compression)

        if metric in ["cov", "ibc2", "ibc3"]:
            state = self._matrices(matrix_type="grm", shape=shape,
                                   compression=compression, options=metric)
        else:
            E.info("%s metric not recognised.  Running with default Fhat1" % metric)
            state = self._matrices(matrix_type="grm", shape=shape,
                                   compression=compression)

        self.statement["matrix"] = state

    def apply_filters(self, filter_type, filter_value):
        '''
        arguments supported by this function.

        * genotype_rate - exclude SNPs with a genotyping rate below this
          value. [float]
        * min_allele_frequency - only include SNPs with cohort/case allele
          frequency above this threshold. [float]
        * max_allele_frequency - include all SNPs with a MAF equal to or below
          this value. [float]
        * exclude_snp - exclude this single variant
        * exclude_snps - text file list of variant IDs to exclude from analysis.
          [file]
        * chromosome - exclude all variants not on the specified chromosome(s).
          [str/list]
        * exclude_chromosome - exclude all variants on the specified
          chromosome(s). [str/list]
        * autosome - exclude all non-place and non-autosomal variants.
          [boolean]
        * pseudo_autosome - include the pseudo-autosomal region of chromosome
          X. [boolean]
        * ignore_indels - remove all indels/multi-character allele coding
          variants. [boolean]
        * snp_bp_range - (from, to) range in bp of variants to include in
          analysis. [tuple]

        '''

        if filter_type == "genotype_rate":
            self._construct_filters(genotype_rate=filter_value)
        elif filter_type == "hwe":
            self._construct_filters(hwe=filter_value)
        elif filter_type == "missingness":
            self._construct_filters(missingness=filter_value)
        elif filter_type == "min_allele_frequency":
            self._construct_filters(min_allele_frequency=filter_value)
        elif filter_type == "max_allele_frequency":
            self._construct_filters(max_allele_frequency=filter_value)
        elif filter_type == "exclude_snp":
            self._construct_filters(exclude_snp=filter_value)
        elif filter_type == "exclude":
            self._construct_filters(exclude=filter_value)
        elif filter_type == "extract":
            self._construct_filters(extract=filter_value)
        elif filter_type == "chromosome":
            self._construct_filters(chromosome=filter_value)
        elif filter_type == "exclude_chromosome":
            self._constuct_filters(exclude_chromosome=filter_value)
        elif filter_type == "autosome":
            self._construct_filters(autosome=filter_value)
        elif filter_type == "pseudo_autosome":
            self._construct_filters(pseudo_autosome=filter_value)
        elif filter_type == "ignore_indels":
            self._construct_filters(ignore_indels=filter_value)
        elif filter_type == "snp_bp_range":
            self._construct_filters(snp_bp_range=filter_value)
        elif filter_type == "conditional_snp":
            self._construct_filters(conditional_snp=filter_value)
        elif filter_type == "keep":
            self._construct_filters(keep=filter_value)
        elif filter_type == "remove":
            self._construct_filters(remove=filter_value)
        elif filter_type == "ignore_indels":
            self._construct_filters(ignore_indels=filter_value)

    def _build_multiple_file_input(self, infiles, file_format):
        '''
        internal function only.  Use it to construct
        the appropriate file input flags
        '''

        statement = None

        if file_format == "oxford":
            statement = " --gen %s --sample %s " % (infiles.gen_file,
                                                    infiles.sample_file)
        elif file_format == "oxford_binary":
            statement = " --bgen %s --sample %s " % (infiles.bgen_file,
                                                     infiles.sample_file)
        elif file_format == "plink":
            statement = " --ped %s --map %s " % (infiles.ped_file,
                                                 infiles.sample_file)
        elif file_format == "plink_binary":
            statement = " --bed %s --bim %s --fam %s " % (infiles.bed_file,
                                                          infiles.bim_file,
                                                          infiles.fam_file)
        elif file_format == "vcf":
            statement = " --vcf %s.vcf.gz " % infiles.vcf_file
        elif file_format == "bcf":
            statement = " --bcf %s " % infiles.vcf_file
        elif file_format == "GRM_binary":
            statement = " --grm-bin %s " % infiles.name
        else:
            raise AttributeError("file format is not defined.  Please "
                                 "define the input file formats when "
                                 "instantiating a FileGroup object")

        return statement

    def _build_single_file_input(self, infiles, file_format):
        '''internal function only. Use it to construct the
        file input flags with --file, --bfile or --data
        '''

        statement = None

        if file_format == "plink":
            statement = " --file %s " % infiles.name
        elif file_format == "plink_binary":
            statement = " --bfile %s " % infiles.name
        elif file_format == "oxford" or file_format == "oxford_binary":
            statement = " --data %s" % infiles.name
        elif file_format == "GRM_plink":
            statement = " --grm.bin  %s " % infiles.name
        elif file_format == "GRM_binary":
            statement = " --grm-bin %s " % infiles.name
        elif file_format == "vcf":
            statement = " --vcf %s.vcf.gz " % infiles.name
        else:
            raise AttributeError("file format is not defined or recognised."
                                 "Please define the input corectly when "
                                 "instantiating a FileGroup object")

        return statement

    def _construct_filters(self, **kwargs):
        '''
        Add filter to each plink run. [data type]

        The filters accepted are defined below.  These are input as keyword
        arguments supported by this function.

        * genotype_rate - exclude SNPs with a genotyping rate below this
          value. [float]
        * missingness - exclude individuals with total genotype missingness
          above this value. [float]
        * hwe - p-value threshold for excluding SNPs deviating from
          Hardy-Weinberg expectations. [float]
        * min_allele_frequency - only include SNPs with cohort/case allele
          frequency above this threshold. [float]
        * max_allele_frequency - include all SNPs with a MAF equal to or below
          this value. [float]
        * mendelian_error - filter out samples/trios exceeding the error
          threshold. [float]
        * keep - keep individuals with matching individual and family IDs.
          [file]
        * remove - remove all individuals with matching individual and family
          IDs. [file]
        * quality_score_file - vcf file with variants and quality scores.  Use
          `qual_score_column` and `var_id_col` to specify which columns
          correspond to the quality score and variant ID columns.
          [file] <int> <int>
        * min_qual_score - alters the lower bound of the quality score
          threshold; default is 0.[int]
        * max_qual_score - sets an upper limit on the quality scores;
          default is Inf. [int]
        * allow_no_sex - prevents phenotypes set to missing if there is no
          gender information. [boolean]
        * enforce_sex - force phenotype missing when using --make-bed, --recode
          or --write-covar. [boolean]
        * subset_filter - filter on a particular subset.  Choices are: cases,
          controls, males, females, founders, nonfounders. [str]
        * extract - text file list of variant IDs to include in analysis,
          ignores all others. [file]
        * exclude - text file list of variant IDs to exclude from analysis.
          [file]
        * chromosome - exclude all variants not on the specified chromosome(s).
          [str/list]
        * exclude_chromosome - exclude all variants on the specified
          chromosome(s). [str/list]
        * autosome - exclude all non-place and non-autosomal variants.
          [boolean]
        * pseudo_autosome - include the pseudo-autosomal region of chromosome
          X. [boolean]
        * ignore_indels - remove all indels/multi-character allele coding
          variants. [boolean]
        * snp_bp_range - (from, to) range in bp of variants to include in
          analysis. [tuple]
        * specific_snp - only load the variant specified. [str]
        * exclude_snp - exclude this single variant
        * window_size - alters behaviour of `specific_snp` and `exclude_snp`
          to include/exclude SNPs within +/- half of this distance (kb) are
          also included. [float]
        * range_resolution - sets the resolution of the (from, to) range.
          Either bp, kb or mb. If set it will take the values from
          `snp_bp_range`. [str/int/float]
        * covariates_file - specify the covariates file with family and
          individual IDs in the first two columns.  Covariates are in the
          (n+2)th column. Only used in conjunction with `covariate_filter`.
          [file]
        * covariate_filter - covariate columns value to filter on.  Can be
          used with non-numeric values to filter out individuals with
          covariate =/= `covariate_filter` value. [str/int/float]
        * covariate_column - column number to apply filtering to if more
          than one covariate in the file. [int]
        '''

        statement = []

        # map of keyword arguments recognised to Plink2 filtering flags
        filter_map = {"genotype_rate": " --geno %s ",
                      "missingness": "--mind %s ",
                      "hwe": " --hwe %s ",
                      "min_allele_frequency": " --maf %s ",
                      "max_allele_frequency": " --max-maf %s ",
                      "mendelian_error": " --me %s ",
                      "keep": " --keep %s ",
                      "remove": " --remove %s ",
                      "quality_score_file": " --qual-scores %s ",
                      "qual_score_column": " %s ",
                      "var_id_col": " %s ",
                      "min_qual_score": " --qual-threshold %s ",
                      "max_qual_score": " --qual-max-threshold %s ",
                      "allow_no_sex": " --allow-no-sex ",
                      "enforce_sex": " --must-have-sex ",
                      "subset_filter": " --filter-%s ",
                      "extract": " --extract %s ",
                      "exclude": " --exclude %s ",
                      "chromosome": " --chr %s ",
                      "exclude_chromosome": " --not-chr %s ",
                      "autosome": " --autosome ",
                      "pseudo_autosome": " --autosome-xy ",
                      "ignore_indels": " --snps-only no-DI ",
                      "snp_id_range": " --from %s --to %s ",
                      "specific_snp": " --snp %s ",
                      "window_size": " --window %s ",
                      "exclude_snp": " --exclude-snp %s ",
                      "snp_bp_range": "--from-bp %s --to-bp %s ",
                      "covariates_file": " --filter %s ",
                      "covariate_filter": " %s ",
                      "covariate_column": " --mfilter %s ",
                      "missing_phenotype": " --prune ",
                      "conditional_snp": " --condition %s ",
                      "haplotype_size": " --blocks-max-kb %s ",
                      "haplotype_frequency": " --blocks-min-maf %s "}

        # compile all filters together, checking for dependencies.
        # use a mapping dictionary to extract the relevant flags and
        # combinations to use.
        filters = []
        filter_dict = {}
        for key, value in kwargs.items():
            filter_dict[key] = value

        # need to check for covariates and qual scores - these
        # are more complex.  Deal with these first and remove
        # from dictionary once complete.
        try:
            assert filter_dict["quality_score_file"]
            assert filter_dict["qual_score_column"]
            assert filter_dict["var_id_col"]

            quals = []
            qual_file = filter_dict["quality_score_file"]
            score_col = filter_dict["qual_score_column"]
            id_col = filter_dict["var_id_col"]

            quals.append(filter_map["quality_score_file"] % qual_file)
            quals.append(filter_map["qual_score_column"] % score_col)
            quals.append(filter_map["var_id_col"] % id_col)

            # remove from dictionary
            filter_dict.pop("qual_score_column", None)
            filter_dict.pop("var_id_col", None)

            filters.append(" ".join(quals))

        except KeyError:
            pass

        try:
            assert filter_dict["covariates_file"]
            assert filter_dict["covariate_filter"]

            covars = []
            covar_file = filter_dict["covariates_file"]
            covar_val = filter_dict["covariate_filter"]
            covars.append(filter_map["covariates_file"] % covar_file)
            covars.append(filter_map["covariate_filter"] % covar_val)

            # check to filter on specific column numnber, default is 3rd file
            # column, i.e. (n+2)th column
            try:
                assert filter_dict["covariate_column"]
                covar_col = filter_dict["covariate_column"]
                covars.append(filter_map["covariate_column"] % covar_col)
                filter_dict.pop("covariate_column", None)
            except KeyError:
                pass

            # remove from dictionary
            filter_dict.pop("covariates_file", None)
            filter_dict.pop("covariate_filter", None)

            filters.append(" ".join(covars))

        except KeyError:
            pass

        # range_resolution and snp_bp_range are used together
        try:
            assert filter_dict["snp_bp_range"]
            flags = filter_map["snp_bp_range"]
            from_pos = filter_dict["snp_bp_range"].split(",")[0]
            to_pos = filter_dict["snp_bp_range"].split(",")[1]
            filters.append(flags % (from_pos, to_pos))

            # remove so they are not duplicated - source of bugs
            filter_dict.pop("snp_bp_range", None)

        except KeyError:
            pass

        for each in filter_dict.keys():
            try:
                assert filter_map[each]
                # check for data type <- behaviour is type dependent
                if type(filter_dict[each]) == bool:
                    filters.append(filter_map[each])
                # handle multiple arguments in string format
                elif len(filter_dict[each].split(",")) > 1:
                    vals = tuple(filter_dict[each].split(","))
                    filters.append(filter_map[each] % vals)
                else:
                    filter_val = filter_dict[each]
                    filters.append(filter_map[each] % filter_val)

            except KeyError:
                E.warn("%s filter not recognised, please see "
                       "documentation for allowed filters" % each)

        self.filters.append(" ".join(filters))
        self.statement["filters"] = " ".join(self.filters)

    def calc_ld(self, ld_statistic, ld_threshold,
                ld_shape="table"):
        '''
        Calculate linkage disequilibrium between all SNP
        pairs.

        Arguments
        ---------
        ld_statistic: string
          The LD statistic to report, either correlation or squared correlation
          of inter-variant allele counts

        ld_threshold: float
          minimum value to report for pair-wise LD

        ld_window: int
          max distance (in Kb) between SNPs for calculating LD

        ld_shape: string
          shape to use for reporting LD, either a table or a matrix.  If a
          matrix then either square, square with diagnonal (square0) or
          triangular.  Square matrices are symmetric.
        '''

        statement = []
        ld_map = {"r": " --r %s dprime ",
                  "r2": "--r2 %s dprime "}

        shape_map = {"table":  "inter-chr gz",
                     "square": "square gz",
                     "square0": "square0 gz",
                     "triangle": "triangle gz"}

        try:
            statement.append(ld_map[ld_statistic] % shape_map[ld_shape])
        except KeyError:
            raise ValueError("%s LD statistic not recognised. Please "
                             "use eithr 'r' or 'r2'" % ld_statistic)

        if type(ld_threshold) == float:
            statement.append(" --ld-window-r2 %0.3f " % ld_threshold)
        else:
            E.warn("threshold type not recognised, setting to default "
                   "value of 0.2")

        self.statement["tasks"] = " ".join(statement)

    def _run_tasks(self, parameter=None, **kwargs):
        '''
        Plink2 is capable of much more than just running basic association
        analyses.

        These include file processing, reformating, filtering, data summaries,
        PCA, clustering, GRM calculation (slow and memory intense), etc.

        multiple tasks can be added by separate calls to this function.
        For instance, adding phenotype and gender information using the
        update_samples task whilst change the file format.


        Tasks
        -----

        * change_format - convert from input format to an alternative format
          after applying filters.
        * change_missing_values - alters the genotype or phenotype missing
          value into the value supplied.
        * update_variants - use this to fill in missing variant IDs, useful
          for data from exome or whole-genome sequencing that have
          non-standard IDs.
        * update_samples - update phenotype and sample information
        * flip_strands - flip the strand for alleles, swaps A for T and
          C for G.
        * flip_scan - use the LD-based scan to check SNPs have not had
          incorrect strand assignment. Particularly useful if cases and
          controls were genotyped separately, or the cohort was genotyped
          in different batches.
        * sort - sort files by individual and/or family IDs
        * merge - merge new filesets with reference fileset.
        * merge_mode - handling of missing values and overwriting values
        * find_duplicates - find and output duplicate variants based on bp position,
          or variant ID.  Useful to output for the --exclude filtering flag.
        * remove_relations - remove one of a pair of individuals with IBS >=
          a threshold.  Recommended minimum is 3rd cousins (IBS >= 0.03125).
        * check_gender - check imputed gender from non-pseudoautosomal X
          chromsome genotypes against self-reported gender
        * estimate_haplotypes - assign SNPs to haplotype blocks and get
          positional information
        '''

        statement = []

        # set up a dictionary of recognised tasks with key word argument
        # values as further dictionaries. Use the parameter argument
        # to pass arguments by value to string formatting

        task_map = {'change_format': {"plink_binary": " --make-bed ",
                                      "plink": " --recode ",
                                      "oxford": " --recode oxford ",
                                      "oxford_binary": " --recode oxford gen-gz ",
                                      "raw": " --recode A tabx "},
                    "change_missing_values": {"genotype": " --missing-genotype %s ",
                                              "phenotype": " --missing-phenotype %s "},
                    "update_variants": {"variant_ids": " --set-missing-var-ids %s ",
                                        "missing_id": " --mising-var-code %s ",
                                        "chromosome": " --update-chr %s ",
                                        "centimorgan": " --update-cm %s ",
                                        "name": " --update-name %s ",
                                        "alleles": " --update-alleles  %s ",
                                        "map": " --update-map %s "},
                    "update_samples": {"sample_ids": " --update-ids %s ",
                                       "parents": " --update-parents %s ",
                                       "gender": " --update-sex %s %s "},
                    "flip_strands": {"all_samples": " --flip %s ",
                                     "subset": " --flip-subset %s "},
                    "flip_scan": {"default": " --flip-scan verbose ",
                                  "window": "--flip-scan --flip-scan-window %s ",
                                  "kb": " --flip-scan  --flip-scan-window-kb %s ",
                                  "threshold": " --flip-scan  --flip-scan-threshold %s "},
                    "sort": {"none": " --indiv-sort %s ",
                             "natural": " --indiv-sort %s ",
                             "ascii": " --indiv-sort %s ",
                             "file": " --indiv-sort %s "},
                    "merge": {"plink": " --merge %s ",
                              "binary_plink": " --bmerge %s "},
                    "merge_mode": {"default": " --merge-mode 1 ",
                                   "orginal_missing": " --merge-mode 2 ",
                                   "new_nonmissing": " --merge-mode 3 ",
                                   "no_overwrite": " --merge-mode 4 ",
                                   "force": " --merge-mode 5 ",
                                   "report_all": " --merge-mode 6 ",
                                   "report_nonmissing": " --merge-mode 7"},
                    "find_duplicates": {"same_ref": " --list-duplicate-vars require-same-ref ",
                                        "id_match": " --list-duplicate-vars ids-only ",
                                        "suppress_first": " --list-duplicate-vars suppress-first"},
                    "remove_relations": {"cutoff": " --rel-cutoff %s "},
                    "check_gender": " --check-sex ",
                    "pca": " --pca %s ",
                    "estimate_haplotypes": " --blocks "}

        for task, value in kwargs.items():
            # check for PCA first as it is not nested in task_map
            if task == "pca":
                try:
                    state = task_map[task] % value
                    statement.append(state)
                except TypeError:
                    statement.append(task_map[task])
                statement.append
            elif task == "check_gender":
                statement.append(task_map[task])
            elif task == "estimate_haplotypes":
                statement.append(task_map[task])
            elif task != "parameter":
                try:
                    # sub_task is a nested dictionary
                    sub_task = task_map[task]
                    try:
                        assert sub_task[value]
                        try:
                            # gender has two string formats
                            if value == "gender":
                                gcol = 1
                                statement.append(sub_task[value] % (parameter,
                                                                    gcol))
                            else:
                                # some tasks do not contain task values for the
                                # parameter argument - catch these with the TypeError
                                # exception
                                statement.append(sub_task[value] % parameter)
                            # the default for parameter is None, check this is appropriate
                            if not parameter:
                                E.warn("Parameter value is set to NoneType. "
                                       "Please check this is an appropriate value "
                                       "to pass for this task")
                            else:
                                pass
                        except TypeError:
                            statement.append(sub_task[value])
                    except KeyError:
                        raise KeyError("No sub task found, see docs for details of "
                                       "recognised tasks")
                except KeyError:
                    raise KeyError("Task not recognised, see docs for details of "
                                   "recognised tasks")
            else:
                pass
        # handle multiple tasks for a single run
        try:
            curr_tasks = self.statement["tasks"]
            new_tasks = " ".join(statement)
            self.statement["tasks"] = " ".join([curr_tasks, new_tasks])
        except KeyError:
            self.statement["tasks"] = " ".join(statement)

    def _output_statistics(self, **kwargs):
        '''
        Summary statistics are written to specific files dictated by the
        type of statistic

        Statistics
        ----------
        * allele_frequency - writes out MAF to `plink`.frq, this can be
          modified with specific keywords.
        * missing_data - generates a report of data missingness, can be subset
          into within family and/or cluster reports
        * hardy_weinberg - calculates all HWE p-values using exact test
          statistics. For case/control studies reports are written for case,
          controls and combined.
        * mendel_errors - generates a Mendelian error report across all trios.
          There are 10 different codes responding to different Mendelian error
          scenarios.
        * inbreeding - calculate observed and expected homozygosity across
          individuals and F statistics.  If the sample size is small then a
          file of MAFs is required. Inbreeding coefficients can also be
          reported on request using inbreeding_coef.
        * gender_checker - checks gender assignment against X chromosome
          genotypes. Gender values can also be imputed based on genotype
          information using gender_impute.
        * wrights_fst - calculate Wright's Fst statistic given a set of
          subpopulations for each autosomal diploid variant.  Used in
          conjunction with the --within flag.
        '''

        stats_map = {"allele_frequency": " --freq %s ",
                     "missing_data": " --missing %s ",
                     "hardy_weinberg": " --hardy midp ",
                     "mendel_errors": " --mendel %s ",
                     "inbreeding": " --het %s ",
                     "inbreeding_coef": " --ibc ",
                     "gender_checker": " --check-sex ",
                     "gender_impute": " --impute-sex ",
                     "wrights_fst": " --fst --within %s ",
                     "case_control_fst": "--fst %s "}

        statement = []
        for key, value in kwargs.tems():
            if value:
                try:
                    assert stats_map[key]
                    statement.append(stats_map[key] % value)
                except KeyError:
                    raise KeyError("statistic not recognised.  Please "
                                   "consult the documentation for allowed "
                                   "options.")
            else:
                try:
                    assert stats_map[key]
                    flag = stats_map[key].rstrip("%s ")
                    statement.append(flag)
                except KeyError:
                    raise KeyError("statistic not recognised.  Please "
                                   "consult the documentation for allowed "
                                   "options.")
        self.statement["stats"] = " ".join(statement)

    def run_association(self, association=None, model=None,
                        run_options=None,
                        permutation=False, n_perms=None,
                        random_seed=None, permutation_options=None,
                        covariates_file=None, covariates=None):
        '''
        Construct a statement for a plink2 association analysis.

        QC filters are constructed from input during instantiation.

        run options include redirecting logging output, using parallelisation,
        defining number of threads to use, etc

        The default association uses the --assoc flag.  Plink will check
        phenotype coding, if it is not case/control it assumes
        it is a continuous trait and uses linear regression.

        Alternative regression models that include covariates can be used,
        i.e. logistic and linear regression.

        key
        ***
        {CC} - applies to case/control analysis only
        {quant} - applies to quantitative trait only
        {CC/quant} - applies to both

        run_options
        -----------
        ``--assoc``:
            * `fisher | fisher-midp` - uses Fisher's exact test to calculate
            association p-values or applies Lancaster's mid-p adjustment. {CC}
            * `counts` - causes --assoc to report allele counts instead of
            frequencies. {CC}
            * `set-test` - implements and tests the significance of variant
            sets.  See documentation below. {CC/quant}
            * `qt-means` - generates a .qassoc.means file reporting trait means
            and standard deviations by genotype. {quant}
            * `lin` - reports the Lin et al (2006) statistic to be reported. If
            multiple testing adjustments and/or permutation is also used, they
            will be based on this statistic. {quant}

        ``--model``:
            * `fisher | fisher-midp | trend-only` - uses Fisher's exact test
            to calculate association p-values or applies Lancaster's mid-p
            adjustment. trend-only forces only a trend test to be performed.
            {CC}
            * `dom | rec | gen | trend` - use the specified test as the basis
            for the model permutation.  If none are defined the result with the
            smallest p-value is reported. {CC}
            * --cell - sets the minimum number of observations per cell in the
            2x3 contingency table.  The default is 0 with the Fisher and
            Fiser-midp test, otherwise 5. {CC}

        ``--linear/logistic``:
            * `set-test` - implements and tests the significance of variant
            sets.  See documentation below. {CC/quant}
            * `hide-covar` - removes the covariate specific sections from the
            results output. {CC/quant
            * `sex | no-x-sex` - `sex` adds sex as covariate to all models,
            whislt `no-x-sex` does not include gender into X-chromosome SNP
            models. {CC/quant}
            * `interaction` - adds in genotype X covariate interaction terms
            into the model. Can only be used with permutation is ``--tests``
            is also specified. {CC/quant}
            * `beta` - reports the beta coefficients instead of the OR in a
            logistic model. {CC}
            * `standard-beta` - standardizes the phenotype and all predictor
            variables to zero mean and unit variance prior to regression
            (separate for each variant analysed). {quant}
            * `intercept` - includes the intercept in the output results.
            {quant}

        model
        -----
        * `recessive` - `recessive` specifies the model assuming the A1 allele
          as recessive. {CC/quant}
        * `dominant` - `dominant` specifies the model assuming the A1 allele is
          dominant. {CC/quant}
        * `genotype` - `genotype` adds an additive effect/dominance deviation
          2df joint test with two genotype variables in the test (coded 0/1/2
          and 0/1/0). {CC/quant}
        * `trend` - forces a trend test to be performed. {CC/quant}
        * `hethom` - `hethom` uses 0/0/1 and 0/1/0 instead of the genotype
          coding. With permutation it will be based on the joint test instead
          of just the additive effects.  This can be overriden using the
          `--tests` flag. {CC/quant}
        * `no-snp` - `no-snp` defines a regression of phenotype on covariates
          without reference to genotype data, except where `--conditon{-list}`
          is specified. If used with permuation, test results will be reported
          for every covariate. {CC/quant}

        permutation
        -----------
        If permutation is True, run an adaptive Monte Carlo permutation test.
        If n_perms is set, this will run a max(T) permutation test with the n
        replications. A random seed will need to be provided.

        * `perm-count` - this alters the permutation output report to include
          counts instead of frequencies

        covariates
        ----------
        These should be provided in a separate file.  Specifying which
        covariates to include can be done as either a comma-separated list
        of covariate names or numbers. These numbers will correspond to the
        (n+2)th covariate file column as per the plink documentation.
        '''

        # model map maps common option effects onto specific syntax
        model_map = {"--logistic": {"recessive": "recssive",
                                    "dominant": "dominant",
                                    "genotype": "genotypic"},
                     "--linear": {"recessive": "recssive",
                                  "dominant": "dominant",
                                  "genotype": "genotypic"},
                     "--model": {"recessive": "rec",
                                 "dominant": "dom",
                                 "genotype": "gen"}}

        statement = []
        # construct analysis flags
        # add model, i.e. additive, recessive, dominant, etc.
        # see docstring for details.  Make sure correct modifier is used
        # with a mapping dictionary

        if association == "logistic":
            statement.append(" --logistic ")
            m_map = model_map["--logistic"]
            if model:
                statement.append(m_map[model])
            else:
                pass

        elif association == "linear":
            statement.append(" --linear ")
            m_map = model_map["--linear"]
            if model:
                statement.append(m_map[model])
            else:
                pass

        elif association == "model":
            statement.append(" --model ")
            m_map = model_map["--model"]
            statement.append(m_map[model])
        else:
            statement.append(" --assoc ")

        # add in run options.  These need to be in their correct
        # format already
        if run_options:
            modifiers = " ".join(run_options)
            statement.append(modifiers)
        else:
            pass

        # permutation should have a random seed set by the user.  Allow
        # this to set it's own seed if one not provided, but report it in
        # the log file
        if permutation:
            try:
                assert random_seed
            except AssertionError:
                rand_seed = random.randint(0, 100000000)
                E.warn("No seed is provided for the permutation test. "
                       "Setting seed to %s. Record this for future "
                       "replicability" % random_seed)
            if n_perms:
                statement.append(" mperm=%i --seed %s " % (n_perms,
                                                           random_seed))
            else:
                statement.append(" perm --seed %s " % (random_seed))
        else:
            pass

        # if using linear or logistic, covariates can be added into the model
        # to adjust for their effects - assumes fixed effects of covariates
        # mixed models are not yet implemented in Plink2.
        if covariates:
            covars = covariates.split(",")
            if len(covars) > 1:
                if type(covars[0]) == str:
                    m_covar = " --covar-name %s " % covariates

                elif type(covars[0]) == int:
                    m_covar = " --covar-number %s " % covariates
                else:
                    # if none are specified then don't adjust the model for any
                    # and log a warning
                    E.warn("Covariate header or numbers are not recognised."
                           "No covariates will be included in the model.  Please"
                           "specifiy them exactly")
                    covariates = None
                    covariates_file = None
            elif len(covars) == 1:
                if type(covars) == str:
                    m_covar = " --covar-name %s " % covariates

                elif type(covars) == int:
                    m_covar = " --covar-number %i " % covariates
                else:
                    # if none are specified then don't adjust the model for any
                    # and log a warning
                    E.warn("Covariate header or numbers are not recognised."
                           "No covariates will be included in the model.  Please"
                           "specifiy them exactly")
                    covariates = None
                    covariates_file = None

        if covariates and covariates_file:
            statement.append(" --covar %s %s " % (covariates_file,
                                                  m_covar))
        elif covariates and not covaries_file:
            E.warn("No covariate file specified.  None included in model.")
        elif covariates_file and not covariates:
            E.warn("No covariates specified to include in the model."
                   "None included")
        else:
            pass

        self.statement["assoc"] = " ".join(statement)

    def PCA(self, n_pcs="20"):
        '''
        Perform PCA analysis on previosly generated GRM, output the number n
        principal componets, default = 20
        '''

        self._run_tasks(pca=n_pcs)

    def _dimension_reduction(self, **kwargs):
        '''
        Use PCA to perform dimensionality reduction on
        input samples.  A PCA  can be calculated using
        a subset of samples which can then be projected on
        to other samples.
        '''

        # FINISH ME!!!!

    def _detect_interactions(self, method=None, modifier=None,
                             set_file=None, set_mode=None,
                             report_threshold=None,
                             sig_threshold=None,
                             covariates_file=None, covariates=None):
        '''
        Detect epistatic interactions between SNPs using either an inaccurate
        scan (fast-epistasis) or a fully saturated linear model

        Methods
        -------
        fast_epistasis - uses an "imprecise but fast" scan of all 3x3 joint genotype
        count tables to test for interactions.  Can be modified to use a likelihood
        ration test `boost` or a joint-effects test `joint-effects`. Default is
        `joint-effects`.

        epistasis - uses a linear model to test for interactions between additive
        effects after main effects.  Logistic regression for case/control and
        linear regression for quantitative traits.

        two_locus - tests a single interaction between two variants using joint genotype
        counts and frequencies.

        adjusted - allows adjustment for covariates in the interaction test, and also adjusts
        for main effects from both the test and target SNP.  Requires and R plugin script.
        '''

        interact_map = {"fast_epistasis": " --fast-epistasis %s ",
                        "epistasis": " --epistasis %s ",
                        "two_locus": " --twolocus %s ",
                        "adjusted": " --R %s "}

        statement = []

        if modifier:
            statement.append(interact_map[method] % modifier)
        else:
            modifier = ""
            statement.append(interact_map[method] % modifier)

        if covariates_file:
            statement.append("--covar %s --covar-name %s " % (covariates_file,
                                                              covariates))
        else:
            pass

        if set_mode and set_file:
            # does not work with two-locus test
            if method == "two_locus" and set_mode:
                E.warn("Two locus test cannot be used in conjunction "
                       "with a set-based test.")
            elif set_mode:
                statement.append(" %s  --set %s " % (set_mode, set_file))
            else:
                pass
        else:
            pass

        # alter reporting of significant interactions and significance
        # level of interactions
        if report_threshold:
            statement.append(" --epi1 %0.3f " % float(report_threshold))
        else:
            pass

        if sig_threshold:
            statement.append(" --epi2 %0.3f " % float(sig_threshold))
        else:
            pass

        self.statement["epistasis"] = " ".join(statement)

    def _matrices(self, matrix_type, shape="triangle", compression=None, options=None):
        '''
        Calculate a number of different distance matrices:
        realised genetic relationship matrix
        relationship covariance matrix
        identity by descent/state matrix
        hamming distance matrix

        * matrix_type - matrix to compute.  Can be either IBS, 1 - IBS,
          Hamming, GRM
        '''

        statement = []
        if matrix_type == "hamming":
            flag = " --distance "
        elif matrix_type == "ibs":
            flag = " --distance ibs "
        elif matrix_type == "genomic":
            flag = " --distance 1-ibs "
        elif matrix_type == "grm":
            flag = " --make-grm-bin "

        if options:
            statement.append(" ".join([flag, shape, compression, options]))
        elif matrix_type == "grm":
            statement.append(flag)
        else:
            statement.append(" ".join([flag, shape, compression]))

        return " ".join(statement)

    def _qc_methods(self, parameter=None, **kwargs):
        ''''
        Perform QC on genotyping data, SNP-wise and sample-wise.
        All arguments are passed as key word arguments, except
        cases detailed in `Parameters` where they are passed with
        the ``parameter`` argument.

        Methods
        -------
        * ld_prune - generate a list of SNPs in linkage equilibrium by
          pruning SNPs on either an LD statistic threshold, i.e. r^2,
          or use a variance inflation factor (VIF) threshold
        * heterozygosity - calculate average heterozygosity from each
          individual across a set of SNPs, threshold on individuals
          with deviation from expected proportions
        * ibd - calculate the genetic relationship of individuals to
          infer relatedness between individuals, threshold on given
          degree of relatedness, e.g. IBD > 0.03125, 3rd cousins
        * genetic_gender - estimate the gender of an individual
          from the X chromosome genotypes - correlate with reported
          gender and output discrepancies
        * ethnicity_pca - perform PCA using a subset of independent
          SNPs to infer genetic ancestry.  Compare and contrast this
          to individuals reported ancestry.  Report discrepancies
          and individuals  greater than a threshold distance away
          from a reference population.
        * homozygosity - identifies sets of runs of homozygosity
          within individuals.  These may be indicative of inbreeding,
          systematic genotyping errors or regions under selection.

        Parameters
        ----------
        Method parameters can also be passed through this function
        as keyword=value pairs.
        * ld_prune:
          `kb` - this modifier changes the window resolution to kb
          rather than bp.
          `r2` - the r^2 threshold above which SNPs are to be removed
          `vif` - the VIF threshold over which SNPs will be removed
          `window` - window size to calculate pair-wise LD over
          `step` - step size to advance window by
        '''

        qc_dict = {"ld_prune": {"R2": " --indep-pairwise %s %s %s ",
                                "VIF": " --indep %s %s %s "},
                   "heterozygosity": {"gz": " --het gz",
                                      "raw": " --het "},
                   "ibd": {"relatives": " --genome gz rel-check ",
                           "full": " --genome gz full ",
                           "norm": " --genome gz "},
                   "genetic_gender": "none",
                   "ethnicity_pca": "none",
                   "homozygosity": {"min_snp": " --homozyg-snp %s ",
                                    "min_kb": " --homozyg-kb %s ",
                                    "default": " --homozyg ",
                                    "density": " --homozyg-density ",
                                    "set_gap": " --homozyg-gap ",
                                    "snp_window": " --homozyg-window-snp %s ",
                                    "het_max": " --homozyg-het %s "}}

        task_dict = {}
        state = []

        # put everything in an accessible dictionary first
        for task, value in kwargs.items():
            task_dict[task] = value

        # LD pruning can be passed multiple parameters,
        # handle this separately
        try:
            sub_task = task_dict["ld_prune"]
            ld_prune_task = qc_dict["ld_prune"]

            try:
                step = task_dict["step"]
            except KeyError:
                raise AttributeError("No step size found, please "
                                     "pass a step size to advance the "
                                     "window by")
            try:
                window = task_dict["window"]
                try:
                    task_dict["kb"]
                    window = "".join([window, "kb"])
                    task_dict.pop("kb", None)
                except KeyError:
                    pass

            except KeyError:
                raise AttributeError("No window size found.  Please input "
                                     "a window size to prune over")
            try:
                threshold = task_dict["threshold"]
            except KeyError:
                raise AttributeError("No threshold value, please input "
                                     "a value to LD prune SNPs on")

            # add in the kb if it is passed as an argument
            state.append(ld_prune_task[sub_task] % (window, step, threshold))

            task_dict.pop("threshold", None)
            task_dict.pop("ld_prune", None)
            task_dict.pop("window", None)
            task_dict.pop("step", None)

        except KeyError:
            pass

        for task, value in task_dict.items():
            try:
                sub_task = qc_dict[task]
                try:
                    state.append(sub_task[value] % parameter)
                except TypeError:
                    state.append(sub_task[value])
            except KeyError:
                raise AttributeError("Task not found, please see "
                                     "documentation for available features")

        self.statement["QC"] = " ".join(state)

    def build_statement(self, infiles, outfile, threads=None,
                        memory="60G", parallel=None):
        '''
        Build statement and execute from components
        '''

        statement = []
        exec_state = self.executable
        # calls to function add to the self.statement dictionary
        try:
            statement.append(self.statement["program"])
        except KeyError:
            raise AttributeError("Input files and format not detected")

        try:
            statement.append(self.statement["QC"])
        except KeyError:
            pass

        try:
            statement.append(self.statement["filters"])
        except KeyError:
            pass

        try:
            statement.append(self.statement["tasks"])
        except KeyError:
            pass

        try:
            statement.append(self.statement["stats"])
        except KeyError:
            pass

        try:
            statement.append(self.statement["assoc"])
        except KeyError:
            pass

        try:
            statement.append(self.statement["matrix"])
        except KeyError:
            pass

        try:
            statement.append(self.statement["epistasis"])
        except KeyError:
            pass

        if threads:
            statement.append(" --threads %i " % threads)
        else:
            pass

        if not memory:
            pass
        elif memory != "60G":
            memory = int(memory.strip("G")) * 1000
            statement.append(" --memory %i " % memory)
        else:
            statement.append(" --memory 60000 ")

        # add output flag
        # outfile needs to be complete path for Plink to save
        # results properly - check if it starts with '/',
        # if so is already a full path
        if not parallel:
            if os.path.isabs(outfile):
                statement.append(" --out %s " % outfile)
            else:
                outpath = "/".join([os.getcwd(), outfile])
                statement.append(" --out %s " % outpath)

            os.system(" ".join(statement))
        else:
            # parallelisation only really applies to GRM calculation
            # at the moment <- need to generalise
            # if parallelisation is used, invoke temp files
            # then agglomerate files

            statements = []

            if os.path.isabs(outfile):
                outpath = outfile
            else:
                outpath = "/".join([os.getcwd(), outfile])

            for i in range(1, parallel+1):
                # copy list, assigning just makes a pointer
                p_state = statement[:]
                p_state.append(" --parallel %i %i " % (i, parallel))
                p_state.append(" --out %s.%i " % (outpath, i))
                statements.append(" ".join(p_state))

            os.system(";".join(statements))


class PlinkDev(Plink2):
    '''
    Run various Plink functions and analysis, including file processing, GRM
    calculation, PCA and other GWA tasks

    Require Plink v1.9_devel to be in the users PATH variable as ``plinkdev`` to
    distinguish it from Plink v1.07 and v1.9.
    Currently uses Nov 11 development build.
    '''

    def __init__(self, files, options=None,
                 settings=None, design=None):
        self.infiles = files
        self.options = options
        self.settings = settings
        self.design = design
        self.executable = "plinkdev"
        self.statement = {}
        self.filters = []


class GWASResults(object):
    '''
    A class for handling the results from a GWA, used for plotting
    and post-analysis QC
    '''

    def __init__(self, assoc_file, **kwargs):
        # if the assoc_file is a list of multiple files,
        # then merge them into a single files
        if type(assoc_file) == list and len(assoc_file) > 1:
            E.info("multiple results files detected")
            self.infiles = assoc_file
            self.infile = None
            self.results = self.parse_genome_wide(assoc_file)
        else:
            E.info("single results file detected")
            self.infile = assoc_file
            self.infiles = None
            # results is a pandas dataframe to operate on
            self.results = self.get_results(assoc_file, **kwargs)

    def parse_genome_wide(self, association_files):
        '''
        Accept a list of results files, merge them together
        and output as a single dataframe

        Will this take a lot of memory??
        '''

        file0 = association_files.pop(0)
        df = self.get_results(file0)

        for afile in association_files:
            _df = self.get_results(afile)
            df = df.append(_df)

        df["CHR"] = df["CHR"].astype(np.int64)
        df.sort_values(by=["CHR", "BP"], inplace=True)

        return df

    def get_results(self, association_file,
                    epistasis=False,
                    file_format="plink"):
        '''
        Parse a GWA or epistasis results file and return the table
        '''

        # use Pandas for now - try something different later
        # SQLite DB maybe?
        # inconsistent number of white spaces between
        # fields means Pandas parsing breaks down
        # fields need to be the correct data type,
        # i.e. BP = int, P = float, SNP = str, etc
        # if the file has already been parsed and processed
        # just assign it instead
        # epistasis results don't have a header

        try:
            peek = pd.read_table(association_file, nrows=5,
                                 sep="\s*", header=0,
                                 index_col=None,
                                 engine='python')
        except StopIteration:
            peek = pd.read_table(association_file, nrows=5,
                                 sep="\t", header=0,
                                 index_col=None)

        if epistasis:
            try:
                results_frame = pd.read_table(association_file,
                                              sep="\s*", header=0,
                                              index_col=None)
            except StopIteration:
                results_frame = pd.read_table(association_file,
                                              sep="\t", header=0,
                                              index_col=None)

            # results from fast epistasis are different to others
            if file_format == "cassi_covar":
                if results_frme.shape[1] == 12:
                    results_frame.columns = ["SNP1", "CHR1", "ID1", "BP1",
                                             "SNP2", "CHR2", "ID2", "BP2",
                                             "OR", "SE", "STAT", "P"]
                elif results_frame.shape[1] == 14:
                    results_frame.columns = ["SNP1", "CHR1", "ID1", "BP1",
                                             "SNP2", "CHR2", "ID2", "BP2",
                                             "OR", "SE", "STAT", "P",
                                             "CASE_RSQ", "CTRL_RSQ"]
                elif results_frame.shape[1] == 16:
                    results_frame.columns = ["SNP1", "CHR1", "ID1", "BP",
                                             "SNP2", "CHR2", "ID2", "BP2",
                                             "OR", "SE", "STAT", "P",
                                             "CASE_RSQ", "CTRL_RSQ",
                                             "CASE_DPRIME" "CTRL_DPRIME"]
                results_frame.loc[:, "BP"] = pd.to_numeric(results_frame["BP"],
                                                           errors="coerce")
            elif file_format == "cassi":
                pass
            elif file_format == "plink":
                if results_frame.shape[1] == 7:
                    results_frame.columns = ["CHR1", "SNP1", "CHR",
                                             "SNP", "OR", "STAT", "P"]
                elif results_frame.shape[1] == 9:
                    results_frame.columns = ["CHR", "SNP", "BP", "A1", "NMISS",
                                             "OR", "SE", "STAT", "P"]
                else:
                    results_frame.columns = ["CHR", "SNP", "BP", "A1", "OR",
                                             "SE", "STAT", "P"]

                results_frame.loc[:, "BP"] = pd.to_numeric(results_frame["BP"],
                                                           errors="coerce")
            results_frame.loc[:, "P"] = pd.to_numeric(results_frame["P"],
                                                      errors="coerce")

            return results_frame

        else:
            try:
                assert peek["log10P"].any()
                results_frame = pd.read_table(association_file,
                                              sep="\t", header=0,
                                              index_col=None,
                                              dtype={"BP": np.int64,
                                                     "NMISS": np.int64})
                return results_frame
            except KeyError:
                pass

        l_count = 0
        E.info("parsing file: %s" % association_file)
        with open(association_file, "r") as ifile:
            for line in ifile:
                # check if spacing is whitespace or tab
                if len(line.split(" ")) > 1:
                    parsed = line.split(" ")
                elif len(line.split("\t")) > 1:
                    parsed = line.split("\t")
                else:
                    raise IOError("file separator not recognised. "
                                  "Must be whitespace or tab")
                # remove multiple blank spaces
                for i in range(parsed.count('')):
                    parsed.remove('')
                # get rid of the newline
                try:
                    parsed.remove('\n')
                except ValueError:
                    parsed = [(px).rstrip("\n") for px in parsed]
                if l_count == 0:
                    header = [iy.upper() for ix, iy in enumerate(parsed)]
                    head_idx = [ix for ix, iy in enumerate(parsed)]
                    map_dict = dict(zip(head_idx, header))
                    res_dict = dict(zip(header, [[] for each in header]))
                    l_count += 1
                else:
                    col_idx = [lx for lx, ly in enumerate(parsed)]
                    col = [ly for lx, ly in enumerate(parsed)]
                    for i in col_idx:
                        res_dict[map_dict[i]].append(col[i])
                    l_count += 1

        # substract one from the index for the header column
        df_idx = range(l_count-1)

        results_frame = pd.DataFrame(res_dict, index=df_idx)
        results_frame.fillna(value=1.0, inplace=True)
        try:
            results_frame = results_frame[results_frame["TEST"] == "ADD"]
        except KeyError:
            pass

        # need to handle NA as strings
        results_frame["P"][results_frame["P"] == "NA"] = 1.0

        results_frame["BP"] = [int(bx) for bx in results_frame["BP"]]
        results_frame["P"] = [np.float64(fx) for fx in results_frame["P"]]
        try:
            results_frame["STAT"][results_frame["STAT"] == "NA"] = 1.0
            results_frame["STAT"] = [np.float64(sx) for sx in results_frame["STAT"]]
        except KeyError:
            try:
                results_frame["CHISQ"][results_frame["CHISQ"] == "NA"] = 1.0
                results_frame["CHISQ"] = [np.float64(sx) for sx in results_frame["CHISQ"]]
            except KeyError:
                try:
                    results_frame["T"][results_frame["T"] == "NA"] = 1.0
                    results_frame["T"] = [np.float64(sx) for sx in results_frame["T"]]
                except KeyError:
                    pass

        try:
            results_frame["F_U"][results_frame["F_U"] == "NA"] = 0.0
            results_frame["F_U"] = [np.float64(ux) for ux in results_frame["F_U"]]
        except KeyError:
            pass

        try:
            results_frame["F_A"][results_frame["F_A"] == "NA"] = 0.0
            results_frame["F_A"] = [np.float64(ax) for ax in results_frame["F_A"]]
        except KeyError:
            pass

        try:
            results_frame["FREQ"][results_frame["FREQ"] == "NA"] = 0.0
            results_frame["FREQ"] = [np.float64(fx) for fx in results_frame["FREQ"]]
        except KeyError:
            pass

        try:
            results_frame["OR"][results_frame["OR"] == "NA"] = 1.0
            results_frame["OR"] = [np.float64(ox) for ox in results_frame["OR"]]
        except KeyError:
            try:
                results_frame["BETA"][results_frame["BETA"] == "NA"] = 1.0
                results_frame["BETA"] = [np.float64(ox) for ox in results_frame["BETA"]]
            except KeyError:
                results_frame["B"][results_frame["B"] == "NA"] = 0.0
                results_frame["B"] = [np.float64(ox) for ox in results_frame["B"]]

        return results_frame

    def plotManhattan(self, save_path, resolution="chromosome",
                      write_merged=True, sig_level=8):
        '''
        Generate a basic manhattan plot of the association results
        Just deal with chromosome-by-chromosome for now.
        '''

        # use the python ggplot plotting package
        # need to calculate -log10P values separately
        self.results["log10P"] = np.log10(self.results["P"])

        # or using rpy2
        py2ri.activate()
        R('''suppressPackageStartupMessages(library(ggplot2))''')
        R('''suppressPackageStartupMessages(library(scales))''')
        R('''suppressPackageStartupMessages(library(qqman))''')
        R('''sink(file="sink.text")''')
        r_df = py2ri.py2ri_pandasdataframe(self.results)
        R.assign("assoc.df", r_df)
        if resolution == "chromosome":
            R('''assoc.df$CHR <- factor(assoc.df$CHR, '''
              '''levels=levels(ordered(unique(assoc.df$CHR))),'''
              '''labels=unique(paste0("chr", assoc.df$CHR)))''')
            R('''nchrom <- length(unique(assoc.df$CHR))''')
            R('''myCols <- rep(c("#ca0020", "#404040"), nchrom)[1:nchrom]''')
            R('''names(myCols) <- sort(unique(assoc.df$CHR))''')
            R('''colScale <- scale_colour_manual(name = "CHR", values=myCols)''')
            R('''bp_indx <- seq_len(dim(assoc.df[1]))''')
            R('''assoc.df$BPI <- bp_indx''')
            R('''p <- ggplot(assoc.df, aes(x=BPI, y=-log10(P), colour=CHR)) + '''
              '''geom_point(size=1) + colScale + '''
              '''geom_hline(yintercept=6, linetype="dashed", colour="blue") + '''
              '''theme_bw() + labs(x="Chromosome position (bp)", '''
              '''y="-log10 P-value") + facet_grid(~CHR, scale="free_x") + '''
              '''theme(axis.text.x = element_text(size=8))''')
            R('''png("%s", res=90, unit="in", height=8, width=12)''' % save_path)
            R('''print(p)''')
            R('''dev.off()''')

        elif resolution == "genome_wide":
            R('''nchroms <- length(unique(assoc.df$CHR))''')
            R('''png("%s", width=720, height=540)''' % save_path)
            R('''p <- manhattan(assoc.df, main="Manhattan plot",'''
              '''ylim=c(0, 50), cex=0.9, suggestiveline=T,'''
              '''genomewideline=-log10(5e-8), chrlabs=c(1:nchroms), '''
              '''col=c("#8B1A1A","#8470FF"))''')
            R('''print(p)''')
            R('''dev.off()''')

        R('''sink(file=NULL)''')

        if write_merged:
            return self.results
        else:
            return False

    def plotQQ(self, save_path, resolution="chromosome"):
        '''
        Generate a QQ-plot of expected vs. observed
        test statistics
        '''

        self.results["log10P"] = np.log(self.results["P"])

        py2ri.activate()
        R('''suppressPackageStartupMessages(library(ggplot2))''')
        R('''suppressPackageStartupMessages(library(scales))''')
        R('''suppressPackageStartupMessages(library(qqman))''')
        r_df = py2ri.py2ri_pandasdataframe(self.results)
        R.assign("assoc.df", r_df)
        R('''png("%s", width=720, height=540)''' % save_path)
        R('''qq(assoc.df$P)''')
        R('''dev.off()''')

    def plotEpistasis(self, save_path, resolution="chromosome"):
        '''
        Generate both manhattan plot of the SNPs tested for
        epistasis with their target SNP, and a QQplot
        of the association test p-values
        '''

        # plot QQplot
        qq_save = "_".join([save_path, "qqplot.png"])
        self.plotQQ(qq_save)

        manhattan_save = "_".join([save_path, "manhattan.png"])
        self.plotManhattan(manhattan_save,
                           resolution=resolution,
                           sig_level=6,
                           write_merged=False)

    def getHits(self, threshold=0.00000005):
        '''
        Pull out regions of association by selecting
        all SNPs with association p-values less than
        a certain threshold.  Defaults is genome-wide
        signifance, p < 5x10-8.
        Then select region +/- 1.5Mb of the index SNP.
        '''

        hits_df = self.results[self.results["P"] <= threshold]
        # find the range of SNPs with 3Mb of each index SNP
        contig_group = hits_df.groupby(["CHR"])

        # there may be multiple independent hits on a given
        # chromosome.  Need to identify independent regions.
        # Independent regions are defined by their statistical
        # independence, not distance.  Just take all SNPs
        # in 3Mb of the lead SNP for each signal
        # this will create overlaps of associatation signals
        for contig, region in contig_group:
            region.index = region["BP"]
            chr_df = self.results[self.results["CHR"] == contig]
            chr_df.index = chr_df["BP"]
            # find independent regions and output consecutively
            # if only a single SNP above threshold then there is
            # only one independent region!!
            if len(region) > 1:
                independents = self.findIndependentRegions(region)
                indi_group = independents.groupby("Group")
            else:
                region["Group"] = 1
                indi_group = region.groupby("Group")

            for group, locus in indi_group:
                # if there is only a single variant should
                # the region be kept?  Likely a false
                # positive
                if min(locus["BP"]) == max(locus["BP"]):
                    pass
                else:
                    try:
                        try:
                            locus.loc[:, "STAT"] = abs(locus["STAT"])
                            locus.sort_values(by="STAT", inplace=True)
                        except KeyError:
                            locus.loc[:, "T"] = abs(locus["T"])
                            locus.sort_values(by="STAT", inplace=True)
                    except KeyError:
                        locus.sort_values(by="CHISQ", inplace=True)

                    index_bp = locus.iloc[0]["BP"]
                    E.info("Lead SNP for regions is: {}".format(locus.iloc[0]["SNP"]))
                    left_end = min(chr_df.loc[chr_df.index >= index_bp - 1500000, "BP"])
                    right_end = max(chr_df.loc[chr_df.index <= index_bp + 1500000, "BP"])

                    range_df = chr_df.loc[left_end: right_end, :]
                    max_stat = max(abs(range_df["STAT"]))

                    yield contig, range_df

    def extractSNPs(self, snp_ids):
        '''
        Extract a specific set of SNP results

        Arguments
        ---------
        snp_ids: list
          a list of SNP IDs to extract from the
          GWAS results

        Returns
        -------
        snp_results: pandasd.Core.DataFrame
        '''

        self.results.index = self.results["SNP"]

        snp_results = self.results.loc[snp_ids]

        return snp_results

    def findIndependentRegions(self, dataframe):
        '''
        Find the number of independent regions on
        a chromsome.  Uses R distance and tree
        cutting functions
        '''

        # mong dataframe into R
        py2ri.activate()
        r_df = py2ri.py2ri_pandasdataframe(dataframe)
        R.assign("rdf", r_df)
        R('''mat <- as.matrix(rdf$BP)''')
        # get distances then cluster, chop tree at 1x10^7bp
        R('''dist.mat <- dist(mat, method="euclidean")''')
        R('''clusts <- hclust(dist.mat, "average")''')
        R('''cut <- cutree(clusts, h=1e6)''')
        R('''out.df <- rdf''')
        R('''out.df$Group <- cut''')

        # need to handle changes in pandas2ri API
        try:
            regions_df = pd.DataFrame(py2ri.ri2py(R["out.df"]))
        except NotImplementedError:
            regions_df = pd.DataFrame(R["out.df"])

        return regions_df

    def mergeFrequencyResults(self, freq_dir, file_regex):
        '''
        Merge GWAS results with frequency information,
        and format for GCTA joint analysis input
        '''

        # create a dummy regex to compare
        # file_regex type against
        test_re = re.compile("A")

        if type(file_regex) == str:
            file_regex = re.compile(file_regex)

        elif type(file_regex) == type(test_re):
            pass

        else:
            raise TypeError("Regex type not recognised.  Must"
                            "be string or re.SRE_Pattern")

        all_files = os.listdir(freq_dir)
        freq_files = [fx for fx in all_files if re.search(file_regex, fx)]

        gwas_df = self.results

        df_container = []
        for freq in freq_files:
            freq_file = os.path.join(freq_dir, freq)
            E.info("Adding information from {}".format(freq_file))
            # files may or may not be tab-delimited
            try:
                _df = pd.read_table(freq_file,
                                    sep="\s*", header=0,
                                    index_col=None,
                                    engine='python')
            except StopIteration:
                _df = pd.read_table(freq_file,
                                    sep="\t", header=0,
                                    index_col=None)

            merge_df = pd.merge(self.results, _df,
                                left_on=["CHR", "SNP"],
                                right_on=["CHR", "SNP"],
                                how='left')
            df_container.append(merge_df)

        count = 0
        for df in df_container:
            if not count:
                gwas_df = df
                count += 1
            else:
                gwas_df = gwas_df.append(df)

        E.info("Calculating Z scores and SEs")
        z_scores = -0.862 + np.sqrt(0.743 - 0.2404 *
                                    np.log(gwas_df.loc[:, "P"]))
        se = np.log(gwas_df.loc[:, "OR"])/z_scores
        gwas_df.loc[:, "Z"] = z_scores
        gwas_df.loc[:, "SE"] = se
        gwas_df.loc[:, "logOR"] = np.log(gwas_df.loc[:, "OR"])

        out_cols = ["SNP", "A1_x", "A2", "MAF", "logOR", "SE", "P", "NMISS"]

        out_df = gwas_df[out_cols]

        # need to remove duplicates, especially those
        # that contain NaN for A2 and MAF

        out_df = out_df.loc[~np.isnan(out_df["MAF"])]

        return out_df


##########################################################
# unbound methods that work on files and data structures #
##########################################################

def plotMapPhenotype(data, coords, coord_id_col, lat_col,
                     long_col, save_path, xvar, var_type,
                     xlabels=None, level=None):
    '''
    Generate a map of the UK, with phenotype data overlaid
    '''

    # merge co-ordinate data with phenotype data
    merged_df = pd.merge(left=coords, right=data, left_on=coord_id_col,
                         right_on=coord_id_col, how='inner')

    # pheno column and set level of categorical variable
    if xlabels and var_type == "categorical":
        # convert to string type as a categorical variable
        # drop NA observations from the merged data frame
        na_mask = pd.isnull(merged_df.loc[:, xvar])
        merged_df = merged_df[~na_mask]

        rvar = merged_df.loc[:, xvar].copy()
        nvar = pd.Series(np.nan_to_num(rvar), dtype=str)
        var = [v for v in set(nvar)]
        var.sort()
        # recode the variables according to the input labels
        xlabs = xlabels.split(",")
        lbls = [str(xlabs[ix]) for ix in range(len(var))]
        for xv in range(len(var)):
            nvar[nvar == var[xv]] = lbls[xv]
        merged_df.loc[:, "cat_var"] = nvar

    else:
        pass

    if level:
        lvar = merged_df.loc[:, "cat_var"].copy()
        mask = lvar.isin([level])
        lvar[mask] = 1
        lvar[~mask] = 0
        lvar = lvar.fillna(0)
        merged_df.loc[:, "dichot_var"] = lvar
    else:
        pass

    # push the df into the R env
    py2ri.activate()
    r_df = py2ri.py2ri_pandasdataframe(merged_df)
    R.assign("pheno.df", r_df)

    # setup the map and plot the points
    R('''suppressPackageStartupMessages(library(maps))''')
    R('''suppressPackageStartupMessages(library(mapdata))''')

    R('''uk_map <- map("worldHires", c("UK", "Isle of Wight",'''
      '''"Ireland", "Isle of Man", "Wales:Anglesey"), '''
      '''xlim=c(-11, 3), ylim=c(50, 60.9), plot=F)''')
    # colour by reference, or a colour for each discrete value
    if level:
        R('''red <- rep("#FF0000", '''
          '''times=length(pheno.df$dichot_var[pheno.df$dichot_var == 1]))''')
        R('''black <- rep("#000000", '''
          '''times=length(pheno.df$dichot_var[pheno.df$dichot_var == 0]))''')

        R('''png("%(save_path)s", width=540, height=540, res=90)''' % locals())
        R('''map(uk_map)''')

        R('''points((-pheno.df[,"%(lat_col)s"])[pheno.df$dichot_var == 1], '''
          '''(-pheno.df[,"%(long_col)s"])[pheno.df$dichot_var == 1], pch=".", col=red)''' % locals())

        R('''points((pheno.df[,"%(long_col)s"])[pheno.df$dichot_var == 0], '''
          '''(pheno.df[,"%(lat_col)s"])[pheno.df$dichot_var == 0], pch=".", col=black)''' % locals())

        R('''legend('topleft', legend=c("not-%(level)s", "%(level)s"),'''
          '''fill=c("#000000", "#FF0000"))''' % locals())
        R('''dev.off()''')
    else:
        R('''png("%(save_path)s", width=540, height=540, res=90)''' % locals())
        R('''map(uk_map)''')

        R('''points(pheno.df[,"%(long_col)s"], pheno.df[,"%(lat_col)s"], pch=".", '''
          '''col=factor(pheno.df$cat_var))''' % locals())

        R('''legend('topleft', legend=unique(pheno.df$cat_var),'''
          '''fill=unique(pheno.df$cat_var))''' % locals())
        R('''dev.off()''')


def plotPhenotype(data, plot_type, x, y=None, group=None,
                  save_path=None, labels=None, xlabels=None,
                  ylabels=None, glabels=None, var_type="continuous"):
    '''
    Generate plots of phenotypes using ggplot
    '''

    # change data format if necessary and convert nan/NA to missing
    if not y and var_type == "categorical":
        var = np.nan_to_num(data.loc[:, x].copy())
        data.loc[:, x] = pd.Series(var, dtype=str)
        if group:
            gvar = np.nan_to_num(data.loc[:, group].copy())
            data.loc[:, group] = pd.Series(gvar, dtype=str)
        else:
            pass

    elif not y and var_type == "integer":
        var = np.nan_to_num(data.loc[:, x].copy())
        data.loc[:, x] = pd.Series(var, dtype=np.int64)
        if group:
            gvar = np.nan_to_num(data.loc[:, group].copy())
            data.loc[:, group] = pd.Series(gvar, dtype=str)
        else:
            pass

    elif not y and var_type == "continuous":
        var = data.loc[:, x].copy()
        data.loc[:, x] = pd.Series(var, dtype=np.float64)
        if group:
            gvar = np.nan_to_num(data.loc[:, group].copy())
            data.loc[:, group] = pd.Series(gvar, dtype=str)
        else:
            pass

    elif y and var_type == "categorical":
        xvar = np.nan_to_num(data.loc[:, x].copy())
        yvar = np.nan_to_num(data.loc[:, y].copy())

        data.loc[:, x] = pd.Series(xvar, dtype=str)
        data.loc[:, y] = pd.Series(yvar, dtype=str)

        if group:
            gvar = np.nan_to_num(data.loc[:, group].copy())
            data.loc[:, group] = pd.Series(gvar, dtype=str)
        else:
            pass

    elif y and var_type == "integer":
        xvar = np.nan_to_num(data.loc[:, x].copy())
        yvar = np.nan_to_num(data.loc[:, y].copy())

        data.loc[:, x] = pd.Series(xvar, dtype=np.int64)
        data.loc[:, y] = pd.Series(yvar, dtype=np.int64)

        if group:
            gvar = np.nan_to_num(data.loc[:, group].copy())
            data.loc[:, group] = pd.Series(gvar, dtype=str)
        else:
            pass

    elif y and var_type == "continuous":
        # NAs and NaNs should be handled by ggplot
        xvar = data.loc[:, x].copy()
        yvar = data.loc[:, y].copy()

        data.loc[:, x] = pd.Series(xvar, dtype=np.float64)
        data.loc[:, y] = pd.Series(yvar, dtype=np.float64)

        if group:
            gvar = np.nan_to_num(data.loc[:, group].copy())
            data.loc[:, group] = pd.Series(gvar, dtype=str)
        else:
            pass

    R('''suppressPackageStartupMessages(library(ggplot2))''')
    # put the pandas dataframe in to R with rpy2
    py2ri.activate()

    r_df = py2ri.py2ri_pandasdataframe(data)
    R.assign("data_f", r_df)

    # plotting parameters, including grouping variables and labels
    # axis labels
    try:
        labs = labels.split(",")
    except AttributeError:
        labs = []

    # if variable labels have been provided then assume they are
    # categorical/factor variables.
    # assume variable labels are input in the correct order
    if xlabels:
        try:
            unique_obs = len(set(data.loc[:, x]))
            xfact = len(xlabels.split(","))
            if xfact == unique_obs:
                R('''lvls <- unique(data_f[,"%(x)s"])''' % locals())
                lbls = ro.StrVector([ri for ri in xlabels.split(",")])
                R.assign("lbls", lbls)
                R('''lvls <- lvls[order(lvls, decreasing=F)]''')
                R('''data_f[,"%(x)s"] <- ordered(data_f[,"%(x)s"], '''
                  '''levels=lvls, labels=lbls)''' % locals())
            else:
                E.warn("the number of labels does not match the "
                       "number of unique observations, labels not "
                       "used.")
        except AttributeError:
            xlabels = None

    else:
        pass

    if glabels:
        unique_obs = len(set(data.loc[:, group]))
        gfact = len(glabels.split(","))
        if gfact == unique_obs:
            R('''lvls <- unique(data_f[, "%(group)s"])''' % locals())
            lbls = ro.StrVector([rg for rg in glabels.split(",")])
            R.assign("lbls", lbls)
            R('''lvls <- lvls[order(lvls, decreasing=F)]''')
            R('''data_f[,"%(group)s"] <- ordered(data_f[,"%(group)s"], '''
              '''levels=lvls, labels=lbls)''' % locals())
        else:
            E.warn("the number of labels does not match the "
                   "number of unique observations, labels not "
                   "used.")

    # start constructing the plot
    # if X and Y are specified, assume Y is a variable to colour
    # observations by, unless group is also set.
    # If Y and group then colour by group and split by Y
    if y:
        R('''p <- ggplot(aes(x=%s, y=%s), data=data_f)''' % (x, y))

        if plot_type == "histogram":
            if group:
                R('''p <- p + geom_histogram(aes(colour=%(group)s)) + '''
                  '''facet_grid(. ~ %(y)s)''' % locals())
            else:
                R('''p <- p + geom_histogram(aes(colour=%(y)s))''' % locals())

        elif plot_type == "barplot":
            if group:
                R('''p <- p + geom_bar(aes(colour=%(group)s)) + '''
                  '''facet_grid(. ~ %(y)s)''' % locals())
            else:
                R('''p <- p + geom_bar(aes(colour=%(y)s))''' % locals())

        elif plot_type == "density":
            if group:
                R('''p <- p + geom_density(aes(colour=%(group)s)) + '''
                  '''facet_grid(. ~ %(y)s)''' % locals())
            else:
                R('''p <- p + geom_density(aes(colour=%(y)s))''' % locals())

        elif plot_type == "boxplot":
            if group:
                R('''p <- p + geom_boxplot(group=%(group)s,'''
                  '''aes(x=factor(%(x)s), y=%(y)s, fill=%(group)s))''' % locals())
            else:
                R('''p <- p + geom_boxplot(aes(colour=%(x)s))''' % locals())

        elif plot_type == "scatter":
            if group:
                R('''p <- p + geom_point(size=1, aes(colour=%(group)s))''' % locals())
            else:
                R('''p <- p + geom_point(size=1)''')

        if len(labs) == 1:
            xlab = labs[0]
            R('''p <- p + labs(x="%s")''' % xlab)
        elif len(labs) == 2:
            xlab = labs[0]
            ylab = labs[1]
            R('''p <- p + labs(x="%(xlab)s", y="%(ylab)s")''' % locals())
        elif len(labs) == 3:
            xlab = labs[0]
            ylab = labs[1]
            title = labs[2]
            R('''p <- p + labs(x="%(xlab)s", y="%(ylab)s", '''
              '''title="%(title)s")''' % locals())
        elif len(labs) == 4:
            xlab = labs[0]
            ylab = labs[1]
            glab = labs[2]
            title = labs[3]
            R('''p <- p + labs(x="%(xlab)s", y="%(ylab)s",'''
              '''title="%(title)s")''' % locals())
            # need to add in guide/legend title

    else:
        R('''p <- ggplot(data=data_f)''')

        if plot_type == "histogram":
            if group:
                R('''p <- p + geom_histogram(aes(%(x)s)) + '''
                  '''facet_grid(. ~ %(group)s)''' % locals())
            else:
                R('''p <- p + geom_histogram(aes(%s))''' % x)

        elif plot_type == "barplot":
            if group:
                R(''' p <- p + geom_bar(aes(%(x)s)) + '''
                  '''facet_grid(. ~ %(group)s)''')
            else:
                R('''p <- p + geom_bar(aes(%s))''' % x)

        elif plot_type == "density":
            if group:
                R('''p <- p + geom_density(aes(%(x)s)) + '''
                  '''facet_grid(. ~ %(group)s)''' % locals())
            else:
                R('''p <- p + geom_density(aes(%s))''' % x)

        elif plot_type == "boxplot":
            if group:
                R('''p <- p + geom_boxplot(aes(y=%(x)s, '''
                  '''x=factor(%(group)s)))''' % locals())
            else:
                raise AttributeError("Y or group variable is missing")

        if len(labs) == 1:
            xlab = labs[0]
            R('''p <- p + labs(x="%s")''' % xlab)
        elif len(labs) == 2:
            xlab = labs[0]
            title = labs[1]
            R('''p <- p + labs(x="%(xlab)s", '''
              '''title="%(title)s")''' % locals())
        elif len(labs) == 3:
            if group:
                xlab = labs[0]
                glab = labs[1]
                title = labs[2]
                R('''p <- p + labs(x="%(glab)s", y="%(xlab)s",'''
                  '''title="%(title)s")''' % locals())
            else:
                E.warn("too many labels provided, assume first is X, "
                       "and second is plot title")
                xlab = labs[0]
                title = labs[1]
                R('''p <- p + labs(x="%(xlab)s", '''
                  '''title="%(title)s")''' % locals())

    # the default theme is bw
    R('''p <- p + theme_bw()''')

    R('''png("%(save_path)s")''' % locals())
    R('''print(p)''')
    R('''dev.off()''')


def parseFlashPCA(pcs_file, fam_file):
    '''
    Parse the principal components file from FlashPCA
    and match with individual identifiers.  This
    assumes the output order of FlashPCA is the same
    as the input order in the .fam file
    '''

    try:
        pc_df = pd.read_table(pcs_file, sep="\s*",
                              header=None, index_col=None)
    except StopIteration:
        pc_df = pd.read_table(pcs_file, sep="\t",
                              header=None, index_col=None)

    # add a header to the pc_df file
    headers = ["PC%i" % (n + 1) for n,
               m in enumerate(pc_df.columns)]

    pc_df.columns = headers

    fam_df = pd.read_table(fam_file, sep="\t",
                           header=None, index_col=None)

    fam_df.columns = ["FID", "IID", "PAR", "MAT", "GENDER",
                      "PHENO"]
    pc_df[["FID", "IID"]] = fam_df.iloc[:, :2]

    return pc_df


def plotPCA(data, nPCs, point_labels, save_path,
            headers, metadata=None, multiplot=False):
    '''
    Plot N principal components from a PCA either as
    a single plot of the first 2 PCs, a grid plot of
    N PCs.

    Arguments
    ---------
    data: string
      PATH to file containing principal components

    nPCs: int
      number of principal components to plot.  If this
      value is > 2, then multiplot will be enabled
      automatically

    point_labels: vector
      a vector of labels of length correpsonding to
      the number of rows in the data file.  These are
      used to colour the points in the plot with relevant
      metadata.  Alternatively, can be the column header
      in the metadata file that corresponds to annotations

    save_path: string
      An absolute PATH to save the plot(s) to

    headers: boolean
      whether the `data` file contains header delineating the
      columns

    metadata: string
      file containing metadata to annotate plot with, includes
      point_labels data

    multiplot: boolean
      If True, generate a grid of scatter plots with successive
      PCs plotted against each other

    Returns
    -------
      None
    '''

    py2ri.activate()

    if metadata:
        meta_df = pd.read_table(metadata, sep="\t", header=0,
                                index_col=None)
    else:
        pass

    labels = meta_df[["FID", "IID", point_labels]]
    merged = pd.merge(data, labels, left_on="FID",
                      right_on="FID", how='inner')

    # TO DO: enable multiplotting of many PCs
    r_df = py2ri.py2ri_pandasdataframe(merged)
    R.assign("pc.df", r_df)
    R('''suppressPackageStartupMessages(library(ggplot2))''')
    R('''pc.df[["%(point_labels)s"]] <- as.factor(pc.df[["%(point_labels)s"]])''' % locals())

    R('''p_pcs <- ggplot(pc.df, aes(x=PC1, y=PC2, colour=%s)) + '''
      '''geom_point(size=1) + theme_bw() + '''
      '''labs(x="PC1", y="PC2", title="PC1 vs. PC2 LD trimmed genotypes")''' % point_labels)
    R('''png("%s")''' % save_path)
    R('''print(p_pcs)''')
    R('''dev.off()''')


def countByVariantAllele(ped_file, map_file):
    '''
    Count the number of individuals carrying the variant allele
    for each SNP.
    Count the number of occurences of each allele with the variant
    allele of each other SNP.

    Requires ped file genotyping to be in format A1(minor)=1, A2=2
    '''

    # parse the ped file - get the variant column headers from
    # the map file - no headers with these files
    # variant order in the map file matters, use an ordered dict
    variants = collections.OrderedDict()
    with open(map_file, "r") as mfile:
        for snp in mfile.readlines():
            attrs = snp.split("\t")
            snpid = attrs[1]
            variants[snpid] = {"chr": attrs[0],
                               "pos": attrs[-1].strip("\n")}

    variant_ids = variants.keys()
    # store genotype matrix as an array
    # rows and columns are variant IDs
    homA1 = np.zeros((len(variant_ids), len(variant_ids)),
                     dtype=np.int64)

    homA2 = np.zeros((len(variant_ids), len(variant_ids)),
                     dtype=np.int64)

    het = np.zeros((len(variant_ids), len(variant_ids)),
                   dtype=np.int64)

    tcount = 0
    with open(ped_file, "r") as pfile:
        for indiv in pfile.readlines():
            indiv = indiv.strip("\n")
            indiv_split = indiv.split("\t")
            fid = indiv_split[0]
            iid = indiv_split[1]
            mid = indiv_split[2]
            pid = indiv_split[3]
            gender = indiv_split[4]
            phen = indiv_split[5]
            genos = indiv_split[6:]
            tcount += 1
            # get genotype counts
            for i in range(len(genos)):
                # missing genotypes are coded '00' in plink format
                if genos[i] == "00":
                    pass
                elif genos[i] == "11":
                    homA1[i, i] += 1
                elif genos[i] == "12":
                    het[i, i] += 1
                else:
                    homA2[i, i] += 1
    allele_counts = ((2 * homA2) + het)/float(2 * tcount)
    mafs = 1 - allele_counts.diagonal()
    maf_df = pd.DataFrame(zip(variant_ids, mafs), columns=["SNP", "MAF"],
                          index=[x for x, y in enumerate(variant_ids)])
    maf_df["A2_HOMS"] = (2 * homA1).diagonal()
    maf_df["A2_HETS"] = het.diagonal()
    maf_df.index = maf_df["SNP"]
    maf_df.drop(["SNP"], axis=1, inplace=True)

    E.info("allele frequencies calculated over %i SNPs and "
           "%i individuals" % (len(genos), tcount))

    return maf_df


def calcMaxAlleleFreqDiff(ped_file, map_file, group_file,
                          test=None, ref=None):
    '''
    Calculate the allele frequency difference between
    two groups of individuals based upon some prior
    assignment.

    Arguments
    ---------
    ped_file: string
      plink text format .ped file - see Plink documentation
      for details (https://www.cog-genomics.org/plink2/input#ped)

    map_file: string
      plink test format .map file - see Plink documentation
      for details (https://www.cog-genomics.org/plink2/input#ped)

    group_file: string
      a file containing grouping information, must be in standard
      Plink format with IID, FID, GROUP as the columns

    test: string
      group label to use as the test case for calculating
      allele frequency differences.  If this isn't set, then
      the first non-ref value encountered will be set as test

    ref: string
      group label to use as the reference case for calculating
      allele frequency differences.  If not set, then the first
      value encountered will be the test.

    Returns
    -------
    freq_diffs: pandas.Core.DataFrame
      dataframe of SNP information and allele frequency difference
      between group labels
    '''

    # group labels need to be of the same type, convert all
    # group values to string
    group_df = pd.read_table(group_file, sep="\t", header=0,
                             index_col=None,
                             converters={"GROUP": str,
                                         "FID": str,
                                         "IID": str})

    group_df["GROUP"] = [str(xg) for xg in group_df["GROUP"]]

    try:
        assert ref
        E.info("Reference label set to %s" % ref)
    except AssertionError:
        ref = set(group_df["GROUP"])[0]
        E.info("Reference label not provided.  Setting "
               "reference label to %s" % ref)

    try:
        assert test
        E.info("Test label set to %s" % test)
    except AssertionError:
        test = [tx for tx in set(group_df["GROUP"]) if not ref][0]
        E.info("Test label not provided, setting test "
               "label to %s." % test)

    # parse the ped file - get the variant column headers from
    # the map file - no headers with these files
    # variant order in the map file matters, use an ordered dict
    variants = collections.OrderedDict()
    with open(map_file, "r") as mfile:
        for snp in mfile.readlines():
            attrs = snp.split("\t")
            snpid = attrs[1]
            variants[snpid] = {"chr": attrs[0],
                               "pos": attrs[-1].strip("\n")}

    variant_ids = variants.keys()
    # store genotype matrix as an array
    # rows and columns are variant IDs
    ref_homA1 = np.zeros((len(variant_ids), len(variant_ids)),
                         dtype=np.int64)

    ref_homA2 = np.zeros((len(variant_ids), len(variant_ids)),
                         dtype=np.int64)

    ref_het = np.zeros((len(variant_ids), len(variant_ids)),
                       dtype=np.int64)

    test_homA1 = np.zeros((len(variant_ids), len(variant_ids)),
                          dtype=np.int64)

    test_homA2 = np.zeros((len(variant_ids), len(variant_ids)),
                          dtype=np.int64)

    test_het = np.zeros((len(variant_ids), len(variant_ids)),
                        dtype=np.int64)
    tcount = 0
    rcount = 0
    ncount = 0
    ref_ids = group_df["IID"][group_df["GROUP"] == ref].values
    test_ids = group_df["IID"][group_df["GROUP"] == test].values
    total = len(group_df)

    with open(ped_file, "r") as pfile:
        for indiv in pfile.readlines():
            indiv = indiv.strip("\n")
            indiv_split = indiv.split("\t")
            fid = indiv_split[0]
            iid = indiv_split[1]
            mid = indiv_split[2]
            pid = indiv_split[3]
            gender = indiv_split[4]
            phen = indiv_split[5]
            genos = indiv_split[6:]

            # check for ref and test conditions
            # ignore individuals in neither camp
            if iid in test_ids:
                tcount += 1
                # get genotype counts
                for i in range(len(genos)):
                    # missing genotypes are coded '00' in plink format
                    if genos[i] == "00":
                        pass
                    elif genos[i] == "11":
                        test_homA1[i, i] += 1
                    elif genos[i] == "12":
                        test_het[i, i] += 1
                    else:
                        test_homA2[i, i] += 1

            elif iid in ref_ids:
                rcount += 1
                # get genotype counts
                for i in range(len(genos)):
                    # missing genotypes are coded '00' in plink format
                    if genos[i] == "00":
                        pass
                    elif genos[i] == "11":
                        ref_homA1[i, i] += 1
                    elif genos[i] == "12":
                        ref_het[i, i] += 1
                    else:
                        ref_homA2[i, i] += 1
            else:
                ncount += 1
                if round((tcount + rcount + ncount)/total, 2) == 0.25:
                    E.info("%i samples counted."
                           "Approximately 25% samples counted" % tcount + rcount + ncount)
                elif round((tcount + rcount + ncount)/total, 2) == 0.50:
                    E.info("%i samples counted."
                           "Approximately 50% samples counted" % tcount + rcount + ncount)
                elif round((tcount + rcount + ncount)/total, 2) == 0.75:
                    E.info("%i samples counted."
                           "Approximately 75% samples counted" % tcount + rcount + ncount)

    E.info("Counted alleles for %i test cases, %i ref cases,"
           " %i neither reference nor test." % (tcount, rcount,
                                                ncount))

    ref_allele_counts = ((2 * ref_homA2) + ref_het)/float(2 * rcount)
    test_allele_counts = ((2 * test_homA2) + test_het)/float(2 * tcount)

    ref_mafs = 1 - ref_allele_counts.diagonal()
    test_mafs = 1 - ref_allele_counts.diagonal()

    ref_maf_df = pd.DataFrame(zip(variant_ids, ref_mafs),
                              columns=["SNP", "ref_MAF"],
                              index=[x for x, y in enumerate(variant_ids)])
    ref_maf_df["ref_A2_HOMS"] = (2 * ref_homA1).diagonal()
    ref_maf_df["ref_A2_HETS"] = ref_het.diagonal()
    ref_maf_df.index = ref_maf_df["SNP"]
    ref_maf_df.drop(["SNP"], axis=1, inplace=True)

    test_maf_df = pd.DataFrame(zip(variant_ids, test_mafs),
                               columns=["SNP", "test_MAF"],
                               index=[x for x, y in enumerate(variant_ids)])
    test_maf_df["test_A2_HOMS"] = (2 * test_homA1).diagonal()
    test_maf_df["test_A2_HETS"] = test_het.diagonal()
    test_maf_df.index = test_maf_df["SNP"]
    test_maf_df.drop(["SNP"], axis=1, inplace=True)

    freq_diffs = pd.merge(ref_maf_df, test_maf_df,
                          left_index=True, right_index=True,
                          how='inner')

    freq_diffs["MAF_diff"] = freq_diffs["ref_MAF"] - freq_diffs["test_MAF"]

    E.info("allele frequencies calculated over %i SNPs and "
           "%i individuals" % (len(genos), tcount + rcount))

    return freq_diffs


def calcPenetrance(ped_file, map_file, mafs=None,
                   subset=None, snpset=None):
    '''
    Calculate the proportion of times an allele is observed
    in the phenotype subset vs it's allele frequency.
    This is the penetrance of the allele
    i.e. if observed in 100% of affected individuals and 0%
    of controls, then penetrance is 100%
    Generates a table of penetrances for each variants/allele
    and a plot of MAF vs # cases carrying the allele

    Generates a heatmap of compound heterozygotes, and homozygotes
    with penetrances.

    Outputs a table of SNPs, homozygote and heterozygote counts
    among subset individuals and proportion of subset individual
    phenotype explained by homozygotes and heterozygotes

    Requires alleles are coded A1(minor)=1, A2=2
    '''
    # check subset is set, if not then throw an error
    # cannot calculate penetrance without a phenotype
    if not subset:
        raise ValueError("Cannot calculate penetrance of alleles "
                         "without a phenotype to subset in")
    else:
        pass

    # parse the ped file - get the variant column headers from
    # the map file - no headers with these files
    # variant order in the map file matters, use an ordered dict
    variants = collections.OrderedDict()
    with open(map_file, "r") as mfile:
        for snp in mfile.readlines():
            attrs = snp.split("\t")
            snpid = attrs[1]
            variants[snpid] = {"chr": attrs[0],
                               "pos": attrs[-1].strip("\n")}

    if snpset:
        with IOTools.openFile(snpset, "r") as sfile:
            snps = sfile.readlines()
            snps = [sx.rstrip("\n") for sx in snps]
            variant_ids = [ks for ks in variants.keys() if ks in snps]
    else:
        variant_ids = variants.keys()

    var_idx = [si for si, sj in enumerate(variant_ids)]
    case_mat = np.zeros((len(variant_ids), len(variant_ids)),
                        dtype=np.float64)

    all_mat = np.zeros((len(variant_ids), len(variant_ids)),
                       dtype=np.float64)

    tcount = 0
    ncases = 0

    # missing phenotype individuals must be ignored, else
    # they will cause the number of individuals explained
    # to be underestimated
    with open(ped_file, "r") as pfile:
        for indiv in pfile.readlines():
            indiv = indiv.strip("\n")
            indiv_split = indiv.split("\t")
            fid = indiv_split[0]
            iid = indiv_split[1]
            mid = indiv_split[2]
            pid = indiv_split[3]
            gender = int(indiv_split[4])
            phen = int(indiv_split[5])
            if phen != -9:
                if subset == "cases":
                    select = phen
                elif subset == "gender":
                    select = gender
                else:
                    select = None
                genos = np.array(indiv_split[6:])
                genos = genos[var_idx]
                tcount += 1

                het = np.zeros(len(genos), dtype=np.float64)
                hom = np.zeros(len(genos), dtype=np.float64)

                for i in range(len(genos)):
                    # missing values are coded '00' in plink format
                    # A2 homs are coded '11' in plink format
                    if genos[i] == "11":
                        hom[i] += 1
                    elif genos[i] == "12":
                        het[i] += 1
                    else:
                        pass

                hom_mat = np.outer(hom, hom)
                het_mat = np.outer(het, het)
                homs = hom_mat.diagonal()
                het_mat[np.diag_indices(len(genos))] = homs

                gen_mat = het_mat

                # separate matrix for subset
                # reference is always level 2 for plink files,
                # either cases or females
                if select == 2:
                    case_mat += gen_mat
                    all_mat += gen_mat
                    ncases += 1
                else:
                    all_mat += gen_mat
            else:
                pass

    E.info("alleles counted over %i SNPs "
           "and %i individuals, of which %i are "
           "in the %s subset" % (len(genos), tcount, ncases, subset))

    penetrance = np.divide(case_mat, all_mat)
    # round for the sake of aesthetics
    penetrance = np.round(penetrance, decimals=5)
    pen_df = pd.DataFrame(penetrance, columns=variant_ids,
                          index=variant_ids)
    pen_df = pen_df.fillna(0.0)

    case_df = pd.DataFrame(case_mat, columns=variant_ids,
                           index=variant_ids)
    all_df = pd.DataFrame(all_mat, columns=variant_ids,
                          index=variant_ids)
    # plot heatmap of penetrances as percentages
    indf = pen_df * 100
    py2ri.activate()

    # only plot penetrances > 0%
    r_pen = py2ri.py2ri_pandasdataframe(indf)
    r_cases = py2ri.py2ri_pandasdataframe(case_df)
    r_all = py2ri.py2ri_pandasdataframe(all_df)
    R.assign("pen.df", r_pen)
    R.assign("case.df", r_cases)
    R.assign("all.df", r_all)
    R('''suppressPackageStartupMessages(library(gplots))''')
    R('''suppressPackageStartupMessages(library(RColorBrewer))''')

    # penetrances
    E.info("plotting penetrance matrix")
    R('''hmcol <- colorRampPalette(brewer.pal(9, "BuGn"))(100)''')
    R('''rowpen <- pen.df[rowSums(pen.df) > 0,]''')
    R('''colpen <- rowpen[,colSums(rowpen) > 0]''')
    R('''png("%s/penetrance-matrix.png", width=720, height=720)''' % os.getcwd())
    R('''heatmap.2(as.matrix(colpen), trace="none", col=hmcol,'''
      '''dendrogram="none", Colv=colnames(colpen), key=FALSE, '''
      '''Rowv=rownames(colpen), margins=c(10,10), cellnote=round(colpen),'''
      '''notecol="white")''')
    R('''dev.off()''')

    E.info("plotting case counts matrix")
    R('''rowcase <- case.df[rowSums(case.df) > 0,]''')
    R('''colcase <- rowcase[,colSums(rowcase) > 0]''')
    R('''png("%s/cases-matrix.png", width=720, height=720)''' % os.getcwd())
    R('''heatmap.2(as.matrix(colcase), trace="none", col=rep("#F0F8FF", 100),'''
      '''dendrogram="none", Colv=colnames(colcase), key=FALSE, '''
      '''colsep=seq(1:length(colnames(colcase))), '''
      '''rowsep=seq(1:length(rownames(colcase))),'''
      '''Rowv=rownames(colcase), margins=c(10,10), cellnote=round(colcase),'''
      '''notecol="black")''')
    R('''dev.off()''')

    E.info("plotting all individuals matrix")
    R('''rowall <- all.df[rownames(colcase),]''')
    R('''colall <- rowall[,colnames(colcase)]''')
    R('''png("%s/all-matrix.png", width=720, height=720)''' % os.getcwd())
    R('''heatmap.2(as.matrix(colall), trace="none", col=rep("#F0F8FF", 100),'''
      '''dendrogram="none", Colv=colnames(colall), key=FALSE, '''
      '''colsep=seq(1:length(colnames(colall))), '''
      '''rowsep=seq(1:length(rownames(colall))), '''
      '''Rowv=rownames(colall), margins=c(10,10), cellnote=round(colall),'''
      '''notecol="black")''')
    R('''dev.off()''')

    # plot MAF vs homozygosity
    maf_df = pd.read_table(mafs, sep="\t", header=0, index_col=0)
    plot_df = pd.DataFrame(columns=["MAF"],
                           index=maf_df.index)
    plot_df["MAF"] = maf_df["MAF"]

    homs = case_mat.diagonal()
    hom_series = pd.Series({x: y for x, y in zip(variant_ids,
                                                 homs)})
    plot_df["explained_by_homozygotes"] = hom_series
    plot_df["SNP"] = plot_df.index
    plot_df.index = [ix for ix, iy in enumerate(plot_df.index)]
    plotPenetrances(plotting_df=plot_df)

    out_df = summaryPenetrance(maf_df=maf_df,
                               case_counts=case_mat,
                               variants=variant_ids,
                               n_cases=ncases,
                               n_total=tcount)
    return out_df, pen_df


def summaryPenetrance(maf_df, case_counts,
                      variants, n_cases, n_total):
    '''
    Summarise genotype counts and proportion of cases explained
    by the observed homozygotes and compound heterozygotes.
    This is a function of the total population size and
    population allele frequency - does this assume 100%
    penetrance of each allele?
    '''

    # homozygous individuals  are on the
    # diagonal of the case_counts array
    homozyg_cases = case_counts.diagonal()
    homozyg_series = pd.Series({x: y for x, y in zip(variants,
                                                     homozyg_cases)})

    # heterozygotes are on the off-diagonal elements
    # get all off diagonal elements by setting diagonals to zero
    # matrix is diagonal symmetric
    np.fill_diagonal(case_counts, 0)

    het_counts = np.sum(case_counts, axis=0)
    het_series = pd.Series({x: y for x, y in zip(variants,
                                                 het_counts)})
    out_df = pd.DataFrame(columns=["homozygote_cases",
                                   "heterozygote_cases"],
                          index=maf_df.index)
    out_df["MAF"] = maf_df["MAF"]
    out_df["homozygote_cases"] = np.round(homozyg_series, 1)
    out_df["expected_cases"] = np.round(((out_df["MAF"] ** 2) * n_total), 3)
    out_df["heterozygote_cases"] = het_series
    out_df["hom_prop_explained"] = np.round(homozyg_series/float(n_cases), 3)
    out_df["het_prop_explained"] = np.round(het_series/float(n_cases), 3)

    return out_df


def plotPenetrances(plotting_df):
    '''
    Plot the proportion of cases/phenotype explained by
    individuals carrying allele vs. population allele frequency.

    Generate final output summary table (should be in separate function)
    '''

    # only need to plot variants with MAF >= 0.01
    low_frq = plotting_df["MAF"] < 0.01
    hi_df = plotting_df[~low_frq]

    # get into R and use ggplot for MAF vs homozygosity amongs cases
    r_plot = py2ri.py2ri_pandasdataframe(hi_df)
    R.assign("hom.df", r_plot)

    R('''suppressPackageStartupMessages(library(ggplot2))''')
    R('''png("%s/penetrance-plot.png", height=720, width=720)''' % os.getcwd())
    R('''pen_p <- ggplot(hom.df, aes(x=explained_by_homozygotes, y=MAF, colour=SNP)) + '''
      '''geom_point(size=4) + theme_bw() + '''
      '''geom_text(aes(label=explained_by_homozygotes),'''
      '''colour="black",vjust=0.5, hjust=0.5) + '''
      '''labs(x="Number of Red haired homozygotes", y="MAF") + '''
      '''theme(axis.title=element_text(size=10, colour="black"))''')
    R('''print(pen_p)''')
    R('''dev.off()''')


def findDuplicateVariants(bim_file, take_last=False):
    '''
    identify variants with duplicate position and reference
    alleles
    '''

    # count the number of lines first to get
    # the necessary array sizes
    E.info("getting number of variants")
    lines = 1
    with open(bim_file, "r") as bfile:
        for line in bfile.readlines():
            lines += 1

    E.info("%i variants found" % lines)

    # setup index arrays
    var_array = np.empty(lines, dtype=object)
    ref_alleles = np.empty(lines, dtype=object)
    pos_array = np.zeros(lines, dtype=np.int64)
    minor_alleles = np.empty(lines, dtype=object)
    idx = 0

    # find duplicates on position
    with open(bim_file, "r") as bfile:
        for line in bfile.readlines():
            line = line.rstrip("\n")
            varline = line.split("\t")
            var = varline[1]
            pos = int(varline[3])
            ref_allele = varline[-1]
            minor_allele = varline[-2]
            var_array[idx] = var
            ref_alleles[idx] = ref_allele
            minor_alleles[idx] = minor_allele
            pos_array[idx] = pos

            idx += 1

    # find duplicates using pandas series
    pos_series = pd.Series(pos_array)
    dup_last = pos_series[pos_series.duplicated(take_last=True)]
    dup_first = pos_series[pos_series.duplicated(take_last=False)]
    var_series = pd.Series(var_array)
    ref_series = pd.Series(ref_alleles)
    alt_series = pd.Series(minor_alleles)

    # a few variants have duplicate IDs - count these as duplicates
    # and add to the exclusion list - these won't be identified
    # based on shared position necessarily - force add them
    ref_first = ref_series[ref_series.duplicated(take_last=False)]
    ref_last = ref_series[ref_series.duplicated(take_last=True)]
    ref_dups = set(ref_first.index).union(ref_last.index)

    # union of take first and take last
    dup_all = set(dup_last.index).union(set(dup_first.index))
    dup_complete = dup_all.union(ref_dups)
    dup_idx = np.array([sx for sx in dup_complete])
    dup_idx.sort()

    # make a dataframe to hold all triallelic and duplicate variants
    dup_dict = {"SNP": var_series[dup_idx],
                "BP": pos_series[dup_idx],
                "REF": ref_series[dup_idx],
                "VAR": alt_series[dup_idx]}
    dup_df = pd.DataFrame(dup_dict)

    # some variants may have more than one ID/entry
    # step through using pandas groupby - group on position
    E.info("looking for duplicates and triallelic variants")
    tri_alleles = []
    dups_alleles = []
    overlap_vars = []
    for names, groups in dup_df.groupby(["BP"]):
        # if there is only one reference allele, indicates a
        # triallelic variant, otherwise its probably a duplicate
        # or overlaping INDEL and SNV
        var_lens = groups["VAR"].apply(len)
        if groups.shape[0] == 1:
            pass
        elif np.mean(var_lens) > 1:
            # probably overlapping variants, exclude, but report
            # separately
            over_vars = groups["SNP"].values.tolist()
            for ovs in over_vars:
                overlap_vars.append(ovs)

        elif len(set(groups["REF"])) == 1:
            tri_vars = groups["SNP"].values.tolist()
            for tri in tri_vars:
                tri_alleles.append(tri)
        else:
            dup_vars = groups["SNP"].values.tolist()
            for dup in dup_vars:
                dups_alleles.append(dup)

    E.info("%i triallelic variants found" % len(tri_alleles))
    E.info("%i duplicate position variants found" % len(dups_alleles))
    E.info("%i overlapping SNVs and INDELs found" % len(overlap_vars))

    return dups_alleles, tri_alleles, overlap_vars


def flagExcessHets(hets_file, plot=True, plot_path=None):
    '''
    Take output from Plink 1.9 --het command
    calculate heterozygosity rate and flag individuals
    with heterozygosity > 3 s.d. from the mean
    value.
    This assumes all individuals are from the same
    population, and thus form a homogenous cluster,
    with only outliers at the extremes.
    Visualise the data, if there are multiple apparent
    clusters then filter for ethnicity/ancestry first
    '''
    if hets_file.endswith("gz"):
        compression = "gzip"
    else:
        compression = None

    het_df = pd.read_table(hets_file, header=0, index_col=None,
                           sep="\t", compression=compression)
    nmiss = pd.Series(het_df.loc[:, "N(NM)"], dtype=np.float64)
    nhoms = het_df.loc[:, "O(HOM)"]
    het_df["het_rate"] = (nmiss - nhoms) / nmiss
    # get mean value and std, set upper and lower thresholds
    mean_het = np.mean(het_df.loc[:, "het_rate"].values)
    sd_het = np.std(het_df.loc[:, "het_rate"].values)

    upper = mean_het + (3 * sd_het)
    lower = mean_het - (3 * sd_het)

    hi_hets = het_df[het_df["het_rate"] > upper]
    lo_hets = het_df[het_df["het_rate"] < lower]

    E.info("%i individuals with high heterozygosity" % len(hi_hets))
    E.info("%i individuals with low heterozygosity" % len(lo_hets))

    hi_hets["exclude"] = "high_heterozygosity"
    lo_hets["exclude"] = "low_heterozygosity"
    all_flags = lo_hets.append(hi_hets)

    if plot:
        E.info("plotting heterozygosity rate distribution")
        py2ri.activate()
        r_df = py2ri.py2ri_pandasdataframe(het_df)
        R.assign("het.df", r_df)
        R('''suppressPackageStartupMessages(library(ggplot2))''')
        R('''p <- ggplot(het.df, aes(het_rate)) + '''
          '''geom_histogram() + '''
          '''labs(title="Distribution of heterozygosity rate") + '''
          '''theme_bw() + '''
          '''geom_vline(xintercept=c(%0.3f, %0.3f), '''
          '''linetype=2, col="#838B83")''' % (lower, upper))
        R('''png("%s/het_rate-hist.png")''' % plot_path)
        R('''print(p)''')
        R('''dev.off()''')

    return all_flags


def flagGender(gender_file, plot=True, plot_path=None):
    '''
    Parse the .sexcheck output report from Plink
    --sex-check and flag gender discordant individuals.

    Arguments
    ---------
    gender_file: string
      the .sexcheck output report file from Plink --sex-check

    plot: boolean
      generate a histogram of F values distributions showing male and
      female clusters, split by reported gender

    plot_path: string
      PATH to save F coefficient histogram

    Returns
    -------
    discords: pandas.Core.DataFrame
      a pandas dataframe of individuals that are gender discordant
    '''

    gender_df = pd.read_table(gender_file, header=0,
                              index_col=None, sep=None)
    genders = lambda x: "male" if x == 1 else "female"
    gender_df["GENDER"] = gender_df["PEDSEX"].apply(genders)

    E.info("checking individuals for discordance")
    discords = gender_df[gender_df["STATUS"] != "OK"]
    discords.drop(labels=["PEDSEX", "SNPSEX", "STATUS", "F",
                          "GENDER"],
                  axis=1, inplace=True)

    E.info("%i individuals with discordant gender" % len(discords))

    if plot:
        E.info("plotting F gender coefficient distributions")
        py2ri.activate()
        r_df = py2ri.py2ri_pandasdataframe(gender_df)
        R.assign("gender.df", r_df)
        R('''suppressPackageStartupMessages(library(ggplot2))''')
        R('''p <- ggplot(gender.df, aes(F, fill=GENDER)) + '''
          '''geom_histogram() + '''
          '''labs(title="F coefficient distributions for gender") + '''
          '''theme_bw() + facet_grid(. ~ GENDER)''')
        R('''png("%s/gender_check-hist.png")''' % plot_path)
        R('''print(p)''')
        R('''dev.off()''')

    else:
        pass

    return discords


def _compare_ibds(ibd_entry, threshold=0.03125):
    '''
    Just for internal use in `flagRelated` function.
    To compare IBD estimates and flag up related
    individuals

    Arguments
    ---------
    ibd_entry: pandas.Core.Series
      a single line entry from an IBD estimates
      file

    threshold: float
      the threshold at which to flag an individual as related

    Returns
    -------
    flag: boolean
      True if related, else false
    '''

    if ibd_entry["PI_HAT"] < threshold:
        return False
    else:
        return True


def flagRelated(ibd_file, chunk_size=None,
                threshold=0.03125, plot=True,
                plotting_path=None):
    '''
    Use IBS estimates to find pairs of related individuals
    above a threshold.

    This will also flag up the number of duplicated/monozygotic
    twin pairs (matrix diagonals).

    Arguments
    ---------
    ibd_file: string
      file containing IBS estimates between pairs from Plink
      or GCTA.

    chunk_size: int
      the file chunk size to read in at a time, should correspond
      to the number of individuals.  If not set, the whole file
      is read in.  Not recommend for large (>2GB) files.

    threshold: float
      IBS threshold, above which individuals will be flagged
      as related. Default is 3rd cousins.

    plot: boolean
      generate a histogram of the distribution of IBS values.
      Default = True

    plotting_path: string
      PATH to plot histogram to

    Returns
    -------
    flagged: pandas.Core.DataFrame
      dataframe of individuals to remove, with the estimated
      relationship to another individual.
    '''

    # need to make this faster
    # sequentially add new IDs only
    related_list = []
    ibds = []
    if ibd_file.endswith("gz"):
        comp = "gzip"
    else:
        pass

    E.info("reading file in chunks of %i lines" % chunk_size)
    if chunk_size:
        # read in and operate on chunks
        df_iter = pd.read_table(ibd_file, header=0, index_col=None,
                                delim_whitespace=True, compression=comp,
                                chunksize=chunk_size)
        count = 0
        for chunk in df_iter:
            count += 1
            entrys = chunk[["FID1", "IID1",
                            "FID2", "IID2",
                            "PI_HAT"]]
            ibds.append(entrys)
            relate_mask = entrys.apply(_compare_ibds, axis=1)
            related = entrys[relate_mask]
            E.info("%i relations found" % len(related))
            related_list.append(related)
    else:
        pass

    df = pd.concat(ibds, axis=0, keys=None)

    if plot:
        # for lots of observations, plot log counts
        E.info("plotting pair-wise IBD distribution")
        py2ri.activate()
        r_df = py2ri.py2ri_pandasdataframe(df)
        R.assign("relate.df", r_df)
        R('''suppressPackageStartupMessages(library(ggplot2))''')
        R('''p <- ggplot(relate.df, aes(PI_HAT+0.5)) + '''
          '''geom_histogram(binwidth=0.01) + '''
          '''labs(title="Proportion of IBD shared distribution") +  '''
          '''theme_bw() + scale_y_log10() + '''
          '''geom_vline(xintercept=%(threshold)f, '''
          '''linetype=4, colour="#838B83")''' % locals())
        R('''png("%s/IBD-hist.png")''' % plotting_path)
        R('''print(p)''')
        R('''dev.off()''')
    else:
        pass

    return related_list


def flagInbred(inbred_file, inbreeding_coefficient,
               ibc_threshold=0.05,
               plot=True, plot_path=None):
    '''
    Use Plink or GCTA's estimate of F, inbreeding coefficient
    to flag individuals that are highly inbred.

    Arguments
    ---------
    inbred_file: string
      file containing estimates of F

    inbreeding_coefficient: string
      coefficient to use to identify inbred individuals.  This name
      should correspond to one of the columns in `inbred_file`.

    ibc_threshold: float
      the threshold above which individuals will be flagged as inbred

    plot: boolean
      generate a histogram of the distribution of F coefficients

    plotting_path: string
      PATH to directoru for plotting F coefficient distribution

    Returns
    -------
    inbreds: padas.Core.DataFrame
      dataframe of inbred individuals to exclude from analysis
    '''

    inbreed_df = pd.read_table(inbred_file, header=0,
                               index_col=None, sep="\t")

    E.info("Identifing individuals with inbreeding coefficient"
           " greater than %0.3f" % ibc_threshold)

    inbreds = inbreed_df[inbreed_df[inbreeding_coefficient] > ibc_threshold]
    inbreds = inbreds[["FID", "IID"]]

    E.info("%i individuals with high inbreeding "
           "coefficient" % len(inbreds))

    if plot:
        E.info("plotting F coefficient distributions")
        py2ri.activate()
        r_df = py2ri.py2ri_pandasdataframe(inbreed_df)
        R.assign("inbreed.df", r_df)
        R('''suppressPackageStartupMessages(library(ggplot2))''')
        R('''p <- ggplot(inbreed.df, aes(%(inbreeding_coefficient)s)) + '''
          '''geom_histogram(binwidth=0.01) + '''
          '''labs(title="Inbreeding coefficient, %(inbreeding_coefficient)s,'''
          '''distribution") + theme_bw() + '''
          '''geom_vline(xintercept=%(ibc_threshold)0.3f, '''
          '''linetype=4, colour="#838B83")''' % locals())
        R('''png("%s/inbreeding-hist.png")''' % plot_path)
        R('''print(p)''')
        R('''dev.off()''')
    else:
        pass

    return inbreds


def mergeQcExclusions(hets_file=None, inbred_file=None,
                      related_file=None, gender_file=None,
                      mask_file=None):
    '''
    Merge sets of excluded individuals into a single file for
    downstream analysis, processing, etc

    Arguments
    ---------
    hets_file: string
      file containing individuals to remove due to excessive or
      reduced heterozygosity

    inbred_file: string
      file of individuals highly related to themselves for
      exclusion

    related_file: string
      file of IDs of individuals pruned due to greater relatedness
      than an arbitrary threshold

    gender_file: string
      individuals with discordant reported vs. genetic gender

    mask_file: string
      individuals to be excluded from analyses, unrelated
      for reasons to QC (i.e. mask out category of individuals)

    Returns
    -------
    exclusions: pandas.Core.DataFrame
      A dataframe of FID and IIDs of the unique set of excluded
      individuals
    '''

    if hets_file:
        hets_df = pd.read_table(hets_file, sep="\t",
                                header=0, index_col=None)
        E.info("%i exclusions due to "
               "heterozygosity deviation" % len(hets_df))
    else:
        hets_df = None
        E.warn("No heterozygosity exclusion file")
    if inbred_file:
        inbred_df = pd.read_table(inbred_file, sep="\t",
                                  header=0, index_col=None)
        E.info("%i exclusions due "
               "to consanguinuity" % len(inbred_df))
    else:
        inbred_df = None
        E.warn("No inbred exclusions")

    if related_file:
        related_df = pd.read_table(related_file, delim_whitespace=True,
                                   header=None, index_col=None)
        related_df.columns = ["FID", "IID"]
        E.info("%i individuals excluded due "
               "to high relatedness" % len(related_df))
    else:
        related_df = None
        E.warn("No individuals excluded on relatedness")

    if gender_file:
        gender_df = pd.read_table(gender_file, sep="\t",
                                  header=0, index_col=None)
        E.info("%i individuals with discordant "
               "gender recorded" % len(gender_df))
    else:
        gender_df = None
        E.warn("No individuals exclued with "
               "discordant gender")

    if mask_file:
        mask_df = pd.read_table(mask_file, sep="\t",
                                header=None, index_col=None)
        E.info("%i individuals to be excluded "
               "for additional reasons" % len(gender_df))
        mask_df.columns = ["FID", "IID"]
    else:
        mask_df = None

    df_list = [hets_df, inbred_df, related_df, gender_df,
               mask_df]

    df_true = [True for x in df_list if x is not False]

    if not all(df_true):
        raise ValueError("no QC files detected - do some QC!!")
    else:
        pass

    # assume all df have FID and IID columns
    real_df = [x for x in df_list if x is not None]
    real_df = [x[["FID", "IID"]] for x in real_df]

    full_df = pd.concat(real_df, keys=None, axis=0)
    exclusions = full_df.drop_duplicates(subset=["FID",
                                                 "IID"],
                                         take_last=True,
                                         inplace=False)

    return exclusions


def selectLdFromTabix(ld_dir, chromosome, snp_pos,
                      ld_threshold=0.01):
    '''
    Select all LD values from a tabix indexed BGZIP
    file of LD.  Assumes Plink format.

    Arguments
    ---------
    ld_dir: string
      path to directory containing LD data

    chromosome: string
      chromosome of SNP to pull out LD values
      assumes chrN format

    snp_pos: int
      bp mapping position of the SNP on the same
      genome build as the LD was calculated

    ld_threshold: float
      minimum LD value to return

    Returns
    -------
    ld_df: pandas.Core.DataFrame
      Pandas dataframe containing LD values over
      target range.
    '''

    tab_dir = [td for td in os.listdir(ld_dir) if re.search(".bgz$", td)]
    contig = int(chromosome.lstrip("chr"))
    start = snp_pos
    end = snp_pos

    tab_query = """
    tabix %(ld_dir)s/%(tab_indx)s %(contig)i:%(start)i-%(end)i |
    awk '{if($7 >= %(ld_threshold)s) print $0}'"""

    tab_indx = [tx for tx in tab_dir if re.search(chromosome,
                                                  tx)][-1]

    E.info("Retrieving LD values at bp: %i" % snp_pos)
    proc = subprocess.Popen(tab_query % locals(),
                            shell=True,
                            stdout=subprocess.PIPE)

    ld_dict = {}
    count = 0
    for line in proc.stdout:
        snp_dict = {}
        parse = line.split("\t")
        snp_dict["CHR_A"] = int(parse[0])
        snp_dict["BP_A"] = int(parse[1])
        snp_dict["SNP_A"] = parse[2]
        snp_dict["CHR_B"] = int(parse[3])
        snp_dict["BP_B"] = int(parse[4])
        snp_dict["SNP_B"] = parse[5]
        snp_dict["R2"] = float(parse[6])
        snp_dict["DP"] = float(parse[7])
        count += 1
        ld_dict[count] = snp_dict

    ld_df = pd.DataFrame(ld_dict).T
    # ld Dataframe may be empty, return
    # empty dataframe
    try:
        ld_df.index = ld_df["SNP_B"]
        ld_df.drop_duplicates(subset="SNP_B",
                              keep="last",
                              inplace=True)
    except KeyError:
        E.info("No SNPs detected in LD "
               "with r^2 > {}".format(ld_threshold))

        ld_df = pd.DataFrame(0.0,
                             index=[snp_pos],
                             columns=["SNP_A",
                                      "DP",
                                      "R2"])
    return ld_df


def selectLdFromDB(database, table_name,
                   index_snp,
                   index_label=None,
                   ld_threshold=None):
    '''
    Select LD values from an SQL table over
    a specific range.  Large regions will consume
    large memory and queries may take several
    minutes to complete.

    Arguments
    ---------
    database: sql.connection
      An SQL database connection to the DB
      containing the LD values

    table_name: string
      The table to query containing LD information

    index_snp: string
      SNP ID to select LD values from the SQL
      database on

    index_label: str
      Column label in SQL database to use as the
      index in the output dataframe

    ld_threshold: float
      minimum LD value to return

    Returns
    -------
    ld_df: pandas.Core.DataFrame
      Pandas dataframe containing LD values over
      target range.
    '''

    # UTF-8 codec struggles to decode ';' in some columns
    database.text_factory = str

    if ld_threshold:
        state = '''
        select SNP_A,SNP_B,R2 FROM %s where %s = "%s" AND
        R2 > %0.3f;
        ''' % (table_name, index_label,
               index_snp, ld_threshold)
    else:
        state = '''
        select SNP_A,SNP_B,R2 FROM %s where %s = "%s";
        ''' % (table_name, index_label, index_snp)

    ld_df = pdsql.read_sql(sql=state, con=database,
                           index_col=index_label)

    return ld_df


def calcLdScores(ld_table, snps,
                 scale=False, metric="R2",
                 snp_list=None):
    '''
    Calculate the LD scores for SNPs across a chromosome,
    stored in a SQL database.

    Arguments
    ---------
    ld_table: pandas.Core.DataFrame
      Pandas dataframe in table format containing LD
      values between SNPs.  Columns are `SNP_A`, `SNP_B`
      and `R2`.

    snps: list
      the snps over which to calculate LD scores

    scale: bool
      Whether to scale LD score by the number of SNPs
      used to calculate the score.  Useful if used
      as a weighting for other SNP scores.

    metric: string
      Use either R^2 or D' as the LD metric

    snp_list: list
      A list of SNP IDs to restrict the
      LD score calculation to
    Returns
    -------
    ld_scores: float
      LD scores for each SNP
    '''

    if len(ld_table) > 0:
        if snp_list:
            try:
                # use np.sum to handle NaN values
                ld_table = ld_table.loc[snp_list]
                if metric == "R2":
                    ld_score = np.sum(ld_table["R2"])
                elif metric == "DP":
                    ld_score = np.sum(ld_table["DP"])
            except KeyError:
                E.warn("None of the SNPs are in LD")
                ld_score = 0
        else:
            if metric == "R2":
                ld_score = np.sum(ld_table["R2"])
            elif metric == "DP":
                ld_score = np.sum(ld_table["DP"])
    else:
        ld_score = 0

    if scale:
        ld_scores = ld_score/len(ld_table)
    else:
        ld_scores = ld_score

    return ld_scores


def calcWeightedEffects(gwas_results, snps, calc_se=True,
                        scale=False):
    '''
    Calculate the standard error weighted effect sizes score
    for each SNP:
    score = sum(ln(OR) * se)

    Arguments
    ---------
    gwas_results: pandas.Core.DataFrame
      A dataframe of the results from a genome_wide association
      study. Assumes SNP IDs are the index column.

    snps: list
      the snps over which to calculate the total weighted
      effect size score.

    calc_se: boolean
      Calculate the standard error from the p-values and
      effect sizes:

      SE = ln(OR)/Z
      Z = -0.862 + sqrt(0.743 - 2.404 * ln(P))

    scale: boolean
      Scale the sum of standard error weighted effect sizes
      by the number of SNPs

    Returns
    -------
    es_score: float
      sum of SE weighted effect sizes
    '''

    # calculate standard error of effect size based on
    # p-value and effect size

    if calc_se:
        # check p-values that = 0 are set to smallest floating
        # point representation instead
        gwas_results["P"][gwas_results["P"] == 0] = np.finfo(np.float64).min
        z_func = lambda x: - 0.862 + sqrt(0.743 - 2.404 * np.log(x))
        gwas_results["Z"] = gwas_results["P"].apply(z_func)
        gwas_results["SE"] = abs(np.log(gwas_results["OR"])/gwas_results["Z"])
    else:
        E.warn("Standard errors have not been calculated, please "
               "make sure they exist in this results table")

    es_score = sum((abs(np.log(gwas_results["OR"])) * gwas_results["SE"]).fillna(0))

    if scale and len(gwas_results):
        return es_score/len(gwas_results)
    else:
        return es_score


def snpPriorityScore(gwas_results, chromosome, ld_dir=None,
                     clean=True, database=None, table_name=None):
    '''
    Generate SNP scores based on the amount of genetic variation
    they capture and the sum of the weighted effect sizes for
    the trait of interest.

    This score can then be integrated with a score based on
    the overlap with functional annotation features
    of interest.

    Arguments
    ---------
    gwas_results: string
      Results from a GWAS, assumed to be in Plink format.

    ld_dir: string
      directory containing tabix index LD files from Plink

    database: string
      Path to an SQL database containing LD values in
      table format

    table_name: string
      Specific table, often referring to a specific
      chromosome, that contains LD values with columns
      SNP_A, SNP_B, BP_A, BP_B and R2.

    chromosome: string
      A chromosome to select from the gwas_results
      file.

    clean: boolean
      Whether the results table has been pre-cleaned to
      remove results not relevant to SNPs. e.g. if
      covariates had been included in the regression
      model these should be removed.

    Returns
    -------
    SNP_scores: pd.Core.DataFrame
      A pandas dataframe of LDscores, weight effect size
      scores and SNP priority score.
    '''

    E.info("Reading association results from %s" % gwas_results)
    gwas_df = pd.read_table(gwas_results, index_col=None,
                            sep="\t", header=0)

    if clean:
        gwas_df = pd.read_table(gwas_results, index_col=None,
                                sep="\t", header=0)

    else:
        gwas_df = pd.read_table(gwas_results, index_col=None,
                                sep="\s*", header=0)
        gwas_df = gwas_df[gwas_df["TEST"] == "ADD"]

    gwas_df.index = gwas_df["SNP"]

    # in order to reduce the computational load it is
    # necessary to break up the SNPs into regions.
    # The logical way would be to use recombination
    # hotspots, however, this will still leave
    # some very large windows
    # Use a moving window over the chromosome of
    # ~250Kb, with 25kb overlap.

    chr_df = gwas_df[gwas_df["CHR"] == int(chromosome)]
    # duplicates cause selection of individual SNPs
    # to break - why are there still duplicates??
    chr_df.drop_duplicates(subset="BP", keep="last",
                           inplace=True)

    priority_list = []
    ld_scores = {}
    es_scores = {}
    priority_scores = {}
    snp_set = chr_df.index

    if database:
        dbh = sql.connect(database)
    else:
        pass

    # iterate over SNPs
    for snp in snp_set:
        if database:
            ld_values = selectLdFromDB(dbh,
                                       table_name=table_name,
                                       index_snp=snp,
                                       index_label="SNP_B")
        elif ld_dir:
            snp_pos = int(chr_df.loc[snp, "BP"])
            ld_values = selectLdFromTabix(ld_dir=ld_dir,
                                          chromosome=chromosome,
                                          snp_pos=snp_pos)

        ldsnps = ld_values.loc[:, "SNP_A"].values
        ldsnps = {sx for sx in ldsnps}

        ldscore = calcLdScores(ld_table=ld_values,
                               snps=ldsnps,
                               scale=False)
        ld_scores[snp] = ldscore

        try:
            gwas_results = chr_df.loc[ldsnps]
            escore = calcWeightedEffects(gwas_results=gwas_results,
                                         snps=ldsnps,
                                         calc_se=True,
                                         scale=True)
        except KeyError:
            gwas_results = chr_df.loc[snp]
            if gwas_results["P"] == 0:
                gwas_results["P"] = np.finfo(np.float64).min
            else:
                pass
            z_func = lambda x: - 0.862 + sqrt(0.743 - 2.404 * np.log(x))
            gwas_results["Z"] = z_func(gwas_results["P"])
            gwas_results["SE"] = abs(np.log(gwas_results["OR"])/gwas_results["Z"])
            escore = gwas_results["SE"] * abs(np.log(gwas_results["OR"]))

        es_scores[snp] = escore

        weight = escore * ldscore
        priority_scores[snp] = weight

    SNP_scores = pd.DataFrame([pd.Series(ld_scores),
                               pd.Series(es_scores),
                               pd.Series(priority_scores)]).T
    SNP_scores.columns = ["LDScore", "WeightEffectSize", "PriorityScore"]
    SNP_scores.sort_values(by="PriorityScore", inplace=True)

    return SNP_scores


def fitPrior(value, distribution, dist_params):
    '''
    Fit a prior probability given a value,
    distribution and distribution parameters.

    You are responsible for defining the appropriate
    distribution and parameters

    Arguments
    ---------
    Value: float
      A value to calculate a prior probability from

    distribution: string
      A distribution from which to calculate a probability.
      Current values are "normal", "t", "gamma",
      "lognormal", "exponential".

    dist_params: tuple
      parameters to define the distribution,
      * normal: (mean, std)
      * t: (df, ncp)
      * gamma: (k, theta)
      * lognormal: (ln(mean), std)
      * exponential: (lambda)

    Returns
    -------
    prior: float
      Prior probability attached to input value
    '''

    # distribution parameters should be passed
    # explicitly
    if distribution == "normal":
        prior = stats.norm(*dist_params).pdf(value)
    elif distribution == "t":
        prior = stats.t(*dist_params).pdf(value)
    elif distribution == "gamma":
        prior = stats.gamma(*dist_params).pdf(value)
    elif distribution == "lognormal":
        prior = stats.lognorm(*dist_params).pdf(value)
    elif distribution == "exponential":
        prior = stats.expon(*dist_params).pdf(value)
    else:
        raise ValueError("Distrubtion %s not "
                         "implemented" % distribution)

    return prior


def calcPriorsOnSnps(snp_list, distribution, params=None):
    '''
    Calculate prior probabilities on SNPs based
    on a predefined value, a distribution and
    parameters to describe the distribution.
    This relies inherently on the correct and
    appropriate distribution to be defined, i.e.
    that it is conjugate to the marginal likelihood
    distribution.
    TO DO: introduce robust Bayesian modelling

    Arguments
    ---------
    snp_list: dict
      SNPs with score/value attached to determine
      the prior probability

    distribution: string
      the distribution from which to draw probabilities

    params: tuple
      parameters to describe the appropriate
      distribution.

    Returns
    -------
    prior_probs: dict
      dictionary of priors for SNPs
    '''

    prior_probs = {}
    # if there is no score for that SNP then use an
    # uninformative or Jeffrey's prior
    for snp in snp_list.keys():
        if snp_list[snp] != 0:
            prior_probs[snp] = fitPrior(value=snp_list[snp],
                                        distribution=distribution,
                                        dist_params=params)
        else:
            prior_probs[snp] = 0.5

    return prior_probs


def estimateDistributionParameters(data,
                                   distribution,
                                   fscale=None,
                                   floc=None,
                                   **kwargs):
    '''
    Use maximum likelihood to estimate the parameters
    of the defined distribution.

    Arguments
    ---------
    data: pd.Series/np.array
      data used to estimate parameters from

    distribution: string
      distribution assumed to underlie the data
      generation process

    fscale: float
      scale parameter of the distribution to fix

    floc: float
      location parameter of the distribution to
      fix

    **kwargs: float
      additional kwargs to pass as fixed parameters

    Returns
    -------
    est_params: tuple
      estimated distribution parameters
    '''

    if distribution == "normal":
        mu, sigma = stats.norm.fit(data)
        est_params = (mu, sigma,)
    elif distribution == "t":
        df, mu, sigma = stats.t.fit(data)
        est_params = (df, mu, sigma,)
    elif distribution == "gamma":
        k, theta, mu = stats.gamma.fit(data)
        est_params = (k, theta, mu,)
    elif distribution == "lognormal":
        exp_mu, sigma, theta = stats.lognorm.fit(data)
        est_params = (exp_mu, sigma, theta,)
    elif distribution == "exponential":
        beta, lambda_x = stats.expon.fit(data)
        est_params = (beta, lambda_x,)
    else:
        raise ValueError("Distrubtion %s not "
                         "implemented" % distribution)

    return est_params


def calculatePicsValues(snp_id, index_log10p, ld_values,
                        priors=None, k=2, analysis_snps=None):
    '''
    Use the PICS method to assign probability to SNPs as
    being causal for association signals at a locus,
    given the strength of their association (log10 P-value),
    and linkage disequilbrium with the lead SNP (smallest
    p-value at the locus).

    This method allows the prioritisation of SNPs, including
    those where there are multiple independent signals.  It
    requires that these independent signals are however
    input as separate SNPs.

    NB::
      Fahr et al initially use k=6.4 based on their observation
      that altering k=[6,8] does not have an appreciable impact
      on the PICS values.  However, when comparing this
      implementation to the PICS webserver, k=2 gives more
      similar values based on the current Phase3 1000 Genomes
      LD values and SNPs.

    Arguments
    ---------
    snp_id: string
      rs ID of the lead SNP from the associated region/
      independent signal of association

    index_log10p: float
      the negative log10(p-value) of association with
      the phenotype of interest

    ld_values: pd.Core.DataFrame
      A pandas dataframe of LD values between the index SNP
      and all other SNPs given an arbitrary threshold.  The
      expected columns are index SNP, SNP of interest, r^2.

    priors: dict
      the prior value to attach to each SNP.  Can be used to
      integrate functional information into the PICS
      calculations.  EXPERIMENTAL

    k: float
      The power to raise the correlation of alleles to.  When
      k=2, this scales the standard deviation of the sample
      distribution for the marignal likelihood by the residual
      LD.  Increasing k downweights the LD difference between
      the index SNP and SNP of interest.

    analysis_snps: list
      The LD query may return SNPs that do not appear in the
      analysis.  This may cause problems downstream, therefore
      output and probability calculation should be limited
      to these variants.

    Returns
    -------
    PICS: pandas.Core.Series
      A pandas series of SNPs and calculated PICS scores.
    '''

    # assume the SNPs of interest are all contained in the
    # ld_values table index

    top_p = stats.norm(index_log10p,
                       sqrt(index_log10p)/2).cdf(index_log10p)

    prob_dict = {}
    prob_dict[snp_id] = top_p

    E.info("calculating scores for %i SNPs" % len(ld_values))
    # If a SNP is in perfect LD with the index SNP this forces
    # the standard deviation to be 0, add a small correction
    # to allow the calculation of marginal likelihood value
    # e.g. 0.0001

    for snp in ld_values.index:
        if snp in analysis_snps:
            try:
                r2 = ld_values.loc[snp]["R2"]
                r = sqrt(r2)
                mu = r2 * index_log10p
                sigma = sqrt(1 - (r ** k)) * (sqrt(index_log10p)/2)
                if sigma == 0:
                    sigma = 0.0001
                else:
                    pass
                # use log likelihoods, these are more numerically
                # stable and avoid the multiplication of very small
                # numbers
                # if priors are not set, force uninformative prior
                # i.e. if not conjugate with likelihood
                likelihood = np.log(stats.norm(mu, sigma).pdf(index_log10p))
                try:
                    prior = np.log(priors[snp])
                except:
                    prior = np.log(1.0)
                    prob_dict[snp] = np.exp(likelihood + prior)
            except KeyError:
                E.warn("SNP %s not found in LD with %s" % (snp,
                                                           snp_id))
        else:
            E.warn("SNP {} is not found in the analysis, "
                   "it will not be included in the output".format(snp))

    # calculate normalized probabilities, where sum of all probs=1
    # use numpy sum to handle NaN values
    sum_probs = np.sum(prob_dict.values())
    pics_dict = {}

    for snp_p in prob_dict.keys():
        pics_dict[snp_p] = prob_dict[snp_p]/sum_probs

    pics_series = pd.Series(pics_dict)
    PICS = pics_series.sort_values(ascending=False)

    return PICS


def getLdValues(database, table_name, index_snp, ld_threshold=0.5):
    '''
    Get all LD values for the index SNP above a given r^2
    threshold

    Arguments
    ---------
    database: sql.connection
      An SQL database connection to the DB
      containing the LD values

    table_name: string
      The table to query containing LD information

    index_snp: string
      SNP ID to select LD values from the SQL
      database on

    ld_threshold: float
      a threshold above which to select LD values
      with the index SNP

    Returns
    -------
    ld_df: pandas.Core.DataFrame
      Pandas dataframe containing LD values over
      target range.
    '''

    E.info("executing SQL query on table: %s" % table_name)
    ld_a = selectLdFromDB(database=database,
                          table_name=table_name,
                          index_snp=index_snp,
                          index_label="SNP_B",
                          ld_threshold=ld_threshold)
    ld_a.columns = ["SNP", "R2"]

    ld_b = selectLdFromDB(database=database,
                          table_name=table_name,
                          index_snp=index_snp,
                          index_label="SNP_A",
                          ld_threshold=ld_threshold)
    ld_b.columns = ["SNP", "R2"]

    ld_df = ld_a.append(ld_b)
    ld_df.index = ld_df["SNP"]
    # drop duplicate indices
    ld_df.drop_duplicates(subset="SNP",
                          keep="last",
                          inplace=True)

    E.info("%i records found matching query" % len(ld_df))

    return ld_df


def PICSscore(gwas_results, chromosome, database=None,
              table_name=None, priors=None, clean=True,
              ld_threshold=0.5, ld_dir=None):
    '''
    Prioritise SNPs based on the conditional probability
    of being the causal SNP at an associated region given
    the strength of association and LD with surrounding SNPs.

    Originally described in::
      Fahr et al Nature 518 (2015) pp337

    The current implementation does not allow the integration
    of a prior probability - this will come in the future.

    Arguments
    ---------
    gwas_results: string
      Results from a GWAS, assumed to be in Plink format.

    ld_dir: string
      directory containing tabix index LD files from Plink

    database: string
      Path to an SQL database containing LD values in
      table format

    table_name: string
      Specific table, often referring to a specific
      chromosome, that contains LD values with columns
      SNP_A, SNP_B, BP_A, BP_B and R2.

    chromosome: string
      A chromosome to select from the gwas_results
      file.

    priors: dict
      the prior value to attach to each SNP.  Can be used to
      integrate functional information into the PICS
      calculations.  EXPERIMENTAL

    clean: boolean
      Whether the results table has been pre-cleaned to
      remove results not relevant to SNPs. e.g. if
      covariates had been included in the regression
      model these should be removed.

    ld_threshold: float
      Threshold above which to select SNPs in LD
      with the lead SNP

    Returns
    -------
    PICS_scores: pd.Core.DataFrame
      A pandas dataframe of PICS scores for SNPs.
    '''

    E.info("Reading association results from %s" % gwas_results)
    if clean:
        gwas_df = pd.read_table(gwas_results, index_col=None,
                                sep="\t", header=0)

    else:
        gwas_df = pd.read_table(gwas_results, index_col=None,
                                sep="\s*", header=0)
        gwas_df = gwas_df[gwas_df["TEST"] == "ADD"]

    gwas_df.index = gwas_df["SNP"]

    E.info("subsetting data on chromosome %s" % chromosome)
    chr_df = gwas_df[gwas_df["CHR"] == int(chromosome)]
    try:
        chr_df.loc[:, "STAT"] = abs(chr_df["STAT"])
        chr_df.sort_values(by="STAT", inplace=True, ascending=False)
    except KeyError:
        chr_df.sort_values(by="CHISQ", inplace=True, ascending=False)

    chr_df.loc[:, "P"][chr_df["P"] == 0] = 1.79769e-308
    chr_df["P"].fillna(1.0)
    chr_df.loc[:, "log10P"] = np.log10(chr_df["P"])

    index_snp = chr_df.iloc[0]["SNP"]

    try:
        indexp = -chr_df.iloc[0]["log10P"]
    except KeyError:
        indexp = -np.log10(chr_df.iloc[0]["P"])

    E.info("index SNP is %s with -log10(p)= %0.3f" % (index_snp,
                                                      indexp))

    if database:
        dbh = sql.connect(database)
        ld_values = getLdValues(database=dbh,
                                table_name=table_name,
                                index_snp=index_snp,
                                ld_threshold=ld_threshold)
    elif ld_dir:
        snp_pos = int(chr_df.loc[index_snp]["BP"])
        ld_values = selectLdFromTabix(ld_dir=ld_dir,
                                      chromosome=chromosome,
                                      snp_pos=snp_pos)

    PICS_scores = calculatePicsValues(snp_id=index_snp,
                                      index_log10p=indexp,
                                      ld_values=ld_values,
                                      priors=priors,
                                      k=2,
                                      analysis_snps=gwas_df.index)
    return PICS_scores


def LdRank(gwas_results, chromosome,
           ld_dir=None, database=None,
           table_name=None, ld_threshold=0.8,
           top_snps=0.01, clean=True):
    '''
    Rank SNPs based on the LD with the lead SNP
    from the association region.  Take the top
    N% SNPs as the SNP set.

    Arguments
    ---------
    gwas_results: string
      Results from a GWAS, assumed to be in Plink format.

    ld_dir: string
      directory containing tabix index LD files from Plink

    database: string
      Path to an SQL database containing LD values in
      table format

    table_name: string
      Specific table, often referring to a specific
      chromosome, that contains LD values with columns
      SNP_A, SNP_B, BP_A, BP_B and R2.

    chromosome: string
      A chromosome to select from the gwas_results
      file.

    ld_threshold: float
      Threshold above which to select SNPs in LD
      with the lead SNP

    top_snps: float
      % SNPs to select, ranked on LD with the lead
      SNP

    Returns
    -------
    '''

    E.info("Reading association results from %s" % gwas_results)
    gwas_df = pd.read_table(gwas_results, index_col=None,
                            sep="\t", header=0)
    if clean:
        gwas_df = pd.read_table(gwas_results, index_col=None,
                                sep="\t", header=0)

    else:
        gwas_df = pd.read_table(gwas_results, index_col=None,
                                sep="\s*", header=0)
        gwas_df = gwas_df[gwas_df["TEST"] == "ADD"]

    gwas_df.index = gwas_df["SNP"]

    E.info("subsetting data on chromosome %s" % chromosome)
    chr_df = gwas_df[gwas_df["CHR"] == int(chromosome)]
    try:
        chr_df.loc[:, "STAT"] = abs(chr_df["STAT"])
        chr_df.sort_values(by="STAT", inplace=True, ascending=False)
    except KeyError:
        chr_df.sort_values(by="CHISQ", inplace=True, ascending=False)

    chr_df.loc[:, "P"][chr_df["P"] == 0] = 1.79769e-308
    chr_df["P"].fillna(1.0)
    chr_df.loc[:, "log10P"] = np.log10(chr_df["P"])

    index_snp = chr_df.iloc[0]["SNP"]

    if database:
        dbh = sql.connect(database)
        ld_values = getLdValues(database=dbh,
                                table_name=table_name,
                                index_snp=index_snp,
                                ld_threshold=ld_threshold)
    elif ld_dir:
        snp_pos = int(chr_df.loc[index_snp]["BP"])
        ld_values = selectLdFromTabix(ld_dir=ld_dir,
                                      chromosome=chromosome,
                                      snp_pos=snp_pos)

    # rank on LD with index SNP
    E.info("sort and rank top %0.3f SNPs in "
           "r2 > %0.3f with SNP %s" % (top_snps,
                                       ld_threshold,
                                       index_snp))
    index_series = pd.DataFrame(
        {"SNP": index_snp,
         "R2": 1.00},
        index=[index_snp])

    if len(ld_values):
        ld_values = ld_values.append(index_series)
    else:
        ld_values = index_series
        ld_values.columns = ["SNP", "R2"]

    ld_values.sort_values(by="R2", inplace=True,
                          ascending=False)
    size = len(ld_values)
    # use the math module ceil function to get
    # smallest integer greater than or equal to
    # the top %n SNPs

    top = int(ceil(size * top_snps))

    top_ld = ld_values.iloc[0:top, ]

    return top_ld


def calcApproxBayesFactor(log_or, standard_error,
                          prior_variance):
    '''
    Calculate the approximate Bayes Factor (ABF) from Wakefield
    Am. J. Hum. Genet.(2015) for a SNP.  The ABF is calculated
    from the effect size (log OR), variance (Standard error ^2)
    and a prior weight on the variance (W).

    Arguments
    ---------
    log_or: float
      The natural logarithm of the odds ratio or the effect
      size estimate on the observed scale.

    standard_error: float
      The standard error estimate on the effect size from
      the appropriate regression model

    prior_variance: float
      A prior variance weighting to apply to the variance for
      calculating the ABF.

    Returns
    -------
    ABF: float
      The calculated Approximate Bayes Factor
    '''

    # the variance on the MLE log OR is the squared standard error
    variance = standard_error ** 2
    _top = sqrt((prior_variance + variance)/variance)
    _exp_left = -((log_or ** 2)/variance)/2.0
    _exp_right = prior_variance/(prior_variance + variance)

    ABF = _top * exp(_exp_left * _exp_right)

    return ABF


def ABFScore(gwas_results, region_size, chromosome,
             prior=None,
             prior_variance=0.04, clean=True):
    '''
    Using approximate Bayes factors calculate the posterior
    association signal for each variant.  Credible intervals
    will be constructed later.

    Arguments
    ---------
    gwas_results: string
      Results from a GWAS, assumed to be in Plink format.

    region_size: int
      The region (in bp) by which to extend around the
      association signal index SNP - taken as the
      fine-mapping region. Region is index bp +/-
      region_size/2

    chromosome: string
      A chromosome to select from the gwas_results
      file.

    prior: float
      Prior probability NOT YET IMPLEMENTED

    prior_variance: float
      The variance prior that weights the standard error

    clean: boolean
      Whether the results table has been pre-cleaned to
      remove results not relevant to SNPs. e.g. if
      covariates had been included in the regression
      model these should be removed.


    Returns
    -------
    out_df: pandas.Core.DataFrame
      All input SNPs in the fine-mapping interval with their
      approximate Bayes Factors and posterior probabilities
    '''

    E.info("Reading association results from %s" % gwas_results)
    try:
        gwas_df = pd.read_table(gwas_results, index_col=None,
                                sep="\s*", header=0)
    except StopIteration:
        gwas_df = pd.read_table(gwas_results, index_col=None,
                                sep="\t", header=0)

    if clean:
        pass
    else:
        gwas_df = gwas_df[gwas_df["TEST"] == "ADD"]

    gwas_df.index = gwas_df["SNP"]

    E.info("subsetting data on chromosome %s" % chromosome)
    chr_df = gwas_df[gwas_df["CHR"] == int(chromosome)]
    try:
        try:
            chr_df.loc[:, "STAT"] = abs(chr_df["STAT"])
            chr_df.sort_values(by="STAT", inplace=True, ascending=False)
        except KeyError:
            chr_df.loc[:, "T"] = abs(chr_df["T"])
            chr_df.sort_values(by="T", inplace=True, ascending=False)
    except KeyError:
        chr_df.sort_values(by="CHISQ", inplace=True, ascending=False)

    # set p = 0 to minimum float value, ~1.79x10-308
    chr_df.loc[:, "P"][chr_df["P"] == 0] = 1.79769e-308
    chr_df["P"].fillna(1.0)
    chr_df.loc[:, "log10P"] = np.log10(chr_df["P"])

    # get the index SNP and calculate standard errors
    # used to calculate the approximate Bayes factor
    E.info("calculating standard errors from association "
           "p-values")
    index_snp = chr_df.iloc[0]["SNP"]
    E.info("The lead SNP is {}".format(index_snp))
    index_bp = chr_df.iloc[0]["BP"]
    z_func = lambda x: - 0.862 + sqrt(0.743 - (2.404 * np.log(x)))
    chr_df["Z"] = abs(chr_df["P"].apply(z_func))
    chr_df["SE"] = np.log(chr_df["OR"])/abs(chr_df["Z"])

    start = index_bp - region_size/2
    end = index_bp + region_size/2
    chr_df.index = chr_df["BP"]

    E.info("Fine mapping region defined as %i - %i "
           "on chromosome %i" % (start, end, int(chromosome)))

    # subsetting on range will create many NaNs due to
    # pandas broadcasting and filling in rows of DF
    sig_df = chr_df.loc[range(start, end+1)]
    sig_df.dropna(axis=0, how='all', inplace=True)
    sig_df.drop_duplicates(subset="SNP", inplace=True)
    sig_df.index = sig_df["SNP"]

    # calculate the approximate bayes factor for
    # each SNP
    E.info("calculating approximate Bayes Factors")
    bayes = {}

    # test overriding the prior on the variance
    # use the standard error on the medina log OR

    med_logor = np.log(np.median(sig_df["OR"]))
    std_logor = np.std(np.log(sig_df["OR"]))
    prior_variance = std_logor/np.sqrt(sig_df.shape[0])

    E.info("The prior variance for this fine-mapping"
           " interval is {}, and the median log OR"
           " is {:f}".format(prior_variance,
                             med_logor))

    for snp in sig_df.index:
        logor = np.log(sig_df.loc[snp]["OR"])
        se = abs(sig_df.loc[snp]["SE"])

        abf = calcApproxBayesFactor(log_or=logor,
                                    standard_error=se,
                                    prior_variance=prior_variance)
        bayes[snp] = abf

    sum_bayes = np.nansum(bayes.values())

    # calculate posterior probabilities as the proportion
    # of bayes factor/ sum all bayes factors
    E.info("calculating posterior probabilities")
    bayes_rank = pd.Series(bayes)
    bayes_rank.sort_values(inplace=True, ascending=False)
    bayes_rank = bayes_rank.fillna(0.0)

    posteriors = bayes_rank/sum_bayes
    posteriors.sort_values(ascending=False,
                           inplace=True)

    # side effect - write all ABFs and Posteriors out to file
    out_df = pd.DataFrame({"Posterior": posteriors,
                           "ApproxBayesFactor": bayes_rank,
                           "SNP": posteriors.index})
    out_df.index = out_df["SNP"]
    out_df.drop(["SNP"], axis=1, inplace=True)
    out_df.sort(["Posterior"], inplace=True, ascending=False)

    index_bayes = out_df.loc[index_snp]["ApproxBayesFactor"]
    index_p = sig_df.loc[index_snp]["log10P"]
    index_or = sig_df.loc[index_snp]["OR"]
    index_se = sig_df.loc[index_snp]["SE"]

    E.info("Bayes factor for lead SNP {} is {}, "
           "p-value {}, OR {} and SE {}".format(index_snp,
                                                index_bayes,
                                                index_p,
                                                index_or,
                                                index_se))

    return out_df


def getSnpIds(snp_set):
    '''
    Parse a text file with SNP IDs,
    one per row. Remove duplicates.

    Arguments
    ---------
    snp_set: string
      file containing SNP IDs

    Returns
    -------
    snp_list: set
      set of unique SNP IDs
    '''

    E.info("Parsing SNP set IDs")
    with IOTools.openFile(snp_set, "r") as sfile:
        snps = [sn.split("\t")[0] for sn in sfile.readlines()]
        snpset = set(snps)
        snp_list = [s.rstrip("\n") for s in snpset]

    return snp_list


def getEigenScores(eigen_dir, bim_file, snp_file):
    '''
    Extract Eigen scores from tabix-index files
    for all SNPs in a provided .bim file

    Arguments
    ---------
    eigen_dir: string
      PATH to directory containing eigen scores, with
      suffix .tab.bgz

    bim_file: string
      plink .bim file containing SNP co-ordinates
      and alleles - assumes minor allele is A2

    snp_file: string
      file containing SNP IDs, one per line

    Returns
    -------
    snp_dict: dict
      SNP eigen scores
    '''

    # setup a generic tabix query to reduce number
    # of operations

    tab_query = """
    tabix %(eigen_dir)s/%(tab_indx)s %(contig)i:%(start)i-%(end)i |
    awk '{if($4 == "%(A1)s") print $0}'
    """

    tab_dir = [td for td in os.listdir(eigen_dir) if re.search(".bgz$", td)]

    snp_list = getSnpIds(snp_file)
    E.info("SNP set of %i SNPs" % len(snp_list))

    snp_dict = {}
    E.info("Parsing SNP co-ordinates")
    # tried straightforward file parsing, took too long
    # as average .bim file contains millions of lines
    # read in chunks in to pandas DataFrame, return
    # a generator
    header = ["CHR", "SNP", "cM", "BP", "A1", "A2"]
    file_iterator = pd.read_table(bim_file, sep="\t",
                                  chunksize=50000,
                                  header=None,
                                  index_col=None,
                                  names=header)

    for dataframe in file_iterator:
        dataframe.index = dataframe["SNP"]
        try:
            snp_frame = dataframe.loc[snp_list]
            # not all SNPs will appear together in a chunk
            # remove NA rows and duplicates
            snp_frame.dropna(axis=0, how='all',
                             inplace=True)
            snp_frame.drop_duplicates(subset="SNP", keep="last",
                                      inplace=True)

            snp_frame.loc[:, "CHR"] = snp_frame["CHR"].astype(np.int64)
            contig = snp_frame["CHR"][0]
            recontig = re.compile("chr%i" % contig)
            tab_indx = [tx for tx in tab_dir if re.search(recontig,
                                                          tx)][-1]

            # redefine float types as int for output
            # prettify and reduce downstream bugs with assumed
            # data types
            snp_frame.loc[:, "BP"] = snp_frame["BP"].astype(np.int64)

            for snp in snp_frame.index:
                # open a process with query, process on the fly
                A1 = snp_frame.loc[snp, "A1"]
                A2 = snp_frame.loc[snp, "A2"]
                start = snp_frame.loc[snp, "BP"]
                end = start

                proc = subprocess.Popen(tab_query % locals(),
                                        shell=True,
                                        stdout=subprocess.PIPE)
                score_line = proc.stdout.readlines()
                if len(score_line):
                    eigen_score = score_line[0].split("\t")[-1].rstrip("\n")
                else:
                    eigen_score = np.nan

                score_dict = {"CHR": contig,
                              "BP": start,
                              "A1": A1,
                              "A2": A2,
                              "SCORE": eigen_score}

                snp_dict[snp] = score_dict
            E.info("Eigen scores found for %i SNPs" % len(snp_dict))
        except KeyError:
            pass

    return snp_dict


def getSNPs(map_file, snp_list):
    '''
    Given a SNP list with GWAS results,
    extract the relevant index

    Arguments
    ---------
    map_file: string
      plink format .map file with SNP positions
      in same order as .ped file

    snp_list: list
      list of SNP rs IDs with GWAS results

    Returns
    -------
    snp_index: dict
      dict of SNP, indices key,value pairs to select
    '''

    # variant order in the map file matters, use an ordered dict
    variants = collections.OrderedDict()
    with open(map_file, "r") as mfile:
        for snp in mfile.readlines():
            attrs = snp.split("\t")
            snpid = attrs[1]
            variants[snpid] = {"chr": attrs[0],
                               "pos": attrs[-1].strip("\n")}

    variant_ids = [vj for vi, vj in enumerate(variants.keys()) if vj in snp_list]
    variant_idx = [i for i, j in enumerate(variants.keys()) if j in snp_list]
    var_idx = dict(zip(variant_ids, variant_idx))

    return var_idx


def flipRiskAlleles(snp_index, snp_results, genos):
    '''
    Given an OR of a SNP on a binary phenotype,
    convert minor alleles to "risk" alleles, i.e.
    where OR > 1, if not, then invert allele

    Arguments
    ---------
    snp_index: list
      list of snp indices with GWAS results

    snp_results: dict
      snp:OR key, value pairs of SNPs and GWAS
      results

    genos: np.ndarray
      array of genotypes in format "11", "12" or
      "22" where 1 = minor allele, 2 = major allele.

    Returns
    -------
    risk_genos: np.ndarray
      Genotypes where all "1" alleles are risk alleles,
      not major alleles.
    '''

    genarray = np.array(genos)

    # find SNP alleles to flip
    flip = []
    for snp in snp_results.keys():
        if snp_results[snp] < 1:
            flip.append(snp_index[snp])
        else:
            pass

    E.info("Flipped alleles: %i" % len(flip))

    # swap alleles for SNPs where the minor (A1) allele
    # is protective
    # use intermediate values to avoid overwriting values
    flip_array = genarray[:, flip]
    np.place(flip_array, flip_array == "22", ["88"])
    np.place(flip_array, flip_array == "11", ["99"])
    np.place(flip_array, flip_array == "88", ["11"])
    np.place(flip_array, flip_array == "99", ["22"])

    genarray[:, flip] = flip_array

    return genarray


def parsePed(ped_file, delim="\t", compound_geno="False"):
    '''
    Parse a plink .ped file into a dataframe

    Arguments
    ---------
    ped_file: string
      Path to a plink .ped file

    delim: string
      delimiter that separates columns
      in ped_file

    compound_geno: boolean
      Whether alleles of genotype
      are separated by a whitespace or not.

    Returns
    -------
    ped_frame: pd.Core.DataFrame
      pandas dataframe representation of
      the ped_file.  Genotypes are presented
      as a numpy array.

    '''

    samples = []

    # parse the ped file, return a dataframe
    with open(ped_file, "r") as pfile:
        for indiv in pfile.readlines():
            ped_dict = {}
            indiv = indiv.strip("\n")
            indiv_split = indiv.split(delim)
            ped_dict["FID"] = indiv_split[0]
            ped_dict["IID"] = indiv_split[1]
            ped_dict["SEX"] = int(indiv_split[4])
            ped_dict["PHEN"] = int(indiv_split[5])
            ped_dict["GENOS"] = np.array(["".join(x.split(" ")) for x in indiv_split[6:]])
            samples.append(ped_dict)

    ped_frame = pd.DataFrame(samples)

    return ped_frame


def countRiskAlleles(ped_frame, snp_index, report, flag):
    '''
    Count the number of risk alleles per individual
    and calculate the probability of the phenotype

    Arguments
    ---------
    ped_frame: pd.Core.DataFrame
      Dataframe of SNP genotypes and phenotype information

    snp_index: list
      list of snp indices denoting which columns of
      ped_frame are the relevant genotypes

    report: string
      either `cases_explained` - the proportion of cases
      explained by risk allele carriage, or
      `probability_phenotype` - the probability (frequency)
      of the binary phenotype amongst all individuals given
      the risk allele carriage

    flag: boolean
      output individuals explained by carriage of 2
      risk alleles

    Returns
    -------
    count_freq: np.ndarray
      cumulative frequency array of #risk alleles
    '''

    # need to include 0 count
    case_freq = np.zeros(shape=(len(snp_index)*2) + 1,
                         dtype=np.float64)
    cntrl_freq = np.zeros(shape=(len(snp_index)*2) + 1,
                          dtype=np.float64)

    # group by phenotype column
    phen_groups = ped_frame.groupby(by="PHEN")
    for name, group in phen_groups:
        genos = group.loc[:, snp_index]

        # convert to 0,1,2 coding for ease of counting
        # treat 00 as missing/NA
        genos.replace({"22": 0,
                       "12": 1,
                       "11": 2,
                       "00": np.nan},
                      inplace=True)
        risk_sums = np.nansum(genos, axis=1)

        for val in risk_sums:
            if name == 1:
                cntrl_freq[val] += 1
            elif name == 2:
                case_freq[val] += 1

    if flag:
        explained = pd.DataFrame(risk_sums)
        explained.index = group["FID"]
        explained["IID"] = explained.index
        explained.columns = ["IID", "riskAlleles"]
        explained = explained[explained["riskAlleles"] == 2.0]
        explained.to_csv("/".join([os.getcwd(), "cases_explained.tsv"]),
                         sep="\t", index_label="FID")
    else:
        pass

    if report == "cases_explained":
        # express as the proportion of cases explained
        cumulative = np.cumsum(case_freq)/np.nansum(case_freq)
        freqs = case_freq/np.nansum(case_freq)
    elif report == "probability_phenotype":
        cumulative = np.cumsum(case_freq + cntrl_freq)/np.nansum(case_freq + cntrl_freq)
        freqs = case_freq/(case_freq + cntrl_freq)

    freqs[np.isnan(freqs)] = 0
    E.info("Individuals with pheno 1: %i" % np.nansum(cntrl_freq))
    E.info("Individuals with pheno 2: %i" % np.nansum(case_freq))

    res_dict = {"freqs": freqs, "cumulative": cumulative,
                "cases": case_freq,
                "controls": cntrl_freq}

    return res_dict


def plotRiskFrequency(bins, frequencies, counts,
                      savepath=None, ytitle=None):
    '''
    Generate a plot of #risk alleles vs.
    P(binary phenotype).

    Arguments
    ---------
    bins: list
      list of histogram bins, i.e. #risk
      alleles

    frequencies: list
      list of frequencies of binary phenotype
      corresponding to #risk allele bins

    counts: list
      counts of cases in each frequency bin. Used to
      scale the point sizes

    Returns
    -------
    None - plot is generated
    '''

    hist_df = pd.DataFrame({"bins": bins,
                            "freq": frequencies,
                            "counts": counts})

    py2ri.activate()
    R('''suppressPackageStartupMessages(library(ggplot2))''')
    R('''suppressPackageStartupMessages(library(scales))''')

    r_df = py2ri.py2ri_pandasdataframe(hist_df)
    R.assign("hist.df", r_df)

    R('''p_hist <- ggplot(hist.df, aes(x=bins, y=freq,'''
      '''size=counts)) + '''
      '''geom_point() + theme_bw() + '''
      '''xlim(c(0, dim(hist.df)[1])) + ylim(c(0, 1)) + '''
      '''labs(x="Number of Risk Alleles", '''
      '''y="%(ytitle)s")''' % locals())

    R('''png("%(savepath)s")''' % locals())
    R('''print(p_hist)''')
    R('''dev.off()''')

    return hist_df


def makeCredibleSet(probs_file, credible_set=0.95, lead_snp_indx=2,
                    filename_sep="_", snp_column=0, probs_column=1):
    '''
    Construct an N% credible set from a list of
    SNPs with posterior probabilities attached.

    If the top SNP has posterior prob >= 80%, then just this SNP
    will be output.

    Otherwise the N% credible set is output.

    In addition to the output credible set, this function
    also outputs several pieces of important information for
    the credible set:
    * The lead SNP
    * The SNP with the highest posterior probability, and whether
      this is also the lead SNP
    * The size of the credible set

    Arguments:
    ----------
    probs_file: string
      Path to a file containing SNP IDs and probabilities. It
      must have these two columns, any others are optional and
      will be ignored

    credible_set: float
      percentage of posterior probability signal to capture in
      the credible set

    lead_snp_indx: int
      0-based index of the lead SNP for the associated region.
      Used in the output file name and summary information

    filename_sep: string
      single character delimiter in the filename that can be
      used to extract the information, i.e. chromosome, position
      and lead SNP.

    snp_column: int
      0-based column number in the input file containing SNP
      IDs.

    probs_column: int
      1-based column number in the input file containing the
      posterior probabilities

    Returns:
    --------
    posterior_set: pandas.Core.DataFrame
      data frame of the N% credible set containing SNP IDs and posterior
      probabilities
    '''

    df = pd.read_table(probs_file,
                       index_col=None, sep="\t",
                       header=None)
    prob_df = df.iloc[:, [snp_column, probs_column]]
    prob_df.columns = ["SNP", "Posterior"]

    # some files may have header, others may not
    if prob_df.iloc[0, 0] == "SNP":
        prob_df = prob_df.iloc[1:, :]
    else:
        pass

    # check probabilities have been properly interpreted as floats
    # prob_df["Posterior"].astype(np.float64)
    prob_df.loc[:, "Posterior"] = pd.to_numeric(prob_df["Posterior"])

    # need to allow for non-rs IDs. check length of split file name
    # expectation is 4, if longer then 2-5 together
    split_name = probs_file.split("/")[-1].split(filename_sep)
    if len(split_name) > 4:
        lead_snp = filename_sep.join(split_name[lead_snp_indx:-1])
    else:
        lead_snp = split_name[lead_snp_indx]
    E.info("Lead SNP is {}".format(lead_snp))

    # sort by posterior signal then create credible set
    prob_df.sort(["Posterior"], inplace=True,
                 ascending=False)

    top_snp = prob_df.iloc[0, 0]
    top_prob = prob_df.iloc[0, 1]

    E.info("Top posterior signal SNP is {} with P = {}".format(top_snp,
                                                               top_prob))
    if top_snp == lead_snp:
        E.info("Lead SNP is the same as top posterior signal SNP")
    else:
        pass

    # often if the top SNP posterior probability is >= 80%
    # the remaining variants have extremely small probs
    # in that case we're only practically interested in the
    # top variant

    if top_prob >= 0.8:
        posterior_set = prob_df[:1]
        posterior_set.index = posterior_set.loc[:, "SNP"]
        posterior_set.drop(["SNP"], inplace=True, axis=1)
        E.info("Size of {}% credible set: 1".format(credible_set * 100))
    else:
        set_indx = []
        prob_set = 0.0
        for ix in range(len(prob_df.index)):
            prob_set += prob_df.iloc[ix, 1]
            set_indx.append(ix)
            if prob_set >= credible_set:
                break
            else:
                continue

        posterior_set = prob_df.iloc[set_indx]
        posterior_set.index = posterior_set.iloc[:, 0]
        posterior_set = pd.DataFrame(posterior_set.iloc[:, 1])
        posterior_set.columns = ["Posterior"]
        E.info("Size of {}% credible set: {}".format(credible_set * 100,
                                                     posterior_set.shape[0]))

    return posterior_set


def summariseResults(file_list):
    '''
    Take a list of files from SNP prioritsation
    and collate into a single table

    Arguments
    ---------
    file_list: list
      list container of input files to collate
      into the results table. File names are expected
      to follow the format:
      <contig>_<position>_<lead_snp>_<method>.tsv

    Returns:
    --------
    summary_table: pandas.Core.DataFrame
      pandas dataframe with columns:
      * lead SNP
      * credible set size
      * top set SNP
      * top set SNP probability
      * chromosome
      * lead SNP position
    '''

    # extract info from file name
    # read in file as temporary dataframe
    # extract info into a dictionary
    # convert dict into a dataframe

    name_re = re.compile(r"(?P<contig>\w{2,5})_(?P<position>\d+)_(?P<snp_id>\w+)_(?P<method>\w+).tsv")

    snp_dicts = []

    for path in file_list:
        filename = re.search(name_re, path.split("/")[-1])
        contig = filename.group("contig")
        position = filename.group("position")
        snp_id = filename.group("snp_id")

        with open(path, "r") as ofile:
            lines = ofile.readlines()
            components = [xl.split("\t") for xl in lines[1:]]
            # snp id is index 0 in first component
            top_snp = components[0][0]
            top_prob = components[0][1].rstrip("\n")
            size = len(components)

        file_dict = {"Lead_SNP": snp_id,
                     "Credible_set_size": size,
                     "Top_set_SNP": top_snp,
                     "Top_set_prob": top_prob,
                     "Chr": contig,
                     "Position": position}
        snp_dicts.append(file_dict)

    summary_table = pd.DataFrame(snp_dicts, index=range(len(snp_dicts)))
    summary_table.index = summary_table["Lead_SNP"]
    summary_table.drop(["Lead_SNP"], axis=1, inplace=True)

    return summary_table
