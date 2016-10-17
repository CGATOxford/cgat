'''CBioPortal.py - Interface with the Sloan-Kettering cBioPortal webservice
========================================================================

The Sloan Kettering cBioPortal webservice provides access to a
database of results of genomics experiments on various cancers. The
database is organised into studies, each study contains a number of
case lists, where each list contains the ids of a set of patients, and
genetic profiles, each of which represents an assay conducted on the
patients in the case list as part of the study.

The main class here is the CBioPortal class representing a connection
to the cBioPortal Database. Query's are represented as methods of the
class. Study ids or names or case lists can be provided to the
constructor to the object, via the setDefaultStudy and
setDefaultCaseList methods or to the indevidual query methods. Where
ever possible the validity of parameters is checked *before* the query
is executed.

Whenever a query requires a genetic profile id or a list of such ids,
but none are given, the list of all profiles for which the
show_in_analysis flag is set will be used.

All of the commands provided in the webservice are implemented here
and as far as possible the name, syntax and paramter names of the
query are identical to the raw commands to the webservice. These
queries are:

* getCancerStudies,
* getCaseLists,
* getProfileData,
* getMutationData,
* getClinicalData,
* getProteinArrayInfo,
* getProteinArrayData,
* getLink,
* getOncoprintHTML.

In addition two new queries are implememented that are not part of the
webservice:

* getPercentAltered and
* getTotalAltered

These emulate the function of the website where the percent of cases
that show any alteration for the gene and profiles given are returned
(getPercentAltered, or the percent of cases that show an alteration in
any of the genes (getTotalAltered) is returned.

examples::

   gene_list = [ "TP53",
   "BCL2",
   "MYC"  ]
   portal = CBioPortal()
   portal.setDefaultStudy(study = "prad_mskcc")
   portal.setDefaultCaseList(case_set_id = "prad_all_complete")
   portal.getPercentAltered(gene_list = gene_list)

or more tersely::

   portal.CBioProtal()
   portal.getPercentAltered(study = "prad_mskcc", case_set_id = "prad_all_complete",
                            gene_list = ["TP53","BCL2","MYC"],
                            genetic_profile_id =["prad_mskcc_mrna"])

Any warnings returned by the query are stored in CBioPortal.last_warnings.

Query's that would give too long an URL are split into smaller querys
and the results combined transparently.

A commandline interface is provided for convenience, syntax::

   python CBioPortal.py [options] command(s)

Reference
---------

'''
from future.moves.urllib.request import urlopen
import re
import optparse
import sys
from CGAT import IOTools as IOTools
from collections import OrderedDict as odict
from CGAT import Experiment as E


class CBioPortal():
    """connect to the cBioPortal Database.

    If no url is specified the default url is used. A list of of valid
    study ids is retrieved from the database. This both confirms that
    the datavase is reachable, and provides cached checking for the
    ids provided. If a study or study name is provided then this is
    set as the defualt study for this session and the details of the
    availible profiles and cases is retrieved.  'Study' is the study
    id. If both study and study_name are specified then the study id
    is used.
    """

    url = "http://www.cbioportal.org/public-portal/webservice.do"

    def __init__(self, url=None, study=None, study_name=None,
                 case_list_id=None):

        if url:
            self.url = url

        self.getCancerStudies()
        self.study = None
        self.case_list = None

        if study:

            if study in self._valid_study_ids:
                self.study = study

                self.profiles = self.getGeneticProfiles(study)
                self.cases = self.getCaseLists(study)
            else:
                raise ValueError("%s is not a valid study id" % study)

        elif study_name:
            if study_name in [x['name'] for x in cancer_studies]:
                study = [x['cancer_study_id']
                         for x in cancer_studies if x['name'] == study_name][0]
                self.study = study
                self.profiles = self.getGeneticProfiles(study)
                self.cases = self.getCaseLists(study)

            else:
                raise ValueError("%s is not a valid study name" % study_name)

        if case_list_id:
            self.setDefaultCaseList(case_list_id)

    def _executeQuery(self, command, args=None):
        """execute the provided command on the database.

        args are specified as a dictionary.  error checking on the
        results is performed here and the returned value is a list of
        lists with each list representing a row of the returned table.
        """
        try:
            argument_string = "&".join(["%s=%s" % (x, args[x]) for x in args])
            command_string = "&".join([command, argument_string])
        except TypeError:
            command_string = command

        query = "%s?cmd=%s" % (self.url, command_string)

        if len(query) > 2048:
            if "gene_list" in args:
                genes = args['gene_list'].split(",")

                if len(genes) < 2:
                    raise ValueError("Request too long")

                args['gene_list'] = ",".join(genes[(len(genes) / 2):])
                query1 = self._executeQuery(command, args)
                warnings = self.last_warnings

                args['gene_list'] = ",".join(genes[:(len(genes) / 2)])
                query2 = self._executeQuery(command, args)
                self.last_warnings = warnings + self.last_warnings

                return query1 + query2

        data = urlopen(query)

        line = data.readline()
        self.last_query = query
        self.last_status = line
        self.last_warnings = []
        self.last_header = [self.last_status]
        return_table = []

        while re.match("^#", line):

            if re.match("^# Warning: (.+)", line):
                self.last_warnings.append(
                    re.match("^# Warning: (.+)", line).groups(1)[0])
                self.last_header.append(line)
                line = data.readline()
                continue

            elif re.match("^#", line):
                self.last_header.append(line)
                line = data.readline()
                continue

        if re.match("^Error: (.+)", line):
            self.last_header.append(line)
            raise CDGSError(re.match("^Error: (.+)", line).groups(1)[0], query)
        line = line.strip()
        headers = line.split("\t")

        for line in data:
            if re.match("^# Warning: (.+)", line):
                self.last_warnings.append(
                    re.match("^# Warning: (.+)", line).groups(1)[0])
                self.last_header.append(line)
                line = data.readline()
                continue
            line = line.strip()
            return_table.append(odict(list(zip(headers, line.split("\t")))))

        return return_table

    def _getStudyId(self, study, study_name):

        if study:
            if study in self._valid_study_ids:
                return study
            else:
                raise ValueError("%s is not a valid study id" % study)

        elif study_name:
            name_lookup = [x['study']
                           for x in self.cancer_studies if x['name'] == study_name]
            if len(name_lookup) == 1:
                return name_lookup[0]
            else:
                raise ValueError(
                    "Cannot find study id for study '%s'" % study_name)

        elif self.study:
            return self.study

        else:
            return None

    def getCancerStudies(self):
        """Fetches the list of cancer studies currently in the database.

        Returns list of dictionaries with three entries
        'cancer_study_id','name' and 'description'.  Also caches this
        data to verify the validity of later calls
        """
        cancer_studies = self._executeQuery("getCancerStudies")
        self.cancer_studies = cancer_studies
        self._valid_study_ids = [x['cancer_study_id'] for x in cancer_studies]

        return self.cancer_studies

    def getGeneticProfiles(self, study=None, study_name=None):
        """Fetches the valid genetic profiles for a particular study.

        study is the study id.  If both study and study_name are
        specified, study is used. If neither study nor study name is
        specified then the default study is used if set, if not a
        value error is raised.
        Returns a list of dictionaries
        """

        study = self._getStudyId(study, study_name)

        if not study:
            raise ValueError("A study must be specified if no default is set")

        genetic_studies = self._executeQuery(command="getGeneticProfiles",
                                             args=dict(cancer_study_id=study))
        return genetic_studies

    def getCaseLists(self, study=None, study_name=None):
        """Retrieves meta-data regarding all case lists stored about a
        specific cancer study.

        For example, a within a particular study, only some cases may
        have sequence data, and another subset of cases may have been
        sequenced and treated with a specific therapeutic
        protocol. Multiple case lists may be associated with each
        cancer study, and this method enables you to retrieve
        meta-data regarding all of these case lists.

        Data is returned as a list of dictionaries with the following
        entries:

         * case_list_id: a unique ID used to identify the case list ID
           in subsequent interface calls. This is a human readable
           ID. For example, "gbm_all" identifies all cases profiles in
           the TCGA GBM study.

         * case_list_name: short name for the case list.

         * case_list_description: short description of the case list.

         * cancer_study_id: cancer study ID tied to this genetic
           profile. Will match the input cancer_study_id.

         * case_ids: space delimited list of all case IDs that make up
           this case list.
        """

        study = self._getStudyId(study, study_name)

        if not study:
            raise ValueError("A study must be specified if no default is set")

        case_lists = self._executeQuery(
            command="getCaseLists",
            args=dict
            ({'cancer_study_id': study}))
        return case_lists

    def getProfileData(self, gene_list, case_set_id=None,
                       genetic_profile_id=None, study=None, study_name=None):
        """Retrieves genomic profile data for one or more genes.

        You can specify one gene and many profiles or one profile and
        many genes.  If you specify no genetic profiles then all
        genetic profiles for the specified or default study are used
        if the case_set_id is from that study otherwise a ValueError
        is raised.

        Return value depends on the parameters. If you specify a
        single genetic profile and multiple genes a list of ordered
        dictionaries with the following entries::

            gene_id: Entrez Gene ID
            common: HUGO Gene Symbol
            entries 3 - N: Data for each case

        If you specify multi genetic profiles and a single gene, a
        list of ordered dictoraries with the following entries is
        returned::

            genetic_profile_id: The Genetic Profile ID.
            alteration_type: The Genetic Alteration Type, e.g. MUTATION, MUTATION_EXTENDED, COPY_NUMBER_ALTERATION, or MRNA_EXPRESSION.
            gene_id: Entrez Gene ID.
            common: HUGO Gene Symbol.
            Columns 5 - N: Data for each case.

        """

        # Do some pre-query checking.

        if len(gene_list) > 1 and len(genetic_profile_id) > 1:
            raise ValueError(
                "%i genes and %i profiles specified\n.Please "
                "specify either one gene many profiles or one "
                "profile many genes" % (
                    len(gene_list), len(genetic_profile_id)))

        study_id = self._getStudyId(study, study_name)
        case_set_id = self._getCaseListId(case_set_id, study_id)

        genetic_profile_id = self._getAndCheckGeneticProfiles(
            genetic_profile_id=genetic_profile_id, study=study)

        gene_list = ",".join(gene_list)
        genetic_profile_id = ",".join(genetic_profile_id)

        profile_data = self._executeQuery(
            command="getProfileData",
            args={'case_set_id': case_set_id,
                  'genetic_profile_id': genetic_profile_id,
                  'gene_list': gene_list})

        return profile_data

    def getMutationData(self, gene_list, genetic_profile_id, case_set_id=None,
                        study=None, study_name=None):
        '''For data of type EXTENDED_MUTATION, you can request the full set of
        annotated extended mutation data.

        This enables you to, for example, determine which sequencing
        center sequenced the mutation, the amino acid change that
        results from the mutation, or gather links to predicted
        functional consequences of the mutation.

        Query Format


            case_set_id= [case set ID] (required)
            genetic_profile_id= [a single genetic profile IDs] (required).
            gene_list= [one or more genes, specified as HUGO Gene Symbols or
                 Entrez Gene IDs](required)

        Response Format

        A list of dictionaries with the following entires

            entrez_gene_id: Entrez Gene ID.
            gene_symbol: HUGO Gene Symbol.
            case_id: Case ID.
            sequencing_center: Sequencer Center responsible for identifying
                this mutation.
                               For example: broad.mit.edu.
            mutation_status: somatic or germline mutation status. all mutations
                         returned will be of type somatic.
            mutation_type: mutation type, such as nonsense, missense, or frameshift_ins.
            validation_status: validation status. Usually valid, invalid, or unknown.
            amino_acid_change: amino acid change resulting from the mutation.

            functional_impact_score: predicted functional impact score,
                 as predicted by: Mutation Assessor.
            xvar_link: Link to the Mutation Assessor web site.
            xvar_link_pdb: Link to the Protein Data Bank (PDB) View within
                       Mutation Assessor web site.
            xvar_link_msa: Link the Multiple Sequence Alignment (MSA) view
                 within the Mutation Assessor web site.
            chr: chromosome where mutation occurs.
            start_position: start position of mutation.
            end_position: end position of mutation.

        If a default study is set then a check will be performed to
        set if the supplied case id is from the specified study. The
        study can be over written using the study and study_name
        parameters

        '''

        study_id = self._getStudyId(study, study_name)

        case_set_id = self._getCaseListId(case_set_id, study_id)

        if study_id:
            if not study_id == self.study:
                profiles = getGeneticProfiles(study_id)
            else:
                profiles = self.profiles

            if genetic_profile_id not in proiles:
                raise ValueError(
                    "%s not a valid genetic profile for study %s" %
                    (genetic_profile_id, gene_id))
        genetic_profile_id = ",".join(genetic_profile_id)
        gene_list = ",".join(gene_list)
        mutation_data = self._executeQuery(
            command="getMutationData",
            args=dict({"case_set_id": case_set_id,
                       "genetic_profile_id": genetic_profile_id,
                       "gene_list": gene_list}))
        return mutation_data

    def getClinicalData(self, case_set_id=None, study=None, study_name=None):
        '''Retrieves overall survival, disease free survival and age at
        diagnosis for specified cases.

        Due to patient privacy restrictions, no other clinical data is
        available.

        Query Format
        ------------

            case_set_id= [case set ID] (required)

        Response Format
        ---------------

        A list of dictionaries with the following entries:

            case_id: Unique Case Identifier.
            overall_survival_months: Overall survival, in months.
            overall_survival_status: Overall survival status, usually
                 indicated as "LIVING" or "DECEASED".
            disease_free_survival_months: Disease free survival, in months.
            disease_free_survival_status: Disease free survival status,
                 usually indicated as "DiseaseFree" or "Recurred/Progressed".
            age_at_diagnosis: Age at diagnosis.

        If a study is specified or a defualt study is set, then the
        case_set_id will be tested to check if it exists for that
        study.

        '''

        study_id = self._getStudyId(study, study_name)

        case_set_id = self._getCaseListId(case_set_id, study_id)

        clincal_data = self._executeQuery(
            command="getClinicalData",
            args=dict({'case_set_id': case_set_id}))
        return clincal_data

    def getProteinArrayInfo(self, protein_array_type=None, gene_list=None,
                            study=None, study_name=None):
        '''Retrieves information on antibodies used by reverse-phase protein
        arrays (RPPA) to measure protein/phosphoprotein levels.

        Query Format
        ------------

            cancer_study_id= [cancer study ID] (required)
            protein_array_type= [protein_level or phosphorylation]
            gene_list= [one or more genes, specified as HUGO Gene
            Symbols or Entrez Gene IDs].

        Response Format
        ---------------

        A list of dictionaries with the following entires:

            ARRAY_ID: The protein array ID.
            ARRAY_TYPE: The protein array antibody type, i.e. protein_level
                or phosphorylation.
            GENE: The targeted gene name (HUGO gene symbol).
            RESIDUE: The targeted resdue(s).

        If no study is specified the default study is used. If that is
        not specified an error is raised.

        '''

        study = self._getStudyId(study, study_name)
        args = dict({'cancer_study_id': study})

        if gene_list:
            args['gene_list'] = ",".join(gene_list)

        if protein_array_type:
            args['protein_array_type'] = protein_array_type

        protein_array_info = self._executeQuery(
            command="getProteinArrayInfo", args=args)

        return protein_array_info

    def getProteinArrayData(self, protein_array_id=None, case_set_id=None,
                            array_info=0, study=None, study_name=None):
        '''Retrieves protein and/or phosphoprotein levels measured by
        reverse-phase protein arrays (RPPA).

        Query Format
        ------------

        case_set_id= [case set ID]
        protein_array_id= [one or more protein array IDs] as list.
        array_info= [1 or 0]. If 1, antibody information will also be exported.

        Response Format 1
        -----------------

        If the parameter of array_info is not specified or it is not
        1, returns a list of dictionaries with the following columns.

        ARRAY_ID: The protein array ID.
        Columns 2 - N: Data for each case.

        Response Format 2
        -----------------

        If the parameter of array_info is 1, you will receive a list
        of ordered dictionaries with the following entires:

        ARRAY_ID: The protein array ID.
        ARRAY_TYPE: The protein array antibody type, i.e. protein_level or
             phosphorylation.
        GENE: The targeted gene name (HUGO gene symbol).
        RESIDUE: The targeted resdue(s).
        Columns 5 - N: Data for each case.

        If the defualt study is set then the case_set_id will be
        check. The default study can be overidden using the study or
        study_name parameters.

        '''

        study_id = self._getStudyId(study, study_name)

        case_set_id = self._getCaseListId(case_set_id, study_id)
        args = dict({"case_set_id": case_set_id,
                     "array_info": array_info})

        if protein_array_id:
            args['protein_array_id'] = ",".join(protein_array_id)

        protein_array_data = self._executeQuery(
            command="getProteinArrayData", args=args)

        return protein_array_data

    def getLink(self, gene_list, study=None, study_name=None,  report="full"):
        '''return a perminant link to the cBioPortal report for the gene_list
            cancer_study_id=[cancer study ID] gene_list=[a comma
            separated list of HUGO gene symbols] (required)
            report=[report to display; can be one of: full (default),
            oncoprint_html]
        '''

        study = self._getStudyId(study, study_name)
        if not study:
            raise ValueError("Study must be specified")
        if report not in ["full", "oncoprint_html"]:
            raise ValueError("%s is not a valid report" % report)

        url = "/".join(self.url.split("/")[:-1])

        gene_list = ",".join(gene_list)
        return "%s/link.do?cancer_study_id=%s&gene_list=%s&report=%s" %\
            (url, study, gene_list, report)

    def getOncoprintHTML(self, gene_list, study=None, study_name=None):
        '''returns the HTML for the oncoprint report for the specified gene
        list and study'''

        url = "/".join(self.url.split("/")[0:-1])
        gene_list = ",".join(gene_list)
        command = "%s/link.do?cancer_study_id=%s&gene_list=%s&report=oncoprint_html" % (
            url, study, gene_list)
        return urlopen(command).read()

    def setDefaultStudy(self, study=None, study_name=None):
        '''sets a new study as the default study. Will check that the study
        id is valid'''
        study = self._getStudyId(study, study_name)
        self.study = study
        self.profiles = self.getGeneticProfiles()
        self.cases = self.getCaseLists()

        all_case_list = [x['case_list_id']
                         for x in self.cases
                         if x['case_list_name'] == "All Tumours"]

        if len(all_case_list) == 1:
            self.setDefaultCaseList(all_case_list)

    def setDefaultCaseList(self, case_set_id, study=None, study_name=None):
        '''set the default case list. If study is not specified the default
        study will be used.

        The study will be used to check that the case_set exists.
        '''

        study = self._getStudyId(study, study_name)
        case_list = self._getCaseListId(case_set_id, study=study)

        if not study == self.study:
            self.setDefaultStudy(study)

        self.case_list = case_list

    def _getCaseListId(self, case_set_id=None, study=None, strict=True):
        ''' checking is only done if study is specified or a default is set '''

        study_id = self._getStudyId(study, None)

        if case_set_id:
            if study_id:
                if not study_id == self.study:
                    case_lists = [x['case_list_id']
                                  for x in self.getCaseLists(study_id)]
                else:
                    case_lists = [x['case_list_id'] for x in self.cases]

                if case_set_id not in case_lists:
                    raise ValueError(
                        "%s is not a valid case list for study %s" % (case_set_id, study_id))
            return case_set_id

        else:
            if self.case_list:
                return self.case_list
            else:
                raise ValueError("No case_set_id provided and no default set")

    def _getAndCheckGeneticProfiles(self, genetic_profile_id=None, study=None):

        study_id = self._getStudyId(study, None)
        if not genetic_profile_id:
            if not study_id:
                raise ValueError(
                    "Either genetic_profile_id or study must be specified")

            if study_id == self.study:
                genetic_profile_id = [x['genetic_profile_id']
                                      for x in self.profiles
                                      if x['show_profile_in_analysis_tab'] == "true"]
            else:
                genetic_profile_id = [x['genetic_profile_id']
                                      for x in self.getGeneticProfiles(study_id)
                                      if x['show_profile_in_analysis_tab'] == "true"]

            return genetic_profile_id

        else:

            if not study_id:
                return genetic_profile_id

            if study_id == self.study:
                genetic_profile_id = [x['genetic_profile_id'] for x in self.profiles
                                      if x['genetic_profile_id'] in genetic_profile_id]
            else:
                genetic_profile_id = [x['genetic_profile_id'] for x in self.getGeneticProfiles(study_id)
                                      if x['genetic_profile_id'] in genetic_profile_id]

            if len(genetic_profile_id) == 0:
                raise ValueError("no valid genetic_profile_ids found")

            return genetic_profile_id

    def getPercentAltered(self, gene_list, study=None, study_name=None,
                          case_set_id=None, genetic_profile_id=None,
                          threshold=2):
        '''Get the percent of cases that have one or more of the specified
        alterations for each gene

        Query Format
        ------------

        study = [cancer_study_id] The study to use.

        study_name = [cancer_study_name] The name of the study to
                     use. If neither this nor study are specified,
                     then the default is used.

        case_set_id = [case_set_id] The case list to use. If not
                      specified, the default case list is used.

        gene_list = [one or more genes, specified as HUGO Gene Symobls
        or ENtrez Gene IDs] (require)

        genetic_profile_id = [one or more genetic profile IDs] If none
        specified all genetic profiles for the specified study are
        used..

        threhold = [z_score_threshold] the numeric threshold at which
        a mrna expression z-score is said to be significant.

        Response Format
        ---------------

        A list of dictionaries with the following entries
        gene_id: The Entrez Gene ID
        common: The Hugo Gene Symbol
        altered_in: The percent of cases in which the gene is altered

        One implementation note is that a guess must be made as to
        wether a returned profile value represents a alteration or
        not. Currently guesses are only made for copy number
        variation, mrna expression and mutionation

        '''

        study = self._getStudyId(study, study_name)
        case_set_id = self._getCaseListId(case_set_id, study)
        genetic_profile_id = self._getAndCheckGeneticProfiles(
            genetic_profile_id, study)

        if study and study == self.study:
            profiles = self.profiles
            case_list = [x['case_ids']
                         for x in self.cases
                         if x['case_list_id'] == case_set_id][0]

        else:
            profiles = self.getGeneticProfiles(study)
            case_list = [x['case_ids'] for x in self.getCaseLists(
                study) if x['case_list_id'] == case_set_id][0]

        return_table = []
        case_list = case_list.split(" ")

        data = []
        warnings = []
        for profile in genetic_profile_id:

            data.append(self.getProfileData(
                gene_list=gene_list, case_set_id=case_set_id,
                genetic_profile_id=[profile]))
            warnings.extend(self.last_warnings)

        # data[profile][gene][case]
        for gene in range(len(data[0])):
            cases_altered = 0.0

            geneProfile = [x[gene] for x in data]

            for case in (set(geneProfile[0]) - set(["gene_id", "common"])):
                if len([geneProfile[x][case] for x in range(len(geneProfile))
                        if self._guessAlteration(
                            geneProfile[x][case],
                            genetic_profile_id[x], profiles)]) > 0:

                    cases_altered += 1

            return_table.append(
                dict({'gene_id': geneProfile[0]['GENE_ID'],
                      'common': geneProfile[0]['COMMON'],
                      'altered_in': cases_altered * 100 / len(case_list)}))

        self.last_warnings = warnings
        return return_table

    def getTotalAltered(self, gene_list, study=None, study_name=None, case_set_id=None, genetic_profile_id=None, threshold=2):
        ''' Calculate the percent of cases in which any one of the specified genes are altered '''

        study = self._getStudyId(study, study_name)
        case_set_id = self._getCaseListId(case_set_id, study)
        genetic_profile_id = self._getAndCheckGeneticProfiles(
            genetic_profile_id, study)

        if study == self.study:
            profiles = self.profiles
            case_list = [x['case_ids']
                         for x in self.cases if x['case_list_id'] == case_set_id][0]
        else:
            profiles = self.getGeneticProfiles(study)
            case_list = [x['case_ids'] for x in self.getCaseLists(
                study) if x['case_list_id'] == case_set_id][0]
        case_list = case_list.split(" ")

        data = []

        # data[gene][profile][case]
        cases_altered = 0.0
        for gene in gene_list:
            data.append(self.getProfileData(
                gene_list=[gene], case_set_id=case_set_id,
                genetic_profile_id=genetic_profile_id))
            # catch special case where only a single
            # genetic_profile_id is passed in, so the query is single
            # gene and single profile and returns genes as columns
            # format.
            if len(genetic_profile_id) == 1 and len(data[0]) == 1:
                data[0][0]['GENETIC_PROFILE_ID'] = genetic_profile_id[0]

        for case_id in set(data[0][0]) - set(
                ["GENETIC_PROFILE_ID", "ALTERATION_TYPE",
                 "GENE_ID", "COMMON"]):

            case_altered = False
            for gene in data:

                altered = len([x for x in gene
                               if self._guessAlteration(x[case_id],
                                                        x['GENETIC_PROFILE_ID'],
                                                        profiles, threshold)])
                if (altered > 0):
                    case_altered = True

            if case_altered is True:
                cases_altered += 1.0

        return cases_altered * 100 / len(case_list)

    def _guessAlteration(self, value, genetic_profile_id, genetic_profiles,
                         threshold=2):

        alteration_type = [x['genetic_alteration_type']
                           for x in genetic_profiles
                           if x['genetic_profile_id'] == genetic_profile_id][0]

        # print alteration_type
        if alteration_type == "COPY_NUMBER_ALTERATION":
            if value == "0" or value == "-1" or value == "1":
                return False
            elif value == "NaN":
                return False
            else:
                return True

        elif alteration_type == "MRNA_EXPRESSION":
            if re.search("[^0-9\.\-]", value):
                # is character string
                return False
            elif re.search("[0-9]+\.[0-9]+", value):
                # is float
                value = float(value)
                if abs(value) > threshold:
                    return True
                else:
                    return False
            else:  # must be int?
                # print value
                if value == "0":
                    return False
                else:
                    return True

        elif alteration_type == "METHYLATION":
            return False

        elif (alteration_type == "MUTATION" or
              alteration_type == "MUTATION_EXTENDED"):
            if value == "NaN":
                return False
            else:
                return True
        else:
            return False


class CDGSError(Exception):

    '''exception that handles errors returned by querys in the database'''

    def __init__(self, error, request):
        self.error = error
        self.request = request

    def __str__(self):
        return "Request %s return error:\n%s" % (self.request, self.error)


def tableToString(intable):

    headers = "\t".join([x for x in intable[0]])

    line_list = []
    for line in intable:

        line_list.append("\t".join([str(line[x]) for x in line]))

    outTable = "\n".join(line_list)
    outTable = headers + "\n" + outTable

    return outTable


def main(argv=None):

    parser = E.OptionParser(
        version="%prog version: $Id$",
        usage=globals()["__doc__"])

    parser.add_option("-o", "--output_file", type="string", default=None,
                      help="[Optional] Filename to output results to. [default=STDOUT]")
    parser.add_option("-u", "--url", type="string", default="http://www.cbioportal.org/public-portal/webservice.do",
                      help="[Optional] Url to the cBioPortal webservice [default=%default]")

    cqueryopts = optparse.OptionGroup(
        parser, "Common parameters", "Common arguments to the query")
    cqueryopts.add_option("-s", "--study_id", dest="study_id", type="string", default=None,
                          help="[Required/OPtional]  cBioPortal ID for study [default=%default].\n This or study_name required for: getGeneticProfiles, getCaseLists, getProteinArrayInfo, getLink,getOncoprintHTML, getPercentAltered, getTotalAltered")
    cqueryopts.add_option("-n", "--study_name", dest="study_name", type="string", default=None,
                          help="[Required/Optional] cBioPortal Name for study [defualt=%default].\n See above for which commands require this.")
    cqueryopts.add_option("-c", "--case_set_id", dest="case_set_id", type="string", default=None,
                          help="[Required for some] cBioPortal case_set_id specifying the case list to use.\nRequired for getProfileData, getMutationData, getClincalData, getProteinArrayData, getPercentAltered, getTotalAltered. Default is case_set_id for case list 'All Tumours' ")
    cqueryopts.add_option("-g", "--gene_list", dest="gene_list", type="string", default=None,
                          help="[Required for some] Comma seperated list of HUGO gene symbols or Entrez gene IDs.\nRequired for getProfileData, getMutationData, getLink, getOncoprintHTML")
    cqueryopts.add_option("-f", "--gene_list_file", dest="gene_list_file", type="string", default=None,
                          help="[Optional] Filename to read in gene_list from")
    cqueryopts.add_option("-p", "--profile_id", dest="profile_id", type="string",
                          help="[Optional] Comma seperated list of cBioPortal genetic_profile_ids. If none are specified then the list of profiles for the study where display in analysis is True is used.")

    squeryopts = optparse.OptionGroup(
        parser, "Query specific parameters", "Arguments specific to a particular query")
    squeryopts.add_option("--protein_array_type", dest="protein_array_type", type="string", default="protein_level",
                          help="[Optional] Either protein_level or phosphorylation [default=%default]")
    squeryopts.add_option("--protein_array_id", dest="protein_array_id", type="string",
                          help="[Required for some] comma seperated list of one or more protein array IDs")
    squeryopts.add_option("--array_info", dest="protein_array_info", type="int",  default=0,
                          help="[Optional] If 1, antibody infomation will also be exported in a getProteinArrayData query [default=%default]")
    squeryopts.add_option("--output-report", dest="report", type="string", default="full",
                          help="[Optional] Report type to display for getLink. Either full or oncoprint_html [default=%default] ")
    squeryopts.add_option("--threshold", dest="threshold", type="int", default=2,
                          help="[Optional] Threshold for deciding if an alteration is significant for continuous metrics [default=%default]")

    parser.add_option_group(cqueryopts)
    parser.add_option_group(squeryopts)

    (options, args) = E.Start(
        parser, add_pipe_options=False, add_output_options=False, argv=argv)

    portal = CBioPortal(url=options.url, study=options.study_id,
                        study_name=options.study_name, case_list_id=options.case_set_id)

    results = []

    if options.gene_list_file:
        infile = IOTools.openFile(options.gene_list_file)
        gene_list = [x.strip() for x in infile]
    elif options.gene_list:
        gene_list = options.gene_list.split(",")

    if options.profile_id:
        profile_id = options.profile_id.split(",")
    else:
        profile_id = None

    if "getCancerStudies" in args:
        results.append(portal.getCancerStudies())

    if "getGeneticProfiles" in args:
        results.append(portal.getGeneticProfiles())

    if "getCaseLists" in args:
        results.append(portal.getCaseLists())

    if "getProfileData" in args:
        results.append(
            portal.getProfileData(gene_list=gene_list,
                                  genetic_profile_id=profile_id))

    if "getMutationData" in args:
        results.append(
            portal.getMutationData(gene_list=gene_list,
                                   genetic_profile_id=profile_id))

    if "getClinicalData" in args:
        results.append(portal.getClinicalData())

    if "getProteinArrayInfo" in args:
        results.append(portal.getProteinArrayInfo(
            gene_list=gene_list,
            protein_array_type=options.protein_array_type))

    if "getProteinArrayData" in args:
        results.append(portal.getProteinArrayData(
            protein_array_id=options.protein_array_id,
            array_info=options.array_info))

    if "getPercentAltered" in args:
        results.append(portal.getPercentAltered(
            gene_list=gene_list, genetic_profile_id=profile_id,
            threshold=options.threshold))

    if "getLink" in args:
        results.append(
            portal.getLink(gene_list=gene_list, report=options.report))

    if "getOncoprintHTML" in args:
        results.append(portal.getOncoprintHTML(gene_list=gene_list))

    if len(results) == 0:
        sys.stderr.write("No recognised query commands provided")
        sys.exit()

    if options.output_file:
        outf = IOTools.openFile(options.output_file, "w")
    else:
        outf = sys.stdout

    for result in results:
        try:
            outf.write(tableToString(result))
        except:
            outf.write(result)

    E.Stop()


if __name__ == "__main__":
    sys.exit(main(sys.argv))
