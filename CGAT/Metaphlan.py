import collections


class ReadMap(object):

    '''
    Each read is mapped to a 
    taxonomic group
    '''

    def __init__(self):
        self.seq_id = None
        self.kingdom = None
        self.phylum = None
        self.c_lass = None
        self.order = None
        self.family = None
        self.genus = None
        self.species = None

    def read(self, seq_id, kingdom, phylum, c_lass, order, family, genus, species):

        self.seq_id, self.kingdom, self.phylum, self.c_lass, self.order, self.family, self.genus, self.species = seq_id, kingdom, phylum, c_lass, order, family, genus, species
        return self


class RelativeAbundance(object):

    '''
    parsing the output from relative
    abundance
    '''

    def __init__(self):
        '''
        each entry will have an associated
        * taxonomic level
        * name within the taxonomic group
        * abundance estimate
        '''
        self.taxon_level = None
        self.taxon = None
        self.abundance = 0

    def read(self, taxon_level, taxon, abundance):
        '''
        read the data
        '''
        self.taxon_level, self.taxon, self.abundance = taxon_level, taxon, abundance
        return self


def relative_abundance_iterator(infile):
    '''
    read in rel_ab outfile from metaphlan 
    and return relative abundances for
    each taxonomic group
    '''
    taxons = collections.defaultdict(list)
    abundances = collections.defaultdict(list)
    for line in infile.readlines():
        data = line[:-1].split("\t")
        abundance = data[1]
        entries = data[0].split("|")
        taxon = entries[len(entries) - 1]
        taxon_names = taxon.split("__")
        if len(taxon_names) == 1:  # unclassified set
            taxons[taxon_names[0]].append(taxon_names[0])
            abundances[taxon_names[0]].append(abundance)
        else:
            taxons[taxon_names[0]].append(taxon_names[1])
            abundances[taxon_names[0]].append(abundance)

    # return the taxonomic group and the specific
    # names within the taxonomic group
    for group in sorted(taxons.keys()):
        for abundance in zip(taxons[group], abundances[group]):
            if group == "k":
                group = "kingdom"
            elif group == "p":
                group = "phylum"
            elif group == "c":
                group = "class"
            elif group == "o":
                group = "order"
            elif group == "f":
                group = "family"
            elif group == "g":
                group = "genus"
            elif group == "s":
                group = "species"
            else:
                group = group
            yield RelativeAbundance().read(group, abundance[0], abundance[1])


def read_map_iterator(infile):
    '''
    iterate over read_map file
    from metaphlan
    '''
    for line in infile.readlines():
        data = line[:-1].split("\t")
        seq_id = data[0]
        data = data[1].split("|")

        # The data follow a hierarchy from more specific to
        # more general taxonomic classes so can iteratively sort them out
        data = [x.split("__")[1] for x in data]
        if len(data) == 7:
            kingdom, phylum, c_lass, order, family, genus, species = data[
                0], data[1], data[2], data[3], data[4], data[5], data[6]
        elif len(data) == 6:
            kingdom, phylum, c_lass, order, family, genus, species = data[
                0], data[1], data[2], data[3], data[4], data[5], "unclassified"
        elif len(data) == 5:
            kingdom, phylum, c_lass, order, family, genus, species = data[0], data[
                1], data[2], data[3], data[4], "unclassified", "unclassified"
        elif len(data) == 4:
            kingdom, phylum, c_lass, order, family, genus, species = data[0], data[
                1], data[2], data[3], "unclassified", "unclassified", "unclassified"
        elif len(data) == 3:
            kingdom, phylum, c_lass, order, family, genus, species = data[0], data[
                1], data[2], "unclassified", "unclassified", "unclassified", "unclassified"
        elif len(data) == 2:
            kingdom, phylum, c_lass, order, family, genus, species = data[0], data[
                1], "unclassified", "unclassified", "unclassified", "unclassified", "unclassified"
        elif len(data) < 2:
            raise ValueError("could not assign taxonomy at the phylum level")
        yield ReadMap().read(seq_id, kingdom, phylum, c_lass, order, family, genus, species)


class Counter(ReadMap):

    '''
    counter class for taxonomic counts
    '''

    def __init__(self):
        self.counts = collections.defaultdict(int)
        self.total = 0

    def count(self, infile):
        '''
        get counts for taxonomic groups i.e. the total of reads 
        that were assigned to a taxonomic group in that
        particular taxonomic group
        '''
        taxonomies = [
            "kingdom", "phylum", "class", "order", "family", "genus", "species"]
        for read in read_map_iterator(infile):
            self.total += 1
            if read.kingdom != "unclassified":
                self.counts["kingdom"] += 1
            if read.phylum != "unclassified":
                self.counts["phylum"] += 1
            if read.c_lass != "unclassified":
                self.counts["class"] += 1
            if read.order != "unclassified":
                self.counts["order"] += 1
            if read.family != "unclassified":
                self.counts["family"] += 1
            if read.genus != "unclassified":
                self.counts["genus"] += 1
            if read.genus != "unclassified":
                self.counts["species"] += 1
        return self.counts

    def total_count(self, infile):
        for read in read_map_iterator(infile):
            self.total += 1
        return self.total

    def count_kingdom(self, infile):
        counts = self.count(infile)
        return x["kingdom"]

    def count_phylum(self, infile):
        counts = self.count(infile)
        return counts["phylum"]

    def count_class(self, infile):
        counts = self.count(infile)
        return counts["class"]

    def count_order(self, infile):
        counts = self.count(infile)
        return counts["order"]

    def count_family(self, infile):
        counts = self.count(infile)
        return counts["family"]

    def count_genus(self, infile):
        counts = self.count(infile)
        return counts["genus"]

    def count_species(self, infile):
        counts = self.count(infile)
        return counts["species"]

    def count_total_species(self, infile):
        species = set()
        for read in read_map_iterator(infile):
            if read.species.find("unclassified") == -1:
                continue
            species.add(read.species)
        return len(species)

    def proportion_with_clade_assignment(self, infile, total_reads):
        return float(total_reads) / self.total_count(infile)
