########################################
########################################
# classes and functions for parsing
# output from lcamapper.sh
########################################
########################################


class LCA(object):
    '''
    lca class describing the taxa associateed with a sequence
    '''

    def __init__(self):

        self.identifier = None
        self.kingdom = None
        self.kingdom_plus = None
        self.phylum = None
        self.phylum_plus = None
        self._class = None
        self._class_plus = None
        self.order = None
        self.order_plus = None
        self.family = None
        self.family_plus = None
        self.genus = None
        self.genus_plus = None
        self.species = None
        self.species_plus = None
        self.subspecies = None
        self.subspecies_plus = None

    def parse(self, line):
        '''
        parse the line
        '''
        data = line.split(";")
        self.identifier = data[0]
        for taxa in data[2:]:
            taxa = taxa.strip()
            # ignore root
            if "root" in taxa:
                continue
            if "[" not in taxa:
                continue
            taxa = taxa.split(" ")
            level, tax = taxa[0], taxa[1:]
            if len(tax) > 1:
                tax = "_".join(tax)
            else:
                tax = tax[0]
            if "Kingdom+" in level:
                self.kingdom_plus = tax
            elif "Kingdom" in level:
                self.kingdon = tax
            elif "Phylum+" in level:
                self.phylum_plus = level
            elif "Phylum" in level:
                self.phylum = tax
            elif "Class+" in level:
                self._class_plus = tax
            elif "Class" in level:
                self._class = tax
            elif "Order+" in level:
                self.order_plus = tax
            elif "Order" in level:
                self.order = tax
            elif "Family+" in level:
                self.family_plus = tax
            elif "Family" in level:
                self.family = tax
            elif "Genus+" in level:
                self.genus_plus = tax
            elif "Genus" in level:
                self.genus = tax
            elif "Species+" in level:
                self.species_plus = tax
            elif "Species" in level:
                self.species = tax
            elif "Subspecies+" in level:
                self.subspecies_plus = tax
            elif "Subspecies" in level:
                self.subspecies = tax

        if not self.kingdom:
            self.kingdom = "NA"
        if not self.kingdom_plus:
            self.kingdom_plus = "NA"
        if not self.phylum:
            self.phylum = "NA"
        if not self.phylum_plus:
            self.phylum_plus = "NA"
        if not self._class:
            self._class = "NA"
        if not self._class_plus:
            self._class_plus = "NA"
        if not self.order:
            self.order = "NA"
        if not self.order_plus:
            self.order_plus = "NA"
        if not self.family:
            self.family = "NA"
        if not self.family_plus:
            self.family_plus = "NA"
        if not self.genus:
            self.genus = "NA"
        if not self.genus_plus:
            self.genus_plus = "NA"
        if not self.species:
            self.species = "NA"
        if not self.species_plus:
            self.species_plus = "NA"
        if not self.subspecies:
            self.subspecies = "NA"
        if not self.subspecies_plus:
            self.subspecies_plus = "NA"

        return self


###############################
###############################
###############################


def iterate(infile):
    '''
    LCA results iterator
    '''
    for line in infile.readlines():
        lca = LCA()
        lca = lca.parse(line)
        yield lca
