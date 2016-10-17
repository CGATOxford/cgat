

class LCA(object):
    '''
    lca class describing the taxa associateed with a sequence
    '''

    def __init__(self):

        self.identifier = None
        self.domain = None
        self.superkingdom_plus = None
        self.superkingdom_plus_plus = None
        self.kingdom = None
        self.kingdom_plus = None
        self.kingdom_plus_plus = None
        self.phylum = None
        self.phylum_plus = None
        self.phylum_plus_plus = None
        self._class = None
        self._class_plus = None
        self._class_plus_plus = None
        self.order = None
        self.order_plus = None
        self.order_plus_plus = None
        self.family = None
        self.family_plus = None
        self.family_plus_plus = None
        self.genus = None
        self.genus_plus = None
        self.genus_plus_plus = None
        self.species = None
        self.species_plus = None
        self.species_plus_plus = None
        self.level = None

    def parse(self, line):
        '''
        parse the line
        '''
        data = line.split(";")
        self.identifier = data[0]
        for taxa in data[2:]:
            taxa = taxa.strip()
            # ignore root
            if "[root" in taxa:
                continue
            if "[" not in taxa:
                continue
            taxa = taxa.split(" ")

            level, tax = taxa[0], taxa[1:]
            level = level.replace("[", "").replace("]", "")

            if len(tax) > 1:
                tax = "_".join(tax)
            else:
                tax = tax[0]

            # Eukaryotes have a lot of sub-taxa e.g. Phylum+++++++.
            # Superkingdom is taken as the Domain and then kingdom
            # and so on. At the moment we are only going to consider
            # up to two sub-levels per global level
            if level == "SuperKingdom":
                self.domain = tax
            elif level == "SuperKingdom+":
                self.superkingdom_plus = tax
            elif level == "SuperKingdom++":
                self.superkingdom_plus_plus = tax
            elif level == "Kingdom":
                self.kingdom = tax
            elif level == "Kingdom+":
                self.kingdom_plus = tax
            elif level == "Kingdom++":
                self.kingdom_plus_plus = tax
            elif level == "Phylum":
                self.phylum = tax
            elif level == "Phylum+":
                self.phylum_plus = tax
            elif level == "Phylum++":
                self.phylum_plus_plus = tax
            elif level == "Class":
                self._class = tax
            elif level == "Class+":
                self._class_plus = tax
            elif level == "Class++":
                self._class_plus_plus = tax
            elif level == "Order":
                self.order = tax
            elif level == "Order+":
                self.order_plus = tax
            elif level == "Order++":
                self.order_plus_plus = tax
            elif level == "Family":
                self.family = tax
            elif level == "Family+":
                self.family_plus = tax
            elif level == "Family++":
                self.family_plus_plus = tax
            elif level == "Genus":
                self.genus = tax
            elif level == "Genus+":
                self.genus_plus = tax
            elif level == "Genus++":
                self.genus_plus_plus = tax
            elif level == "Species":
                self.species = tax
            elif level == "Species+":
                self.species_plus = tax
            elif level == "Species++":
                self.species_plus_plus = tax

        # make NA id doesn't exist in taxonomy
        if not self.domain:
            self.domain = "NA"
        if not self.superkingdom_plus:
            self.superkingdom_plus = "NA"
        if not self.superkingdom_plus_plus:
            self.superkingdom_plus_plus = "NA"
        if not self.kingdom:
            self.kingdom = "NA"
        if not self.kingdom_plus:
            self.kingdom_plus = "NA"
        if not self.kingdom_plus_plus:
            self.kingdom_plus_plus = "NA"
        if not self.phylum:
            self.phylum = "NA"
        if not self.phylum_plus:
            self.phylum_plus = "NA"
        if not self.phylum_plus_plus:
            self.phylum_plus_plus = "NA"
        if not self._class:
            self._class = "NA"
        if not self._class_plus:
            self._class_plus = "NA"
        if not self._class_plus_plus:
            self._class_plus_plus = "NA"
        if not self.order:
            self.order = "NA"
        if not self.order_plus:
            self.order_plus = "NA"
        if not self.order_plus_plus:
            self.order_plus_plus = "NA"
        if not self.family:
            self.family = "NA"
        if not self.family_plus:
            self.family_plus = "NA"
        if not self.family_plus_plus:
            self.family_plus_plus = "NA"
        if not self.genus:
            self.genus = "NA"
        if not self.genus_plus:
            self.genus_plus = "NA"
        if not self.genus_plus_plus:
            self.genus_plus_plus = "NA"
        if not self.species:
            self.species = "NA"
        if not self.species_plus:
            self.species_plus = "NA"
        if not self.species_plus_plus:
            self.species_plus_plus = "NA"

        return self


def iterate(infile):
    '''
    LCA results iterator
    '''
    for line in infile.readlines():
        lca = LCA()
        lca = lca.parse(line)
        yield lca
