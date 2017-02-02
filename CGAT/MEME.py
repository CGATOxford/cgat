import re


class MemeMotif:

    primary_id = None
    motif_line = None
    motif_len = None
    letter_probability_line = None
    letter_probability_matrix = None
    evalue = None
    properties = None
    alphabet = 'ACGT'

    def __init__(self, motif_line):

        try:
            self.primary_id = re.match("MOTIF\s+(\S+)", motif_line).groups()[0]
        except AttributeError:
            raise ValueError('Motif line: %s does not match re: "MOTIF\s(\S+)"'
                             % motif_line)

        try:
            self.secondary_id = re.match(
                "MOTIF\s+\S+\s+(\S+)", motif_line).groups()[0]
        except AttributeError:
            self.secondary_id = None

        self.motif_line = motif_line
        self.properties = {}
        self.letter_probability_matrix = []

    def parse_probability_line(self, line):
        
        line = line.strip()
        properties = re.findall("(\S+)=\s?(\S+)", line)

        for key, value in properties:
  
            try:
                value = int(value)
            except ValueError:
                try:
                    value = float(value)
                except ValueError:
                    pass

            self.properties[key] = value
        
        if "E" in self.properties:
            self.evalue = self.properties["E"]

        if "nsites" in self.properties:
            self.nsites = self.properties["nsites"]

        if "w" in self.properties:
            self.motif_len = self.properties["w"]
             
        self.letter_probability_line = line

    def consensus(self):
        ''' Returns a sequence defined as most probable letter at each base
        No ambiguity will be used '''

        probs = [map(float, line) for line in self.letter_probability_matrix]
        consensus = [self.alphabet[line.index(max(line))] for line in probs]
        return ''.join(consensus)
            
    def __str__(self):

        letter_probability_matrix = "\n".join(
            [" ".join(line) for line in self.letter_probability_matrix])

        return self.motif_line + "\n\n" + self.letter_probability_line + \
            "\n" + letter_probability_matrix + "\n"


class MotifList:

    def __init__(self, motifs=None):
        
        if motifs is None:
            self.motifs = []
        else:
            self.motifs = motifs

    def __getitem__(self, key):

        try:
            return self.motifs[key]
        except TypeError or IndexError:

            for motif in self.motifs:
                if motif.primary_id == key:
                    return motif

            raise KeyError("%s not in list" % key)

    def sort(self, key, *args, **kwargs):

        self.motifs.sort(key=lambda motif: motif.__dict__[key],
                         *args, **kwargs)

    def append(self, motif):
        
        self.motifs.append(motif)

    def remove(self, key_or_motif):

        try:
            self.motifs.remove(key_or_motif)
        except ValueError:
            self.motifs.remove(self[key_or_motif])

    def take(self, key_or_index):
        '''Removes a motif from the list and returns it'''
        
        motif = self[key_or_index]
        self.remove(motif)
        return motif

    def __iter__(self):

        return iter(self.motifs)

    def __len__(self):
        
        return len(self.motifs)

    def keys(self):

        return [motif.primary_id for motif in self.motifs]

    def extend(self, other):

        try:
            self.motifs += other.motifs
        except AttributeError:
            self.motifs += other


class MemeMotifFile(MotifList):

    header = None
    motif_line_re = re.compile("^MOTIF")
    matrix_line_re = re.compile(
        "letter-probability matrix: .+w=\s?([0-9]+).+E=\s?(.+)")

    def __init__(self, buffer_or_Motiffile):
        
        MotifList.__init__(self)

        try:
            self.header = buffer_or_Motiffile.header
        except AttributeError:

            self.header = []
            for line in buffer_or_Motiffile:
                line.strip()
                if self.motif_line_re.match(line):
                    current_motif = MemeMotif(line)
                    break
                else:
                    self.header.append(line)
                    
            self.header = '\n'.join(self.header)

            for line in buffer_or_Motiffile:
                line.strip()

                if self.motif_line_re.match(line):
                    self.motifs.append(current_motif)
                    current_motif = MemeMotif(line)

                if self.matrix_line_re.match(line):
                    current_motif.parse_probability_line(line)
                    for i in range(current_motif.motif_len):
                        line = buffer_or_Motiffile.next()
                        line.strip()
                        current_motif.letter_probability_matrix.append(
                            line.split())
            else:
                self.motifs.append(current_motif)

    def __str__(self):
        
        return self.header + "\n" + '\n'.join(map(str, self.motifs))


class MotifCluster(MotifList):

    def __init__(self, seed):

        MotifList.__init__(self)

        assert len(self.motifs) == 0
        self.seed = seed
        self.motifs.append(seed)
        assert len(self.motifs) == 1
