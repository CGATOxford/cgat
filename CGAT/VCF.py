'''
VCF.py - Tools for working with VCF files
=========================================

The parser for VCF files is very simplistic.

.. note::
   Another way to access the information in :term:`vcf` formatted
   files is through pysam_.

The Variant Call Format (:term:`vcf`) is described
at http://www.1000genomes.org/wiki/doku.php?id=1000_genomes:analysis:vcf4.0

Reference
---------

'''
import sys


class VCFEntry:
    """A VCF Entry"""

    def __init__(self, data, samples):

        assert len(data) == len(samples) + 9
        self.contig, self.pos, self.id, self.ref, self.alt, self.qual, \
            self.filter, self.info, self.format = \
            data[:9]

        self.genotypes = dict(list(zip(samples, data[9:])))
        self.order = samples

    def __str__(self):
        return "\t".join(map(str, (
            self.contig, self.pos, self.id, self.ref, self.alt, self.qual,
            self.filter, self.info, self.format,
            "\t".join([self.genotypes[x] for x in self.order]))))


class VCFFile:
    """A VCF File"""

    def __init__(self, infile):

        self.infile = infile
        self.format = {}
        self.info = {}
        self.fileformat = None

        while 1:
            line = self.infile.readline()

            if line.startswith("##"):
                self.addMetaFromLine(line[2:-1])
                continue
            elif line.startswith("#CHROM"):
                self.samples = line[:-1].split("\t")[9:]
                continue
            elif line.startswith("#"):
                continue
            break
        self.line = line

    def writeHeader(self, outfile, order=None):
        outfile.write("##fileformat=%s\n" % self.fileformat)
        for key, values in self.format.items():
            outfile.write("##FORMAT=%s,%s\n" % (key, ",".join(values)))
        for key, values in self.info.items():
            outfile.write("##INFO=%s,%s\n" % (key, ",".join(values)))
        outfile.write(
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t")
        if order:
            assert len(order) == len(self.samples), \
                "number of samples do not match: %i != %i" % (
                    len(order), len(self.samples))
            outfile.write("\t".join(order))
        else:
            outfile.write("\t".join(self.samples))
        outfile.write("\n")

    def __iter__(self):
        return self

    def addMetaFromLine(self, line):

        key, value = line.split("=", 1)
        if key == "INFO":
            data = value.split(",")
            self.info[data[0]] = data[1:]
        elif key == "FORMAT":
            data = value.split(",")
            self.format[data[0]] = data[1:]
        elif key == "fileformat":
            self.fileformat = value

    def __next__(self):

        data = self.line[:-1].split("\t")
        self.line = self.infile.readline()
        if not self.line:
            raise StopIteration
        return VCFEntry(data, self.samples)

    def next(self):
        return self.__next__()

if __name__ == "__main__":

    inf = VCFFile(sys.stdin)

    for x in inf:
        print(str(x))
