"""
RateEstimation.py - utilities for computing rate estimates for codon models.
=============================================================================

:Tags: Python

"""

import Bio.Data.CodonTable
try:
    from XGram.Generator.Prebuilt import Codons
    import XGram.Exceptions
except ImportError:
    pass


def evaluateCodonPair(codon1, codon2):
    """evaluate differences between codon pair."""

    codon_table = Bio.Data.CodonTable.standard_dna_table

    changes = []
    for x in range(0, 3):
        if codon1[x] != codon2[x]:
            changes.append((x, codon1[x], codon2[x]))

    if len(changes) != 1:
        return False, False, False

    change = (changes[0][1], changes[0][2])

    is_synonymous = codon_table.forward_table[
        codon1] == codon_table.forward_table[codon2]
    is_transition = change in (("A", "G"), ("G", "A"), ("T", "C"), ("C", "T"))

    return True, is_synonymous, is_transition


def countSubstitutions(pi, Q):
    """count substitituions given a matrix Q and frequencies pi."""

    rI, rV, rS, rN = 0, 0, 0, 0

    codon_table = Bio.Data.CodonTable.standard_dna_table
    codons = list(codon_table.forward_table.keys())

    for codon_i in codons:
        for codon_j in codons:
            if codon_i == codon_j:
                continue

            is_single, is_synonymous, is_transition = evaluateCodonPair(
                codon_i, codon_j)

            v = Q[codon_i][codon_j] * pi[codon_i]

            if is_transition:
                rI += v
            else:
                rV += v

            if is_synonymous:
                rS += v
            else:
                rN += v

    return rI, rV, rS, rN


def initializeQMatrix(codons):
    """get an initialized Q matrix."""

    Q = {}
    for codon_i in codons:
        Q[codon_i] = {}
        for codon_j in codons:
            Q[codon_i][codon_j] = 0.0

    return Q


def getQMatrix(pi, Rsi, Rsv, Rni, Rnv):
    """build a q matrix.

    Diagonal elements are set to the negative of the row sums.
    The matrix is normalized such that trace of the matrix is -1.
    """

    codons = list(Bio.Data.CodonTable.standard_dna_table.forward_table.keys())

    Q = initializeQMatrix(codons)

    trace = 0.0
    for codon_i in codons:
        row_sum = 0.0
        for codon_j in codons:
            if codon_i == codon_j:
                continue

            is_single, is_synonymous, is_transition = evaluateCodonPair(
                codon_i, codon_j)

            if not is_single:
                continue

            if is_synonymous:
                if is_transition:
                    v = Rsi
                else:
                    v = Rsv
            else:
                if is_transition:
                    v = Rni
                else:
                    v = Rnv

            v *= pi[codon_j]
            Q[codon_i][codon_j] = v
            row_sum += v

        Q[codon_i][codon_i] = -row_sum
        trace += pi[codon_i] * row_sum

    for codon_i in codons:
        for codon_j in codons:
            Q[codon_i][codon_j] /= trace

    return Q, trace


def getRateMatrix(trained_model, terminals=None):
    """return a rate matrix from an xrate grammar.

    terminals: return rate matrix and frequencies for these
        terminals. If none are given, a dictionaries of
        matrices and frequencies are returned.


    """

    if terminals:
        xt = (terminals,)
    else:
        xt = trained_model.mGrammar.getTerminals()
        pis, matrices = {}, {}

    for tt in xt:
        # retrieve the terminal frequencies for all codons
        xpi = trained_model.evaluateTerminalFrequencies()[tt]

        pi = {}
        for codon, f in list(xpi.items()):
            pi["".join(codon).upper()] = f

        # retrieve the rate matrix for all codons
        xmatrix = trained_model.evaluateRateMatrix()[tt]
        matrix = {}
        for codon1, v in list(xmatrix.items()):
            x = {}
            for codon2, f in list(xmatrix[codon1].items()):
                x["".join(codon2).upper()] = v
            matrix["".join(codon1).upper()] = x

        if terminals:
            return pi, matrix
        else:
            pis[tt] = pi
            matrices[tt] = matrix

    return pis, matrices


def setFrequencies(model, mali, prefix=""):
    """set frequencies in a model according to those observed in data.

    prefix: prefix for rate parameters.

    Frequencies are labelled:
    pa0, pc0, ..., pa1, pc1, ..., pa2, pc2, ...
    """

    try:
        frequencies = Codons.getFrequenciesPerCodonPosition(
            [x.mString for x in list(mali.values())])
    except XGram.Exceptions.UsageError:
        return

    # build a dummy grammar to insert frequencies
    dummy_grammar = XGram.Model.Grammar()
    for x in range(0, 3):
        params = []
        for a in ('A', 'C', 'G', 'T'):
            params.append(
                ("%sp%s%i" % (prefix, a.lower(), x), frequencies[x][a]))
        dummy_grammar.addVariable(params)

    model.mGrammar.copyParameters(dummy_grammar,
                                  ignore_missing=True)


def getDistanceGTR(pi, matrix):
    """obtain distance from a GTR model.
    see Felsenstein 1994, pp 209
    """
    alphabet = list(pi.keys())
    Q = initializeQMatrix(alphabet)

    # fill matrix AtD and compute trace
    trace = 0
    for i in alphabet:
        row_sum = 0
        for j in alphabet:
            if i == j:
                continue
            v = pi[i] * matrix[i][j]
            Q[i][j] = v
            row_sum += v

        Q[i][i] = -row_sum
        trace += row_sum

    for i in alphabet:
        for j in alphabet:
            Q[i][j] /= trace

    return Q, trace
