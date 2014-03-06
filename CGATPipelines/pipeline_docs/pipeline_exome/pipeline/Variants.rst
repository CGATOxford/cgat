== == == ==
Variants
== == == ==

Sequence Variants within Target Regions
== == == == == == == == == == == == == == == == == == == =

Variants
----------

The following table presents an overview of the variants detected within the target regions across all of the:
    term:
        `tracks`.

.. report:
    :
        VarStats.VariantSummary
    :
        render:
            table
    :
        slices:
            count, snp_count, indel_count, nalt_1, nalt_2, nalt_3, nalt_4, nalt_5

    Variants Summary

.. report:
    :
        VarStats.SnpSummary
    :
        render:
            table
    :
        slices:
            A_C, A_G, A_T, C_A, C_G, C_T, G_A, G_C, G_T, T_A, T_C, T_G

    SNV Profile

.. report:
    :
        VarStats.IndelSummary
    :
        render:
            table
    :
        slices:
            indel_length, indel_count

    Indel Summary

.. report:
    :
        VarStats.IndelSummary
    :
        render:
            line - plot
    :
        slices:
            indel_length, indel_count

    Indel Summary

.. report:
    :
        VarStats.SharedSummary
    :
        render:
            table
    :
        slices:
            no_samples, var_count

    Shared variants Summary
