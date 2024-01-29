# disambiguateR
Genotype List String Disambiguation Using Human Leukocyte Antigen Allele Frequencies and
    a Generic String Canonicalization Algorithm
WARNING: This package is in development. Both code and documentation are in flux and are not
yet suitable for use outside the development team. Stay tuned!
    This package provides an algorithm to disambiguate Genotype List
    (GL) string data based on human leukocyte antigen (HLA) allele frequency data,
    deleted allele data, geographic region, and International Immunogenetics
    information system (IMGT) version. Arguments must
    be provided to define the behavior of the function with respect to the desired
    output, including complete or partial disambiguation, whether allele names
    should be standardized to IMGT version data or HLA accession numbers should be
    output in GL string format. Additionally, a generic function to canonicalize
    GL strings is provided, so that two non-equal GL strings that describe the same
    possible genotypes will be converted to the same, reasonably compact GL string.
    Updating and management of source data used for disambiguation, including allele
    frequency, allele name history, and deleted allele names, is made possible using
    functions that fetch data from preconfigured sources.

