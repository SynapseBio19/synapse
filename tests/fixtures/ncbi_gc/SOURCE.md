# NCBI genetic code table (gc.prt)

- Source: https://ftp.ncbi.nlm.nih.gov/entrez/misc/data/gc.prt
- Retrieved: 2026-11-19, version 4.6 (header of the file), 16,444 bytes
- SHA-256: 2aecbdd09ecc5e0d233c105090a593737298ef93d88f6bc0c1bef26df6051a04
- Stored verbatim. NCBI data are in the public domain (US government work).

`core/seqtools.py` copies the 27 tables from this file;
`tests/test_genetic_codes_ncbi.py` checks every table against it.
