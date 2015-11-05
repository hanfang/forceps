forceps
==============
forceps is a tool for designing RNA antisense oligonucleotides. It is built with the Smith–Waterman algorithm for end gap free local sequence alignment, while allowing both Wobble pairing and Watson-Crick pairing. 

forceps reports the number of edit distance (ED), mismatches (MM), insertions (INS), deletions (DEL), back-to-back (B2B) ED, and ED near either end of the query sequence.

## forceps report format:

| chr | start | end | query | strand | ED | MM | INS | DEL | B2B | near_edge | Wobble |

forceps is currently under development.

Contact: Han Fang hanfang.cshl@gmail.com

https://github.com/hanfang/forceps
