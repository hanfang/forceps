marlin
==============
marlin is a tool for designing RNA antisense oligonucleotides. It is built with the Smith–Waterman algorithm for end gap free local sequence alignment, while allowing both Wobble pairing and Watson-Crick pairing. 

marlin reports the number of edit distance (ED), mismatches (MM), insertions (INS), deletions (DEL), back-to-back (B2B) ED, and ED near either end of the query sequence.

marlin report format:

chr	start	end	query	strand	edit_distance	mis_match	insertion	deletion	back_to_back	near_edge	Wobble

marlin is currently under development.

Contact: Han Fang hanfang.cshl@gmail.com

https://github.com/hanfang/marlin
