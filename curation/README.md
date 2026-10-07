# Curated elements

Full-length consensi of elements curated by hand (README "Curating a family",
"Naming curated elements"). The curated-families table that maps library
families to these elements (classify.curated_families) is built from blastn
of these consensi against the library: `family_groups.py members` per element,
then `family_groups.py resolve` across related elements.

| element | class | built from | notes |
|---|---|---|---|
| Gypsy-N1_Esto | LTR/Gypsy, non-autonomous | Esto_rnd-1_family-332 extended (18 copies) | was "EstoLTR332". ~8.1 kb, LTRs ~1.7 kb, no ORFs; GC-rich tandem repeats replace gag-pol. Shares PBS, PPT, LTR termini/U5, 5' leader and the post-pol 3' region with Gypsy-1_Esto. Inserts in (TA)n. |
| Gypsy-1_Esto | LTR/Gypsy, autonomous | Esto_rnd-5_family-3054 extended (14 copies) | ~11.8 kb, LTRs ~1.6 kb; gag (681 aa) + pol (1272 aa; PR, RT, RH, INT gmr1-like). Same PBS (tRNA primer) and PPT as Gypsy-N1_Esto. Inserts in (TA)n. |
