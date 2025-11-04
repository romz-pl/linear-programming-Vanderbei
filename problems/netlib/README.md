# Netlib LP Test Problems

The [Netlib](https://www.netlib.org/lp/data/) mathematical software repository maintains a comprehensive collection of linear programming test problems.

The original readme file for Netlib repository can be found at [readme](https://www.netlib.org/lp/data/readme).

## Overview

This collection contains standard linear programming (LP) test problems that are used to benchmark optimization solvers. The problems are stored in MPS-standard input files.

## Problem Collection

### Summary Statistics

The collection includes **109** standard problems with the following characteristics:
- **Row counts**: 28 to 16,676
- **Column counts**: 32 to 22,275
- **Nonzero counts**: 88 to 151,120


## Problem Sources and Contributors

### Major Contributors
- **Michael Saunders** (Stanford SOL): 13 problems including ADLITTLE, AFIRO, PILOT, STAIR
- **Bob Fourer**: 44 problems including SCFXM series, SHIP series, GROW series
- **Nick Gould** (Harwell): 9 problems including BLEND, BOEING1, BOEING2
- **John Tomlin** (Ketron): 10 problems including BNL1, BNL2, CYCLE
- **Istvan Maros**: MAROS, MAROS-R7, MODSZK1
- **Additional contributors**: Linus Schrage, Bob Bixby, Gus Gassmann

## References and Documentation

### Key Publications
- Lustig, I.J. **"An Analysis of an Available Set of Linear Programming Test Problems"**, *Computers & Operations Research*, vol. 16, no. 2, pp. 173-184 (1989).
- Carolan et al. **"An Empirical Evaluation of the KORBX Algorithms for Military Airlift Applications"**, *Operations Research*, vol. 38, no. 2, pp. 240-248 (1990).
- Ho, J.K. and Loute, E. **"A Set of Staircase Linear Programming Test Problems"**, *Mathematical Programming*, vol. 20, pp. 245-250 (1981).

## Problem Summary Table

| Name                                   |    Rows |    Cols |     Nonzeros |       Bytes |   BR |        Optimal Value |
|----------------------------------------|--------:|--------:|-------------:|------------:|------|---------------------:|
| [25FV47](25fv47.mps)                   |     822 |    1571 |        11127 |       70477 |      |     5.5018458883E+03 |
| [80BAU3B](80bau3b.mps)                 |    2263 |    9799 |        29063 |      298952 |    B |     9.8723216072E+05 |
| [ADLITTLE](adlittle.mps)               |      57 |      97 |          465 |        3690 |      |     2.2549496316E+05 |
| [AFIRO](afiro.mps)                     |      28 |      32 |           88 |         794 |      |    -4.6475314286E+02 |
| [AGG](agg.mps)                         |     489 |     163 |         2541 |       21865 |      |    -3.5991767287E+07 |
| [AGG2](agg2.mps)                       |     517 |     302 |         4515 |       32552 |      |    -2.0239252356E+07 |
| [AGG3](agg3.mps)                       |     517 |     302 |         4531 |       32570 |      |     1.0312115935E+07 |
| [BANDM](bandm.mps)                     |     306 |     472 |         2659 |       19460 |      |    -1.5862801845E+02 |
| [BEACONFD](beaconfd.mps)               |     174 |     262 |         3476 |       17475 |      |     3.3592485807E+04 |
| [BLEND](blend.mps)                     |      75 |      83 |          521 |        3227 |      |    -3.0812149846E+01 |
| [BNL1](bnl1.mps)                       |     644 |    1175 |         6129 |       42473 |      |     1.9776292856E+03 |
| [BNL2](bnl2.mps)                       |    2325 |    3489 |        16124 |      127145 |      |     1.8112365404E+03 |
| [BOEING1](boeing1.mps)                 |     351 |     384 |         3865 |       25315 |   BR |    -3.3521356751E+02 |
| [BOEING2](boeing2.mps)                 |     167 |     143 |         1339 |        8761 |   BR |    -3.1501872802E+02 |
| [BORE3D](bore3d.mps)                   |     234 |     315 |         1525 |       13160 |    B |     1.3730803942E+03 |
| [BRANDY](brandy.mps)                   |     221 |     249 |         2150 |       14028 |      |     1.5185098965E+03 |
| [CAPRI](capri.mps)                     |     272 |     353 |         1786 |       15267 |    B |     2.6900129138E+03 |
| [CYCLE](cycle.mps)                     |    1904 |    2857 |        21322 |      166648 |    B |    -5.2263930249E+00 |
| [CZPROB](czprob.mps)                   |     930 |    3523 |        14173 |       92202 |    B |     2.1851966989E+06 |
| [D2Q06C](d2q06c.mps)                   |    2172 |    5167 |        35674 |      258038 |      |     1.2278423615E+05 |
| [D6CUBE](d6cube.mps)                   |     416 |    6184 |        43888 |      167633 |    B |     3.1549166667E+02 |
| [DEGEN2](degen2.mps)                   |     445 |     534 |         4449 |       24657 |      |    -1.4351780000E+03 |
| [DEGEN3](degen3.mps)                   |    1504 |    1818 |        26230 |      130252 |      |    -9.8729400000E+02 |
| [DFL001](dfl001.mps)                   |    6072 |   12230 |        41873 |      353192 |    B |       1.12664E+07 ** |
| [E226](e226.mps)                       |     224 |     282 |         2767 |       17749 |      |    -1.8751929066E+01 |
| [ETAMACRO](etamacro.mps)               |     401 |     688 |         2489 |       21915 |    B |    -7.5571521774E+02 |
| [FFFFF800](fffff800.mps)               |     525 |     854 |         6235 |       39637 |      |     5.5567961165E+05 |
| [FINNIS](finnis.mps)                   |     498 |     614 |         2714 |       23847 |    B |     1.7279096547E+05 |
| [FIT1D](fit1d.mps)                     |      25 |    1026 |        14430 |       51734 |    B |    -9.1463780924E+03 |
| [FIT1P](fit1p.mps)                     |     628 |    1677 |        10894 |       65116 |    B |     9.1463780924E+03 |
| [FIT2D](fit2d.mps)                     |      26 |   10500 |       138018 |      482330 |    B |    -6.8464293294E+04 |
| [FIT2P](fit2p.mps)                     |    3001 |   13525 |        60784 |      439794 |    B |     6.8464293232E+04 |
| [FORPLAN](forplan.mps)                 |     162 |     421 |         4916 |       25100 |   BR |    -6.6421873953E+02 |
| [GANGES](ganges.mps)                   |    1310 |    1681 |         7021 |       60191 |    B |    -1.0958636356E+05 |
| [GFRD-PNC](gfrd-pnc.mps)               |     617 |    1092 |         3467 |       24476 |    B |     6.9022359995E+06 |
| [GREENBEA](greenbea.mps)               |    2393 |    5405 |        31499 |      235711 |    B |    -7.2462405908E+07 |
| [GREENBEB](greenbeb.mps)               |    2393 |    5405 |        31499 |      235739 |    B |    -4.3021476065E+06 |
| [GROW15](grow15.mps)                   |     301 |     645 |         5665 |       35041 |    B |    -1.0687094129E+08 |
| [GROW22](grow22.mps)                   |     441 |     946 |         8318 |       50789 |    B |    -1.6083433648E+08 |
| [GROW7](grow7.mps)                     |     141 |     301 |         2633 |       17043 |    B |    -4.7787811815E+07 |
| [ISRAEL](israel.mps)                   |     175 |     142 |         2358 |       12109 |      |    -8.9664482186E+05 |
| [KB2](kb2.mps)                         |      44 |      41 |          291 |        2526 |    B |    -1.7499001299E+03 |
| [LOTFI](lotfi.mps)                     |     154 |     308 |         1086 |        6718 |      |    -2.5264706062E+01 |
| [MAROS](maros.mps)                     |     847 |    1443 |        10006 |       65906 |    B |    -5.8063743701E+04 |
| [MAROS-R7](maros-r7.mps)               |    3137 |    9408 |       151120 |     4812587 |      |     1.4971851665E+06 |
| [MODSZK1](modszk1.mps)                 |     688 |    1620 |         4158 |       40908 |    B |     3.2061972906E+02 |
| [NESM](nesm.mps)                       |     663 |    2923 |        13988 |      117828 |   BR |     1.4076073035E+07 |
| [PEROLD](perold.mps)                   |     626 |    1376 |         6026 |       47486 |    B |    -9.3807580773E+03 |
| [PILOT](pilot.mps)                     |    1442 |    3652 |        43220 |      278593 |    B |    -5.5740430007E+02 |
| [PILOT.JA](pilot.ja.mps)               |     941 |    1988 |        14706 |       97258 |    B |    -6.1131344111E+03 |
| [PILOT.WE](pilot.we.mps)               |     723 |    2789 |         9218 |       79972 |    B |    -2.7201027439E+06 |
| [PILOT4](pilot4.mps)                   |     411 |    1000 |         5145 |       40936 |    B |    -2.5811392641E+03 |
| [PILOT87](pilot87.mps)                 |    2031 |    4883 |        73804 |      514192 |    B |     3.0171072827E+02 |
| [PILOTNOV](pilotnov.mps)               |     976 |    2172 |        13129 |       89779 |    B |    -4.4972761882E+03 |
| [QAP8](qap8.mps)                       |     913 |    1632 |         8304 | (see NOTES) |      |     2.0350000000E+02 |
| [QAP12](qap12.mps)                     |    3193 |    8856 |        44244 | (see NOTES) |      |     5.2289435056E+02 |
| [QAP15](qap15.mps)                     |    6331 |   22275 |       110700 | (see NOTES) |      |     1.0409940410E+03 |
| [RECIPE](recipe.mps)                   |      92 |     180 |          752 |        6210 |    B |    -2.6661600000E+02 |
| [SC105](sc105.mps)                     |     106 |     103 |          281 |        3307 |      |    -5.2202061212E+01 |
| [SC205](sc205.mps)                     |     206 |     203 |          552 |        6380 |      |    -5.2202061212E+01 |
| [SC50A](sc50a.mps)                     |      51 |      48 |          131 |        1615 |      |    -6.4575077059E+01 |
| [SC50B](sc50b.mps)                     |      51 |      48 |          119 |        1567 |      |    -7.0000000000E+01 |
| [SCAGR25](scagr25.mps)                 |     472 |     500 |         2029 |       17406 |      |    -1.4753433061E+07 |
| [SCAGR7](scagr7.mps)                   |     130 |     140 |          553 |        4953 |      |    -2.3313892548E+06 |
| [SCFXM1](scfxm1.mps)                   |     331 |     457 |         2612 |       19078 |      |     1.8416759028E+04 |
| [SCFXM2](scfxm2.mps)                   |     661 |     914 |         5229 |       37079 |      |     3.6660261565E+04 |
| [SCFXM3](scfxm3.mps)                   |     991 |    1371 |         7846 |       53828 |      |     5.4901254550E+04 |
| [SCORPION](scorpion.mps)               |     389 |     358 |         1708 |       12186 |      |     1.8781248227E+03 |
| [SCRS8](scrs8.mps)                     |     491 |    1169 |         4029 |       36760 |      |     9.0429998619E+02 |
| [SCSD1](scsd1.mps)                     |      78 |     760 |         3148 |       17852 |      |     8.6666666743E+00 |
| [SCSD6](scsd6.mps)                     |     148 |    1350 |         5666 |       32161 |      |     5.0500000078E+01 |
| [SCSD8](scsd8.mps)                     |     398 |    2750 |        11334 |       65888 |      |     9.0499999993E+02 |
| [SCTAP1](sctap1.mps)                   |     301 |     480 |         2052 |       14970 |      |     1.4122500000E+03 |
| [SCTAP2](sctap2.mps)                   |    1091 |    1880 |         8124 |       57479 |      |     1.7248071429E+03 |
| [SCTAP3](sctap3.mps)                   |    1481 |    2480 |        10734 |       78688 |      |     1.4240000000E+03 |
| [SEBA](seba.mps)                       |     516 |    1028 |         4874 |       38627 |   BR |     1.5711600000E+04 |
| [SHARE1B](share1b.mps)                 |     118 |     225 |         1182 |        8380 |      |    -7.6589318579E+04 |
| [SHARE2B](share2b.mps)                 |      97 |      79 |          730 |        4795 |      |    -4.1573224074E+02 |
| [SHELL](shell.mps)                     |     537 |    1775 |         4900 |       38049 |    B |     1.2088253460E+09 |
| [SHIP04L](ship04l.mps)                 |     403 |    2118 |         8450 |       57203 |      |     1.7933245380E+06 |
| [SHIP04S](ship04s.mps)                 |     403 |    1458 |         5810 |       41257 |      |     1.7987147004E+06 |
| [SHIP08L](ship08l.mps)                 |     779 |    4283 |        17085 |      117083 |      |     1.9090552114E+06 |
| [SHIP08S](ship08s.mps)                 |     779 |    2387 |         9501 |       70093 |      |     1.9200982105E+06 |
| [SHIP12L](ship12l.mps)                 |    1152 |    5427 |        21597 |      146753 |      |     1.4701879193E+06 |
| [SHIP12S](ship12s.mps)                 |    1152 |    2763 |        10941 |       82527 |      |     1.4892361344E+06 |
| [SIERRA](sierra.mps)                   |    1228 |    2036 |         9252 |       76627 |    B |     1.5394362184E+07 |
| [STAIR](stair.mps)                     |     357 |     467 |         3857 |       27405 |    B |    -2.5126695119E+02 |
| [STANDATA](standata.mps)               |     360 |    1075 |         3038 |       26135 |    B |     1.2576995000E+03 |
| [STANDGUB](standgub.mps)               |     362 |    1184 |         3147 |       27836 |    B |          (see NOTES) |
| [STANDMPS](standmps.mps)               |     468 |    1075 |         3686 |       29839 |    B |     1.4060175000E+03 |
| [STOCFOR1](stocfor1.mps)               |     118 |     111 |          474 |        4247 |      |    -4.1131976219E+04 |
| [STOCFOR2](stocfor2.mps)               |    2158 |    2031 |         9492 |       79845 |      |    -3.9024408538E+04 |
| [STOCFOR3](stocfor3.mps)               |   16676 |   15695 |        74004 | (see NOTES) |      |    -3.9976661576E+04 |
| [TRUSS](truss.mps)                     |    1001 |    8806 |        36642 | (see NOTES) |      |     4.5881584719E+05 |
| [TUFF](tuff.mps)                       |     334 |     587 |         4523 |       29439 |    B |     2.9214776509E-01 |
| [VTP.BASE](vtp.base.mps)               |     199 |     203 |          914 |        8175 |    B |     1.2983146246E+05 |
| [WOOD1P](wood1p.mps)                   |     245 |    2594 |        70216 |      328905 |      |     1.4429024116E+00 |
| [WOODW](woodw.mps)                     |    1099 |    8405 |        37478 |      240063 |      |     1.3044763331E+00 |

## Notes

- BR: Indicates bounded and ranged constraints
- B: Indicates bounded constraints
- **: Special notation for DFL001
- Some entries reference additional notes for specific values (Bytes or Optimal Value columns)
