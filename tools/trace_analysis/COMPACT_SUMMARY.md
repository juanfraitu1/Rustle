# Compact divergence summary per locus

Per locus: StringTie final output, kill reasons, and notable patterns.

## 44669742-44715973(+)

### StringTie kill-reason histogram

| reason | count |
|---|---|
| BAD_MM_NEG | 153 |
| GOOD | 140 |
| EXISTING_BUNDLE | 140 |
| RETAINED_INTRON | 43 |
| NO_MAXC | 27 |
| fwd_fail | 27 |
| BAD_NO_STRAND | 20 |
| zero_flux | 17 |
| NO_MAXP | 17 |
| back_fail | 17 |
| DEMOTED_RIGHT | 5 |
| INCLUDED | 5 |
| INCOMPLETE_PATH | 4 |

### FINAL decisions: 26 KEEP, 58 DROP

#### Kept

- pred[0] 44409602-44509287 cov=174.36 - exons=13
- pred[1] 44409602-44509287 cov=4.76 - exons=12
- pred[2] 44409602-44509287 cov=3.14 - exons=12
- pred[3] 44409602-44509287 cov=2.94 - exons=14
- pred[10] 44509333-44566790 cov=6.24 + exons=12
- pred[13] 44566655-44607934 cov=100.72 - exons=12
- pred[14] 44566655-44607934 cov=56.30 - exons=13
- pred[15] 44566655-44607934 cov=55.85 - exons=12
- pred[16] 44566655-44607934 cov=47.93 - exons=13
- pred[17] 44566655-44608798 cov=3.74 - exons=13
- pred[18] 44566655-44607934 cov=3.24 - exons=13
- pred[35] 44607586-44643633 cov=11.41 + exons=12
- pred[37] 44608321-44643633 cov=351.29 + exons=11
- pred[38] 44608321-44643633 cov=8.12 + exons=10
- pred[41] 44646789-44670832 cov=182.24 - exons=11
- pred[42] 44646789-44670832 cov=19.06 - exons=11
- pred[44] 44646789-44670832 cov=4.92 - exons=10
- pred[45] 44646789-44670832 cov=3.48 - exons=12
- pred[46] 44646789-44670832 cov=3.40 - exons=11
- pred[47] 44646789-44670832 cov=3.13 - exons=12
- pred[48] 44646789-44670832 cov=2.91 - exons=7
- pred[58] 44669742-44715973 cov=21.55 + exons=11
- pred[78] 44670958-44715973 cov=1.17 + exons=10
- pred[80] 44750688-44781194 cov=3.42 + exons=5
- pred[81] 44758672-44781194 cov=6.20 + exons=6
- pred[83] 45034136-45079269 cov=2.67 + exons=5

#### Dropped (sample first 10)

- pred[4] 44409602-44471342 cov=0.97 - exons=11
- pred[5] 44409602-44509287 cov=0.92 - exons=12
- pred[6] 44409602-44509287 cov=0.88 - exons=11
- pred[7] 44409602-44509287 cov=0.82 - exons=11
- pred[8] 44444091-44509287 cov=1.50 - exons=8
- pred[9] 44450177-44509287 cov=0.93 - exons=7
- pred[11] 44564480-44566790 cov=5.30 + exons=1
- pred[12] 44565698-44607934 cov=2.58 - exons=11
- pred[19] 44566655-44607934 cov=2.49 - exons=11
- pred[20] 44566655-44607934 cov=1.92 - exons=14

## 17433800-17456446(+)

### StringTie kill-reason histogram

| reason | count |
|---|---|
| BAD_MM_NEG | 98 |
| GOOD | 93 |
| EXISTING_BUNDLE | 91 |
| RETAINED_INTRON | 38 |
| zero_flux | 24 |
| NO_MAXP | 24 |
| back_fail | 24 |
| fwd_fail | 22 |
| BAD_NO_STRAND | 21 |
| NO_MAXC | 20 |
| DEMOTED_RIGHT | 6 |
| DEMOTED_LEFT | 2 |
| EXCLUDE_NO_COV | 2 |
| SINGLETON_ARTIFACT | 1 |
| SECONDARY_SUPPLEMENTARY | 1 |
| INCOMPLETE_PATH | 1 |
| INCLUDED | 1 |

### FINAL decisions: 17 KEEP, 41 DROP

#### Kept

- pred[0] 17190254-17433495 cov=11.91 - exons=19
- pred[1] 17190254-17433819 cov=8.32 - exons=20
- pred[2] 17190254-17371471 cov=3.02 - exons=11
- pred[4] 17191767-17433819 cov=127.01 - exons=19
- pred[5] 17191767-17431270 cov=1.97 - exons=19
- pred[24] 17433800-17451950 cov=49.38 + exons=7
- pred[25] 17433800-17451950 cov=6.02 + exons=8
- pred[26] 17433800-17451950 cov=2.09 + exons=8
- pred[27] 17433800-17451950 cov=2.05 + exons=8
- pred[29] 17451304-17465073 cov=10.40 - exons=14
- pred[31] 17452179-17465073 cov=719.16 - exons=15
- pred[37] 17452251-17456446 cov=72.65 + exons=2
- pred[43] 17465144-17510532 cov=39.07 + exons=23
- pred[44] 17465144-17510532 cov=3.58 + exons=24
- pred[47] 17465144-17510532 cov=1.67 + exons=24
- pred[48] 17466198-17510532 cov=14.56 + exons=23
- pred[49] 17466198-17510532 cov=3.88 + exons=24

#### Dropped (sample first 10)

- pred[3] 17190254-17433495 cov=1.04 - exons=18
- pred[6] 17191767-17433819 cov=0.91 - exons=20
- pred[7] 17194325-17433819 cov=0.89 - exons=20
- pred[8] 17215401-17433819 cov=1.82 - exons=14
- pred[9] 17215401-17433495 cov=0.51 - exons=14
- pred[10] 17301999-17433819 cov=8.90 - exons=13
- pred[11] 17301999-17433819 cov=2.80 - exons=14
- pred[12] 17301999-17433495 cov=0.93 - exons=13
- pred[13] 17301999-17433819 cov=0.88 - exons=14
- pred[14] 17323597-17433819 cov=11.94 - exons=12

## 110888244-110930143(+)

### StringTie kill-reason histogram

| reason | count |
|---|---|
| GOOD | 23 |
| BAD_MM_NEG | 23 |
| EXISTING_BUNDLE | 16 |
| back_fail | 16 |
| NO_MAXP | 9 |
| EXCLUDE_NO_COV | 7 |
| RETAINED_INTRON | 7 |
| BAD_NO_STRAND | 4 |
| SECONDARY_SUPPLEMENTARY | 1 |
| DEMOTED_LEFT | 1 |
| INCLUDED | 1 |

### FINAL decisions: 13 KEEP, 12 DROP

#### Kept

- pred[0] 110888244-110930143 cov=12.84 + exons=16
- pred[1] 110888244-110930143 cov=2.58 + exons=15
- pred[2] 110888244-110930143 cov=2.26 + exons=17
- pred[9] 110903774-110930143 cov=10.00 + exons=13
- pred[10] 110903774-110930143 cov=4.54 + exons=11
- pred[11] 110903955-110930143 cov=2.06 + exons=13
- pred[13] 110904260-110930143 cov=35.21 + exons=13
- pred[14] 110904260-110930143 cov=30.59 + exons=11
- pred[15] 110904260-110930143 cov=14.26 + exons=12
- pred[16] 110904260-110930143 cov=8.40 + exons=11
- pred[17] 110904260-110930143 cov=6.97 + exons=14
- pred[18] 110904260-110930143 cov=5.00 + exons=14
- pred[19] 110904260-110930143 cov=2.86 + exons=12

#### Dropped (sample first 10)

- pred[3] 110888244-110930143 cov=1.30 + exons=16
- pred[4] 110888244-110930143 cov=1.26 + exons=17
- pred[5] 110888244-110930143 cov=1.16 + exons=15
- pred[6] 110888244-110901897 cov=0.94 + exons=2
- pred[7] 110888244-110890490 cov=0.94 + exons=2
- pred[8] 110888642-110930143 cov=1.19 + exons=15
- pred[12] 110903955-110930143 cov=1.09 + exons=13
- pred[20] 110904260-110930143 cov=2.85 + exons=11
- pred[21] 110904329-110929497 cov=1.05 + exons=9
- pred[22] 110919128-110930143 cov=2.20 + exons=5

## 57118036-57129805(+)

### StringTie kill-reason histogram

| reason | count |
|---|---|
| GOOD | 39 |
| EXISTING_BUNDLE | 39 |
| BAD_MM_NEG | 9 |
| zero_flux | 7 |
| RETAINED_INTRON | 4 |
| SECONDARY_SUPPLEMENTARY | 3 |
| SINGLETON_ARTIFACT | 3 |
| INCLUDED | 2 |
| BAD_NO_STRAND | 1 |
| DEMOTED_RIGHT | 1 |
| NO_REACH | 1 |
| back_fail | 1 |
| NO_MAXC | 1 |
| fwd_fail | 1 |
| READTHR_PAIRWISE | 1 |

### FINAL decisions: 9 KEEP, 7 DROP

#### Kept

- pred[0] 57112612-57118037 cov=2.48 - exons=2
- pred[1] 57112612-57118037 cov=1.30 - exons=2
- pred[2] 57118036-57129805 cov=66.07 + exons=7
- pred[3] 57118036-57129805 cov=27.94 + exons=7
- pred[4] 57118036-57129805 cov=4.99 + exons=6
- pred[5] 57118036-57129805 cov=4.80 + exons=7
- pred[6] 57118036-57129805 cov=2.86 + exons=6
- pred[9] 57118036-57129805 cov=1.92 + exons=7
- pred[12] 57133274-57238547 cov=6.98 - exons=30

#### Dropped (sample first 10)

- pred[7] 57118036-57129805 cov=2.66 + exons=5
- pred[8] 57118036-57129805 cov=1.95 + exons=5
- pred[10] 57118036-57129805 cov=1.81 + exons=4
- pred[11] 57118036-57129805 cov=0.94 + exons=5
- pred[13] 57148337-57238547 cov=1.24 - exons=24
- pred[14] 57152483-57238547 cov=1.29 - exons=23
- pred[15] 57168991-57238547 cov=0.99 - exons=9

## 20532687-20568095(-)

### StringTie kill-reason histogram

| reason | count |
|---|---|
| GOOD | 257 |
| EXISTING_BUNDLE | 243 |
| BAD_MM_NEG | 215 |
| RETAINED_INTRON | 44 |
| fwd_fail | 40 |
| zero_flux | 39 |
| NO_MAXP | 33 |
| back_fail | 33 |
| NO_MAXC | 31 |
| BAD_NO_STRAND | 29 |
| NO_CHILDPAT_REACH | 9 |
| DEMOTED_RIGHT | 7 |
| DEMOTED_LEFT | 6 |
| INCOMPLETE_PATH | 4 |
| INCLUDED | 4 |
| READTHR_PAIRWISE | 3 |
| STRAND_MISMATCH | 1 |
| INCLUDED_REVERSE | 1 |

### FINAL decisions: 55 KEEP, 66 DROP

#### Kept

- pred[0] 20128504-20153110 cov=35.02 + exons=22
- pred[2] 20161425-20164285 cov=1.94 - exons=5
- pred[3] 20204465-20223147 cov=35.33 + exons=2
- pred[5] 20232547-20239074 cov=7.96 - exons=4
- pred[6] 20232547-20243244 cov=5.79 - exons=5
- pred[7] 20232547-20243244 cov=4.66 - exons=6
- pred[8] 20232547-20239074 cov=3.69 - exons=5
- pred[9] 20232547-20239074 cov=2.10 - exons=5
- pred[10] 20232547-20239074 cov=1.47 - exons=4
- pred[11] 20233343-20239074 cov=2.49 - exons=4
- pred[12] 20233343-20239074 cov=1.57 - exons=3
- pred[14] 20238071-20241114 cov=3.24 + exons=2
- pred[15] 20243230-20266946 cov=2.00 + exons=22
- pred[17] 20264902-20305865 cov=1.80 + exons=18
- pred[19] 20303936-20309518 cov=3.89 - exons=2
- pred[20] 20324219-20335718 cov=81.07 - exons=2
- pred[21] 20324219-20335718 cov=5.15 - exons=2
- pred[22] 20331572-20335718 cov=2.67 - exons=2
- pred[23] 20335679-20354913 cov=2.75 + exons=2
- pred[25] 20429749-20533661 cov=281.86 + exons=18
- pred[26] 20429749-20533661 cov=74.77 + exons=19
- pred[27] 20429749-20533661 cov=72.45 + exons=19
- pred[28] 20429749-20533661 cov=23.76 + exons=20
- pred[29] 20429749-20533661 cov=17.32 + exons=19
- pred[30] 20429749-20533661 cov=12.56 + exons=20
- pred[31] 20429749-20533661 cov=5.40 + exons=20
- pred[44] 20532687-20568095 cov=246.55 - exons=17
- pred[72] 20580926-20588070 cov=541.36 - exons=6
- pred[76] 20590817-20617650 cov=3.28 - exons=12
- pred[77] 20590817-20617650 cov=1.50 - exons=11
- pred[78] 20591174-20617650 cov=90.68 - exons=12
- pred[79] 20591174-20617650 cov=10.16 - exons=11
- pred[81] 20591174-20617650 cov=2.06 - exons=11
- pred[82] 20591174-20617650 cov=2.04 - exons=11
- pred[85] 20619868-20635436 cov=270.22 - exons=7
- pred[86] 20619868-20635436 cov=21.20 - exons=6
- pred[87] 20619868-20635436 cov=12.59 - exons=7
- pred[88] 20619868-20635436 cov=7.63 - exons=8
- pred[90] 20619868-20635436 cov=3.47 - exons=5
- pred[95] 20622915-20635436 cov=3.80 - exons=7
- pred[101] 20634427-20639429 cov=33.59 + exons=3
- pred[102] 20634427-20639429 cov=29.64 + exons=2
- pred[106] 20635318-20639429 cov=4.53 + exons=3
- pred[107] 20635318-20639429 cov=2.08 + exons=3
- pred[108] 20648002-20649293 cov=5.50 - exons=2
- pred[109] 20649910-20733035 cov=1.40 - exons=17
- pred[110] 20649910-20708598 cov=1.15 - exons=12
- pred[111] 20660267-20708677 cov=1.43 + exons=14
- pred[112] 20660267-20708677 cov=1.29 + exons=14
- pred[113] 20660267-20708677 cov=1.20 + exons=16
- pred[114] 20662731-20733035 cov=2.59 - exons=18
- pred[115] 20662731-20733035 cov=2.36 - exons=15
- pred[116] 20662731-20683510 cov=1.14 - exons=7
- pred[117] 20663392-20708677 cov=1.20 + exons=14
- pred[118] 20681421-20708677 cov=1.32 + exons=11

#### Dropped (sample first 10)

- pred[1] 20131980-20137533 cov=1.38 + exons=7
- pred[4] 20232349-20239074 cov=1.38 - exons=4
- pred[13] 20233343-20239074 cov=0.85 - exons=3
- pred[16] 20253815-20266946 cov=1.01 + exons=16
- pred[18] 20298782-20305865 cov=3.74 + exons=4
- pred[24] 20335679-20360050 cov=0.99 + exons=2
- pred[32] 20429749-20533661 cov=4.34 + exons=18
- pred[33] 20429749-20533661 cov=3.80 + exons=17
- pred[34] 20429749-20533661 cov=2.84 + exons=18
- pred[35] 20429749-20486893 cov=2.16 + exons=5

## 40800807-40875796(-)

### StringTie kill-reason histogram

| reason | count |
|---|---|
| BAD_MM_NEG | 167 |
| GOOD | 132 |
| EXISTING_BUNDLE | 130 |
| fwd_fail | 45 |
| NO_MAXC | 37 |
| BAD_NO_STRAND | 21 |
| RETAINED_INTRON | 16 |
| zero_flux | 12 |
| EXCLUDE_NO_COV | 8 |
| back_fail | 8 |
| NO_MAXP | 7 |
| READTHR_PAIRWISE | 4 |
| INCOMPLETE_PATH | 3 |
| DEMOTED_LEFT | 2 |
| DEMOTED_RIGHT | 2 |
| NO_REACH | 1 |
| INCLUDED | 1 |

### FINAL decisions: 29 KEEP, 36 DROP

#### Kept

- pred[0] 40800807-40829994 cov=8.91 - exons=7
- pred[1] 40800807-40829737 cov=5.94 - exons=7
- pred[2] 40800807-40825304 cov=1.05 - exons=5
- pred[4] 40830314-40875063 cov=377.18 - exons=18
- pred[5] 40830314-40875063 cov=165.83 - exons=19
- pred[6] 40830314-40875063 cov=132.80 - exons=18
- pred[7] 40830314-40875063 cov=56.94 - exons=19
- pred[8] 40830314-40875796 cov=46.70 - exons=18
- pred[28] 40875306-40880088 cov=139.31 + exons=3
- pred[29] 40875306-40893796 cov=41.13 + exons=5
- pred[30] 40875306-40880088 cov=19.94 + exons=3
- pred[33] 40875540-40880088 cov=10.21 + exons=3
- pred[34] 40875540-40880088 cov=5.02 + exons=3
- pred[35] 40875540-40880088 cov=3.01 + exons=3
- pred[38] 40912981-40964963 cov=4.75 + exons=12
- pred[39] 40912981-40964963 cov=2.60 + exons=12
- pred[42] 40967471-41010547 cov=42.30 - exons=12
- pred[43] 40967471-41010547 cov=3.83 - exons=11
- pred[45] 40967471-40978735 cov=1.76 - exons=5
- pred[46] 40967471-41010547 cov=1.46 - exons=13
- pred[47] 40967471-41010547 cov=1.41 - exons=10
- pred[48] 40967471-41010547 cov=1.39 - exons=9
- pred[50] 41013194-41030428 cov=7.34 - exons=4
- pred[51] 41037715-41060793 cov=87.98 - exons=9
- pred[52] 41037715-41060793 cov=1.48 - exons=8
- pred[54] 41070990-41082101 cov=305.53 - exons=9
- pred[59] 41118350-41155377 cov=5.89 - exons=15
- pred[61] 41119154-41155377 cov=260.03 - exons=16
- pred[62] 41119154-41155377 cov=15.55 - exons=15

#### Dropped (sample first 10)

- pred[3] 40804630-40829737 cov=0.95 - exons=5
- pred[9] 40830314-40875063 cov=4.81 - exons=17
- pred[10] 40830314-40875063 cov=4.78 - exons=17
- pred[11] 40830314-40875063 cov=3.65 - exons=18
- pred[12] 40830314-40875063 cov=2.55 - exons=18
- pred[13] 40830314-40875063 cov=2.42 - exons=20
- pred[14] 40830314-40874805 cov=2.35 - exons=18
- pred[15] 40830314-40875063 cov=2.22 - exons=19
- pred[16] 40830314-40875063 cov=1.93 - exons=17
- pred[17] 40830314-40875063 cov=1.91 - exons=19

