# root-hair-patterning
Mills (2025) "A theoretical investigation provides an integrated understanding of the complex regulatory network controlling Arabidopsis root hair patterning"

# run-model
to use the model developed in Mills (2025)

# parameter-space-analysis
the code used to find the successful parameter sets in Mills (2025)

# paper-figures
the code used to produce the figures in Mills (2025)

# SuccessData.mat
a string containing the model description, parameter set number, signal strength, and noise in the initial conditions for each successful parameter set. This contains the columns:
nonlinearity | CPC movement type | wer repressor | EGL3 diffusion | parameter set number | signal strength | noise
which relate to the notation used in the code throughout:
nonlinearity [linear, myb, complex, cpc, mybcpc, mybcomplex, cpccomplex, all] - can be none; "linear", a single nonlinearity; "myb"; "complex"; or "cpc", a double nonlinearity; "mybcpc"; "mybcomplex"; or "cpccomplex", or all three; "all".
movement [DM, Diff] - (of CPC), eiter by both diffusion and directed movement, "DM", or diffusion only "Diff".
werR [complex, protein, CI] - the mechanism of wer repression which can be via the CPC complex, "complex", CPC protein, "protein", or competitive inhibition, "CI".
DEGL3 [Y, N] - can EGL3 diffuse "Y" or "N".
parameter set number - the index of parspace.mat giving the successful parameter set.
signal - strength of cortical signal.
noise - maximum (C0max) of the interaval of random initial conditions.
There are 129 combinations of successful models and parameter sets (rows), the 10 successful models highlighted in Mills (2025) appear from rows 116-125.

# parspace.mat
a matrix of the 20,000 random sets of the 48 parameter values. 
parameters in the matrix are in the order such that:
column | parameter
1      | k11
2      | k12
3      | k13
4      | k14
5      | k15
6      | k16
7      | k17
8      | k18
9      | k21
10     | k22
11     | k23
12     | k24
13     | DG
14     | DE
15     | DC
16     | k3 (D~C)
17     | bg
18     | be
19     | bw
20     | dG
21     | dE
22     | dC
23     | dW
24     | dM
25     | dg
26     | de
27     | dc
28     | dw
29     | dm
30     | qG
31     | qE
32     | qC
33     | qW
34     | qM
35     | pc
36     | pm
37     | rg
38     | re
39     | rw11
40     | rw12
41     | rw21
42     | rw22
43     | n1
44     | n2
45     | n3
46     | n4
47     | KM
48     | KC
