Alignment Steps

1. Place in expresult\experiment\ folder the following matlab files.

(Alignuser.m)
Alignuser2.m

Readoneframe.m  This creates a bmp of one frame (a good one usually in the
middle of a run) so one can specific region of interest pixels for
analysis with MSPAINT.  Use curser to identify pixel boarders and enter
these into alignuser.  Keep region of interest within ccd part.

Specify file, pixels and folder in Alignuser.m

Run alignuser.m and pray.

Alignuser.m generates ***.2.txt, and ***.3.txt

2. Place \goodstim\ subfolder in expresult\experiment\ folder.
OIAshift2goodstim3.m
OIAshift2goodstim2.m

Copy ***.2.txt, and ***.3.txt files to goodstim folder

Alignuser.m generated ***.2.txt, and ***.3.txt
Copy ***.2.txt, and ***.3.txt filenames into OIAshift2goodstim2.m
and OIAshift2goodstim3.m, repectively.

After running “OIAshift2goodstim3.m and OIAshift2goodstim2.m” you will get
“***goodstim.2.txt” and “***goodstim.3.txt”

3.

Combine these two goodstim, save results in “goodstim.txt” copy to
results\run folder
Copy “***.2.txt” into Excel

4. , and copy first 7 columns into file “shiftinput.2.txt” save at
result\run folder.

5. run suninuser. remember flagalign=41