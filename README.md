Code written to solve density classification problem in C.

Usage: ./project <file.txt

File Format:
Automaton Sequence Size
Update Rule Number
Initial Automaton Sequence
Times to Run Simulation
Cell Position Index, Starting Sequence to Count (1)
Cell Position Index, Starting Sequence to Count (2)

Example File:
31
150
**..........*..*..*..........*.
10
0,5
30,10

Output Format:
Stage 0:
Size
Update Rule
Initial Sequence

Stage 1:
Sequences Over Time
On Off Cell Count

Stage 2:
Sequences After Using Rule 184
Sequences After Using Rule 232
On Off Cell Count
Density Classifcation at End of Simulation

Example Output:
==STAGE 0============================
SIZE: 31
RULE: 30
-------------------------------------
 000 001 010 011 100 101 110 111
  0   1   1   1   1   0   0   0 
-------------------------------------
   0: **..........*..*..*..........*.
==STAGE 1============================
   0: **..........*..*..*..........*.
   1: *.*........*********........**.
   2: *.**......**........*......**..
   3: *.*.*....**.*......***....**.**
   4: ..*.**..**..**....**..*..**..*.
   5: .**.*.***.***.*..**.******.****
   6: .*..*.*...*...****..*......*...
   7: *****.**.***.**...****....***..
   8: *.....*..*...*.*.**...*..**..**
   9: .*...******.**.*.*.*.*****.***.
  10: ***.**......*..*.*.*.*.....*..*
#ON=3 #OFF=3 CELL#0 START@5
==STAGE 2============================
RULE: 184; STEPS: 14.
-------------------------------------
  10: ***.**......*..*.*.*.*.....*..*
  11: **.**.*......*..*.*.*.*.....*.*
  12: *.**.*.*......*..*.*.*.*.....**
  13: .**.*.*.*......*..*.*.*.*....**
  14: **.*.*.*.*......*..*.*.*.*...*.
  15: *.*.*.*.*.*......*..*.*.*.*...*
  16: .*.*.*.*.*.*......*..*.*.*.*..*
  17: *.*.*.*.*.*.*......*..*.*.*.*..
  18: .*.*.*.*.*.*.*......*..*.*.*.*.
  19: ..*.*.*.*.*.*.*......*..*.*.*.*
  20: *..*.*.*.*.*.*.*......*..*.*.*.
  21: .*..*.*.*.*.*.*.*......*..*.*.*
  22: *.*..*.*.*.*.*.*.*......*..*.*.
  23: .*.*..*.*.*.*.*.*.*......*..*.*
  24: *.*.*..*.*.*.*.*.*.*......*..*.
-------------------------------------
RULE: 232; STEPS: 15.
-------------------------------------
  24: *.*.*..*.*.*.*.*.*.*......*..*.
  25: .*.*....*.*.*.*.*.*...........*
  26: *.*......*.*.*.*.*.............
  27: .*........*.*.*.*..............
  28: ...........*.*.*...............
  29: ............*.*................
  30: .............*.................
  31: ...............................
  32: ...............................
  33: ...............................
  34: ...............................
  35: ...............................
  36: ...............................
  37: ...............................
  38: ...............................
  39: ...............................
-------------------------------------
#ON=10 #OFF=20 CELL#30 START@10
-------------------------------------
  10: ***.**......*..*.*.*.*.....*..*
AT T=10: #ON/#CELLS < 1/2
==THE END============================


