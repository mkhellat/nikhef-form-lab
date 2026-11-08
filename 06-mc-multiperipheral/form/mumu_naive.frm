* Naive ("standard fashion") squared amplitude for gamma*gamma* -> mu- mu+
* in FORM (report Section 3 / 5). Same N / Nbar split as the Python
* scripts, evaluated with Trace4.
*
* Does not by itself prove equality with the paper's M_ij form; that needs
* a transcribed (III.11). This file is the FORM reimplementation of the
* naive route, for independent checking and as input to later FORM work.
*
* Run: form mumu_naive.frm

Off statistics;
Vectors e1,e2,k1,k2,p6,p7;
Symbols m,D1,D2;

Local N =
  g_(1,e1,k1,e2)/D1 + m*g_(1,e1,e2)/D1
  + g_(1,e2,k2,e1)/D2 + m*g_(1,e2,e1)/D2;

Local Nbar =
  g_(1,e2,k1,e1)/D1 + m*g_(1,e2,e1)/D1
  + g_(1,e1,k2,e2)/D2 + m*g_(1,e1,e2)/D2;

Local S = (g_(1,p6)+gi_(1)*m) * N * (g_(1,p7)-gi_(1)*m) * Nbar;
Trace4,1;
id p6.p6 = m^2;
id p7.p7 = m^2;
.sort

Print +s;
.end
