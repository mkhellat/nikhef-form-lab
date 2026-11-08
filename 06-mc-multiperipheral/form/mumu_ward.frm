* Squared Ward identity for gamma*gamma* -> mu- mu+
* (report/mumu_cancellation.tex Section 5).
*
* Computes
*   Sg = Tr[(p6slash+m) Ng (p7slash-m) Nbar]
* with e1 -> q1 only on Ng (diagram1: q1=p6-k1; diagram2: q1=p7+k2),
* Nbar kept physical, then clears D1^2 D2^2, substitutes
*   D_i = k_i.k_i - m^2,  p7 = p6-k1-k2,  p6.p6 = m^2.
* Pass criterion: Sg = 0.
*
* Run:  form mumu_ward.frm
* Or:   make -C form check

Off statistics;
Vectors e1,e2,k1,k2,p6,p7;
Symbols m,D1,D2;

Local Ng =
  (g_(1,p6,k1,e2) - g_(1,k1,k1,e2))/D1
  + m*(g_(1,p6,e2) - g_(1,k1,e2))/D1
  + (g_(1,e2,k2,p7) + g_(1,e2,k2,k2))/D2
  + m*(g_(1,e2,p7) + g_(1,e2,k2))/D2;

Local Nbar =
  g_(1,e2,k1,e1)/D1 + m*g_(1,e2,e1)/D1
  + g_(1,e1,k2,e2)/D2 + m*g_(1,e1,e2)/D2;

Local Sg = (g_(1,p6)+gi_(1)*m) * Ng * (g_(1,p7)-gi_(1)*m) * Nbar;
Trace4,1;
id p6.p6 = m^2;
id p7.p7 = m^2;
.sort

Multiply D1^2*D2^2;
id D1 = k1.k1 - m^2;
id D2 = k2.k2 - m^2;
.sort

* Momentum conservation: q1+q2=p6+p7 with q1=p6-k1, q2=p6-k2
* implies p7 = p6 - k1 - k2.
id p7 = p6 - k1 - k2;
id p6.p6 = m^2;
.sort

Print;
.end
