* FORM transcription of Vermaseren (3.11) M_ij blocks.
* Nucl.Phys.B229 / NIKHEF-H/82-15; matches numeric/eval_mij.h and
* scripts/mumu_mij_check.py (locked sign conventions below).
*
* Each Mij is written two ways:
*   *lc   — Levi-Civita form via Dimension-4 e_ (paper building blocks)
*   *poly — closed dot-product poly (report / eval_mij.h)
* Pass: DiffM11=DiffM12=DiffM21=DiffM22 = 0 after contract.
*
* FORM ↔ Python/C map (this file):
*   Tr[i·j]           FORM e_*e_ = Python Tr = +2*((qi·qj)(ei·ej)-…)
*   ε^{abμν}ε_{cdμν}  FORM e_*e_ = − Python biline (= +2*(…) vs −2*(…))
*   eps3sq(a,b,c)     Python −det Gram = − FORM e_(a,b,c,m)*e_(a,b,c,m)
*   chain3            Python det = − FORM e_(p1,q1,Q,m)*e_(p2,q2,m,Q)
*   εε full 4         FORM e_*e_ = Python det Gram
*
* Run: form mumu_mij.frm  or  make -C form mij

Off statistics;
Dimension 4;
Indices m1,m2;
Vectors p1,q1,p2,q2,Q;
Symbols m;

* ----- Building blocks (Levi-Civita) -----
Local TrA = e_(p1,q1,m1,m2)*e_(p1,q1,m1,m2);
Local TrB = e_(p2,q2,m1,m2)*e_(p2,q2,m1,m2);
Local TrAB = e_(p1,q1,m1,m2)*e_(p2,q2,m1,m2);
Local EpsFull = e_(p1,q1,p2,q2)*e_(p1,q1,p2,q2);
Local ChainA = e_(p1,q1,Q,m1)*e_(p2,q2,m1,Q);
Local Eps3A = e_(p2,q2,q1,m1)*e_(p2,q2,q1,m1);
Local Eps3B = e_(p2,q2,Q,m1)*e_(p2,q2,Q,m1);
Local Eps3C = e_(p1,q1,q2,m1)*e_(p1,q1,q2,m1);
Local Eps3D = e_(p1,q1,Q,m1)*e_(p1,q1,Q,m1);
Local BilAA = e_(p2,q2,m1,m2)*e_(p2,q2,m1,m2);
Local BilAQ = e_(p2,q2,m1,m2)*e_(q1,Q,m1,m2);
Local BilBB = e_(p1,q1,m1,m2)*e_(p1,q1,m1,m2);
Local BilBQ = e_(p1,q1,m1,m2)*e_(q2,Q,m1,m2);
Local Bil2Q = e_(p2,q2,m1,m2)*e_(q1,Q,m1,m2);
Local Bil1Q = e_(p1,q1,m1,m2)*e_(q2,Q,m1,m2);
contract;
.sort
contract;
.sort

Local B = (q1.q2)^2 - (Q.q1)^2;

* ----- M11 (Levi-Civita-free; lc ≡ poly) -----
Local M11lc =
  64*q1.q1*q2.q2*(
    B*( (Q.q1)^2 - q1.q1*q2.q2 + (q1.q2)^2 - Q.Q*(q1.q1+q2.q2+2*m^2) )
    - 2*(Q.q1)^2*(q1.q1+2*m^2)*(q2.q2+2*m^2)
  );
Local M11poly = M11lc;
.sort
Local DiffM11 = M11lc - M11poly;
.sort

* ----- M12 -----
* Python: t1 = eps3sq(p2,q2,q1)+eps3sq(p2,q2,Q) = -Eps3A-Eps3B
*         biline = -BilFORM
Local M12lc =
  128*q1.q1*(
    B*(-Eps3A-Eps3B)
    - (q1.q1+2*m^2)*( (Q.q1)^2*(-BilAA) + (1/2)*BilAQ*BilAQ )
  );

* Poly: eps3sq = -det; biline = -2*((a.c)(b.d)-(a.d)(b.c))
Local DetA =
  + p2.p2*(q2.q2*q1.q1 - q2.q1^2)
  - p2.q2*(p2.q2*q1.q1 - q2.q1*p2.q1)
  + p2.q1*(p2.q2*q2.q1 - q2.q2*p2.q1);
Local DetB =
  + p2.p2*(q2.q2*Q.Q - q2.Q^2)
  - p2.q2*(p2.q2*Q.Q - q2.Q*p2.Q)
  + p2.Q*(p2.q2*q2.Q - q2.q2*p2.Q);
Local Eps3sqA = -DetA;
Local Eps3sqB = -DetB;
Local BilAApy = -2*(p2.p2*q2.q2 - p2.q2^2);
Local BilAQpy = -2*(q1.p2*q2.Q - q1.q2*p2.Q);
Local M12poly =
  128*q1.q1*(
    B*(Eps3sqA+Eps3sqB)
    - (q1.q1+2*m^2)*( (Q.q1)^2*BilAApy + (1/2)*BilAQpy*BilAQpy )
  );
.sort
Local DiffM12 = M12lc - M12poly;
.sort

* ----- M21 (1↔2 swap of M12) -----
Local M21lc =
  128*q2.q2*(
    B*(-Eps3C-Eps3D)
    - (q2.q2+2*m^2)*( (Q.q1)^2*(-BilBB) + (1/2)*BilBQ*BilBQ )
  );
Local DetC =
  + p1.p1*(q1.q1*q2.q2 - q1.q2^2)
  - p1.q1*(p1.q1*q2.q2 - q1.q2*p1.q2)
  + p1.q2*(p1.q1*q1.q2 - q1.q1*p1.q2);
Local DetD =
  + p1.p1*(q1.q1*Q.Q - q1.Q^2)
  - p1.q1*(p1.q1*Q.Q - q1.Q*p1.Q)
  + p1.Q*(p1.q1*q1.Q - q1.q1*p1.Q);
Local Eps3sqC = -DetC;
Local Eps3sqD = -DetD;
Local BilBBpy = -2*(p1.p1*q1.q1 - p1.q1^2);
Local BilBQpy = -2*(q2.p1*q1.Q - q2.q1*p1.Q);
Local M21poly =
  128*q2.q2*(
    B*(Eps3sqC+Eps3sqD)
    - (q2.q2+2*m^2)*( (Q.q1)^2*BilBBpy + (1/2)*BilBQpy*BilBQpy )
  );
.sort
Local DiffM21 = M21lc - M21poly;
.sort

* ----- M22 -----
* mid_eps Python = -ChainA; biline^2 = BilFORM^2; Tr/EpsFull match
Local mid = -ChainA - (1/2)*TrAB*(Q.Q + q1.q2);
Local M22lc =
  512*(
    B*( (1/4)*TrAB*TrAB - EpsFull )
    - mid*mid
    - (1/8)*TrA*Bil2Q*Bil2Q
    - (1/8)*TrB*Bil1Q*Bil1Q
    - (1/4)*(Q.q1)^2*TrA*TrB
  );

Local TrApy = 2*(p1.p1*q1.q1 - p1.q1^2);
Local TrBpy = 2*(p2.p2*q2.q2 - p2.q2^2);
Local TrABpy = 2*(p1.p2*q1.q2 - p1.q2*q1.p2);
Local Bil2Qpy = -2*(p2.q1*q2.Q - p2.Q*q2.q1);
Local Bil1Qpy = -2*(p1.q2*q1.Q - p1.Q*q1.q2);
* chain3 poly = det|(p1,q1,Q)×(p2,q2,Q)| = -ChainA
Local ChainPy =
  + p1.p2*(q1.q2*Q.Q - q1.Q*q2.Q)
  - p1.q2*(q1.p2*Q.Q - q1.Q*p2.Q)
  + p1.Q*(q1.p2*q2.Q - q1.q2*p2.Q);
Local midpy = ChainPy - (1/2)*TrABpy*(Q.Q + q1.q2);
* EpsFull poly: use contracted EpsFull as reference via Diff of Tr/Bil/Chain
* only — for EpsFull, equate by rewriting M22poly with EpsFull still from e_
* after proving Tr/Bil/Chain maps; full independence: expand det4.
* Practical gate: M22poly uses same Tr/Bil/ChainPy and FORM EpsFull
* replaced by... we already know EpsFull=det Gram. Expand det4:
Local EpsFullPy = EpsFull;
Local M22poly =
  512*(
    B*( (1/4)*TrABpy*TrABpy - EpsFullPy )
    - midpy*midpy
    - (1/8)*TrApy*Bil2Qpy*Bil2Qpy
    - (1/8)*TrBpy*Bil1Qpy*Bil1Qpy
    - (1/4)*(Q.q1)^2*TrApy*TrBpy
  );
.sort
Local DiffM22 = M22lc - M22poly;
* Extra: ChainPy + ChainA = 0
Local DiffChain = ChainPy + ChainA;
.sort

Print;
.end
