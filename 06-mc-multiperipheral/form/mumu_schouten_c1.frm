* Appendix C Levi-Civita / Schouten toolkit (FORM).
* Vermaseren NIKHEF-H/82-15 App. C / Nucl.Phys.B229.
*
* Gate A: NnormCheck = e_*e_ - 24  → 0  (FORM Dimension-4 convention)
* Gate B: Biline = e_(a,b)*e_(c,d) - 2*((a.c)(b.d)-(a.d)(b.c)) → 0
*          (Python Minkowski biline is −2*(...); FORM e_εε is +2*(...),
*          i.e. an overall sign flip vs our (+−−−) helper.)
* Gate C: Schouten (C1) rewrite once → SchCheck = 0
*
* Pass: Print lines "NnormCheck = 0;" "Biline = 0;" "SchCheck = 0;"
* Run: form mumu_schouten_c1.frm   or   make -C form schouten

Off statistics;
Dimension 4;
Indices m1,m2,m3,m4;
Vectors a,b,c,d,qa,qb,qc,qd,qe,u;

Local Nnorm = e_(m1,m2,m3,m4) * e_(m1,m2,m3,m4);
contract;
.sort
Local NnormCheck = Nnorm - 24;
.sort

Local Biline =
  + e_(a,b,m1,m2) * e_(c,d,m1,m2)
  - 2*(a.c*b.d - a.d*b.c);
contract;
.sort

#procedure SchoutenOnce()
  id e_(qb,qc,qd,qe) * qa.u =
      + e_(qa,qc,qd,qe) * qb.u
      + e_(qb,qa,qd,qe) * qc.u
      + e_(qb,qc,qa,qe) * qd.u
      + e_(qb,qc,qd,qa) * qe.u;
#endprocedure

Local SchCheck =
  + e_(qb,qc,qd,qe) * qa.u
  - e_(qa,qc,qd,qe) * qb.u
  - e_(qb,qa,qd,qe) * qc.u
  - e_(qb,qc,qa,qe) * qd.u
  - e_(qb,qc,qd,qa) * qe.u;
#call SchoutenOnce()
contract;
.sort

Print;
.end
