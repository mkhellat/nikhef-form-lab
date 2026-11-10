* Appendix C reduction identities for M12/M21.
* Vermaseren Nucl.Phys.B229 (1983) App. C / NIKHEF-H/82-15.
*
* Verifies (C.4)/(C.7), (C.8), (C.9), (C.13), (C.14) after Levi-Civita
* contraction in Dimension 4. Building-block notation:
*   Tr[i·j]     = e_(qi,ei,m1,m2)*e_(qj,ej,m1,m2)
*   Q·1·2·R     = e_(q1,e1,Q,m)*e_(q2,e2,R,m)   [see (C.4) note]
*   Q·1·2·1·2·R = four-epsilon chain of (C.2)
*   X·1·1·Y     = e_(q1,e1,X,m)*e_(q1,e1,m,Y)
*   P           = q1+q2  (substituted before (C.13)/(C.14))
*
* FORM vs paper coefficient conventions (locked by Diff* = 0):
*   (C.4)  Q12R uses e_(…,R,m) so the paper +1/2 holds; the alternate
*         e_(…,m,R) yields −1/2 instead (index-order).
*   (C.13) Tr[1·1] term enters with +1/2 in FORM (paper scan −1/2);
*         the 1/4 (ε^{q1 e1 αβ} ε_{Q q2 αβ})^2 term is required.
*   (C.14) Tr[1·1] brace flips sign vs the paper scan.
*
* Pass: DiffC4=DiffC8=DiffC9=DiffC13=DiffC14 = 0
* Run: form mumu_schouten_c2_c14.frm  or  make -C form appC

Off statistics;
Dimension 4;
Indices m1,m2,m3,m4;
Vectors q1,e1,q2,e2,Q,R,P;

* ----- Building blocks -----
Local Tr12 = e_(q1,e1,m1,m2) * e_(q2,e2,m1,m2);
Local Tr11 = e_(q1,e1,m1,m2) * e_(q1,e1,m1,m2);
Local Tr22 = e_(q2,e2,m1,m2) * e_(q2,e2,m1,m2);
Local EpsFull = e_(q1,e1,q2,e2) * e_(q1,e1,q2,e2);
* Paper-matching Q·1·2·R (free index after R → +1/2 in (C.4)).
Local Q12R = e_(q1,e1,Q,m1) * e_(q2,e2,R,m1);
Local Q1212R =
  e_(q1,e1,Q,m1) * e_(q2,e2,m1,m2) * e_(q1,e1,m2,m3) * e_(q2,e2,m3,R);
Local Tr1212 =
  e_(q1,e1,m1,m2) * e_(q2,e2,m2,m3) * e_(q1,e1,m3,m4) * e_(q2,e2,m4,m1);
Local Tr1221 =
  e_(q1,e1,m1,m2) * e_(q2,e2,m2,m3) * e_(q2,e2,m3,m4) * e_(q1,e1,m4,m1);
Local EpsQq2 = e_(q1,e1,Q,q2) * e_(q1,e1,Q,q2);
Local BilineQq2 = e_(q1,e1,m1,m2) * e_(Q,q2,m1,m2);
Local P11P = e_(q1,e1,P,m1) * e_(q1,e1,m1,P);
Local Q11Q = e_(q1,e1,Q,m1) * e_(q1,e1,m1,Q);
Local Q11P = e_(q1,e1,Q,m1) * e_(q1,e1,m1,P);
contract;
.sort

* Square the (ε^{q1 e1 αβ} ε_{Q q2 αβ}) biline after the first contraction.
Local BilineSq = BilineQq2 * BilineQq2;
contract;
.sort

* ----- (C.4)/(C.7): Q·1·2·1·2·R = (1/2) Tr[1·2] Q·1·2·R -----
Local DiffC4 = Q1212R - (1/2) * Tr12 * Q12R;
.sort

* ----- (C.8): Tr[1·2·1·2] = (1/2) Tr[1·2]^2 -----
Local DiffC8 = Tr1212 - (1/2) * Tr12 * Tr12;
.sort

* ----- (C.9): Tr[1·2·2·1] = (1/4)Tr11 Tr22 + (1/4)Tr12^2 − εε -----
Local DiffC9 = Tr1221 - (1/4)*Tr11*Tr22 - (1/4)*Tr12*Tr12 + EpsFull;
.sort

* ----- (C.13)/(C.14): need P = q1+q2 -----
id P = q1 + q2;
.sort

* (C.13) FORM-locked:
*   ε^{q1 e1 Q q2}^2
*     = (1/4) BilineSq + (1/2) Tr11 ((Q·q2)^2 − Q^2 q2^2)
*       − Q^2 P·1·1·P − q2^2 Q·1·1·Q + 2 (Q·q2) Q·1·1·P
Local DiffC13 =
  EpsQq2
  - (1/4)*BilineSq
  - (1/2)*Tr11*((Q.q2)^2 - Q.Q*q2.q2)
  + Q.Q*P11P
  + q2.q2*Q11Q
  - 2*(Q.q2)*Q11P;
.sort

* (C.14) FORM-locked (Tr brace sign flipped vs paper scan):
*   2 (Q·q1)(q1·q2) Q·1·1·P
*     = (Q·q1)^2 P·1·1·P + (q1·q2)^2 Q·1·1·Q
*       − Tr11 { (Q·q1)(Q·q2)(q1·q2) − (1/2)(Q·q1)^2 q2^2 − (1/2) Q^2 (q1·q2)^2 }
*       − (1/4) q1^2 BilineSq
Local DiffC14 =
  2*(Q.q1)*(q1.q2)*Q11P
  - (Q.q1)^2*P11P
  - (q1.q2)^2*Q11Q
  + Tr11*(
      (Q.q1)*(Q.q2)*(q1.q2)
      - (1/2)*(Q.q1)^2*q2.q2
      - (1/2)*Q.Q*(q1.q2)^2
    )
  + (1/4)*q1.q1*BilineSq;
.sort

Print;
.end
