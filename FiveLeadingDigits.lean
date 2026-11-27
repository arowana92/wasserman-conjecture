import Lean

namespace Wasserman.FiveLeadingDigits

abbrev Row := Nat × Nat × Nat
private def total (r : Row) := r.1+r.2.1+r.2.2
private def upper (r : Row) (j : Nat) : Row :=
  (r.1 / 3^j, r.2.1 / 3^j, r.2.2 / 3^j)
private def addDigit (r d : Row) : Row :=
  (3*r.1+d.1, 3*r.2.1+d.2.1, 3*r.2.2+d.2.2)
private def digits : List Row :=
  [(0,0,0),(0,0,1),(0,0,2),(0,1,0),(0,1,1),(0,2,0),(1,0,0),(1,0,1),(1,1,0),(2,0,0)]
private def caps (r : Row) : Bool :=
  let n := total r + 1
  0 < total r &&
  703687441776640*r.1 < n*422212465065984 &&
  703687441776640*r.2.1 < n*274727834419201 &&
  703687441776640*r.2.2 < n*6747142291457
private def possible : Nat → Row → Bool
  | 0, _ => true
  | n+1, r => digits.any (fun d => caps (addDigit r d) && possible n (addDigit r d))

set_option maxRecDepth 20000
set_option maxHeartbeats 0

private theorem noPath : possible 28 (0,0,0) = false := by decide

private theorem step (r : Row) (j : Nat)
    (fits : total (upper r j) = total r / 3^j)
    (next : total (upper r (j+1)) = total r / 3^(j+1)) :
    ∃ d ∈ digits, upper r j = addDigit (upper r (j+1)) d := by
  let u := upper r j
  have div3 : (upper r (j+1)) = (u.1/3,u.2.1/3,u.2.2/3) := by
    simp [u,upper,Nat.pow_succ,Nat.div_div_eq_div_mul]
  have rel : total (upper r (j+1)) = total u/3 := by
    rw [next, show total u = total r/3^j from fits, Nat.div_div_eq_div_mul, Nat.pow_succ]
  rw [div3] at rel
  change u.1/3 + u.2.1/3 + u.2.2/3 = (u.1+u.2.1+u.2.2)/3 at rel
  refine ⟨(u.1%3,u.2.1%3,u.2.2%3), ?_, ?_⟩
  · simp [digits]; omega
  · rw [div3]; change u = _
    rcases u with ⟨a,b,c⟩
    apply Prod.ext
    · simp [addDigit]; omega
    apply Prod.ext <;> simp [addDigit] <;> omega

private theorem rangeOf (x : Nat) (pos : x>0) : ∃ j, 3^j ≤ x ∧ x < 3^(j+1) := by
  induction x using Nat.strongRecOn with
  | ind x ih =>
    by_cases small : x < 3
    · exact ⟨0,by simp; omega,by simpa using small⟩
    obtain ⟨j,h₁,h₂⟩ := ih (x/3) (Nat.div_lt_self pos (by decide)) (by omega)
    exact ⟨j+1,by rw [Nat.pow_succ]; omega,by rw [Nat.pow_succ]; omega⟩

private theorem trace (r : Row) (fit : ∀ j, total (upper r j) = total r / 3^j)
    (h j : Nat) (bound : ∀ t, j ≤ t → t < j+h → caps (upper r t) = true) :
    possible h (upper r (j+h)) = true := by
  induction h generalizing j with
  | zero => rfl
  | succ h ih =>
    obtain ⟨d,mem,id⟩ := step r (j+h) (fit _) (fit _)
    simp only [possible,List.any_eq_true]
    refine ⟨d,mem,?_⟩
    rw [show j+(h+1)=(j+h)+1 by omega, ← id, bound (j+h) (by omega) (by omega)]
    simpa using ih j (by intro t ht ht'; apply bound t ht; omega)

private theorem oneCap (v n x N m A B : Nat) (_mpos : m>0)
    (hx : v*m ≤ x) (hn : N < (n+1)*m) (ratio : A*x ≤ B*N) (bpos : B>0) :
    A*v < (n+1)*B := by
  have lo := Nat.mul_le_mul_left A hx
  have hi := Nat.mul_lt_mul_of_pos_left hn bpos
  have both : (A*v)*m < ((n+1)*B)*m := by
    calc (A*v)*m = A*(v*m) := Nat.mul_assoc ..
         _ ≤ A*x := lo
         _ ≤ B*N := ratio
         _ < B*((n+1)*m) := hi
         _ = ((n+1)*B)*m := by ac_rfl
  exact Nat.lt_of_mul_lt_mul_right both

private def A := 6442450944
private def B := 4192014075
private def C := 102953221
private def row (w t : Nat) : Row := (A*w,B*w,C*w+t)

/-- For all large positive scalars, either choice of the final unit gives a
    strictly positive ternary factorial layer. -/
theorem largeScalar (w t : Nat) (hw : 65536 ≤ w) (ht : t ≤ 1) :
    ∃ j, total (upper (row w t) j) < total (row w t) / 3^j := by
  let r := row w t
  let N := total r
  have pos : 3^28 ≤ N := by dsimp [N,r,row,total,A,B,C]; omega
  obtain ⟨s,hs,hs'⟩ := rangeOf N (by omega)
  have s28 : 28 ≤ s := by
    by_cases hh : 28 ≤ s
    · exact hh
    have := Nat.pow_le_pow_right (by decide : 0<3) (show s+1≤28 by omega)
    omega
  apply Classical.byContradiction
  intro none
  have fail (j : Nat) : N/3^j ≤ total (upper r j) :=
    Nat.le_of_not_gt (fun h => none ⟨j,h⟩)
  have fit (j : Nat) : total (upper r j) = N/3^j := by
    have h₁ := Nat.div_mul_le_self r.1 (3^j)
    have h₂ := Nat.div_mul_le_self r.2.1 (3^j)
    have h₃ := Nat.div_mul_le_self r.2.2 (3^j)
    have le : total (upper r j)*3^j ≤ N := by
      dsimp [total,upper,N]; rw [Nat.add_mul,Nat.add_mul]; omega
    have posm : 0<3^j := Nat.pow_pos (by decide)
    have other := (Nat.le_div_iff_mul_le posm).2 le
    exact Nat.le_antisymm other (fail j)
  have bound (j : Nat) (hj : j≤s) : caps (upper r j) = true := by
    let m := 3^j
    have hm : m>0 := Nat.pow_pos (by decide)
    have hn : N < (total (upper r j)+1)*m := by
      rw [fit]; change N < (N/m+1)*m
      simpa only [Nat.mul_comm m] using Nat.lt_mul_div_succ N hm
    have small : m ≤ N := Nat.le_trans (Nat.pow_le_pow_right (by decide) hj) hs
    have hp : 0 < total (upper r j) := by rw [fit]; exact Nat.div_pos small hm
    have cap1 := oneCap _ _ _ _ _ 703687441776640 422212465065984 hm
       (Nat.div_mul_le_self r.1 m) hn (by dsimp [r,N,total,row,A,B,C]; omega) (by decide)
    have cap2 := oneCap _ _ _ _ _ 703687441776640 274727834419201 hm
       (Nat.div_mul_le_self r.2.1 m) hn (by dsimp [r,N,total,row,A,B,C]; omega) (by decide)
    have cap3 := oneCap _ _ _ _ _ 703687441776640 6747142291457 hm
       (Nat.div_mul_le_self r.2.2 m) hn (by dsimp [r,N,total,row,A,B,C]; omega) (by decide)
    unfold caps
    simp only [Bool.and_eq_true]
    exact ⟨⟨⟨decide_eq_true hp,decide_eq_true cap1⟩,decide_eq_true cap2⟩,decide_eq_true cap3⟩
  have path := trace r fit 28 (s+1-28) (by intro j lo hi; apply bound; omega)
  have id : upper r (s+1-28+28) = (0,0,0) := by
    rw [show s+1-28+28=s+1 by omega]
    dsimp [upper]; dsimp [N,total] at hs'
    simp [Nat.div_eq_of_lt (show r.1 < 3^(s+1) by omega),
          Nat.div_eq_of_lt (show r.2.1 < 3^(s+1) by omega),
          Nat.div_eq_of_lt (show r.2.2 < 3^(s+1) by omega)]
  rw [id,noPath] at path
  contradiction

private def remSum (w t m : Nat) : Nat :=
  (A*w)%m+(B*w)%m+(C*w+t)%m

private def prefixes (t : Nat) : Nat → List Nat
  | 0 => [0]
  | j+1 =>
    let m := 3^j
    ((prefixes t j).flatMap (fun r => ([0,1,2] : List Nat).map (fun b => r+m*b))).filter
      (fun r => r < 65536 && remSum r t (m*3) < m*3)

private theorem lowEq (w t m : Nat) : remSum (w%m) t m = remSum w t m := by
  unfold remSum
  simp only [Nat.mul_mod_mod, Nat.add_mod (C*(w%m)), Nat.add_mod (C*w), Nat.mul_mod_mod]

private theorem keepPrefix (w t : Nat) (small : w<65536)
    (clear : ∀ j, remSum w t (3^j)<3^j) (j : Nat) : w%3^j ∈ prefixes t j := by
  induction j with
  | zero => simp [prefixes]; omega
  | succ j ih =>
    let m := 3^j
    let d := w/m%3
    have mem : d ∈ ([0,1,2] : List Nat) := by simp [d]; omega
    have id : w%(m*3) = w%m+m*d := Nat.mod_mul
    simp only [prefixes,Nat.pow_succ,List.mem_filter,List.mem_flatMap,List.mem_map,
               Bool.and_eq_true,decide_eq_true_iff]
    refine ⟨⟨w%m, ih, d, mem, id.symm⟩, ?_, ?_⟩
    · exact Nat.lt_of_le_of_lt (Nat.mod_le ..) small
    · rw [lowEq]; simpa only [Nat.pow_succ] using clear (j+1)

private theorem endPrefix : prefixes 0 12 = [0] ∧ prefixes 1 12 = [0] := by decide

/-- Both the homogeneous and the affine trinomial have a ternary carry at
    every positive scalar. -/
theorem everyScalar (w t : Nat) (pos : 0 < w) (ht : t ≤ 1) :
    ∃ j, total (upper (row w t) j) < total (row w t) / 3^j := by
  by_cases big : 65536 ≤ w
  · exact largeScalar w t big ht
  have bad : ∃ j, 3^j ≤ remSum w t (3^j) := by
    apply Classical.byContradiction
    intro no
    have kept := keepPrefix w t (by omega) (fun j => Nat.lt_of_not_ge (fun h => no ⟨j,h⟩)) 12
    have terminal : prefixes t 12 = [0] := by
      rcases (show t=0 ∨ t=1 by omega) with h|h
      · simpa [h] using endPrefix.1
      · simpa [h] using endPrefix.2
    rw [terminal,Nat.mod_eq_of_lt (show w<3^12 by omega)] at kept
    simp at kept
    omega
  obtain ⟨j,overflow⟩ := bad
  refine ⟨j,?_⟩
  let m := 3^j
  have hm : 0<m := Nat.pow_pos (by decide)
  apply (Nat.le_div_iff_mul_le hm).2
  change ((A*w/m+B*w/m+(C*w+t)/m)+1)*m ≤ A*w+B*w+(C*w+t)
  simp only [Nat.add_mul,Nat.one_mul]
  have a := Nat.div_add_mod (A*w) m
  have b := Nat.div_add_mod (B*w) m
  have c := Nat.div_add_mod (C*w+t) m
  simp only [Nat.mul_comm m] at a b c
  change m ≤ (A*w)%m+(B*w)%m+(C*w+t)%m at overflow
  omega

#print axioms everyScalar
end Wasserman.FiveLeadingDigits
