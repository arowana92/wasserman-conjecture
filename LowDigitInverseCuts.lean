import Lean

namespace Wasserman.LowDigitInverseCuts

/-- Lift a modular identity when the positive side is strictly inside the
modulus and the other side can reach its upper endpoint. -/
theorem liftPositive {a b m : Nat} (ha : 0 < a) (ham : a < m)
    (hb : b ≤ m) (hc : a % m = b % m) : a = b := by
  rw [Nat.mod_eq_of_lt ham] at hc
  by_cases h : b = m
  · subst b
    simp at hc
    omega
  · rwa [Nat.mod_eq_of_lt (by omega : b < m)] at hc

/-- At a small predecessor residue the equal fractions are exact. -/
theorem predecessor {x P d e A : Nat} (hx : x < P) (hxpos : 0 < x)
    (hd : 0 < d) (he : e ≤ d) (hA : P*d ≤ A)
    (hc : (x*d) % A = (P*e) % A) : x*d = P*e :=
  liftPositive (Nat.mul_pos hxpos hd)
    (Nat.lt_of_lt_of_le (Nat.mul_lt_mul_of_pos_right hx hd) hA)
    (Nat.le_trans (Nat.mul_le_mul_left _ he) hA) hc

/-- Converting residues to proportions for each part of the companion. -/
theorem companion {g U d e A : Nat} (hg : g < U) (hgpos : 0 < g)
    (hd : 0 < d) (he : e ≤ d) (hA : U*d < A)
    (hc : (g*d) % A = (U*e) % A) : g*d = U*e ∧ 0 < e := by
  have H := predecessor hg hgpos hd he (Nat.le_of_lt hA) hc
  refine ⟨H, ?_⟩
  cases e with
  | zero =>
    rw [Nat.mul_zero] at H
    have hpos := Nat.mul_pos hgpos hd
    omega
  | succ k => omega

/-- With one low digit the first-row cut cannot have positive residue
strictly below that digit. -/
theorem oneDigit {P R x : Nat} (hx : x < P) (hr : 0 < R)
    (e : Nat) (he : e ≤ 1) (hn : 0 < (R-1)*x)
    (hh : (R-1)*x = R*P*e) : False := by
  have hlt : (R-1)*x < R*P :=
    Nat.lt_of_le_of_lt (Nat.mul_le_mul_right x (by omega : R-1 ≤ R))
      (Nat.mul_lt_mul_of_pos_left hx hr)
  cases e with
  | zero => simp at hh; omega
  | succ k =>
    have hk : k = 0 := by omega
    subst k
    simp at hh
    omega

end Wasserman.LowDigitInverseCuts
