import Lean

namespace Wasserman.TwoDigitWeightCuts

/-- Multiples of the same scalar lift to an exact rational ratio
when the two cross products lie in one residue interval. -/
theorem lift_ratio (B u x y a b : Int)
    (hx : B ∣ u*x-a) (hy : B ∣ u*y-b)
    (low : -B < x*b-y*a) (high : x*b-y*a < B) : x*b = y*a := by
  obtain ⟨c, hc⟩ := hx
  obtain ⟨d, hd⟩ := hy
  have eqn : x*b-y*a = B*(y*c-x*d) := by
    have h₁ := congrArg (fun z : Int => y*z) hc
    have h₂ := congrArg (fun z : Int => x*z) hd
    grind
  have small : (x*b-y*a).natAbs < B.natAbs := by omega
  have zero := Int.eq_zero_of_dvd_of_natAbs_lt_natAbs
    (show B ∣ x*b-y*a from ⟨_, eqn⟩) small
  omega

/-- The exact fraction for a distinguished part P*z*W. -/
theorem cut_fraction (B W D P T z e : Int)
    (weight : B ∣ W*D-(P*T-1)) (part : B ∣ W*z-e)
    (d0 : 0 ≤ D*e) (d1 : D*e < B)
    (z0 : 0 ≤ z*(P*T-1)) (z1 : z*(P*T-1) < B) :
    D*e = z*(P*T-1) :=
  lift_ratio B W D z (P*T-1) e weight part (by omega) (by omega)

/-- If D*G divides x*s*t but D exceeds (1+r)*x*(s+t),
then s+t exceeds 4*(1+r)*G, even at an equality in AM-GM. -/
theorem split_sum (D G x r s t : Nat) (dg : 0 < G)
    (positive : 0 < x*s*t) (divides : D*G ∣ x*s*t)
    (size : (1+r)*x*(s+t) < D) : 4*(1+r)*G < s+t := by
  have h := Nat.le_of_dvd positive divides
  have strict := Nat.mul_lt_mul_of_pos_right size dg
  have cancelling : x*((1+r)*(s+t)*G) < x*(s*t) := by
    calc
      _ = ((1+r)*x*(s+t))*G := by ac_rfl
      _ < D*G := strict
      _ ≤ x*s*t := h
      _ = _ := by ac_rfl
  have cut := Nat.lt_of_mul_lt_mul_left cancelling
  have amgm : 4*s*t ≤ (s+t)*(s+t) := by
    by_cases order : s ≤ t
    · obtain ⟨z, hz⟩ := Nat.exists_eq_add_of_le order
      rw [hz]
      grind
    · obtain ⟨z, hz⟩ := Nat.exists_eq_add_of_le (show t ≤ s by omega)
      rw [hz]
      grind
  have scaled := Nat.mul_lt_mul_of_pos_left cut (by decide : 0 < 4)
  apply Nat.lt_of_mul_lt_mul_right (a := s+t)
  calc
    _ = 4*((1+r)*(s+t)*G) := by ac_rfl
    _ < 4*(s*t) := scaled
    _ ≤ (s+t)*(s+t) := by simpa [Nat.mul_assoc] using amgm

end Wasserman.TwoDigitWeightCuts
