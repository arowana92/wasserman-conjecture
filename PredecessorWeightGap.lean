import Lean

namespace Wasserman.PredecessorWeightGap

/-- A part P*z*W cannot receive the single low digit. The two input
congruences express D*W = N-1 and the residue of that part after division by P. -/
theorem no_single_digit_cut (B D P W z : Int)
    (positive : 0 < z) (power : 1 < P)
    (small : P*z < D) (gap : D ≤ B)
    (weight : B ∣ D*W-(P-1)) (cut : B ∣ z*W-1) : False := by
  obtain ⟨a, ha⟩ := weight
  obtain ⟨b, hb⟩ := cut
  have difference : D-z*(P-1) = B*(z*a-D*b) := by
    have h₁ := congrArg (fun u : Int => z*u) ha
    have h₂ := congrArg (fun u : Int => D*u) hb
    grind
  have hz : 0 < z*(P-1) := Int.mul_pos positive (by omega)
  have id : z*(P-1) = P*z-z := by grind
  have hsmall : (D-z*(P-1)).natAbs < B.natAbs := by omega
  have zero := Int.eq_zero_of_dvd_of_natAbs_lt_natAbs
    (show B ∣ D-z*(P-1) from ⟨_, difference⟩) hsmall
  omega

/-- If the residue sum is one, two specified parts cannot both contribute. -/
theorem one_zero (a b c : Nat) (sum : a+b+c = 1) : a = 0 ∨ b = 0 := by
  omega

end Wasserman.PredecessorWeightGap
