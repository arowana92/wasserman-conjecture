import Lean

namespace Wasserman

/-- Necessary factor and part inequalities bound Q independently of the two exponents. -/
theorem higherSuccessorSize (Q e f T U d y : Nat)
    (factors : Q = e * f * T) (companion : 27 * T^2 ≤ U^3)
    (exceptional : e ≤ d) (unshifted : U*f ≤ y) (scale : U ≤ y) (small : y < d) :
    27 * Q^2 < d^5 := by
  have productBound : e * (U*f) ≤ d*d :=
    Nat.mul_le_mul exceptional (Nat.le_trans unshifted (Nat.le_of_lt small))
  have bound : 27 * Q^2 ≤ (d*d)^2 * U := calc
    27 * Q^2 = (e*f)^2 * (27*T^2) := by rw [factors]; grind
    _ ≤ (e*f)^2 * U^3 := Nat.mul_le_mul_left _ companion
    _ = (e*(U*f))^2 * U := by grind
    _ ≤ (d*d)^2 * U := Nat.mul_le_mul (Nat.pow_le_pow_left productBound 2) (Nat.le_refl U)
  have positive : 0 < d := by omega
  calc
    27 * Q^2 ≤ (d*d)^2 * U := bound
    _ < (d*d)^2 * d := Nat.mul_lt_mul_of_pos_left (by omega)
      (Nat.pow_pos (Nat.mul_pos positive positive))
    _ = d^5 := by grind

private theorem fourProduct (x y : Nat) : 4*x*y ≤ (x+y)^2 := by
  rcases Nat.le_total x y with h | h
  · obtain ⟨t, ht⟩ := Nat.exists_eq_add_of_le h
    rw [ht]
    have expand : (x+(x+t))^2 = 4*x*(x+t)+t*t := by grind
    omega
  · obtain ⟨t, ht⟩ := Nat.exists_eq_add_of_le h
    rw [ht]
    have expand : ((y+t)+y)^2 = 4*(y+t)*y+t*t := by grind
    omega

/-- When both small parts are unshifted the same conditions give the quadratic bound. -/
theorem higherSuccessorShiftedLarge (Q e f T U d x y : Nat)
    (factors : Q = e * f * T) (companion : 27 * T^2 ≤ U^3)
    (scale : 3 ≤ U) (first : U*e ≤ x) (second : U*f ≤ y) (small : x+y ≤ d) :
    36 * Q ≤ d^2 := by
  have square : (9*Q)^2 ≤ (x*y)^2 := calc
    (9*Q)^2 = (27*Q^2)*3 := by grind
    _ ≤ (27*Q^2)*U := Nat.mul_le_mul_left _ scale
    _ = ((e*f)^2 * U) * (27*T^2) := by rw [factors]; grind
    _ ≤ ((e*f)^2 * U) * U^3 := Nat.mul_le_mul_left _ companion
    _ = ((U*e)*(U*f))^2 := by grind
    _ ≤ (x*y)^2 := Nat.pow_le_pow_left (Nat.mul_le_mul first second) 2
  have linear : 9*Q ≤ x*y := (Nat.pow_le_pow_iff_left (by decide : 2 ≠ 0)).mp square
  calc
    36*Q = 4*(9*Q) := by grind
    _ ≤ 4*(x*y) := Nat.mul_le_mul_left 4 linear
    _ = 4*x*y := by grind
    _ ≤ (x+y)^2 := fourProduct x y
    _ ≤ d^2 := Nat.pow_le_pow_left small 2

end Wasserman
