import Lean

namespace Wasserman.ExponentCubeTransfer

/-- Shifting the 5-power changes the inverse by the old quotient times the
    shift divided by the product of the other two powers. -/
theorem shiftInverse (D W x t s : Int) (h : x * W = 1 + D * t) :
    (x + t*s) * W = 1 + (D + W*s) * t := by
  grind

/-- The shifted inverse remains the representative strictly within the new range. -/
theorem shiftRange (D W x t s : Nat) (positive : 0 < x) (within : x < D)
    (quotient : t ≤ W) : 0 < x + t*s ∧ x + t*s < D + W*s := by
  have bound := Nat.mul_le_mul_right s quotient
  omega

/-- A further congruent shift of the two remaining powers preserves the
    quotient through the prescribed layer. The integral witness k can be signed. -/
theorem shiftProduct (D W x t s m k : Int) (h : x * W = 1 + D * t) :
    (x+t*s) * (W + (D+W*s)*m*k) =
      1 + (D+W*s) * (t + m*((x+t*s)*k)) := by
  grind

end Wasserman.ExponentCubeTransfer
