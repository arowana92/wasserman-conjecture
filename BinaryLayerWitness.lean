import Lean

set_option exponentiation.threshold 5000
set_option maxRecDepth 16000

namespace Wasserman

/-- This congruence fixes the three predecessor groups for every exponent. -/
theorem binary_predecessor_residue (k : Nat) :
    (4 * 8 ^ (98 + 4760*k)) % 65025 = 9436 := by
  rw [Nat.pow_add, Nat.pow_mul, ← Nat.mul_assoc, Nat.mul_mod]
  have start : (4 * 8^98) % 65025 = 9436 := by decide
  have period : 8^4760 % 65025 = 1 := by decide
  rw [start, Nat.pow_mod (8^4760) k 65025, period, Nat.one_pow]

/-- The exceptional residue at each fixed group and the two base-seven
    digits already overflowing after the common factor seven is removed. -/
theorem binary_witness_groups (r : Nat) (pred : r % 255 = 37) (wt : r % 7 = 1) :
    (238*r) % 15 = 1 ∧ (210*r) % 17 = 1 ∧
    34*r % 7 + 30*r % 7 ≥ 7 := by
  omega

end Wasserman
