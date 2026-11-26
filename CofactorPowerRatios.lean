import Lean

namespace Wasserman

/-- The quadratic remainder in the first-order power expansion is integral. -/
theorem principalUnitExpansion (p t n : Nat) :
    ∃ k : Nat, (1+p*t)^n = 1+p*t*n+p*p*k := by
  induction n with
  | zero => exact ⟨0, by simp⟩
  | succ n ih =>
    obtain ⟨k, hk⟩ := ih
    refine ⟨t*t*n+k+k*p*t, ?_⟩
    rw [Nat.pow_succ, hk]
    grind

set_option maxRecDepth 10000
set_option maxHeartbeats 1000000

/-- On the two allowed residues at 29, the mod-7 ratio test permits
only the identity exponent modulo 28. -/
theorem cofactor986_ratio (r : Fin 28)
    (seven : r.val % 7 = 1 ∨ r.val % 7 = 6)
    (twentynine : Nat.powMod 10 r.val 29 = 10 ∨ Nat.powMod 10 r.val 29 = 21 ∨
                   Nat.powMod 21 r.val 29 = 10 ∨ Nat.powMod 21 r.val 29 = 21) :
    r.val = 1 := by
  have table : ∀ z : Fin 28, (z.val % 7 = 1 ∨ z.val % 7 = 6) →
      (Nat.powMod 10 z.val 29 = 10 ∨ Nat.powMod 10 z.val 29 = 21 ∨
       Nat.powMod 21 z.val 29 = 10 ∨ Nat.powMod 21 z.val 29 = 21) → z.val = 1 := by decide
  exact table r seven twentynine

end Wasserman
