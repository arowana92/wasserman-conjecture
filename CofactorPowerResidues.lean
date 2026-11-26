import Lean

namespace Wasserman

set_option maxRecDepth 12000
set_option maxHeartbeats 4000000
set_option exponentiation.threshold 1000

private theorem pow_reduce (x e m : Nat) :
    x^e % m = Nat.powMod (x%m) e m := by
  rw [Nat.powMod_def, ← Nat.pow_mod]

/-- Neither residue of a cube-eligible scale at denominator 986 is a seventh power. -/
theorem denominator986_seventh (x : Nat) :
    x^7 % 29 ≠ 10 ∧ x^7 % 29 ≠ 21 := by
  have table : ∀ z : Fin 29, Nat.powMod z.val 7 29 ≠ 10 ∧
      Nat.powMod z.val 7 29 ≠ 21 := by decide
  rw [pow_reduce]
  exact table ⟨x%29, Nat.mod_lt _ (by decide)⟩

/-- A seventeenth power of residue 16 mod 17 has just one lift mod 289. -/
theorem denominator986_seventeenth (x : Nat) (h : x^17 % 17 = 16) :
    x^17 % 289 = 288 := by
  have table : ∀ z : Fin 289, Nat.powMod (z.val % 17) 17 17 = 16 →
      Nat.powMod z.val 17 289 = 288 := by decide
  rw [pow_reduce] at h ⊢
  apply table ⟨x%289, Nat.mod_lt _ (by decide)⟩
  change Nat.powMod ((x%289)%17) 17 17 = 16
  have congr : (x%289)%17 = x%17 := by omega
  rw [congr]; exact h

/-- The lift at 29 needed for the second cube-eligible row. -/
theorem denominator986_power841 (x : Nat) (h : x^841 % 29 = 21) :
    x^841 % 24389 = 9388 := by
  have small : ∀ z : Fin 29, Nat.powMod z.val 841 29 = z.val := by decide
  have low : x%29 = 21 := by
    rw [pow_reduce] at h
    exact (small ⟨x%29, Nat.mod_lt _ (by decide)⟩).symm.trans h
  have table : ∀ k : Fin 841, Nat.powMod (21+29*k.val) 841 24389 = 9388 := by decide
  rw [pow_reduce]
  let k : Fin 841 := ⟨(x%24389)/29, by omega⟩
  have eqn : x%24389 = 21+29*k.val := by dsimp [k]; omega
  rw [eqn]
  exact table k

/-- The first cube-eligible row already carries at a twenty-ninth power. -/
theorem denominator986_power29 (x : Nat) (h : x^29 % 29 = 10) :
    x^29 % 841 = 416 := by
  have table : ∀ z : Fin 841, Nat.powMod (z.val % 29) 29 29 = 10 →
      Nat.powMod z.val 29 841 = 416 := by decide
  rw [pow_reduce] at h ⊢
  apply table ⟨x%841, Nat.mod_lt _ (by decide)⟩
  change Nat.powMod ((x%841)%29) 29 29 = 10
  have congr : (x%841)%29 = x%29 := by omega
  rw [congr]; exact h

end Wasserman
