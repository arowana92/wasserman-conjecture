import Lean

namespace Wasserman.FiveExponentGaps

/-- Kummer's carry inequalities at the first four powers of five. -/
def fits (x : Nat) : Bool :=
  [5,25,125,625].all (fun m => (3*x)%m + (2*x+1)%m < m)

set_option maxRecDepth 12000
set_option maxHeartbeats 0
set_option exponentiation.threshold 1000

private theorem finiteGap :
    (List.range 625).all (fun x =>
      !(x%5 == 2 && fits x && fits ((16*x)%625))) = true := by decide

private theorem reduced (t : Nat) : fits t = fits (t%625) := by
  have digit (m : Nat) (h : m = 5 ∨ m = 25 ∨ m = 125 ∨ m = 625) :
      (3*t)%m = (3*(t%625))%m ∧ (2*t+1)%m = (2*(t%625)+1)%m := by
    rcases h with h | h | h | h <;> subst m <;> constructor <;> omega
  have h₁ := digit 5 (by decide)
  have h₂ := digit 25 (by decide)
  have h₃ := digit 125 (by decide)
  have h₄ := digit 625 (by decide)
  unfold fits
  simp only [List.all_cons,List.all_nil]
  rw [h₁.1,h₁.2,h₂.1,h₂.2,h₃.1,h₃.2,h₄.1,h₄.2]

private theorem block (k : Nat) : 2^(500*k+4)%625 = 16 := by
  induction k with
  | zero => decide
  | succ k ih =>
    rw [show 500*(k+1)+4 = (500*k+4)+500 by omega, Nat.pow_add, Nat.mul_mod, ih]

/-- No input in the necessary unit class passes all four layers both before
and after any of the infinitely many dilations with exponent 4 modulo 500. -/
theorem incompatible (t k : Nat) (unitClass : t%5 = 2) :
    ¬ (fits t = true ∧ fits (2^(500*k+4)*t) = true) := by
  rintro ⟨first, second⟩
  have mem : t%625 ∈ List.range 625 := List.mem_range.mpr (Nat.mod_lt _ (by decide))
  have finite := List.all_eq_true.mp finiteGap (t%625) mem
  rw [reduced t] at first
  rw [reduced, Nat.mul_mod, block] at second
  have unitRed : (t%625)%5 = 2 := by omega
  rw [unitRed, first, second] at finite
  contradiction

end Wasserman.FiveExponentGaps
