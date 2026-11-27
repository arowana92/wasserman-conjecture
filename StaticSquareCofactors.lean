import Lean

namespace Wasserman.StaticSquareCofactors

theorem scale_options (u : Nat) (lo : 3 ≤ u) (hi : u < 10)
    (odd : u % 2 ≠ 0) (five : u % 5 ≠ 0) : u = 3 ∨ u = 7 ∨ u = 9 := by
  omega

theorem scale_three (w : Nat) (assigned : 3 ∣ 10*w+1)
    (full : ¬ 9 ∣ 10*w+1) : w % 9 = 2 ∨ w % 9 = 5 := by
  omega

theorem scale_seven (x w : Nat) (size : 0 < x ∧ x < 10) (part : 7 ∣ x)
    (assigned : 7 ∣ 10*w+1) (five : (x*w)%5 = 1) :
    w%7 = 2 ∧ w%5 = 3 := by
  have eq : x = 7 := by omega
  subst x
  omega

private theorem pow_reduce (x e m : Nat) : x^e % m = Nat.powMod (x%m) e m := by
  rw [Nat.powMod_def, ← Nat.pow_mod]

theorem square (t : Nat) : (t^2 % 9 ≠ 2 ∧ t^2 % 9 ≠ 5) ∧ t^2 % 5 ≠ 3 := by
  have h₉ : ∀ z : Fin 9, Nat.powMod z.val 2 9 ≠ 2 ∧ Nat.powMod z.val 2 9 ≠ 5 := by decide
  have h₅ : ∀ z : Fin 5, Nat.powMod z.val 2 5 ≠ 3 := by decide
  simp only [pow_reduce]
  exact ⟨h₉ ⟨_, Nat.mod_lt t (by decide)⟩, h₅ ⟨_, Nat.mod_lt t (by decide)⟩⟩

theorem cube (t : Nat) : (t^3 % 9 ≠ 2 ∧ t^3 % 9 ≠ 5) ∧ t^3 % 7 ≠ 2 := by
  have h₉ : ∀ z : Fin 9, Nat.powMod z.val 3 9 ≠ 2 ∧ Nat.powMod z.val 3 9 ≠ 5 := by decide
  have h₇ : ∀ z : Fin 7, Nat.powMod z.val 3 7 ≠ 2 := by decide
  simp only [pow_reduce]
  exact ⟨h₉ ⟨_, Nat.mod_lt t (by decide)⟩, h₇ ⟨_, Nat.mod_lt t (by decide)⟩⟩

/-- Independent residues for the static prime power and the possible row divisor. -/
def NoSquareScale (h b : Nat) : Prop :=
  ∀ x : Fin (h*b), 3 ≤ x.val →
    (∃ z : Fin b, (x.val*z.val^2)%b = 1) →
    ∀ y : Fin (max 1 x.val), Nat.gcd x.val (h*b*y.val^2+1) ≤ 2

instance (h b : Nat) : Decidable (NoSquareScale h b) := by
  unfold NoSquareScale
  infer_instance

set_option maxRecDepth 32000 in
set_option maxHeartbeats 4000000 in
theorem square_tables : NoSquareScale 2 3 ∧ NoSquareScale 2 5 ∧ NoSquareScale 2 9 ∧
    NoSquareScale 2 11 ∧ NoSquareScale 2 29 := by decide

set_option maxRecDepth 32000 in
set_option maxHeartbeats 4000000 in
theorem prime_power_tables : NoSquareScale 1 2 ∧ NoSquareScale 1 3 ∧ NoSquareScale 1 4 ∧
    NoSquareScale 1 5 ∧ NoSquareScale 1 8 ∧ NoSquareScale 1 9 ∧ NoSquareScale 1 13 ∧
    NoSquareScale 1 16 ∧ NoSquareScale 1 25 ∧ NoSquareScale 1 37 := by decide

end Wasserman.StaticSquareCofactors
