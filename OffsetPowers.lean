import Lean

namespace Wasserman

private theorem power_difference (w x y : Int) (base : w ∣ x-y) (n : Nat) :
    w ∣ x^n-y^n := by
  induction n with
  | zero =>
    change w ∣ (1 : Int)-1
    exact Int.dvd_zero w
  | succ n ih =>
    rw [Int.pow_succ, Int.pow_succ]
    have identity : x^n*x-y^n*y = (x^n-y^n)*x + y^n*(x-y) := by
      rw [Int.sub_mul, Int.mul_sub]
      omega
    rw [identity]
    exact Int.dvd_add (Int.dvd_trans ih (Int.dvd_mul_right _ _))
      (Int.dvd_trans base (Int.dvd_mul_left _ _))

/-- A relation between positive powers gives a signed integer divisible by
    every common divisor of the two shifted bases. -/
theorem shifted_power_elimination (w Q P b s left right : Int) (h n : Nat)
    (weight : w ∣ Q+b) (middle : w ∣ P+s)
    (same : left * Q^h = right * P^n) :
    w ∣ left*(-b)^h-right*(-s)^n := by
  have first := power_difference w Q (-b) (by simpa using weight) h
  have second := power_difference w P (-s) (by simpa using middle) n
  have ident : left*(-b)^h-right*(-s)^n =
      right*(P^n-(-s)^n)-left*(Q^h-(-b)^h) := by
    rw [Int.mul_sub, Int.mul_sub, same]
    omega
  rw [ident]
  exact Int.dvd_sub (Int.dvd_trans second (Int.dvd_mul_left _ _))
    (Int.dvd_trans first (Int.dvd_mul_left _ _))

/-- An exact shifted binary row whose scale forces F+1 divisible by four
    already overflows a layer in the two parts RF and 2RF. -/
theorem shifted_row_layer (R F : Nat) (one : 4 ∣ F+1) :
    4*R ≤ (R*F)%(4*R)+(2*R*F)%(4*R) := by
  have residue : F%4 = 3 := by omega
  have doubled : (2*F)%4 = 2 := by omega
  rw [Nat.mul_comm 4 R, show 2*R*F = R*(2*F) by ac_rfl,
    Nat.mul_mod_mul_left, Nat.mul_mod_mul_left, residue, doubled]
  omega

end Wasserman
