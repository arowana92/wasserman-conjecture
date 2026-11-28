import Lean

namespace Wasserman.ShiftedQuadraticRoots

/-- Recover a monic quadratic from the two distinct negative roots. -/
theorem monic_roots (b c r d : Int) (distinct : r ≠ d)
    (first : r^2-b*r+c = 0) (second : d^2-b*d+c = 0) :
    b = r+d ∧ c = r*d := by
  have prod : (r-d)*(r+d-b) = 0 := by grind
  have sum := (Int.mul_eq_zero.mp prod).resolve_left (by omega)
  constructor <;> grind

/-- A sum reaching one entire digit unit can be carry-free at that
    modulus only if the other summand is zero. -/
private theorem one_unit (x y w : Nat) (hx : 0 < x) (total : x+y = w)
    (free : x%w+y%w < w) : y = 0 := by
  apply Classical.byContradiction
  intro hy
  have h₁ : x < w := by omega
  have h₂ : y < w := by omega
  rw [Nat.mod_eq_of_lt h₁, Nat.mod_eq_of_lt h₂] at free
  omega

/-- In the remaining polynomial (k-r)(X+d), put E=k-r. A line
    alpha*(X+r) cannot leave a positive part while passing the two
    Kummer layers E*d and E*X. Both moduli are powers of the base. -/
theorem two_units (X E d r alpha beta gamma : Nat)
    (ha : 0 < alpha) (hr : 0 < r)
    (linear : alpha+gamma = E) (constant : alpha*r+beta = E*d)
    (disjoint : E*d < X) (alignment : E*d ∣ X)
    (low : (alpha*X+alpha*r)%(E*d)+(gamma*X+beta)%(E*d) < E*d)
    (high : (alpha*X+alpha*r)%(E*X)+(gamma*X+beta)%(E*X) < E*X) :
    beta = 0 ∧ gamma = 0 := by
  have left : E*d ∣ alpha*X := Nat.dvd_trans alignment (Nat.dvd_mul_left X alpha)
  have right : E*d ∣ gamma*X := Nat.dvd_trans alignment (Nat.dvd_mul_left X gamma)
  simp [Nat.add_mod, Nat.mod_eq_zero_of_dvd left, Nat.mod_eq_zero_of_dvd right] at low
  have hbeta := one_unit (alpha*r) beta (E*d) (Nat.mul_pos ha hr) constant low
  refine ⟨hbeta, ?_⟩
  apply Classical.byContradiction
  intro hgamma
  have la : alpha+1 ≤ E := by omega
  have lg : gamma+1 ≤ E := by omega
  have boundA := Nat.mul_le_mul_right X la
  have boundG := Nat.mul_le_mul_right X lg
  rw [Nat.add_mul, Nat.one_mul] at boundA boundG
  have hA : alpha*X+alpha*r < E*X := by omega
  have hG : gamma*X+beta < E*X := by omega
  rw [Nat.mod_eq_of_lt hA, Nat.mod_eq_of_lt hG] at high
  have total := congrArg (fun t => t*X) linear
  rw [Nat.add_mul] at total
  omega

end Wasserman.ShiftedQuadraticRoots
