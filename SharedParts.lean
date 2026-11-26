import Lean

namespace Wasserman

/-- At a carry-free layer dividing the sum, both parts have zero residue. -/
theorem divides_both_at_layer (x y g : Nat) (total : g ∣ x+y)
    (free : x%g + y%g < g) : g ∣ x ∧ g ∣ y := by
  have h := Nat.mod_eq_zero_of_dvd total
  rw [Nat.add_mod, Nat.mod_eq_of_lt free] at h
  exact ⟨Nat.dvd_of_mod_eq_zero (by omega), Nat.dvd_of_mod_eq_zero (by omega)⟩

/-- CRT-matched parts agree when their range is shorter than the product
of the common layer modulus and the CRT modulus. The smaller part is named x. -/
theorem shared_cut_rigidity (x y X Y t g : Nat)
    (ordered : x ≤ X) (matched : t ∣ X-x) (coprime : g.Coprime t)
    (total : x+y = X+Y) (multiple : g ∣ x+y)
    (free₁ : x%g+y%g < g) (free₂ : X%g+Y%g < g)
    (width : X-x < g*t) : x = X ∧ y = Y := by
  have hx := (divides_both_at_layer x y g multiple free₁).1
  rw [total] at multiple
  have hX := (divides_both_at_layer X Y g multiple free₂).1
  have bound := coprime.mul_dvd_of_dvd_of_dvd (Nat.dvd_sub hX hx) matched
  by_cases h : X-x = 0
  · omega
  · have := Nat.le_of_dvd (by omega : 0 < X-x) bound
    omega

/-- For a common multiple g coprime to PQ, the shared parts have size
at least gP and gQ when they are respectively multiples of P and Q. -/
theorem shared_scale_budget (x y g P Q : Nat)
    (positive₁ : 0 < x) (positive₂ : 0 < y)
    (coprime₁ : g.Coprime P) (coprime₂ : g.Coprime Q)
    (gx : g ∣ x) (gy : g ∣ y) (px : P ∣ x) (qy : Q ∣ y) :
    g*(P+Q) ≤ x+y := by
  rw [Nat.mul_add]
  exact Nat.add_le_add
    (Nat.le_of_dvd positive₁ (coprime₁.mul_dvd_of_dvd_of_dvd gx px))
    (Nat.le_of_dvd positive₂ (coprime₂.mul_dvd_of_dvd_of_dvd gy qy))

end Wasserman
