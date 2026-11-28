import Lean

namespace Wasserman.InverseWeight

/-- A higher-column excess persists after including all lower digits. -/
theorem column_excess (place x y n : Nat) (h : n / place < x / place + y / place) :
    n < x + y := by
  apply Nat.lt_of_not_ge
  intro hn
  have upper := Nat.div_le_div_right (c := place) hn
  have lower : x / place + y / place ≤ (x + y) / place := by
    by_cases hz : place = 0
    · simp [hz]
    · apply (Nat.le_div_iff_mul_le (Nat.pos_of_ne_zero hz)).2
      rw [Nat.add_mul]
      exact Nat.add_le_add (Nat.div_mul_le_self x place) (Nat.div_mul_le_self y place)
  omega

/-- The two required binary layers for the CRT row in the class 5 modulo 8. -/
theorem example_layers (w : Nat) (hw : w % 8 = 5) :
    (35599*w+1) % 8 = 4 ∧
    (4404*w)/8 + (2716*w)/8 + (28479*w+1)/8 + 1 ≤ (35599*w+1)/8 ∧
    (4404*w)/16 + (2716*w)/16 + (28479*w+1)/16 + 1 ≤ (35599*w+1)/16 := by
  omega

end Wasserman.InverseWeight
