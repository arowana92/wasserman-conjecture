import Lean

namespace Wasserman.ShiftedPowerCarries

/-- Distinct carries in two possibly overlapping stretches bound the retained exponent. -/
theorem retained_exponent (v i j carries r : Nat) (hij : i ≤ j)
    (hr : r = v - carries)
    (h : (v-i) + (i + min i v - j) ≤ carries) :
    r ≤ min i (j-i) := by
  omega

/-- The tied binary case loses at most one carry from the second stretch. -/
theorem equal_binary (v i carries : Nat)
    (h : (v-i) + min i v ≤ carries + 1) : v-carries ≤ 1 := by
  omega

/-- Both exponent bounds putting the gcd product below the required size. -/
theorem height_bounds (a n : Nat) (h : 90*n + 9 < a) (hn : 0 < n) :
    a + 59*n + 26 < 2*a ∧ 3*a + 20*n + 20 < 4*a := by
  omega

end Wasserman.ShiftedPowerCarries
