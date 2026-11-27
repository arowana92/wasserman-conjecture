import Lean

namespace Wasserman

/-- A predecessor divisor at which the pure part is one and the low
singleton part is zero divides a number independent of both free powers. -/
theorem singleton_elimination (w k l q low : Int)
    (hpure : w ∣ k*q-1) (hlow : w ∣ l*q+low) : w ∣ l+low*k := by
  obtain ⟨a, ha⟩ := hpure
  obtain ⟨b, hb⟩ := hlow
  refine ⟨k*b-l*a, ?_⟩
  calc
    l+low*k = k*(l*q+low)-l*(k*q-1) := by grind
    _ = w*(k*b-l*a) := by rw [ha, hb]; grind

/-- Numerical implication of the split-singleton divisibilities. Here u is
the weight factor with its base-prime part removed, and s the predecessor factor. -/
theorem singleton_cutoff (low m k l u s : Nat)
    (hpos : 1 ≤ low) (hsum : l+k ≤ m) (hu : u ≤ k)
    (hs : u*s ≤ k*(l+low*k)) : (low*u)^3*s^2 ≤ (low*m)^5 := by
  have hl : l ≤ low*l := calc
    l = 1*l := by simp
    _ ≤ low*l := Nat.mul_le_mul_right l hpos
  have hedge : l+low*k ≤ low*m := calc
    l+low*k ≤ low*l+low*k := Nat.add_le_add_right hl _
    _ = low*(l+k) := by grind
    _ ≤ low*m := Nat.mul_le_mul_left _ hsum
  have hk : k ≤ m := by omega
  have hus := Nat.le_trans hs (Nat.mul_le_mul hk hedge)
  calc
    (low*u)^3*s^2 = low^3*u*(u*s)^2 := by grind
    _ ≤ low^3*m*(m*(low*m))^2 :=
      Nat.mul_le_mul (Nat.mul_le_mul_left _ (Nat.le_trans hu hk))
        (Nat.pow_le_pow_left hus 2)
    _ = (low*m)^5 := by grind

end Wasserman
