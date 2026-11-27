import Lean

namespace Wasserman.BlockRatios

/-- A common evaluation divisor divides the determinant of two affine forms. -/
theorem divides_affine_minor (u p i j k l : Int)
    (hx : u ∣ i*p+j) (hy : u ∣ k*p+l) : u ∣ i*l-k*j := by
  obtain ⟨s, hs⟩ := hx
  obtain ⟨t, ht⟩ := hy
  refine ⟨i*t-k*s, ?_⟩
  grind

/-- Two blocks with the same quotient magnify their cross difference by the base. -/
theorem block_step (b q x y x' y' d f : Int)
    (hx : b*x = q*d+x') (hy : b*y = q*f+y') :
    x'*f-y'*d = b*(x*f-y*d) := by
  grind

/-- A submask inverse in two blocks gives a small signed determinant modulo the total. -/
theorem inverse_minor (base high low upper lower w k : Int)
    (hinv : w*(upper*base+lower) = 1+k*(high*base+low)) :
    high*base+low ∣ w*(high*lower-low*upper)-high := by
  refine ⟨high*k-w*upper, ?_⟩
  grind

/-- The two independent determinant bounds contradict the endpoint size inequality. -/
theorem successor_size (d u s e : Nat) (h : 27*d^2 ≤ u^3*s^2)
    (hu : u ≤ e) (hs : s ≤ e^2) : 27*d^2 ≤ e^7 := by
  have hprod := Nat.mul_le_mul (Nat.pow_le_pow_left hu 3) (Nat.pow_le_pow_left hs 2)
  calc
    27*d^2 ≤ u^3*s^2 := h
    _ ≤ e^3*(e^2)^2 := hprod
    _ = e^7 := by grind

end Wasserman.BlockRatios
