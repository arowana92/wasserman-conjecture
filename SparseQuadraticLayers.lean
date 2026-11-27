import Lean

namespace Wasserman.SparseQuadraticLayers

private theorem residue (a m k r : Int) (h : a = r + m*k)
    (lo : 0 ≤ r) (hi : r < m) : a % m = r := by
  rw [h, Int.add_mul_emod_self_left, Int.emod_eq_of_lt lo hi]

/-- Exact no-carry test after removing the bottom block of the repeated-part row. -/
theorem middle_layer (b w m k z t : Int) (hw : w = m*k)
    (hb : b = t + m*z) (pos : 0 < t) (lt : t < m) :
    (b^2*w-2*b) % m + 2*((b^2*w-b) % m) < m ↔ 3*m < 4*t := by
  have h₁ : (b^2*w-b) % m = m-t := by
    apply residue _ _ (b^2*k-z-1) <;> first | grind | omega
  rw [h₁]
  by_cases half : 2*t ≤ m
  · have h₂ : (b^2*w-2*b) % m = m-2*t := by
      apply residue _ _ (b^2*k-2*z-1) <;> first | grind | omega
    rw [h₂]
    omega
  · have h₂ : (b^2*w-2*b) % m = 2*m-2*t := by
      apply residue _ _ (b^2*k-2*z-2) <;> first | grind | omega
    rw [h₂]
    omega

/-- The quadratic weight splits into two fixed-position blocks and its bottom unit. -/
theorem separated_blocks (d q p h l : Int)
    (hp : 3*p = d*q+1) (hh : 3*h = d^2-1) (hl : 3*l = q+4*d) :
    p*(d*q+3) = 1+l*q+h*q^2 := by
  grind

end Wasserman.SparseQuadraticLayers
