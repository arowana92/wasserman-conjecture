import Lean

namespace Wasserman

private theorem shifted_residue (n n₀ c c₀ w v q : Int)
    (hw : 0 < w) (hn₀ : 0 ≤ n₀ ∧ n₀ < w) (hc₀ : 0 ≤ c₀ ∧ c₀ < w)
    (hn : n = n₀ + (2*w)*v) (hc : c = c₀+w+(2*w)*q) :
    n % (2*w) < c % (2*w) := by
  rw [hn, hc, Int.add_mul_emod_self_left, Int.add_mul_emod_self_left]
  rw [Int.emod_eq_of_lt (by omega : 0 ≤ n₀) (by omega : n₀ < 2*w),
      Int.emod_eq_of_lt (by omega : 0 ≤ c₀+w) (by omega : c₀+w < 2*w)]
  omega

/-- A cut at a lift of one of the four representatives violates the necessary
low-residue inequality whenever a side of its representative is small. -/
theorem adjacent_lift_carry (n n₀ a a₀ d r w v z : Int)
    (hw : 0 < w) (hn₀ : 0 ≤ n₀ ∧ n₀ < w)
    (hn : n = n₀ + 4*w*v) (hd : d = r + 2*w*(2*z+1))
    (cut : 2*a = n-d-1) (rep : 2*a₀ = n₀-r-1)
    (offset : Int)
    (small : (0 ≤ a₀+offset ∧ a₀+offset < w) ∨
      (0 ≤ n₀-(a₀+offset) ∧ n₀-(a₀+offset) < w)) :
    n % (2*w) < (a+offset) % (2*w) ∨
      n % (2*w) < (n-(a+offset)) % (2*w) := by
  have hn' : n = n₀+(2*w)*(2*v) := by grind
  rcases small with h | h
  · left
    apply shifted_residue n n₀ _ (a₀+offset) w (2*v) (v-z-1) hw hn₀ h hn'
    grind
  · right
    apply shifted_residue n n₀ _ (n₀-(a₀+offset)) w (2*v) (v+z) hw hn₀ h hn'
    grind

/-- The cubic controlling the boundary of the five-sixths estimate. -/
theorem adjacent_boundary_cubic (δ e j l g i h z : Int)
    (gdef : g = 4*l+δ) (idef : i = 4*g*j+4*δ*j+e)
    (edef : e*l = 1+δ^2*j) (hdef : h^2 = i*(i-4*g*j)-4*j)
    (zdef : z = 4*δ*j+e) :
    δ*e*h^2 = δ*z^3+(z-e)*(4*z-e) := by
  rw [hdef, idef, gdef, zdef]
  grind

end Wasserman
