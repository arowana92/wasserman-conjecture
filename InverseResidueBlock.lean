import Lean

namespace Wasserman

/-- Dominance of digit residues fails if r is a short negative multiple of t
modulo a block, with t itself in the initial part of the block. -/
theorem inverse_wrap (r t B u : Int)
    (hu : 0 < u) (ht : 0 < t % B)
    (hsize : (u + 1) * (t % B) < B)
    (hcong : B ∣ r + u*t) : t % B < r % B := by
  let a := t % B
  have ha : 0 < a := ht
  have hb : u*a + a < B := by simpa [a, Int.add_mul] using hsize
  have hpos : 0 < u*a := Int.mul_pos hu ha
  have hrepr : a + B*(t/B) = t := Int.emod_add_mul_ediv t B
  obtain ⟨k, hk⟩ := hcong
  have hx : r - (B-u*a) = B*(k-1-u*(t/B)) := by
    have hz : u*a + B*(u*(t/B)) = u*t := calc
      _ = u*(a+B*(t/B)) := by rw [Int.mul_add]; ac_rfl
      _ = _ := by rw [hrepr]
    rw [Int.mul_sub, Int.mul_sub, Int.mul_one]
    omega
  have hmod := Int.emod_eq_emod_iff_emod_sub_eq_zero.mpr
    (Int.dvd_iff_emod_eq_zero.mp (show B ∣ r-(B-u*a) from ⟨_,hx⟩))
  rw [Int.emod_eq_of_lt (by omega : 0 ≤ B-u*a) (by omega)] at hmod
  change a < r % B
  omega

/-- The two residue bounds for a divisor of H+1 cannot both hold. -/
theorem inverse_complement (B b u : Nat)
    (hb : 0 < b) (hu : 0 < u) (ub : u*b < B)
    (upper : B-u*b ≤ b-1)
    (lower : u*b ≤ B-b) : False := by
  have b_le : b ≤ u*b := by simpa using Nat.mul_le_mul_right b hu
  omega

end Wasserman
