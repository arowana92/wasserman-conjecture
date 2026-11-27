import Lean

namespace Wasserman.InverseScale

private theorem residue (q z k w : Int) (hz : 0 < z) (hk : 0 < k)
    (small : k*z < w) (cong : w ∣ q+k) : (q*z)%w = w-k*z := by
  obtain ⟨t, ht⟩ := cong
  have h : w ∣ q*z-(w-k*z) := by
    refine ⟨t*z-1, ?_⟩
    have id := congrArg (fun v : Int => v*z) ht
    change (q+k)*z = (w*t)*z at id
    rw [Int.add_mul, Int.mul_assoc] at id
    rw [Int.mul_sub, Int.mul_one]
    omega
  have eq := Int.emod_eq_emod_iff_emod_sub_eq_zero.mpr (Int.dvd_iff_emod_eq_zero.mp h)
  have pos := Int.mul_pos hk hz
  rwa [Int.emod_eq_of_lt (by omega : 0 ≤ w-k*z) (by omega)] at eq

/-- Two positive multiples of a short negative residue force a carry at w. -/
theorem two_residues (q k x y w : Int) (hx : 0 < x) (hy : 0 < y)
    (hk : 0 < k) (small : k*(x+y) < w) (cong : w ∣ q+k) :
    w < (q*x)%w + (q*y)%w := by
  rw [Int.mul_add] at small
  have px := Int.mul_pos hk hx
  have py := Int.mul_pos hk hy
  rw [residue q x k w hx hk (by omega) cong, residue q y k w hy hk (by omega) cong]
  omega

end Wasserman.InverseScale
