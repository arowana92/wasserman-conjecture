import Lean

namespace Wasserman.MixedRowMinors

/-- Combine the weight divisibilities of the two cross products and the
predecessor divisibility at the exceptional positions. Here x is the
product of the minors, y the product of the distinguished entries. -/
theorem retained_cell (n m N M x y s r : Int)
    (coprime : Int.gcd m n = 1) (hn : n ∣ N) (hm : m ∣ M)
    (hx : n ∣ x) (hy : n ∣ y) (hp : m ∣ x * y) :
    n^2 * m ∣ x * (s * y - N * M * r) := by
  obtain ⟨a, ha⟩ := hx
  obtain ⟨b, hb⟩ := hy
  obtain ⟨c, hc⟩ := hn
  obtain ⟨d, hd⟩ := hm
  have hxy : x*y = n^2*(a*b) := by rw [ha, hb, Int.pow_succ]; simp only [Int.pow_one]; ac_rfl
  rw [hxy] at hp
  have hquot := Int.dvd_gcd_mul_iff_dvd_mul.mpr hp
  have hsq : Int.gcd m (n^2) = 1 := Int.gcd_pow_right_of_gcd_eq_one coprime
  rw [hsq] at hquot
  simp at hquot
  obtain ⟨k, hk⟩ := hquot
  refine ⟨s*k - a*c*d*r, ?_⟩
  calc
    x * (s*y-N*M*r) = s*(x*y)-x*N*M*r := by rw [Int.mul_sub]; congr 1 <;> ac_rfl
    _ = n^2*m*(s*k) - n^2*m*(a*c*d*r) := by
      rw [hxy, hk, ha, hc, hd, Int.pow_succ]; simp only [Int.pow_one]; congr 1 <;> ac_rfl
    _ = _ := (Int.mul_sub ..).symm

/-- A nonzero pair of minors in the strict integer window fixes the cell. -/
theorem short_cell (modulus left right cell : Int)
    (divides : modulus ∣ left * right * cell)
    (small : (left*right*cell).natAbs < modulus.natAbs)
    (hl : left ≠ 0) (hr : right ≠ 0) : cell = 0 :=
  (Int.mul_eq_zero.mp (Int.eq_zero_of_dvd_of_natAbs_lt_natAbs divides small)).resolve_left
    (fun h => (Int.mul_eq_zero.mp h).elim hl hr)

end Wasserman.MixedRowMinors
