import Lean

namespace Wasserman

/-- Recover the complementary divisor from the congruence at the common
    primitive weight. The analytic size estimate supplies k < n. -/
theorem recover_complement (H n g k : Int) (hk : 0 ≤ k) (hlt : k < n)
    (heq : H = (g*n-1)*k) : (-H) % n = k := by
  have factor : -H = (-(g*k))*n + k := by
    rw [heq, Int.sub_mul, Int.one_mul]
    rw [Int.neg_mul, Int.mul_right_comm g k n]
    omega
  rw [factor, Int.mul_add_emod_self_right, Int.emod_eq_of_lt hk hlt]

/-- Divisors whose product exceeds their common ambient integer must
    have a common factor greater than one. -/
theorem large_endpoint_divisors (m a b : Nat) (hm : 0 < m)
    (ha : a ∣ m) (hb : b ∣ m) (large : m < a*b) : 1 < a.gcd b := by
  have lbound := Nat.le_of_dvd hm (Nat.lcm_dvd ha hb)
  have formula := Nat.gcd_mul_lcm a b
  apply Nat.lt_of_not_ge
  intro small
  have h := Nat.mul_le_mul_right (a.lcm b) small
  simp only [Nat.one_mul] at h
  omega

end Wasserman
