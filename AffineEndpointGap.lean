import Lean

namespace Wasserman.AffineEndpointGap

/-- The companion product bound after replacing the retained factors by
    a linear form and its predecessor lcm. -/
theorem gap (M d j u s e l : Nat) (hu : u ≤ e) (hs : s ≤ l) (hj : 1 ≤ j)
    (bound : 27*M^2 ≤ j^2*d^3*u^3*s^2) :
    27*M^2 ≤ (d*j)^3*e^3*l^2 := by
  have hj' : j^2 ≤ j^3 := by
    calc j^2 = j^2*1 := by simp
         _ ≤ j^2*j := Nat.mul_le_mul_left _ hj
         _ = j^3 := by simp [Nat.pow_succ]
  calc
    27*M^2 ≤ j^2*d^3*u^3*s^2 := bound
    _ ≤ j^3*d^3*e^3*l^2 :=
      Nat.mul_le_mul (Nat.mul_le_mul (Nat.mul_le_mul_right _ hj')
        (Nat.pow_le_pow_left hu 3)) (Nat.pow_le_pow_left hs 2)
    _ = (d*j)^3*e^3*l^2 := by simp [Nat.mul_pow]; ac_rfl

/-- The signed divisor equation used in counting near the two-ninths bound. -/
theorem divisor_pair (a w i j k x y z c : Int)
    (eq : a*i*j*k = w*(x*j*k+y*i*k+z*i*j)+c) :
    ((a*i-w*x)*j-w*i*y)*((a*i-w*x)*k-w*i*z) =
      (a*i-w*x)*c + w^2*i^2*y*z := by
  have rearr : (a*i-w*x)*j*k-w*i*y*k-w*i*z*j = c := by grind
  rw [← rearr]
  grind

end Wasserman.AffineEndpointGap
