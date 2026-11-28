import Lean

namespace Wasserman.ConsecutiveCuts

/-- Eliminate N from the consecutive-cut equations. -/
theorem integral_curve (a m N q : Int)
    (hprod : a*(a+1) = m*N) (hpred : q*(N-1) = m*(m-2)) :
    q*a*(a+1) = m^3-2*m^2+q*m := by
  calc
    _ = q*(a*(a+1)) := by ac_rfl
    _ = m*(q*(N-1))+q*m := by rw [hprod]; grind
    _ = _ := by rw [hpred]; grind

/-- Map the cut to an integral short Weierstrass equation. -/
theorem short_model (d t q a m : Int) (hsq : q = d*t^2)
    (h : q*a*(a+1) = m^3-2*m^2+q*m) :
    let X := 36*d*m-24*d
    let Y := 108*d^2*t*(2*a+1)
    let c4 := d^2*(64-48*q)
    let c6 := d^3*(512-792*q)
    Y^2 = X^3-27*c4*X-54*c6 := by
  dsimp only
  have aux : 46656*d^3*(q*a*(a+1)) =
      46656*d^3*(m^3-2*m^2+q*m) := by rw [h]
  rw [hsq] at aux ⊢
  grind only

theorem discriminant (d q : Int) :
    (d^2*(64-48*q))^3-(d^3*(512-792*q))^2 =
      -1728*d^6*q*(64*q^2+107*q-128) := by grind

/-- Difference of squares when m and q have the same squarefree part. -/
theorem square_gap (a m q d t v : Int) (hd : 0 < d)
    (hm : m = d*v^2) (hq : q = d*t^2)
    (h : q*(2*a+1)^2 = 4*m^3-8*m^2+4*q*m+q) :
    (t*(2*a+1))^2-(2*v*(m-1))^2 = 4*(q-1)*v^2+t^2 := by
  apply Int.eq_of_mul_eq_mul_left (show d ≠ 0 by omega)
  rw [hm, hq] at h ⊢
  grind only

end Wasserman.ConsecutiveCuts
