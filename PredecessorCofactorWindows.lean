import Lean

namespace Wasserman.PredecessorCofactorWindows

/-- Extract the signed integer t directly from the quotient of the forced part. -/
theorem signedLift (P T h d x e z n : Int) (hP : 0 < P)
    (low : x*d = 1+P*e)
    (part : x*(h*(P*T)+d) = 1+P*z+P*(P*T)*n) :
    z = e+T*(x*h-P*n) ∧ d*(x*h-P*n) = h+P*(h*e-d*n) := by
  have cancel : P*(z-(e+T*(x*h-P*n))) = 0 := by grind
  have first := Int.mul_eq_zero.mp cancel
  constructor <;> grind

/-- Bounding both ends avoids a choice of sign for the CRT representative. -/
theorem signedBounds (d T t z e limit : Int)
    (hd : 0 < d) (hz : z = e+T*t)
    (ze : 0 ≤ z ∧ z ≤ limit) (ee : 0 ≤ e ∧ e ≤ limit) :
    -(d*limit) ≤ T*(d*t) ∧ T*(d*t) ≤ d*limit := by
  have hlo : -limit ≤ T*t := by omega
  have hhi : T*t ≤ limit := by omega
  have a := Int.mul_le_mul_of_nonneg_left hlo (Int.le_of_lt hd)
  have b := Int.mul_le_mul_of_nonneg_left hhi (Int.le_of_lt hd)
  constructor <;> grind

/-- Small tails can only divide the signed central residue. -/
theorem isolate (d t c P j : Int) (hP : 0 < P)
    (hv : d*t = c+j*P) (lo : c-P < d*t) (hi : d*t < c+P) : d ∣ c := by
  have hj : j = 0 := by
    by_cases a : j < 0
    · have b := Int.mul_le_mul_of_nonneg_right (show j ≤ -1 by omega) (Int.le_of_lt hP)
      grind
    · by_cases b : 0 < j
      · have q := Int.mul_le_mul_of_nonneg_right (show 1 ≤ j by omega) (Int.le_of_lt hP)
        grind
      · omega
  exact ⟨t, by grind⟩

/-- In the branch where a different CRT part takes the leading digit,
    the small complementary factor gives the square-root gap. -/
theorem gap (d B delta A k m T : Nat) (large : delta*B ≤ d)
    (identity : d*(delta*B) = m*(A*B)+T+d*A*k) :
    (m+delta*k)*(A*B)+T ≤ d*d := by
  have left := Nat.mul_le_mul_left d large
  have right := Nat.mul_le_mul_right (A*k) large
  grind

/-- At the boundary P = S^2, the multiplier divides N exactly as it divides
    z-1, because N-(z-1) has the following integral factor. -/
theorem boundary (S z : Int) :
    (1+S^2*((1+z*S^2)*S^2+S+1))-(z-1) =
      (S^2-S+1)*(z*(S^4+S^3-S-1)+S^2+2*S+2) := by grind

end Wasserman.PredecessorCofactorWindows
