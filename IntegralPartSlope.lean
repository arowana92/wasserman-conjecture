import Lean

namespace Wasserman.IntegralPartSlope

/-- Compose the primitive-sum, retained-predecessor and row-product bounds. -/
theorem direction_bound (M N u v H K D J : Nat)
    (weight : u*v ≤ D*N) (predecessor : M ≤ (J*H)*K)
    (rowProduct : 27*K^2 ≤ v^3) :
    27*M^2*u^3 ≤ (D*N)^3*(J*H)^2 := by
  calc
    27*M^2*u^3 ≤ 27*((J*H)*K)^2*u^3 :=
      Nat.mul_le_mul_right _ (Nat.mul_le_mul_left _ (Nat.pow_le_pow_left predecessor 2))
    _ = (27*K^2)*u^3*(J*H)^2 := by rw [Nat.mul_pow]; ac_rfl
    _ ≤ v^3*u^3*(J*H)^2 :=
      Nat.mul_le_mul_right _ (Nat.mul_le_mul_right _ rowProduct)
    _ = (u*v)^3*(J*H)^2 := by rw [Nat.mul_pow u v 3]; ac_rfl
    _ ≤ (D*N)^3*(J*H)^2 := Nat.mul_le_mul_right _ (Nat.pow_le_pow_left weight 3)

/-- The real minimum behind the sharp prime-power constants. -/
theorem coprime_factor_bound (d e : Int) (hd : 0 ≤ d) (he : 0 ≤ e) :
    27*d*e^2 ≤ 4*(d+e)^3 := by
  have identity : 4*(d+e)^3 = 27*d*e^2+(e-2*d)^2*(4*e+d) := by grind
  have h := Int.mul_nonneg (Int.sq_nonneg (e-2*d)) (show 0 ≤ 4*e+d by omega)
  omega

/-- Squared form of the estimate using the product of the two divisors. -/
theorem factor_product_bound (d e : Nat) : 64*(d*e)^3 ≤ (d+e)^6 := by
  have diff := Int.sq_nonneg ((d:Int)-(e:Int))
  have h : 4*d*e ≤ (d+e)^2 := by
    have z : 4*(d:Int)*(e:Int) ≤ ((d:Int)+(e:Int))^2 := by grind
    exact_mod_cast z
  have cubed := Nat.pow_le_pow_left h 3
  grind

end Wasserman.IntegralPartSlope
