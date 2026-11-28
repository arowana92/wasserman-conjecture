import Lean

namespace Wasserman.CrossCutRatios

/-- Elimination of the weight and the other ratio term. -/
theorem eliminate (c J m x y k h L N t : Int)
    (hN : h*N = c*k*x*y)
    (hM : h*(N-1) = J*m*t*(c*y-m)) (hL : L = J*m*t-k*x) :
    c*L*y + h = J*m^2*t := by
  grind

/-- The two discarded terms in the fifth-power bound are positive. -/
theorem exact_deficit (c J m x y k h L N t : Nat)
    (hN : h*N = c*k*x*y) (hL : k*x+L = J*m*t)
    (he : c*L*y+h = J*m^2*t) :
    h*L*N + k*x*h + L*(J*m^2*t) = J^2*m^3*t^2 := by
  calc
    _ = (k*x)*(c*L*y+h) + L*(J*m^2*t) := by
      rw [show h*L*N = (h*N)*L by ac_rfl, hN]
      grind
    _ = (k*x+L)*(J*m^2*t) := by rw [he]; grind
    _ = _ := by rw [hL]; grind

/-- The sharp maximum of m^3 t^2 at fixed m+t. -/
theorem fifth_maximum (m t : Nat) : 3125*m^3*t^2 ≤ 108*(m+t)^5 := by
  have ident (a b : Int) : 108*(a+b)^5 - 3125*a^3*b^2 =
      (3*b-2*a)^2*(27*a^3+216*a^2*b+76*a*b^2+12*b^3) := by grind
  have h := ident (m : Int) (t : Int)
  have hp : 0 ≤ ((3*(t : Int)-2*(m : Int))^2) := Int.sq_nonneg _
  have hq : 0 ≤ (27*(m : Int)^3+216*(m : Int)^2*t+76*(m : Int)*t^2+12*(t : Int)^3) := by
    exact Int.natCast_nonneg (27*m^3+216*m^2*t+76*m*t^2+12*t^3)
  have hn := Int.mul_nonneg hp hq
  have hi : 3125*(m : Int)^3*t^2 ≤ 108*((m : Int)+t)^5 := by omega
  exact Int.ofNat_le.mp hi

/-- The strict cross-cut bound for the integer parameters obtained from the cuts. -/
theorem fifth_bound (c J m x y k h L N t : Nat)
    (hsum : m+t = c*x) (hk : 0 < k) (hx : 0 < x) (hh : 0 < h)
    (hN : h*N = c*k*x*y) (hL : k*x+L = J*m*t)
    (he : c*L*y+h = J*m^2*t) :
    3125*(h*L*N) < 108*J^2*(c*x)^5 := by
  have deficit := exact_deficit c J m x y k h L N t hN hL he
  have pos := Nat.mul_pos (Nat.mul_pos hk hx) hh
  have strict : h*L*N < J^2*m^3*t^2 := by omega
  calc
    _ < 3125*(J^2*m^3*t^2) := Nat.mul_lt_mul_of_pos_left strict (by decide)
    _ = J^2*(3125*m^3*t^2) := by ac_rfl
    _ ≤ J^2*(108*(m+t)^5) := Nat.mul_le_mul_left _ (fifth_maximum m t)
    _ = _ := by rw [hsum]; ac_rfl

end Wasserman.CrossCutRatios
