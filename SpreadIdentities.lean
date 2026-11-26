import Lean

namespace Wasserman

private theorem cross_term (g b c : Int) :
    (g*b-1)*c-(g*c+1)*b = -(b+c) := by
  rw [Int.sub_mul, Int.add_mul, Int.one_mul, Int.one_mul]
  have same : g*b*c = g*c*b := by ac_rfl
  omega

/-- The two integer quotients from the outer parts of the sorted row
    determine its step and the quadratic equation for its primitive gaps. -/
theorem spread_quadratic (r s h g P B C : Int)
    (span : h = r+s) (nonzero : h ≠ 0)
    (left : B*r*h = P*(g*(r+h)-1))
    (right : C*s*h = P*(g*(s+h)+1)) :
    B*r+C*s = 3*P*g ∧ B*r*(r+2*s)-C*s*(2*r+s) = -(3*P) := by
  have total : (r+h)+(s+h) = 3*h := by omega
  constructor
  · apply Int.eq_of_mul_eq_mul_right nonzero
    calc
      (B*r+C*s)*h = P*(g*(r+h)-1)+P*(g*(s+h)+1) := by
        rw [Int.add_mul, left, right]
      _ = P*(g*((r+h)+(s+h))) := by
        rw [← Int.mul_add, Int.mul_add g (r+h) (s+h)]
        congr 1
        omega
      _ = 3*P*g*h := by rw [total]; ac_rfl
  · have first : r+2*s = s+h := by omega
    have second : 2*r+s = r+h := by omega
    rw [first, second]
    apply Int.eq_of_mul_eq_mul_right nonzero
    calc
      (B*r*(s+h)-C*s*(r+h))*h = (B*r*h)*(s+h)-(C*s*h)*(r+h) := by
        rw [Int.sub_mul]
        congr 1 <;> ac_rfl
      _ = P*((g*(r+h)-1)*(s+h)-(g*(s+h)+1)*(r+h)) := by
        rw [left, right, Int.mul_sub P ((g*(r+h)-1)*(s+h)) ((g*(s+h)+1)*(r+h))]
        simp only [Int.mul_assoc]
      _ = -(3*P)*h := by
        rw [cross_term, total, Int.mul_neg, Int.neg_mul]
        congr 1
        ac_rfl

end Wasserman
