import Lean

namespace Wasserman.SevenQuadraticCofactor

/-- The two occupied parts above the prime-power shift force a strict residue bound. -/
theorem strict_layer (t u D E z w k l K : Int)
    (ht : 0 < t) (hu : 0 < u) (hD : 0 < D) (hE : 0 < E)
    (budget : t*(z+w) < D)
    (projection : (t*(k+l)-K)*E ≤ t^2*(z+w)-D*(t+u)) :
    t*(k+l) < K := by
  have bound := Int.mul_lt_mul_of_pos_left budget ht
  have pos := Int.mul_pos hD hu
  have negative : t^2*(z+w)-D*(t+u) < 0 := by grind
  by_cases nonneg : 0 ≤ t*(k+l)-K
  · have scaled := Int.mul_nonneg nonneg (Int.le_of_lt hE)
    omega
  · omega

/-- Exact endpoints for each accepted binary word of the common scale. -/
theorem input_square (P : Int) :
    (10*P+6)^2 = 5*((2*P+1)*(10*P+7))+1 := by grind

/-- Appending an input bit and doubling the dyadic block preserves the square invariant. -/
theorem root_lift (K n s e b c : Int) (h : K*n = s^2+e) :
    (2*K)*(2*n+b) = (2*s+c)^2+(4*e+2*K*b-c*(4*s+c)) := by grind

/-- The independent progression search doubles the index of its current grid point. -/
theorem progression_lift (K n y e y₀ b : Int) (h : K*n = y^2+e) :
    (2*K)*(2*n+b) = (2*y-y₀)^2+(4*e+y₀*(4*y-y₀)+2*K*b) := by grind

/-- Distance between consecutive values of the unrounded quadratic on the input class. -/
theorem parameter_step (P D : Int) :
    (2*(P+D)+1)*(10*(P+D)+7) - (2*P+1)*(10*P+7) =
      D*(40*P+24+20*D) := by grind

end Wasserman.SevenQuadraticCofactor
