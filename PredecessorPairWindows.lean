import Lean

namespace Wasserman.PredecessorPairWindows

/-- Eliminate P from a weight and a linear quotient of one of its parts. -/
theorem norm_divisor (U P D h a b : Int) (weight : U ∣ D*P*(P+h)+1)
    (part : U ∣ a*P+b) : U ∣ a*a+D*b*(b-a*h) := by
  obtain ⟨v, hv⟩ := weight
  obtain ⟨w, hw⟩ := part
  refine ⟨a*a*v-D*(a*P+a*h-b)*w, ?_⟩
  have h₁ := congrArg (fun z : Int => a*a*z) hv
  have h₂ := congrArg (fun z : Int => D*(a*P+a*h-b)*z) hw
  grind

/-- The determinant is an inverse residue at the second predecessor power. -/
theorem inverse_residue (P Q H m h l c : Int) (layer : Q = m*H+h)
    (position : Q ∣ P*(l*H+c)-1) : Q ∣ P*(m*c-l*h)-m := by
  obtain ⟨s, hs⟩ := position
  refine ⟨m*s-P*l, ?_⟩
  have h₁ := congrArg (fun z : Int => m*z) hs
  subst Q
  grind

/-- A double zero in the sum and simple zeros in the two parts. -/
theorem last_shape (P : Int) (m Q : Int) (hm : m = 2+4*P+4*P^2)
    (hq : Q = m*P+1) :
    P*(2*Q-m) = 4*P^3*(2*P+1) ∧
    P*(2*Q-m)+Q*(2*P+1) = (2*P+1)^2*(4*P^2+1) := by
  subst m
  subst Q
  grind

end Wasserman.PredecessorPairWindows
