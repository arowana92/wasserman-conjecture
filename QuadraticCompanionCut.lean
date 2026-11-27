import Lean

namespace Wasserman.QuadraticCompanionCut

/-- The exceptional congruence in the smallest group of predecessor powers. -/
theorem group_projection (d S n Z t u B W : Int)
    (weight : t*W = S*d*Z-u) (slot : d ∣ Z*n-B*W) :
    d ∣ t*n*Z+B*u := by
  obtain ⟨c,hc⟩ := slot
  refine ⟨t*c+B*S*Z, ?_⟩
  have h₁ := congrArg (fun x : Int => t*x) hc
  have h₂ := congrArg (fun x : Int => B*x) weight
  grind

/-- Changes of the two cut parts, after clearing the common denominator.
The coefficients of s are the three derivatives in the note. -/
theorem finite_changes (P s t u L n : Int) :
    let H := t*u*P+t+u
    let T := 1+t*P
    let E := L-t*n
    let H' := H+t*u*s
    let T' := T+t*s
    n*(T'*H'-T*H) = s*n*(2*t*H-t^2)+s^2*n*t^2*u ∧
    L*((P+s)*H'-P*H)-n*(T'*H'-T*H) =
      s*(E*(2*H-(t+u))-t*u*n)+s^2*E*t*u := by
  dsimp
  constructor <;> grind

/-- The upper binary block of a quadratic agrees with its tangent. -/
theorem tangent_congruence (a b c x H y : Int) :
    H^2 ∣ (a*(x+H*y)^2+b*(x+H*y)+c) -
             (a*x^2+b*x+c+H*y*(2*a*x+b)) := by
  refine ⟨a*y^2, ?_⟩
  grind

end Wasserman.QuadraticCompanionCut
