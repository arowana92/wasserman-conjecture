import Lean

namespace Wasserman.PredecessorDigitSplits

/-- Eliminate P between a quadratic weight and the linear equation of a split part. -/
theorem quadratic_divisor (P h M V b x : Int)
    (weight : V ∣ 1+h*P+M*P*P) (part : V ∣ b*P+x) :
    V ∣ b*b-h*x*b+M*x*x := by
  obtain ⟨u, hu⟩ := weight
  obtain ⟨i, hi⟩ := part
  refine ⟨b*b*u-i*(h*b+M*(b*P-x)), ?_⟩
  have h₁ := congrArg (fun z : Int => b*b*z) hu
  have h₂ := congrArg (fun z : Int => z*(h*b+M*(b*P-x))) hi
  grind

/-- The two unshifted parts force their combined low residue below the root. -/
theorem root_budget (P W t i j : Nat) (hp : 1 ≤ P)
    (room : W*(P*i+P*j) < W*(1+t*P)) : i+j ≤ t := by
  have small := Nat.lt_of_mul_lt_mul_left room
  by_cases h : i+j ≤ t
  · exact h
  have big := Nat.mul_le_mul_left P (show t+1 ≤ i+j by omega)
  rw [Nat.mul_add, Nat.mul_one, Nat.mul_add] at big
  rw [Nat.mul_comm t P] at small
  omega

end Wasserman.PredecessorDigitSplits
