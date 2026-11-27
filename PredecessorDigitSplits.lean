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

/-- The uniform error margin after using h >= t+1 and M >= t. -/
theorem margin_identity (t : Int) :
    t*((t+1)^2+4*t)-16*t^2+32*(t+1) = t*((t-5)^2+8)+32 := by
  grind

/-- Necessary small-root data. The outputs retain the ordered roots, the two
low residues, the scale denominator, and D modulo the product of the roots. -/
def smallRoots (h : Nat) : List (List Nat) := Id.run do
  let mut result := []
  for t in [2:h] do
    let u := h-t
    let m := t*u
    for i in [1:t/2+1] do
      for j in [i:t-i+1] do
        for b in [1:i+1] do
          if t != u && (i % b == 0) && (j % b == 0) && Nat.gcd b t == 1 &&
              (i != j || b % 2 != 0) && b*u != 1 then
            for d in [0:m] do
              let k := if d*d % m == 0 then m else d*d % m
              if h % Nat.gcd d m == 0 && (2*t < b*k || (2*t == b*k && (i+j)*h < t*t)) then
                result := result ++ [[t,u,i,j,b,d]]
  return result

set_option maxRecDepth 8000 in
set_option maxHeartbeats 2000000 in
theorem small_cases :
    ([1,2,3,4,5,6,7,8].map smallRoots) =
      [[],[],[],[],[],[],[[2,5,1,1,1,3],[2,5,1,1,1,7]],[]] := by decide +kernel

end Wasserman.PredecessorDigitSplits
