import Lean

namespace Wasserman.StaticPowerTransfer

/-- Two local bounds on the same primitive ratio and one cross-prime
congruence force both local multipliers to be exactly the same. -/
theorem glue (D X Y u v : Nat) (xpos : 0 < X) (ypos : 0 < Y)
    (positive : 1 < u*X+v*Y) (cong : D ∣ u*X+v*Y-1)
    (upper₁ : u*(X+Y) ≤ D+1) (upper₂ : v*(X+Y) ≤ D+1) :
    u = v ∧ u*(X+Y) = D+1 := by
  have enough := Nat.le_of_dvd (by omega : 0 < u*X+v*Y-1) cong
  rw [Nat.mul_add] at upper₁ upper₂
  have eq : u = v := by
    by_cases h : u < v
    · have lt := Nat.mul_lt_mul_of_pos_right h xpos
      omega
    by_cases h' : v < u
    · have lt := Nat.mul_lt_mul_of_pos_right h' ypos
      omega
    omega
  subst v
  rw [Nat.mul_add]
  omega

/-- A grouped scale divides a bounded integer after replacing W by one. -/
theorem scale_divisor (u a b j e : Int) (weight : u ∣ a*b+1)
    (part : u ∣ j*b+e) : u ∣ a*e-j := by
  obtain ⟨c, hc⟩ := weight
  obtain ⟨d, hd⟩ := part
  refine ⟨a*d-j*c, ?_⟩
  have h₁ := congrArg (fun z : Int => a*z) hd
  have h₂ := congrArg (fun z : Int => j*z) hc
  grind

end Wasserman.StaticPowerTransfer
