import Lean

namespace Wasserman

/-- With b = B*s, the positive dividend b*d*d*y - z has U as a factor.
    The additive form keeps the entire identity in the natural numbers. -/
theorem scaled_bin_dividend (b d t u v z y : Nat)
    (shape : u = z*t + d*y) (total : v*u = b*d*t + 1)
    (companions : 1 < v*(d*y)) :
    ∃ k : Nat, 0 < k ∧ b*d*d*y = z + k*u := by
  have expanded : v*z*t + v*(d*y) = b*d*t + 1 := by
    rw [← total, shape, Nat.mul_add, Nat.mul_assoc]
  have lt : v*z < b*d := by
    by_cases h : b*d ≤ v*z
    · have := Nat.mul_le_mul_right t h
      omega
    · omega
  let k := b*d - v*z
  have factor : b*d = v*z+k := by dsimp [k]; omega
  have first : (b*d)*u = (b*d)*(z*t) + (b*d*d*y) := by
    rw [shape, Nat.mul_add]
    ac_rfl
  have second : (b*d)*u = (b*d)*(z*t) + z + k*u := calc
    (b*d)*u = (v*z+k)*u := congrArg (fun x => x*u) factor
    _ = z*(v*u) + k*u := by rw [Nat.add_mul]; ac_rfl
    _ = z*(b*d*t) + z + k*u := by rw [total, Nat.mul_add, Nat.mul_one]
    _ = (b*d)*(z*t) + z + k*u := by ac_rfl
  refine ⟨k, ?_, ?_⟩
  · dsimp [k]; omega
  · omega

/-- If one bin is small, normalized positivity in its position forces
    the modulus s*u to be large compared with the cofactor s*d1*d2*d3. -/
theorem cofactor_from_small_bin (s u d₁ d₂ d₃ bound : Nat)
    (h₁ : d₂*d₃ ≤ u) (h₂ : d₁*d₃ ≤ u) (h₃ : d₁*d₂ ≤ u)
    (small : d₁ ≤ bound ∨ d₂ ≤ bound ∨ d₃ ≤ bound) :
    s*(d₁*d₂*d₃) ≤ s*u*bound := by
  have lift (x y : Nat) (hx : x ≤ u) (hy : y ≤ bound) :
      s*(x*y) ≤ s*u*bound := by
    rw [Nat.mul_assoc]
    exact Nat.mul_le_mul_left s (Nat.mul_le_mul hx hy)
  rcases small with ha | hb | hc
  · exact (show s*(d₁*d₂*d₃) = s*((d₂*d₃)*d₁) by ac_rfl) ▸ lift _ _ h₁ ha
  · exact (show s*(d₁*d₂*d₃) = s*((d₁*d₃)*d₂) by ac_rfl) ▸ lift _ _ h₂ hb
  · exact lift _ _ h₃ hc

end Wasserman
