import Lean

namespace Wasserman.NestedLeadingLayers

/-- The cell index z records floor(M W / H), with H the leading base-r power of W. -/
def cover (x t k r c : Nat) : Bool := Id.run do
  let d := x+t+k
  let m := Nat.lcm (Nat.lcm x t) (Nat.lcm (c*k+1) d)
  return (List.range ((r-1)*m)).all fun i =>
    let z := m+i
    (List.range 4).any fun j =>
      let w := m*r^j
      decide (x*z/w + t*z/w + (c*k+1)*z/(c*w) < d*z/w)

set_option maxRecDepth 16000 in
set_option maxHeartbeats 4000000 in
theorem tables :
    cover 3 1 1 2 4 ∧
    cover 3 3 1 2 1 ∧ cover 3 2 2 2 1 ∧ cover 5 1 1 2 6 ∧
    cover 3 4 1 2 1 ∧ cover 5 2 1 2 1 ∧
    cover 4 4 1 2 1 ∧ cover 5 3 1 3 2 ∧ cover 5 2 2 2 4 ∧
    cover 3 6 1 3 2 ∧ cover 3 5 2 3 1 ∧ cover 3 3 4 3 1 ∧
    cover 6 3 1 3 2 ∧ cover 6 2 2 2 2 ∧ cover 7 2 1 5 9 := by decide

theorem seventeen (q : Nat) (low : q % 5 = 4) :
    25 ≤ (119*q)%25 + (16*q+1)%25 + (35*q)%25 := by omega

/-- Passing to the integral cell index is exact for each of the four floor functions. -/
theorem grid_floor (a b M W H : Nat) (ha : 0 < a) (hm : 0 < M) (h : a ∣ M) :
    a*W/(H*b) = a*(M*W/H)/(M*b) := by
  obtain ⟨l, rfl⟩ := h
  have hl : 0 < l := by
    cases l with
    | zero => simp at hm
    | succ n => exact Nat.zero_lt_succ n
  rw [Nat.mul_assoc a l b, Nat.mul_div_mul_left _ _ ha, Nat.div_div_eq_div_mul]
  have hnum : a*l*W = l*(a*W) := by ac_rfl
  have hden : H*(l*b) = l*(H*b) := by ac_rfl
  rw [hnum, hden, Nat.mul_div_mul_left _ _ hl]

end Wasserman.NestedLeadingLayers
