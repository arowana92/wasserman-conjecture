import Lean

namespace Wasserman.NestedPowerResidues

/-- Finite residue obstructions for all three possible nonempty ten-cofactor shapes. -/
def Inner (u x : Nat) : Prop :=
  ∀ k : Fin (10-x), 0 < k.val → ∀ z : Fin (5*u), (k.val*z.val^2+1) % (5*u) ≠ 0

def Outer (u x : Nat) : Prop :=
  ∀ k : Fin (10-x), 0 < k.val → ∀ z : Fin (5*u), Nat.gcd z.val (5*u) = 1 →
    (x*z.val^4+k.val)%5 ≠ 0 ∨ (10*z.val^4)%u ≠ k.val%u

instance (u x : Nat) : Decidable (Inner u x) := by unfold Inner; infer_instance
instance (u x : Nat) : Decidable (Outer u x) := by unfold Outer; infer_instance

set_option maxRecDepth 12000 in
set_option maxHeartbeats 1000000 in
theorem inner_squares : Inner 3 3 ∧ Inner 3 6 ∧ Inner 7 7 := by decide

set_option maxRecDepth 12000 in
set_option maxHeartbeats 1000000 in
theorem outer_fourths : Outer 3 3 ∧ Outer 3 6 ∧ Outer 7 7 := by decide

end Wasserman.NestedPowerResidues
