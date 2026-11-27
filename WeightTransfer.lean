import Lean

namespace Wasserman

/-- The high digits of a shifted part at the enlarged layer give its
    original residue. Dropping the bottom digit gives a lower bound. -/
theorem lift_shifted_layer (Q m a b c : Nat) (hq : 1 < Q)
    (overflow : m ≤ a % m + b % m + c % m) :
    Q*m ≤ (Q*a) % (Q*m) + (Q*b) % (Q*m) + (Q*c+1) % (Q*m) := by
  have bound := Nat.div_mul_le_self ((Q*c+1) % (Q*m)) Q
  rw [Nat.mod_mul_right_div_self, Nat.mul_add_div (by omega : 0 < Q),
      Nat.div_eq_of_lt hq, Nat.add_zero, Nat.mul_comm (c % m) Q] at bound
  rw [Nat.mul_mod_mul_left, Nat.mul_mod_mul_left]
  have scaled := Nat.mul_le_mul_left Q overflow
  rw [Nat.mul_add, Nat.mul_add] at scaled
  exact Nat.le_trans scaled (Nat.add_le_add_left bound _)

/-- Any choice of the part to decrease still overflows the layer eight.
    The statement is symmetric up to naming that part c. -/
theorem transfer_binary_obstruction (Q a b c : Nat) (hq : 1 < Q)
    (sum : 8 < a % 8 + b % 8 + c % 8) (hc : 0 < c) :
    Q*8 ≤ (Q*a) % (Q*8) + (Q*b) % (Q*8) + (Q*(c-1)+1) % (Q*8) := by
  apply lift_shifted_layer Q 8 a b (c-1) hq
  omega

/-- The complete list of sorted low-residue triples left by the small-weight
    computation satisfies the strict inequality used in the obstruction. -/
theorem transfer_final_residues :
    [(2,5,5), (2,7,7), (3,6,7), (4,5,7), (4,6,6), (5,5,6)].all
        (fun (v : Nat × Nat × Nat) => decide (8 < v.1 + v.2.1 + v.2.2)) = true := by
  decide

end Wasserman
