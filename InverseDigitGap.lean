import Lean

namespace Wasserman

/-- A block of zero digits below a nonzero digit of Q obstructs a
    summand z defined by P*z+1=u*Q. Here B is the place of the forced
    digit of z and P is the length of the shift, expressed as a power. -/
theorem inverse_digit_gap (p P B H D u z : Nat)
    (pP : p ∣ P) (hp : 0 < p) (shift : 0 < P)
    (low : 0 < D) (gap : D < B) (unit : 0 < u) (small : u ≤ P)
    (digit : 0 < (u*H)%p) (eqn : P*z+1 = u*(P*B*H+D)) :
    (P*B*H+D)%(p*B) < z%(p*B) := by
  let top := u*H*B
  have eq' : P*z+1 = P*top + u*D := by
    rw [eqn, Nat.mul_add]
    congr 1
    dsimp [top]
    ac_rfl
  have udpos := Nat.mul_pos unit low
  have udlt : u*D < P*B := Nat.lt_of_le_of_lt (Nat.mul_le_mul_right D small)
    (Nat.mul_lt_mul_of_pos_left gap shift)
  have lower : top ≤ z := by
    apply Nat.le_of_mul_le_mul_left (c := P)
    omega
    exact shift
  have upper : z < top+B := by
    apply Nat.lt_of_mul_lt_mul_left (a := P)
    rw [Nat.mul_add]
    omega
  let tail := z-top
  have zform : z = top+tail := by dsimp [tail]; omega
  have tail_lt : tail < B := by omega
  let d := (u*H)%p
  have dlt : d < p := Nat.mod_lt _ hp
  have chunk : d*B+tail < p*B := calc
    _ < d*B+B := Nat.add_lt_add_left tail_lt _
    _ = (d+1)*B := (Nat.succ_mul ..).symm
    _ ≤ _ := Nat.mul_le_mul_right _ (by omega)
  have split : top = p*B*((u*H)/p) + d*B := by
    change u*H*B = _
    conv => lhs; rw [← Nat.div_add_mod (u*H) p]
    rw [Nat.add_mul]
    congr 1 <;> ac_rfl
  have zmod : z%(p*B) = d*B+tail := by
    rw [zform, split, Nat.add_assoc, Nat.mul_add_mod, Nat.mod_eq_of_lt chunk]
  obtain ⟨j, hj⟩ := pP
  have weight : P*B*H = p*B*(j*H) := by rw [hj]; ac_rfl
  have B_le : B ≤ p*B := by simpa using Nat.mul_le_mul_right B hp
  have bottom : B ≤ d*B := by simpa using Nat.mul_le_mul_right B digit
  rw [weight, Nat.mul_add_mod, Nat.mod_eq_of_lt (Nat.lt_of_lt_of_le gap B_le), zmod]
  omega

/-- At a layer dividing a*h, the last part a*(h-1) forces a carry
    as soon as the residue of h reaches the positive cofactor a. -/
theorem inverse_cofactor_carry (h a m : Nat) (positive : 0 < h) (ha : 0 < a)
    (hm : 0 < m) (layer : m ∣ a*h) (large : a ≤ h%m) :
    m ≤ h%m + (a*(h-1))%m := by
  have bound := Nat.mod_lt h hm
  have small : a < m := by omega
  have id : a*(h-1)+a = a*h := by
    rw [← Nat.mul_succ, show (h-1).succ = h by omega]
  have zero : (a*(h-1)+a)%m = 0 := by rw [id]; exact Nat.mod_eq_zero_of_dvd layer
  by_cases carry : m ≤ h%m + (a*(h-1))%m
  · exact carry
  · have stays : (a*(h-1))%m+a < m := by omega
    rw [Nat.add_mod, Nat.mod_eq_of_lt small, Nat.mod_eq_of_lt stays] at zero
    omega

end Wasserman
