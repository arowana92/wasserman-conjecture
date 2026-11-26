import Lean

namespace Wasserman

private theorem negative_residue (x y m : Int)
    (hy : 0 < y) (hym : y < m) (h : m ∣ x + y) :
    x % m = m - y := by
  obtain ⟨v, hv⟩ := h
  have hx : x = -y + v * m := by
    rw [Int.mul_comm] at hv
    omega
  rw [hx, Int.add_mul_emod_self, Int.neg_emod, Int.emod_eq_of_lt] <;> omega

private theorem carry_from_deficit (a b q m u : Int)
    (ha : 0 < a) (hb : 0 < b) (hu : 0 < u)
    (size : (a+b)*u ≤ m) (congr : m ∣ q+u) :
    m ≤ (a*q) % m + (b*q) % m := by
  obtain ⟨v,hv⟩ := congr
  have mult (x : Int) : m ∣ x*q + x*u := by
    refine ⟨x*v, ?_⟩
    rw [← Int.mul_add, hv]
    ac_rfl
  have apos := Int.mul_pos ha hu
  have bpos := Int.mul_pos hb hu
  rw [Int.add_mul] at size
  rw [negative_residue (a*q) (a*u) m apos (by omega) (mult a),
      negative_residue (b*q) (b*u) m bpos (by omega) (mult b)]
  omega

/-- If v lies t below a multiple of d, the indicated gap forces a carry at m. -/
theorem carry_at_orbit_gap (a b c d q m v s t : Int)
    (ha : 0 < a) (hb : 0 < b) (hc : 0 ≤ c) (hd : 0 < d) (hm : 0 < m)
    (scaled : d*q = v*m+c) (nearby : d*s = v+t)
    (gap : (a+b)*t ≤ d) (offset : c < t*m) :
    m ≤ (a*q) % m + (b*q) % m := by
  have deficit : d*(s*m-q) = t*m-c := by
    rw [Int.mul_sub, ← Int.mul_assoc, nearby, Int.add_mul, scaled]
    omega
  have pos : 0 < s*m-q := by
    apply Int.lt_of_mul_lt_mul_left (a := d) _ (Int.le_of_lt hd)
    rw [Int.mul_zero, deficit]
    omega
  have size := Int.mul_le_mul_of_nonneg_right gap (Int.le_of_lt hm)
  have wc : 0 ≤ (a+b)*c := Int.mul_nonneg (by omega) hc
  have bound : d*((a+b)*(s*m-q)) ≤ d*m := calc
    d*((a+b)*(s*m-q)) = (a+b)*(d*(s*m-q)) := by ac_rfl
    _ = ((a+b)*t)*m-(a+b)*c := by rw [deficit, Int.mul_sub, Int.mul_assoc]
    _ ≤ d*m := by omega
  apply carry_from_deficit a b q m (s*m-q) ha hb pos (Int.le_of_mul_le_mul_left bound hd)
  refine ⟨s, ?_⟩
  rw [Int.mul_comm m s]
  omega

/-- If the modulus divides n and is 1 modulo d, the two residues force a carry. -/
theorem carry_at_order_layer (a b c d q m n : Int)
    (ha : 0 < a) (hb : 0 < b) (hc : 0 < c) (hd : 0 < d) (hm : 1 < m)
    (size : (a + b) * c ≤ d) (scaled : d * q = n + c)
    (layer : m ∣ n) (order : d ∣ m - 1) :
    m ≤ (a * q) % m + (b * q) % m := by
  obtain ⟨s, hs⟩ := order
  have spos : 0 < s := by
    by_cases h : 0 < s
    · exact h
    · have := Int.mul_nonpos_of_nonneg_of_nonpos (Int.le_of_lt hd) (by omega : s ≤ 0)
      omega
  obtain ⟨t, ht⟩ := layer
  have ident : q + c * s = (q - t * s) * m := by
    have eq1 : q * (d * s) = n*s + c*s := by
      calc
        q * (d * s) = (d*q)*s := by ac_rfl
        _ = n*s + c*s := by rw [scaled, Int.add_mul]
    rw [← hs, Int.mul_sub, Int.mul_one] at eq1
    have eq2 : n*s = (t*s)*m := by rw [ht]; ac_rfl
    rw [eq2] at eq1
    rw [Int.sub_mul]
    omega
  have bound := Int.mul_le_mul_of_nonneg_right size (Int.le_of_lt spos)
  rw [Int.mul_assoc, ← hs] at bound
  apply carry_from_deficit a b q m (c*s) ha hb (Int.mul_pos hc spos) (by omega)
  exact ⟨q-t*s, ident.trans (Int.mul_comm _ _)⟩

/-- Carry freedom through h digits rules out every order exponent from 1 to h. -/
theorem no_small_order (a b c d q : Int) (p h : Nat)
    (hp : 1 < p) (ha : 0 < a) (hb : 0 < b) (hc : 0 < c) (hd : 0 < d)
    (size : (a+b)*c ≤ d) (scaled : d*q = (p^h : Nat) + c)
    (free : ∀ j : Nat, 0 < j → j ≤ h →
      (a*q) % (p^j : Nat) + (b*q) % (p^j : Nat) < (p^j : Nat)) :
    ∀ j : Nat, 0 < j → j ≤ h → ¬ d ∣ (p^j : Nat) - (1 : Int) := by
  intro j hj jh hdvd
  have small : p^1 ≤ p^j := Int.pow_le_pow_of_le_right (by omega) hj
  simp at small
  have hm : (1 : Int) < (p^j : Nat) := by omega
  have hpow : (↑(p^j) : Int) ∣ (↑(p^h) : Int) :=
    Int.ofNat_dvd.mpr ((Nat.pow_dvd_pow_iff_le_right hp).mpr jh)
  have overflow := carry_at_order_layer a b c d q (p^j : Nat) (p^h : Nat)
    ha hb hc hd hm size scaled hpow hdvd
  have strict := free j hj jh
  omega

/-- Necessary residues for the factor-five and factor-seven rows clash at matched inputs. -/
theorem five_seven_clash (t v u z : Nat)
    (same_three : t%3 = v%3) (same_four : u%4 = z%4) (odd : u%2 = 1)
    (five_scale : (5*t+1)%3 = 0) (five_cut : (3*u)%4+(2*u)%4 < 4)
    (seven_shape : (7*v+1)%3 = 0 ∨ (5*z)%4+(2*z)%4 < 4) : False := by
  omega

/-- Two necessary binary layers after removing the common power of two. -/
abbrev lowCut (u a b : Nat) : Prop :=
  (a*u)%4+(b*u)%4 < 4 ∧ (a*u)%8+(b*u)%8 < 8

abbrev fivePossible (t u : Nat) : Prop :=
  (3*t)%5 = 1 ∧ (5*t+1)%3 = 0 ∧ lowCut u 3 2

abbrev sevenPossible (t u : Nat) : Prop :=
  ((3*t)%7 = 1 ∧ (7*t+1)%3 = 0 ∧ lowCut u 3 4) ∨
  ((6*t)%7 = 1 ∧ (7*t+1)%3 = 0 ∧ lowCut u 6 1) ∨
  ((5*t)%7 = 1 ∧ (7*t+1)%5 = 0 ∧ lowCut u 5 2)

abbrev ninePossible (t u : Nat) : Prop :=
  ((5*t)%9 = 1 ∧ (9*t+1)%5 = 0 ∧ lowCut u 5 4) ∨
  ((7*t)%9 = 1 ∧ (9*t+1)%7 = 0 ∧ lowCut u 7 2)

abbrev compatibleAtMostOne (t u : Nat) : Prop :=
  ¬(fivePossible t u ∧ sevenPossible t u) ∧
  ¬(fivePossible t u ∧ ninePossible t u) ∧
  ¬(sevenPossible t u ∧ ninePossible t u)

set_option maxRecDepth 16000 in
set_option maxHeartbeats 2000000 in
/-- The admissible residue pairs for the three families are pairwise disjoint. -/
theorem three_families_table :
    (List.range 315).all (fun t => [1,3,5,7].all (fun u =>
      decide (compatibleAtMostOne t u))) = true := by
  decide

end Wasserman
