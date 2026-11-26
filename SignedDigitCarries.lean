import Lean

namespace Wasserman.SignedDigits

/-- Repeated division drops the specified number of base-b digits. -/
def dropDigits (b value : Int) : Nat → Int
  | 0 => value
  | i+1 => dropDigits b value i / b

def digit (b value : Int) (i : Nat) : Int := dropDigits b value i % b

/-- One Euclidean division normalizes a signed coefficient and its incoming borrow. -/
theorem normalize (b c h : Int) (base : 4 ≤ b)
    (lo : 2-b < 2*c) (hi : 2*c < b) (borrow : -1 ≤ h ∧ h ≤ 0) :
    (-1 ≤ (c+h)/b ∧ (c+h)/b ≤ 0) ∧
    (c < 0 → b < 2*((c+h)%b)) ∧ (0 < c → 2*((c+h)%b) < b) := by
  by_cases neg : c+h < 0
  · have rem := Int.emod_eq_of_lt (a := c+h+b) (b := b) (by omega) (by omega)
    have quot := Int.ediv_eq_zero_of_lt (a := c+h+b) (b := b) (by omega) (by omega)
    have adj := Int.add_mul_ediv_right (c+h) 1 (c := b) (by omega)
    have modAdj := Int.add_mul_emod_self (a := c+h) (b := (1 : Int)) (c := b)
    simp only [Int.one_mul] at adj modAdj
    rw [quot] at adj
    rw [rem] at modAdj
    omega
  · rw [Int.ediv_eq_zero_of_lt (by omega) (by omega),
        Int.emod_eq_of_lt (by omega) (by omega)]
    omega

/-- tail i represents the signed polynomial starting at coefficient i.
    Its actual quotient differs from that evaluation by at most one borrow. -/
theorem propagated (b : Int) (c tail : Nat → Int) (base : 4 ≤ b)
    (size : ∀ i, 2-b < 2*c i ∧ 2*c i < b)
    (split : ∀ i, tail i = c i + tail (i+1)*b) (i : Nat) :
    ∃ h : Int, (-1 ≤ h ∧ h ≤ 0) ∧ dropDigits b (tail 0) i = tail i+h := by
  induction i with
  | zero => exact ⟨0, by decide, by simp [dropDigits]⟩
  | succ i ih =>
    obtain ⟨h, hh, eqn⟩ := ih
    refine ⟨(c i+h)/b, (normalize b (c i) h base (size i).1 (size i).2 hh).1, ?_⟩
    simp only [dropDigits, eqn, split i]
    rw [show c i + tail (i+1)*b + h = (c i+h) + tail (i+1)*b by omega,
        Int.add_mul_ediv_right _ _ (by omega)]
    omega

theorem digit_sign (b : Int) (c tail : Nat → Int) (base : 4 ≤ b)
    (size : ∀ i, 2-b < 2*c i ∧ 2*c i < b)
    (split : ∀ i, tail i = c i + tail (i+1)*b) (i : Nat) :
    (c i < 0 → b < 2*digit b (tail 0) i) ∧
    (0 < c i → 2*digit b (tail 0) i < b) := by
  obtain ⟨h, hh, eqn⟩ := propagated b c tail base size split i
  unfold digit
  rw [eqn, split i,
      show c i + tail (i+1)*b + h = (c i+h) + tail (i+1)*b by omega,
      Int.add_mul_emod_self]
  exact (normalize b (c i) h base (size i).1 (size i).2 hh).2

/-- If each part has a negative position and the total permits at most
    one such position, the ordinary base-b addition must carry. -/
theorem two_negative_positions (b x y z : Int) (limit : Nat) (base : 4 ≤ b)
    (a f t : Nat → Int)
    (negA : ∃ i < limit, a i < 0) (negF : ∃ i < limit, f i < 0)
    (uniquePos : ∀ i < limit, ∀ j < limit, t i ≤ 0 → t j ≤ 0 → i = j)
    (highA : ∀ i, a i < 0 → b < 2*digit b x i)
    (highF : ∀ i, f i < 0 → b < 2*digit b y i)
    (lowT : ∀ i, 0 < t i → 2*digit b z i < b) :
    ¬ ∀ i, digit b x i + digit b y i = digit b z i := by
  intro free
  obtain ⟨i, ib, ai⟩ := negA
  obtain ⟨j, jb, fj⟩ := negF
  have hx := highA i ai
  have hy := highF j fj
  have nneg (v : Int) (k : Nat) : 0 ≤ digit b v k := Int.emod_nonneg _ (by omega)
  have bounded (v : Int) (k : Nat) : digit b v k < b := Int.emod_lt_of_pos _ (by omega)
  have ti : t i ≤ 0 := by
    by_cases pos : 0 < t i
    · have := lowT i pos
      have := free i
      have := nneg y i
      omega
    omega
  have tj : t j ≤ 0 := by
    by_cases pos : 0 < t j
    · have := lowT j pos
      have := free j
      have := nneg x j
      omega
    omega
  have eqn := uniquePos i ib j jb ti tj
  subst j
  have := free i
  have := bounded z i
  omega

end Wasserman.SignedDigits
