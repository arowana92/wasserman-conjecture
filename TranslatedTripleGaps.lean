import Lean

namespace Wasserman.TranslatedTripleGaps

def gap (a b n : Int) : Prop := b - a = n ∨ a - b = n

def twoPairs (a b c d n : Int) : Prop :=
  (gap a b n ∧ gap c d n) ∨ (gap a c n ∧ gap b d n) ∨
  (gap a d n ∧ gap b c n)

def fourOfFive (a b c d e n : Int) : Prop :=
  twoPairs b c d e n ∨ twoPairs a c d e n ∨ twoPairs a b d e n ∨
  twoPairs a b c e n ∨ twoPairs a b c d n

def distinct (a b c d e : Int) : Prop :=
  a ≠ b ∧ a ≠ c ∧ a ≠ d ∧ a ≠ e ∧ b ≠ c ∧ b ≠ d ∧ b ≠ e ∧
  c ≠ d ∧ c ≠ e ∧ d ≠ e

def candidate (a b c d e n : Int) : Prop :=
  distinct a b c d e ∧ fourOfFive a b c d e n

set_option maxHeartbeats 2000000 in
/-- After a single binary merger in a translated triple, two equal gaps
    large enough for the companion inequality force the arithmetic progression. -/
theorem one_merge (i j k n : Int) (hi : 0 < i) (hj : i < j)
    (large : 2*(j+k)+3 < 3*n)
    (hc : (k = i ∧ candidate 0 (i+1) (2*i) j (j+i) n) ∨
          (k = j ∧ candidate 0 i (j+1) (i+j) (2*j) n) ∨
          (k = j-i ∧ candidate 0 i (j-i) (j+1) (2*j-i) n)) :
    j = 2*i ∧ k = j ∧ n = 3*i := by
  simp only [candidate, distinct, fourOfFive, twoPairs, gap] at hc
  rcases hc with ⟨h, _hd, hp⟩ | ⟨h, _hd, hp⟩ | ⟨h, _hd, hp⟩ <;>
    rcases hp with h1 | h1 | h1 | h1 | h1 <;>
    rcases h1 with ⟨h1,h2⟩ | ⟨h1,h2⟩ | ⟨h1,h2⟩ <;>
    rcases h1 with h1 | h1 <;> rcases h2 with h2 | h2 <;> omega

end Wasserman.TranslatedTripleGaps
