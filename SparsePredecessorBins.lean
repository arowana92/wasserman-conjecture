import Lean

namespace Wasserman.SparsePredecessorBins

/-- Two nonempty groups give coprime odd indices at least 3 and 5. -/
theorem two_groups (a b : Int) (ha : 3 ≤ a) (hb : 5 ≤ b) :
    5 * (a+b-2) ≤ 2*a*b := by
  have h : 0 ≤ (a-3)*(b-5) := Int.mul_nonneg (by omega) (by omega)
  grind

/-- Three different coprime odd indices are separated by at least two. -/
theorem three_groups (a b c : Int) (ha : 3 ≤ a)
    (hb : a+2 ≤ b) (hc : b+2 ≤ c) :
    35*(a+b+c-3) ≤ 12*b*c := by
  let s := b-5
  let t := c-b-2
  have hs : 0 ≤ s := by dsimp [s]; omega
  have ht : 0 ≤ t := by dsimp [t]; omega
  have hst := Int.mul_nonneg hs ht
  have hss := Int.mul_nonneg hs hs
  have h : 12*b*c-35*(a+b+c-3) =
      12*(s*s)+12*(s*t)+39*s+25*t+35*(b-2-a) := by
    dsimp [s,t]
    grind
  omega

end Wasserman.SparsePredecessorBins
