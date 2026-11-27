import Lean

namespace Wasserman.PartGcdGap

/-- Every common divisor of a and the predecessor divides this small resultant. -/
theorem resultant (a M r n d e h S : Int)
    (hn : n = a + M*r) (hh : h = d*a-e*M)
    (ha : S ∣ a) (hp : S ∣ n-1) : S ∣ e+h*r := by
  obtain ⟨A,hA⟩ := ha
  obtain ⟨B,hB⟩ := hp
  refine ⟨(d*r+e)*A-e*B, ?_⟩
  grind

/-- The finite search identity after clearing denominators. -/
theorem scale_identity (m k U d e r g T h0 beta : Int)
    (h : k*(U*m*(e+d*r)+beta) = T*(e*g*m+h0)) :
    m*(k*U*(e+d*r)-e*g*T) = T*h0-k*beta := by
  grind

/-- A vanishing scale remainder forces the exceptional ratio. -/
theorem zero_remainder (k U d e r g T h beta : Int)
    (hk : k ≠ 0) (hu : U ≠ 0) (hd : d ≠ 0)
    (h1 : k*U*(e+d*r) = e*g*T)
    (h2 : T*h = k*U*beta) (hb : g*beta = h-d) :
    e+h*r = 0 := by
  have hz : k*U*d*(e+h*r) = 0 := by
    calc
      k*U*d*(e+h*r) = h*(k*U*(e+d*r)) - k*U*e*(h-d) := by grind
      _ = e*g*(T*h) - k*U*e*(g*beta) := by rw [h1, hb]; grind
      _ = 0 := by rw [h2]; grind
  rcases Int.mul_eq_zero.mp hz with h | h
  · rcases Int.mul_eq_zero.mp h with h | h
    · exact False.elim ((Int.mul_ne_zero hk hu) h)
    · exact False.elim (hd h)
  · exact h

end Wasserman.PartGcdGap
