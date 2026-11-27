import Lean

namespace Wasserman.FortyFourPowerRemainders

theorem heavy (k H T : Nat) (lower : 15*k ≤ H) (upper : H < 22*k)
    (parts : 15*k ≤ T) (period : 959512576 ∣ H+1-15*k) :
    15*959512576 ≤ 7*T := by
  have h := Nat.le_of_dvd (by omega : 0 < H+1-15*k) period
  omega

theorem light (k H : Nat) (lower : 2*k ≤ H) (upper : H < 3*k)
    (small : k < 959512576) (period : 959512576 ∣ 3*k-1-H) :
    H+1 = 3*k := by
  by_cases eq : H+1 = 3*k
  · exact eq
  · have hh := Nat.le_of_dvd (by omega : 0 < 3*k-1-H) period
    omega

private theorem crt_lower (k s w c d a b : Int) (pos : 0 < k)
    (per : 0 < s*w) (bound : d ≤ s*w) (bez : a*s+b*w = 1)
    (hs : s ∣ k) (hw : w ∣ k-c) (ds : s ∣ d) (dw : w ∣ d-c) :
    d ≤ k := by
  obtain ⟨u, hu⟩ := Int.dvd_sub hs ds
  have hwd := Int.dvd_sub hw dw
  have id : (k-c)-(d-c) = k-d := by omega
  rw [id] at hwd
  obtain ⟨v, hv⟩ := hwd
  have h₁ := congrArg (fun z : Int => z*(k-d)) bez
  have h₂ := congrArg (fun z : Int => a*s*z) hv
  have h₃ := congrArg (fun z : Int => b*w*z) hu
  have eqn : k-d = (s*w)*(a*v+b*u) := by grind
  by_cases neg : a*v+b*u ≤ -1
  · have le := Int.mul_le_mul_of_nonneg_left neg (by omega : 0 ≤ s*w)
    rw [Int.mul_neg_one] at le
    omega
  · have le := Int.mul_nonneg (by omega : 0 ≤ s*w) (by omega : 0 ≤ a*v+b*u)
    omega

theorem at_three (k : Int) (pos : 0 < k) (hs : 15625 ∣ k)
    (res : 59049 ∣ k-1 ∨ 59049 ∣ k-4) : 144906250 ≤ k := by
  rcases res with h | h
  · have := crt_lower k 15625 59049 1 497546875 31843 (-8426)
        pos (by decide) (by decide) (by decide) hs h (by decide) (by decide)
    omega
  · exact crt_lower k 15625 59049 4 144906250 31843 (-8426)
        pos (by decide) (by decide) (by decide) hs h (by decide) (by decide)

theorem at_five (k : Int) (pos : 0 < k) (hs : 59049 ∣ k)
    (res : 78125 ∣ k-1 ∨ 78125 ∣ k-2 ∨ 78125 ∣ k-6 ∨ 78125 ∣ k-7) :
    207734382 ≤ k := by
  rcases res with h | h | h | h
  · have := crt_lower k 59049 78125 1 1347734376 22824 (-17251)
        pos (by decide) (by decide) (by decide) hs h (by decide) (by decide)
    omega
  · have := crt_lower k 59049 78125 2 2695468752 22824 (-17251)
        pos (by decide) (by decide) (by decide) hs h (by decide) (by decide)
    omega
  · have := crt_lower k 59049 78125 6 3473203131 22824 (-17251)
        pos (by decide) (by decide) (by decide) hs h (by decide) (by decide)
    omega
  · exact crt_lower k 59049 78125 7 207734382 22824 (-17251)
        pos (by decide) (by decide) (by decide) hs h (by decide) (by decide)

theorem other_base (k : Int) (pos : 0 < k) (hs : 59049 ∣ k) (hw : 15625 ∣ k) :
    922640625 ≤ k := by
  exact crt_lower k 59049 15625 0 922640625 7199 (-27206)
      pos (by decide) (by decide) (by decide) hs (by simpa using hw) (by decide) (by decide)

end Wasserman.FortyFourPowerRemainders
