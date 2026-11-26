import Lean

namespace Wasserman.CutWindows

/-- The strict product bound and the cross-row divisibilities force one entire
vector of discrepancies to vanish. This works for any index type. -/
theorem zero_row {ι κ : Type} (N : Int) (x : ι → Int) (y : κ → Int)
    (divides : ∀ i j, N ∣ x i * y j)
    (small : ∀ i j, (x i * y j).natAbs < N.natAbs) :
    (∀ i, x i = 0) ∨ (∀ j, y j = 0) := by
  classical
  by_cases hx : ∀ i, x i = 0
  · exact Or.inl hx
  · right
    intro j
    apply Classical.byContradiction
    intro hy
    apply hx
    intro i
    have prod := Int.eq_zero_of_dvd_of_natAbs_lt_natAbs (divides i j) (small i j)
    exact (Int.mul_eq_zero.mp prod).resolve_right hy

/-- In the second interval, two different powers at a gap other than d force
both residues to be exact. The possible shifts are zero and d. -/
theorem two_shifted_multiples (Q R d z x y : Int) (hd : 0 < d)
    (hQ : d < Q) (hR : d < R) (distinct : Q ≠ R)
    (gap₁ : Q-R ≠ d) (gap₂ : R-Q ≠ d)
    (shift₁ : x = 0 ∨ x = d) (shift₂ : y = 0 ∨ y = d)
    (div₁ : Q ∣ z-x) (div₂ : R ∣ z-y)
    (bound₁ : -(2*Q) < z-x ∧ z-x < 2*Q)
    (bound₂ : -(2*R) < z-y ∧ z-y < 2*R) :
    z = x ∧ z = y := by
  obtain ⟨m, hm⟩ := div₁
  obtain ⟨n, hn⟩ := div₂
  have qm : Q*(-2) < Q*m ∧ Q*m < Q*2 := by omega
  have rn : R*(-2) < R*n ∧ R*n < R*2 := by omega
  have qpos : 0 < Q := by omega
  have rpos : 0 < R := by omega
  have hm₁ := (Int.mul_lt_mul_left qpos).mp qm.1
  have hm₂ := (Int.mul_lt_mul_left qpos).mp qm.2
  have hn₁ := (Int.mul_lt_mul_left rpos).mp rn.1
  have hn₂ := (Int.mul_lt_mul_left rpos).mp rn.2
  have mc : m = -1 ∨ m = 0 ∨ m = 1 := by omega
  have nc : n = -1 ∨ n = 0 ∨ n = 1 := by omega
  rcases mc with rfl | rfl | rfl <;> rcases nc with rfl | rfl | rfl <;>
    simp at hm hn <;> omega

end Wasserman.CutWindows
