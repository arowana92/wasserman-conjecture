import Lean

namespace Wasserman.Ternary247

abbrev R := Nat × Nat × Nat

private def row (s : R) (w : Nat) : R :=
  (114*w+s.1, 117*w+s.2.1, 16*w+s.2.2)
private def quot (s : R) : R := (s.1/3,s.2.1/3,s.2.2/3)
private def layerSum (s : R) (m : Nat) : Nat := s.1%m+s.2.1%m+s.2.2%m
private def clear (s : R) : Prop := ∀ j : Nat, layerSum s (3^j) < 3^j

private theorem clearQuot (s : R) (h : clear s) : clear (quot s) := by
  intro j
  have h' := h (j+1)
  rw [Nat.pow_succ, Nat.mul_comm (3^j) 3] at h'
  unfold layerSum at h' ⊢
  simp only [quot, ← Nat.mod_mul_right_div_self]
  omega

private def trap : List R := [
  (38,39,5),(76,78,11),(76,78,10),(50,52,7),(63,65,9),(101,104,14),
  (25,26,3),(21,21,3),(7,7,1),(45,46,6),(83,85,11),(78,80,11),
  (15,15,2),(53,54,7),(64,65,9),(5,5,0),(43,44,6),(81,83,11),
  (93,96,13),(65,66,9),(31,32,4),(69,71,9),(107,110,15),(21,22,3),(23,23,3)]
private def initial : List R := [(0,0,1),(0,0,0)]

private theorem trapped : ∀ s ∈ trap, ∃ j : Fin 6, 3^j.val ≤ layerSum s (3^j.val) := by decide

private theorem closed : ∀ s ∈ trap, ∀ d : Fin 3,
    3 ≤ layerSum (row s d.val) 3 ∨ quot (row s d.val) ∈ trap := by decide

private theorem enter : ∀ s ∈ initial, ∀ d : Fin 3,
    d.val ≠ 0 → quot (row s d.val) ∈ trap := by decide

private theorem waitZero : ∀ s ∈ initial, quot (row s 0) ∈ initial := by decide

private theorem low (s : R) (n : Nat) : layerSum (row s n) 3 = layerSum (row s (n%3)) 3 := by
  unfold layerSum row
  dsimp
  simp only [Nat.add_mod _ s.1, Nat.add_mod _ s.2.1, Nat.add_mod _ s.2.2, Nat.mul_mod_mod]

private theorem high (s : R) (n : Nat) : quot (row s n) = row (quot (row s (n%3))) (n/3) := by
  apply Prod.ext
  · simp [quot, row]; omega
  apply Prod.ext <;> simp [quot, row] <;> omega

private theorem noClear (n : Nat) : ∀ s : R,
    s ∈ trap ∨ (s ∈ initial ∧ 0 < n) → ¬clear (row s n) := by
  induction n using Nat.strongRecOn with
  | ind n ih =>
    intro s active ok
    by_cases zero : n = 0
    · subst n
      have mem : s ∈ trap := active.elim id (fun h => by omega)
      obtain ⟨j, bad⟩ := trapped s mem
      have good := ok j.val
      have eqn : row s 0 = s := by rcases s with ⟨a,b,c⟩; simp [row]
      rw [eqn] at good
      omega
    · let d : Fin 3 := ⟨n%3, Nat.mod_lt _ (by decide)⟩
      have digit : layerSum (row s d.val) 3 < 3 := by
        have one := ok 1
        simpa only [Nat.pow_one, low s n] using one
      have next := clearQuot (row s n) ok
      rw [high] at next
      apply ih (n/3) (Nat.div_lt_self (by omega) (by decide)) _ _ next
      rcases active with hs | ⟨hs, _⟩
      · exact Or.inl ((closed s hs d).resolve_left (Nat.not_le_of_gt digit))
      · by_cases dz : d.val = 0
        · change n%3 = 0 at dz
          exact Or.inr ⟨by rw [dz]; exact waitZero s hs, by omega⟩
        · exact Or.inl (enter s hs d dz)

/-- Every positive scalar gives a ternary carry at a finite Legendre layer. -/
theorem forcedCarry (w : Nat) (pos : 0 < w) :
    ∃ j : Nat, 3^j ≤ (114*w)%(3^j) + (117*w)%(3^j) + (16*w+1)%(3^j) := by
  apply Classical.byContradiction
  intro fail
  apply noClear w (0,0,1) (Or.inr ⟨by decide, pos⟩)
  intro j
  have hj : ¬(3^j ≤ (114*w)%(3^j) + (117*w)%(3^j) + (16*w+1)%(3^j)) :=
    fun h => fail ⟨j,h⟩
  simpa [layerSum,row] using (Nat.lt_of_not_ge hj)

end Wasserman.Ternary247
