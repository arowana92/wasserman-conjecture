import Lean

namespace Wasserman.ThreeUnitLayers

/-- After clearing a common denominator, the low digits contribute
either zero or the distinguished power, and either zero or one. -/
def atPosition (d n X P : Int) : Prop :=
  ∃ e f : Int, (e = 0 ∨ e = d) ∧ (f = 0 ∨ f = d) ∧ X = (e - n)*P + f

private theorem commonSlope (c d X P Q R : Int)
    (pq : P ≠ Q) (pr : P ≠ R) (qr : Q ≠ R)
    (hp : X = c*P ∨ X = c*P + d)
    (hq : X = c*Q ∨ X = c*Q + d)
    (hr : X = c*R ∨ X = c*R + d) :
    c = 0 ∧ (X = 0 ∨ X = d) := by
  have hz : c = 0 := by
    apply Classical.byContradiction
    intro hc
    have hpq := mt (Int.eq_of_mul_eq_mul_left hc (b := P) (c := Q)) pq
    have hpr := mt (Int.eq_of_mul_eq_mul_left hc (b := P) (c := R)) pr
    have hqr := mt (Int.eq_of_mul_eq_mul_left hc (b := Q) (c := R)) qr
    omega
  exact ⟨hz, by simpa [hz] using hp⟩

private theorem locate (d n X P : Int)
    (hd : 0 < d) (hn : 0 ≤ n) (hnd : n ≤ d) (large : 2 ≤ P)
    (fits : atPosition d n X P) :
    (X ≤ d-2*n → X = (-n)*P ∨ X = (-n)*P + d) ∧
    (d-2*n < X → X = (d-n)*P ∨ X = (d-n)*P + d) := by
  obtain ⟨e, f, he, hf, h⟩ := fits
  have lo := Int.mul_le_mul_of_nonpos_left (show -n ≤ 0 by omega) large
  have hi := Int.mul_le_mul_of_nonneg_left large (show 0 ≤ d-n by omega)
  rcases he with h₀ | h₁ <;> subst e <;> rcases hf with h₀ | h₁ <;> subst f <;>
    (try simp only [Int.zero_sub, Int.add_zero] at h) <;> constructor <;> intro hx <;> omega

/-- Three distinct positions with a common rational high-part proportion
force that proportion to be zero or one and the residual offset to be zero or one. -/
theorem threePositions (d n X P Q R : Int)
    (hd : 0 < d) (hn : 0 ≤ n) (hnd : n ≤ d)
    (p2 : 2 ≤ P) (q2 : 2 ≤ Q) (r2 : 2 ≤ R)
    (pq : P ≠ Q) (pr : P ≠ R) (qr : Q ≠ R)
    (hp : atPosition d n X P) (hq : atPosition d n X Q) (hr : atPosition d n X R) :
    (n = 0 ∨ n = d) ∧ (X = 0 ∨ X = d) := by
  have a := locate d n X P hd hn hnd p2 hp
  have b := locate d n X Q hd hn hnd q2 hq
  have c := locate d n X R hd hn hnd r2 hr
  by_cases low : X ≤ d-2*n
  · have h := commonSlope (-n) d X P Q R pq pr qr (a.1 low) (b.1 low) (c.1 low)
    omega
  · have high : d-2*n < X := by omega
    have h := commonSlope (d-n) d X P Q R pq pr qr (a.2 high) (b.2 high) (c.2 high)
    omega

end Wasserman.ThreeUnitLayers
