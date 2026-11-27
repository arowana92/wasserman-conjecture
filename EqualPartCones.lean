import Lean

namespace Wasserman.EqualPartCones

/-- At a largest cyclic rotation x, the representative x+y-z is too
small to be the sum of two copies of x and a third nonnegative block. -/
theorem maximumObstruction (M x y z w : Int)
    (hy : 0 < y) (hz : 0 < z) (byx : y ≤ x) (bzx : z ≤ x)
    (lo : 2*x ≤ w) (hi : w ≤ M) :
    ¬ M ∣ w - (x+y-z) := by
  rintro ⟨k,h⟩
  have pos : 0 < w-(x+y-z) := by omega
  have small : w-(x+y-z) < M := by omega
  by_cases hk : k ≤ 0
  · have : M*k ≤ 0 := Int.mul_nonpos_of_nonneg_of_nonpos (by omega) hk
    omega
  · have : M*1 ≤ M*k := Int.mul_le_mul_of_nonneg_left (by omega) (by omega)
    omega

/-- Produce a proposed rank table at elaboration time. Its correctness is
checked by the kernel from the resulting single integer literal below. -/
private def makeRanks : Nat := Id.run do
  let mut ranks := Array.replicate 9099 0
  for _ in [:72] do
    let old := ranks
    for s in [1:9099] do
      let a := old[s/3]!
      let b := if s%3 == 1 then old[6066+s/3]! else 0
      ranks := ranks.set! s (1 + max a b)
  let mut packed := 0
  for s in [:9099] do
    packed := packed + (ranks[s]! <<< (7*s))
  return packed

open Lean Elab Term in
elab "proposedConeRanks" : term => return mkNatLit makeRanks

private def table : Nat := proposedConeRanks
def rank (s : Nat) := (table >>> (7*s)) % 128

set_option maxRecDepth 100000
set_option maxHeartbeats 0

private theorem checked : (List.range 9099).all (fun s =>
    rank s ≤ 72 && (s == 0 || (rank (s/3) < rank s &&
      (s%3 != 1 || rank (6066+s/3) < rank s)))) = true := by decide

theorem rank_le (s : Nat) (bound : s < 9099) : rank s ≤ 72 := by
  have h := List.all_eq_true.mp checked s (List.mem_range.mpr bound)
  simpa using (Bool.and_eq_true_iff.mp h).1

theorem decrease (s d : Nat) (bound : s < 9099) (nz : 0 < s)
    (allowed : 2*d + (18197*d+s)%3 < 3) :
    rank ((18197*d+s)/3) < rank s := by
  have h := List.all_eq_true.mp checked s (List.mem_range.mpr bound)
  have h0 : (s == 0) = false := by simp; omega
  simp [h0] at h
  have hd : d = 0 ∨ d = 1 := by omega
  rcases hd with rfl | rfl
  · simpa using h.2.1
  · have hr : s % 3 = 1 := by omega
    have step : (18197*1+s)/3 = 6066+s/3 := by omega
    rw [step]
    simpa [hr] using h.2.2

end Wasserman.EqualPartCones
