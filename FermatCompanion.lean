import Lean

namespace Wasserman.FermatCompanion

abbrev Row := Int × Int × Int

private def live : List Row := (List.range 16).map fun i =>
  let w := (2 : Int) ^ (i+1)
  (16384/w, 32768/w, 16129/w)

private def digits : List Int := [0, 1]
private def units : List Row := [(0,0,0), (1,0,0), (0,1,0), (0,0,1)]

private def bounds (s : Row) : List (Int × Int) :=
  [(0,1), (s.1,16384), (s.2.1,32768), (s.2.2,16129),
   (s.1+s.2.1+s.2.2,65281)]

-- Each lower endpoint is strictly smaller than each upper endpoint.
private def feasible (s : Row) : Bool :=
  (bounds s).all fun l => (bounds s).all fun u =>
    decide (l.1 * u.2 < (u.1+1) * l.2)

private def previous (t : Row) (d : Int) (u : Row) : Row :=
  (2*t.1+u.1-d*16384, 2*t.2.1+u.2.1-d*32768, 2*t.2.2+u.2.2-d*16129)

private theorem closed : ∀ t ∈ live, ∀ d ∈ digits, ∀ u ∈ units,
    feasible (previous t d u) = true → previous t d u ∈ live := by decide

/-- A finite legal reading of the remaining input, padded with zeros to state zero. -/
inductive Finishes : Row → Prop where
  | halt : Finishes (0,0,0)
  | step (t : Row) (d : Int) (u : Row) (hd : d ∈ digits) (hu : u ∈ units)
      (ok : feasible (previous t d u) = true) (rest : Finishes t) :
      Finishes (previous t d u)

private theorem mustBeLive (s : Row) (h : Finishes s) : s ∈ live := by
  induction h with
  | halt => decide
  | step t d u hd hu ok _ ih => exact closed t ih d hd u hu ok

/-- No finite completion is carry-free after the required 14-bit prefix. -/
theorem prefixCannotFinish : ¬Finishes (257,514,253) := by
  intro h
  exact (by decide : (257,514,253) ∉ live) (mustBeLive _ h)

end Wasserman.FermatCompanion
