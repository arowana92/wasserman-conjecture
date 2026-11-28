import Lean

namespace Wasserman.ResidualPowerClasses

/-- A unit times an arbitrary square can give the specified inverse residue. -/
def squareInverse (c target m : Nat) : Bool :=
  (List.range m).any fun y => (c*y*y) % m == target

/-- Necessary data for a grouped part x*2^a*3^b. A and B are exponent parities.
The prime ell divides both the multiplier and the assigned weight factor. -/
def parityCandidates (D : Nat) : List (List Nat) := Id.run do
  let mut result := []
  for x in [1:D] do
    for ell in [5,7,11,13,17] do
      if x % ell == 0 then
        for A in [0:2] do
          for B in [0:2] do
            let P := 2^A
            let Q := 3^B
            let binary := (x*Q) &&& ((D-x)*Q) &&& 7
            let ternary := (x*P) % 3 + ((D-x)*P) % 3
            if binary == 0 && ternary < 3 &&
                squareInverse (x*P*Q) 1 D &&
                squareInverse (D*P*Q) (ell-1) ell then
              result := [D,x,ell,A,B] :: result
  return result

set_option maxRecDepth 8192 in
set_option maxHeartbeats 2000000 in
theorem unrestrictedCofactors :
    ([5,7,11,19].flatMap parityCandidates) = [] := by decide +kernel

/-- The only positive unordered possible normalized row taking 3 with sum 7
has the two required multiples of three equal to 3. -/
theorem companionAtSeven (i j k : Nat) (hi : 0 < i) (hj : 0 < j) (hk : 0 < k)
    (hs : i+j+k = 7) (di : 3 ∣ i) (dj : 3 ∣ j) :
    i = 3 ∧ j = 3 ∧ k = 1 := by omega

end Wasserman.ResidualPowerClasses
