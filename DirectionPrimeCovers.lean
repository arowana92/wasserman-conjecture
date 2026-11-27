import Lean

set_option maxRecDepth 16000
set_option maxHeartbeats 4000000

namespace Wasserman.DirectionCovers

/-- Each bin still has a common divisor from the corresponding position. -/
inductive Linked : List Nat → List Nat → Prop
  | nil : Linked [] []
  | cons {p b : Nat} {ps bs : List Nat} :
      1 < p → p ∣ b → Linked ps bs → Linked (p :: ps) (b :: bs)

def moves (n : Nat) : List Nat → List (List Nat)
  | [] => []
  | b :: bs => (Nat.gcd b n :: bs) :: (moves n bs).map (b :: ·)

def canPlace : List Nat → List Nat → Bool
  | [], _ => true
  | n :: ns, bins => (moves n bins).any fun next =>
      next.all (· != 1) && canPlace ns next

theorem bins_valid {ds bs : List Nat} (h : Linked ds bs) : bs.all (· != 1) = true := by
  induction h with
  | nil => rfl
  | @cons p b ps tail hp hd h ih =>
    have hb : b ≠ 1 := by
      intro he
      subst b
      have := Nat.le_of_dvd (by decide : 0 < 1) hd
      omega
    simp [List.all_cons, hb, ih]

theorem move_available {ds bs : List Nat} {n : Nat} (h : Linked ds bs)
    (covered : ∃ p ∈ ds, p ∣ n) : ∃ next ∈ moves n bs, Linked ds next := by
  induction h with
  | nil => simp at covered
  | @cons p b ps tail hp hd h ih =>
    obtain ⟨q, hq, hqn⟩ := covered
    rcases List.mem_cons.mp hq with first | rest
    · subst q
      exact ⟨_, List.mem_cons_self, .cons hp (Nat.dvd_gcd hd hqn) h⟩
    · obtain ⟨next, hn, hl⟩ := ih ⟨q, rest, hqn⟩
      refine ⟨b :: next, List.mem_cons_of_mem _ ?_, .cons hp hd hl⟩
      exact List.mem_map.mpr ⟨next, hn, rfl⟩

theorem cover_sound (ns ds bs : List Nat) (h : Linked ds bs)
    (covered : ∀ n ∈ ns, ∃ p ∈ ds, p ∣ n) : canPlace ns bs = true := by
  induction ns generalizing bs with
  | nil => rfl
  | cons n ns ih =>
    obtain ⟨next, hn, hl⟩ := move_available h (covered n List.mem_cons_self)
    apply List.any_eq_true.mpr
    exact ⟨next, hn, Bool.and_eq_true_iff.mpr
      ⟨bins_valid hl, ih next hl (fun m hm => covered m (List.mem_cons_of_mem _ hm))⟩⟩

theorem empty_bins (ds : List Nat) (h : ∀ p ∈ ds, 1 < p) :
    Linked ds (List.replicate ds.length 0) := by
  induction ds with
  | nil => exact .nil
  | cons p ps ih =>
    exact .cons (h p List.mem_cons_self) (Nat.dvd_zero p)
      (ih (fun q hq => h q (List.mem_cons_of_mem _ hq)))

abbrev Row := Nat × Nat × Nat
structure Witness where
  n : Nat
  twos : Nat
  x : Row
  y : Row

def rowSum (r : Row) := r.1 + r.2.1 + r.2.2
def primitive (r : Row) := Nat.gcd (Nat.gcd r.1 r.2.1) r.2.2 == 1
def profile (r : Row) := Nat.gcd (r.1*r.2.1) (Nat.gcd (r.1*r.2.2) (r.2.1*r.2.2))
def oddRow (m : Nat) (r : Row) :=
  (m*r.1 &&& m*r.2.1) == 0 && ((m*r.1 ||| m*r.2.1) &&& m*r.2.2) == 0
def oddEndpoint (w : Witness) := w.n*(w.n-1)/2^w.twos
def valid (w : Witness) :=
  let u := rowSum w.x
  let v := rowSum w.y
  primitive w.x && primitive w.y && Nat.gcd u v == 1 && w.n % (u*v) == 0 &&
    Nat.lcm (profile w.x) (profile w.y) % (w.n-1) == 0 &&
    2^w.twos * oddEndpoint w == w.n*(w.n-1) && oddEndpoint w % 2 == 1 &&
    max u v <= 100547

def parityWitnesses : List Witness := [
  ⟨4019768513, 6, (643,41152,58752), (2852,3151,33976)⟩,
  ⟨1681782257, 4, (1361,10912,21776), (7021,9086,33286)⟩,
  ⟨12648641, 6, (1,64,128), (12615,13395,39527)⟩,
  ⟨29290496, 12, (100,641,6410), (703,988,2405)⟩,
  ⟨536010749, 2, (665,2716,3686), (10387,23912,41548)⟩,
  ⟨7972256513, 8, (13061,27904,55808), (25973,26334,30074)⟩,
  ⟨5871960064, 16, (1442,6169,81988), (6183,27522,31831)⟩]

def otherWitnesses : List Witness := [
  ⟨3879389, 2, (212,464,1537), (491,631,631)⟩,
  ⟨5157437, 2, (316,456,1501), (551,859,859)⟩,
  ⟨3688253, 2, (3,743,743), (68,1168,1241)⟩,
  ⟨1403273, 3, (95,131,393), (104,824,1339)⟩,
  ⟨2385197, 2, (44,604,1661), (315,359,359)⟩,
  ⟨2440853, 2, (184,344,989), (375,617,617)⟩]

theorem witness_data :
    (parityWitnesses ++ otherWitnesses).all valid = true ∧
    (parityWitnesses.all fun w =>
      oddRow (w.n / rowSum w.x) w.x || oddRow (w.n / rowSum w.y) w.y) = true := by decide

theorem cover_lower (ws : List Witness) (k : Nat)
    (table : ∀ j : Fin k, canPlace (ws.map oddEndpoint) (List.replicate j.val 0) = false)
    (ds : List Nat) (h : ∀ p ∈ ds, 1 < p)
    (hit : ∀ m ∈ ws.map oddEndpoint, ∃ p ∈ ds, p ∣ m) : k ≤ ds.length := by
  by_cases enough : k ≤ ds.length
  · exact enough
  · have fail := table ⟨ds.length, by omega⟩
    have ok := cover_sound _ ds _ (empty_bins ds h) hit
    simp_all

theorem five_for_parity (ds : List Nat) (h : ∀ p ∈ ds, 1 < p)
    (hit : ∀ m ∈ parityWitnesses.map oddEndpoint, ∃ p ∈ ds, p ∣ m) : 5 ≤ ds.length :=
  cover_lower parityWitnesses 5 (by decide) ds h hit

theorem six_without_two (ds : List Nat) (h : ∀ p ∈ ds, 1 < p)
    (hit : ∀ m ∈ otherWitnesses.map oddEndpoint, ∃ p ∈ ds, p ∣ m) : 6 ≤ ds.length :=
  cover_lower otherWitnesses 6 (by decide) ds h hit

end Wasserman.DirectionCovers
