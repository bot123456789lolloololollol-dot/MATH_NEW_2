/*
 * solver.js — Exact maximum size of a Sidon set in Z_n (strict convention).
 *
 * Definition (strict / "sums of any two elements" = OEIS A260999):
 *   A subset A of Z_n is Sidon iff all sums a+b (a,b in A, unordered, doubles
 *   included) are distinct in Z_n.  Equivalent formulation used by the search:
 *   all C(k,2) circular distances min((b-a) mod n, (a-b) mod n) are distinct,
 *   AND (for even n) no pair lies at distance n/2.  [See SPEC.md, Lemma 1/2.]
 *
 * Search: DFS over sets containing 0 as minimum element, elements chosen in
 * increasing order; bitmask-free boolean array of used "slots" (circular
 * distance values 1..floor(n/2)).  Deterministic: fixed candidate order.
 *
 * Completeness argument (see SPEC.md): every Sidon set has a translate that
 * contains 0 and every element of Z_n\{0} has a representative in 1..n-1;
 * the DFS enumerates ALL subsets {0} ∪ S with S ⊆ {1..n-1}, |S|=k-1,
 * in increasing order, testing exactly the defining condition on insertion.
 * Hence for each k the search is a complete decision procedure for
 * "exists k-element Sidon set in Z_n".
 *
 * Usage:
 *   node solver.js --from A --to B [--deadline-per-level-ms 240000]
 *                  [--out FILE.json] [--only n1,n2,...]
 * Output JSON per n: { n, ub, levels:[{k,status,nodes,ms,witness?}],
 *                      k_max (or null), witness, certified }
 * certified == true iff descending levels resolved (SAT or UNSAT, no timeout)
 * down to and including the first SAT level.
 */
'use strict';

function parseArgs() {
  const args = process.argv.slice(2);
  const opt = { from: 2, to: 64, deadlineMs: 240000, out: null, only: null };
  for (let i = 0; i < args.length; i++) {
    if (args[i] === '--from') opt.from = parseInt(args[++i]);
    else if (args[i] === '--to') opt.to = parseInt(args[++i]);
    else if (args[i] === '--deadline-per-level-ms') opt.deadlineMs = parseInt(args[++i]);
    else if (args[i] === '--out') opt.out = args[++i];
    else if (args[i] === '--only') opt.only = args[++i].split(',').map(s => parseInt(s));
  }
  return opt;
}

// Largest k with C(k,2) <= floor((n-1)/2).  Slot count is floor((n-1)/2):
// odd n -> distances 1..(n-1)/2 ; even n -> distances 1..n/2-1 (distance n/2 forbidden).
function upperBoundK(n) {
  const slots = Math.floor((n - 1) / 2);
  let k = 1;
  while (((k + 1) * k) / 2 <= slots) k++;  // advance while C(k+1,2) still fits
  return k;                                 // max k with C(k,2) <= slots
}

/*
 * Decision procedure: exists k-element Sidon set in Z_n ?
 * Returns { status: 'SAT'|'UNSAT'|'TIMEOUT', witness?, nodes }
 */
function solve(n, k, deadlineMs) {
  const halfN = Math.floor(n / 2);           // slot index of the self-inverse distance (even n)
  const used = new Uint8Array(halfN + 1);    // used[slot] for slot values 1..floor(n/2)
  const elems = new Int32Array(k);           // sorted chosen elements; elems[0] = 0
  elems[0] = 0;
  let nodes = 0;
  const deadline = Date.now() + deadlineMs;
  let timedOut = false;
  const witnessOut = [];

  function dfs(t, last) {
    // t elements already placed; need k-t more, each > last.
    nodes++;
    if (nodes % 4096 === 0 && Date.now() > deadline) { timedOut = true; throw 'timeout'; }
    if (t === k) {
      for (let i = 0; i < k; i++) witnessOut.push(elems[i]);
      return true;
    }
    // remaining r elements must fit into {last+1 .. n-1}
    if (n - 1 - last < k - t) return false; // completeness-safe counting prune
    for (let x = last + 1; x <= n - 1; x++) {
      // validity: circular distances x-a must be unused and not the n/2 slot
      let ok = true;
      const addedSlots = []; // slots we would mark for this candidate
      for (let i = 0; i < t && ok; i++) {
        const d = (x - elems[i]) % n;          // in 1..n-1 since 0<=elems[i]<x<n... elems[0]=0<x
        let dd;
        if (d > halfN) dd = n - d; else dd = d; // circular fold
        if (dd === halfN && n % 2 === 0) { ok = false; break; } // strict: forbid n/2 pairs
        if (used[dd]) { ok = false; break; }
        used[dd] = 1;
        addedSlots.push(dd);
      }
      if (ok) {
        elems[t] = x;
        if (dfs(t + 1, x)) return true;
      }
      for (const s of addedSlots) used[s] = 0;
      if (timedOut) throw 'timeout';
    }
    return false;
  }

  try {
    const ok = dfs(1, 0);
    return { status: ok ? 'SAT' : 'UNSAT', witness: ok ? witnessOut.slice() : null, nodes };
  } catch (e) {
    if (e === 'timeout') return { status: 'TIMEOUT', witness: null, nodes };
    throw e;
  }
}

function main() {
  const opt = parseArgs();
  const ns = opt.only ? opt.only : [];
  if (!opt.only) for (let n = opt.from; n <= opt.to; n++) ns.push(n);
  const results = [];
  for (const n of ns) {
    if (n < 2) continue;
    const t0 = Date.now();
    const rec = { n, ub: upperBoundK(n), levels: [], k_max: null, witness: null, certified: false };
    for (let k = rec.ub; k >= 1; k--) {
      const t1 = Date.now();
      const r = solve(n, k, opt.deadlineMs);
      const lvl = { k, status: r.status, nodes: r.nodes, ms: Date.now() - t1 };
      if (r.status === 'SAT') { lvl.witness = r.witness; rec.k_max = k; rec.witness = r.witness; }
      rec.levels.push(lvl);
      if (r.status === 'SAT' || r.status === 'TIMEOUT') break;
    }
    // certified iff every level above/at the answer resolved without timeout
    rec.certified = rec.k_max !== null &&
      rec.levels.every(l => l.status !== 'TIMEOUT');
    rec.total_ms = Date.now() - t0;
    results.push(rec);
    const lv = rec.levels.map(l => `${l.k}:${l.status}[${l.nodes}n/${l.ms}ms]`).join(' ');
    console.log(`n=${n} ub=${rec.ub} k_max=${rec.k_max} certified=${rec.certified} | ${lv}`);
    if (opt.out) require('fs').writeFileSync(opt.out, JSON.stringify(results, null, 1));
  }
  if (opt.out) require('fs').writeFileSync(opt.out, JSON.stringify(results, null, 1));
  console.log('DONE ' + results.length + ' cases');
}

main();
