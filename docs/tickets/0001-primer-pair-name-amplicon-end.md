# Primer pair name overstates amplicon end by the reverse-primer length

**Status**: Done
**Created**: 2026-06-15
**Fixed**: 2026-06-15
**Scope**: primer-design, pipeline naming

## Summary

The primer pair name `{short_tid}:{amplicon_start}-{amplicon_end}` had its
`amplicon_end` inflated by `reverse_len - 1`, so the span implied by the name
did not match the reported product size. Example, before fix:
`95388:646-794` with `Size 130` (a 130 bp amplicon cannot span 646–794).
After fix: `95388:646-775`.

## Root cause

`pipeline.py` computed `amplicon_end` as
`transcript_offset + reverse_start + reverse_len`, assuming Primer3's
`PRIMER_RIGHT` is the interior (3') end of the reverse primer.

In fact Primer3's `PRIMER_RIGHT` **is already the 5' (outermost) base of the
reverse primer** — the last template base the amplicon covers. Adding
`reverse_len` counts the reverse primer twice.

Empirically verified with `primer3-py 2.3.0`:

```
LEFT=[70,20]  RIGHT=[182,20]  PRIMER_PAIR_PRODUCT_SIZE=113
  right − left + 1 = 182 − 70 + 1 = 113   ✓ (product size)
  right + len − left = 132               ✗ (old name span)
```

So `product_size == reverse_start - forward_start + 1`, and the correct name
end is `transcript_offset + reverse_start + 1`.

## Fix

- `src/brimer_plast/pipeline.py`: extracted `make_pair_name()` which computes
  `amplicon_end = transcript_offset + reverse_start + 1`. Replaced the inline
  (buggy) computation with a call to it.
- `tests/test_pipeline.py`: added `TestMakePairName` asserting the name span
  equals `reverse_start - forward_start + 1` (product size).

## Verification

- Unit: `TestMakePairName` (3 tests) pass.
- End-to-end: re-ran Brimer-PLAST for `HLA-DRA` on
  `mEonSpe_hap2_dv_hom_wMT.fa` + the `Re.Hiller2A_*_withMT_sorted.gtf`
  annotation; all 10 filtered names now have `end - start + 1 == Size`.
  Verification run: `hla-dra_run2.log`.

## Reference

The archived decision ADR `docs/adr/0002-primer-pair-naming.md` documented
the buggy formula; a pointer to this ticket was added at its top without
altering its archived contents.