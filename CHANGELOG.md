# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

## 1.0.0 - 2026-10-20
### Breaking Changes
- None
### Added
### Fixed
- Neighbor-gap handling: when two ranges sit next to each other on the same 
  strand, the previous version was meant to split the gap between them at 
  the midpoint. That test compared char* pointers and could never be true, 
  so it never ran.  This version gives the whole gap to whichever copy 
  extends into it first, and limits the neighbor to where that copy stopped. 
  Extensions can now run longer than before (e.g. ce10-fam2 under Refiner 
  parameters extends 193 bp left where the previous version stopped at 131).
- Right-extension line wrapping in -cons: The previous version wrapped the 
  right-extension sequence at 80 columns, but measured from the start of 
  the left extension.  Lines in that block came out at uneven lengths. This 
  version wraps the block cleanly at 80 columns, and since Refiner strips 
  whitespace, its results are unaffected.
- anchor_range in -outtsv: The previous version left out the subsequence 
  offset and didn't swap the coordinates for minus-strand ranges, so the 
  TSV disagreed with its own stdout report.  This version writes the same 
  1-based anchor coordinates the stdout report shows.
### Removed
- None
### Known Issues
- When neighboring copies on the minus strand extend toward each other, the 
  second can reuse 2 bp the first already took.  Refiner merges overlapping 
  ranges, so in tandem arrays seeded on the minus strand several copies can 
  be merged into one sequence.  See TODO.md.


