# Plain Language Summary

**What we studied:** Liver-chip devices are small plastic wells containing living liver cells, used to predict how fast the human liver will break down a new drug. You add a known amount of drug, then sample the surrounding fluid over hours or days to see how its concentration falls.

**The problem:** That concentration doesn't fall for just one reason. Four things happen at once, all shrinking the measured amount simultaneously:

1. Drug moves out of the fluid and into the cells (and some drifts back out)
2. Once inside the cells, some of it gets broken down by the liver's metabolism — this is the number researchers actually want
3. Some drug sticks permanently to the plastic walls of the device
4. The fluid itself slowly evaporates over a multi-day experiment, concentrating whatever drug is left even as some disappears

**What we found:** If you only watch the fluid's concentration fall and don't measure the other three processes separately, you cannot mathematically tell them apart. Watching that one curve gives you three independent pieces of information, but there are four unknown quantities tangled up inside it — like trying to solve for four unknowns with only three equations. A number that looks like "the metabolism rate" could actually be masking a completely different split between metabolism, plastic-sticking, and cell uptake. We proved this precisely: we identified the exact family of alternative explanations that all fit the same data equally well, and showed two very different combinations of numbers producing curves that agree to seven decimal places — genuinely indistinguishable from the data alone.

**The fix:** Run a separate control experiment with no cells in the well, just to measure how much drug sticks to the plastic on its own. Once that number is known independently, the ambiguity disappears and the remaining numbers become solvable. Some labs already run this control as routine good practice; our result shows it isn't just good practice — it's mathematically required for the reported numbers to mean what people assume they mean.

**A further complication:** Even after removing that ambiguity in principle, we found that realistic experiments — with limited data points and normal measurement noise — can still fail to pin down the metabolism rate reliably. This surprised us partway through the work: a "quick check" numerical method gave a confident-looking answer that turned out to be a calculation artifact, not real information in the data. Only a more careful, computationally expensive check revealed the true picture. We report this as a caution for the field: a method can look precise and still be wrong, especially near this kind of ambiguity.

**What this means in practice:** When reading a liver-chip clearance study, check whether the authors measured drug-sticking-to-plastic independently, and whether they verified — rather than assumed — that their reported number was actually well-determined by the data collected. We propose a short checklist of what such studies should report so readers can make that judgment for themselves.
