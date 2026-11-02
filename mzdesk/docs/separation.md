# Separation

Give a compound and a mobile-phase pH. Get the charge at that pH, logD and suitable separation modes (pKa from MolGpKa).

## Results

### pH behaviour

- **logD at pH**: logD = logP (Crippen) + log10(fraction of fully uncharged molecules). Ions are assumed not to partition and zwitterions count as charged, so logD is a lower bound at low or high pH.
- **Dominant charge state**: Most probable net charge at this pH, summed over the predicted sites (a zwitterion is net 0).
- **Fraction neutral**: Share of the molecules with no charge at all; a zwitterion counts as charged.

### logD and net charge vs pH

- Solid line, left axis: logD. Dashed line, right axis: net charge.
- Dotted line: the chosen pH.

### Ionizable sites on the structure

- Coral: acid site. Teal: base site. Violet: both.
- Label: atom and pKa.

### Separation mode

- Simple rules on logD, net charge, TPSA, MW and stereocentres (thresholds: mzdesk/platform_fit/separation_rules.py).
- Green: good. Amber: possible. Grey: poor or not needed. Hover or click a mode for its reason.

### Ionizable sites

- Site: element and atom number, counted from 1 as in the SMILES.
- Acid: loses H. Base: takes H.
- Base sites on amide, sulfonamide and pyrrole-type N are ignored (see Notes).

### Run notes

- Sites are treated as independent (Henderson–Hasselbalch per site), so successive ionizations of one group, such as phosphate pKa1 / pKa2, are not modelled.
- Base sites on amide, sulfonamide and pyrrole-type aromatic N are dropped as not basic. Other weakly basic sites (aniline N, for example) still count.
- Predicted pKa values carry errors of about one unit, sometimes more (validation/molgpka).
- For zwitterions the logD is likely strongly underestimated: they partition more than this model assumes.
- Ion-pairing or a low pH can help retention of a permanently charged compound.
- A site is more than 99 % one form outside pKa ± 2. When these windows cover the whole pH range, no pH keeps every site in one form and retention depends on pH however it is set.
- Near a site's pKa a small pH shift changes the charge and the retention, so buffer well or move away from these pKa values.
